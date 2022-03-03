package main

import (
	"bytes"
	"fmt"
	"os"
	"regexp"
	"testing"
	"text/template"

	"gonum.org/v1/gonum/mat"
)

func charComp(got, want string) {
	for c := range got {
		if len(got) <= c {
			fmt.Println("got too short")
			return
		} else if len(want) <= c {
			fmt.Println("want too short")
			return
		}
		if got[c] != want[c] {
			fmt.Printf("got\n%q, wanted\n%q\n",
				got[:c+1], want[:c+1])
			return
		}
	}
}

func TestWriteParams(t *testing.T) {
	var b bytes.Buffer
	p := LoadParams(LoadConfig("testfiles/test.in").Params)
	WriteParams(&b, p)
	got := b.String()
	want := `USS            H    -11.246958000000
ZS             H      1.268641000000
BETAS          H     -8.352984000000
GSS            H     14.448686000000
USS            C    -51.089653000000
UPP            C    -39.937920000000
ZS             C      2.047558000000
ZP             C      1.702841000000
BETAS          C    -15.385236000000
BETAP          C     -7.471929000000
GSS            C     13.335519000000
GPP            C     10.778326000000
GSP            C     11.528134000000
GP2            C      9.486212000000
HSP            C      0.717322000000

`
	if got != want {
		charComp(got, want)
		t.Errorf("got\n%#v, wanted\n%#v\n", got, want)
	}
}

func TestWriteCom(t *testing.T) {
	var b bytes.Buffer
	p := LoadParams("testfiles/opt.out")
	os.MkdirAll("tmparam", 0744)
	defer os.RemoveAll("tmparam")
	WriteCom(&b, []string{"C", "C", "C", "H", "H"},
		[]float64{
			0.0000000000, 0.0000000000, -1.6794733900,
			0.0000000000, 1.2524327590, 0.6959098120,
			0.0000000000, -1.2524327590, 0.6959098120,
			0.0000000000, 3.0146272390, 1.7138963510,
			0.0000000000, -3.0146272390, 1.7138963510,
		}, p)
	got := b.String()
	want := `XYZ A0 scfcrt=1.D-21 aux(precision=14) external=tmparam/1598109253 1SCF charge=0 PM6
blank line
blank line
C      0.000000000000      0.000000000000     -1.679473390000
C      0.000000000000      1.252432759000      0.695909812000
C      0.000000000000     -1.252432759000      0.695909812000
H      0.000000000000      3.014627239000      1.713896351000
H      0.000000000000     -3.014627239000      1.713896351000


`
	random := regexp.MustCompile(`(external=)/?.*tmparam/[^ ]+`)
	got = random.ReplaceAllString(got, "$1")
	want = random.ReplaceAllString(want, "$1")
	if got != want {
		charComp(got, want)
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestMain(t *testing.T) {
	if !testing.Short() {
		t.Skip()
	}
	gauss := LoadEnergies("testfiles/gauss.nrg.out")
	mopac := LoadEnergies("testfiles/mopac.nrg.out")
	labels := []string{"C", "C", "C", "H", "H"}
	geoms := LoadGeoms("testfiles/file07")
	os.RemoveAll("inp")
	os.Mkdir("inp", 0755)
	defer os.RemoveAll("inp")
	os.MkdirAll("tmparam", 0744)
	defer os.RemoveAll("tmparam")
	params := LoadParams("testfiles/opt.out")
	got := mat.NewDense(len(geoms), 1, nil)
	jobs := SEnergy(labels, geoms, params, 0, None)
	tmp := PBS_TEMPLATE
	var err error
	PBS_TEMPLATE, err = template.ParseFiles("testfiles/local.tmpl")
	if err != nil {
		panic(err)
	}
	cmd := SUBMIT_CMD
	SUBMIT_CMD = "bash"
	defer func() {
		PBS_TEMPLATE = tmp
		SUBMIT_CMD = cmd
	}()
	RunJobs(jobs, got)
	if norm, fail := vecNorm(got, gauss, 7e-3); fail {
		t.Errorf("gaussian mismatch with norm %.8e\n", norm)
		vecDiff(got, gauss)
	}
	if norm, fail := vecNorm(got, mopac, 1e-10); fail {
		t.Errorf("mopac mismatch with norm %.8e\n", norm)
		vecDiff(got, mopac)
	}
}

// return the norm of the difference between got and want and whether or not it
// is greater than eps
func vecNorm(got, want *mat.Dense, eps float64) (float64, bool) {
	var diff mat.Dense
	diff.Sub(got, want)
	norm := mat.Norm(&diff, 2)
	return norm, norm > eps
}

// print the difference between column vectors got and want
func vecDiff(got, want *mat.Dense) {
	l, _ := got.Dims()
	fmt.Printf("\n%20s%20s%20s\n", "Got", "Want", "Diff")
	for i := 0; i < l; i++ {
		g := got.At(i, 0)
		w := want.At(i, 0)
		fmt.Printf("%20.12f%20.12f%20.12f\n",
			g, w, g-w,
		)
	}
}
