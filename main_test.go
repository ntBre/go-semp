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

func test_takedown() {
	takedown()
	os.RemoveAll("inp")
}

func TestMain(t *testing.T) {
	if !testing.Short() {
		t.Skip()
	}
	gauss := LoadEnergies("testfiles/gauss.nrg.out")
	mopac := LoadEnergies("testfiles/mopac.nrg.out")
	labels := []string{"C", "C", "C", "H", "H"}
	geoms := LoadGeoms("testfiles/file07")
	setup()
	defer test_takedown()
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

func compMat(a, b *mat.Dense, eps float64) bool {
	var diff mat.Dense
	diff.Sub(a, b)
	return mat.Norm(&diff, 2) < eps
}

func TestNumJac(t *testing.T) {
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
		test_takedown()
	}()
	setup()
	got := NumJac(
		[]string{"C", "C", "C", "H", "H"},
		LoadGeoms("testfiles/three07"),
		LoadParams(LoadConfig("testfiles/test.in").Params),
	)
	want := mat.NewDense(3, 15, []float64{
		-0.01497188, 0.20599558, 0.04540870,
		0.00756935, -0.07643601, 0.09140744,
		0.60444017, 1.26628479, 0.06613604,
		0.14448405, -0.06807515, 0.00109442,
		-0.05038643, 0.25973557, -0.04459850,
		-0.01499188, 0.20528582, 0.04541566,
		0.00754239, -0.07667791, 0.09167065,
		0.60234381, 1.26036888, 0.06667187,
		0.14533904, -0.06778789, -0.00037060,
		-0.05114822, 0.26220733, -0.04444806,
		-0.01500130, 0.20545525, 0.04541131,
		0.00754616, -0.07681023, 0.09181168,
		0.60351258, 1.26152403, 0.06675347,
		0.14539861, -0.06789877, -0.00041641,
		-0.05106358, 0.26227226, -0.04426819,
	})
	if !compMat(got, want, 1e-7) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
