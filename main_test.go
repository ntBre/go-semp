package main

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"regexp"
	"testing"
	"text/template"

	"gonum.org/v1/gonum/mat"
)

func TestWriteParams(t *testing.T) {
	var b bytes.Buffer
	p := LoadConfig("testfiles/test.in").Params
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

func TestEnergies(t *testing.T) {
	if !testing.Short() {
		t.Skip()
	}
	gauss := LoadEnergies("testfiles/gauss.nrg.out")
	mopac := LoadEnergies("testfiles/mopac.nrg.out")
	labels := []string{"C", "C", "C", "H", "H"}
	geoms := LoadGeoms("testfiles/file07")
	setup()
	defer test_takedown()
	params := LoadConfig("testfiles/test.in").Params
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

func TestWork(t *testing.T) {
	if !testing.Short() {
		t.Skip()
	}
	setup()
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
	got := work("testfiles/work.in")
	want := []float64{447.1502, 42.6167, 0.9004}
	if !compFloat(got, want, 1e-4) {
		t.Errorf("got %v, wanted %v\n", got, want)
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

var TESTJAC *mat.Dense

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
		LoadConfig("testfiles/test.in").Params,
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
	eps := 1e-7
	if !compMat(got, want, eps) {
		g := got.RawMatrix().Data
		w := want.RawMatrix().Data
		for i := range g {
			if diff := g[i] - w[i]; math.Abs(diff) > eps {
				fmt.Printf("%5d%20.12g\n", i, diff)
			}
		}
		t.Errorf("got %v, wanted %v\n", got, want)
	}
	TESTJAC = got
}

var (
	TESTSTEP   *mat.Dense
	TESTPARAMS []Param
)

func TestLevMar(t *testing.T) {
	TESTPARAMS, _, TESTSTEP = LevMar(TESTJAC,
		LoadEnergies("testfiles/three.dat"),
		LoadEnergies("testfiles/three.nrg.dat"),
		LoadConfig("testfiles/test.in").Params,
		1.0, 1e-8,
	)
	want := []Param{
		{"USS", "H", -10.367883047437},
		{"ZS", "H", 1.215268220252},
		{"BETAS", "H", -8.739427051094},
		{"GSS", "H", 12.724944577073},
		{"USS", "C", -51.004152272936},
		{"UPP", "C", -40.021839956047},
		{"ZS", "C", 2.041117221309},
		{"ZP", "C", 1.694803991853},
		{"BETAS", "C", -15.546321646344},
		{"BETAP", "C", -7.578605930451},
		{"GSS", "C", 13.419037030316},
		{"GPP", "C", 11.232884059306},
		{"GSP", "C", 12.127694620927},
		{"GP2", "C", 9.420331106671},
		{"HSP", "C", 1.662217856082},
	}
	if !compParams(TESTPARAMS, want, 1e-12) {
		t.Errorf("got %v, wanted %v\n", TESTPARAMS, want)
	}
}

func TestBroyden(t *testing.T) {
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
	geoms := LoadGeoms("testfiles/three07")
	atoms := []string{"C", "C", "C", "H", "H"}
	params := LoadConfig("testfiles/test.in").Params
	nrg := mat.NewDense(len(geoms), 1, nil)
	jobs := SEnergy(atoms, geoms, params, 0, None)
	RunJobs(jobs, nrg)
	seOld := Relative(nrg)
	jobs = SEnergy(atoms, geoms, TESTPARAMS, 0, None)
	nrg.Zero()
	RunJobs(jobs, nrg)
	seNew := Relative(nrg)
	got := Broyden(TESTJAC, TESTSTEP, seOld, seNew)
	r, c := TESTJAC.Dims()
	want := mat.NewDense(r, c, []float64{
		0.06197344, 0.20132387, 0.01158339,
		-0.14330951, -0.06895214, 0.08406194,
		0.60387641, 1.26558131, 0.05203623,
		0.13514664, -0.06076483, 0.04088183,
		0.00209303, 0.25396902, 0.03810811,
		0.06192401, 0.20061590, 0.01160329,
		-0.14327877, -0.06919690, 0.08432796,
		0.60178026, 1.25966567, 0.05257746,
		0.13600519, -0.06048036, 0.03940159,
		0.00131118, 0.25644298, 0.03822692,
		0.06188915, 0.20078687, 0.01161012,
		-0.14322513, -0.06933170, 0.08447141,
		0.60294922, 1.26082106, 0.05266372,
		0.13606785, -0.06059366, 0.03934263,
		0.00137847, 0.25650983, 0.03837945,
	})
	if !compMat(got, want, 1e-7) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
