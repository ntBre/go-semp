package main

import (
	"fmt"
	"math"
	"reflect"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestLoadConfig(t *testing.T) {
	got := LoadConfig("testfiles/test.in")
	want := Config{
		GeomFile:   "file07",
		EnergyFile: "rel.dat",
		Atoms:      []string{"C", "C", "C", "H", "H"},
		Params: []Param{
			{"USS", "H", -11.246958000000},
			{"ZS", "H", 1.268641000000},
			{"BETAS", "H", -8.352984000000},
			{"GSS", "H", 14.448686000000},
			{"USS", "C", -51.089653000000},
			{"UPP", "C", -39.937920000000},
			{"ZS", "C", 2.047558000000},
			{"ZP", "C", 1.702841000000},
			{"BETAS", "C", -15.385236000000},
			{"BETAP", "C", -7.471929000000},
			{"GSS", "C", 13.335519000000},
			{"GPP", "C", 10.778326000000},
			{"GSP", "C", 11.528134000000},
			{"GP2", "C", 9.486212000000},
			{"HSP", "C", 0.717322000000},
		},
		MaxIt:  250,
		Lambda: 1e-8,
	}
	if !compParams(got.Params, want.Params, 1e-12) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
	if got.GeomFile != want.GeomFile {
		t.Errorf("GeomFile: got %v, wanted %v\n", got, want)
	}
	if got.EnergyFile != want.EnergyFile {
		t.Errorf("EnergyFile: got %v, wanted %v\n", got, want)
	}
	if !reflect.DeepEqual(got.Atoms, want.Atoms) {
		t.Errorf("Atoms: got %v, wanted %v\n", got, want)
	}
	if got.MaxIt != want.MaxIt {
		t.Errorf("MaxIt: got %v, wanted %v\n", got, want)
	}
	if got.Lambda != want.Lambda {
		t.Errorf("Lambda: got %v, wanted %v\n", got, want)
	}
}

func TestLoadGeoms(t *testing.T) {
	got := LoadGeoms("testfiles/three07")
	want := [][]float64{
		{
			0.0000000000, 0.0000000000, -1.6794733900,
			0.0000000000, 1.2524327590, 0.6959098120,
			0.0000000000, -1.2524327590, 0.6959098120,
			0.0000000000, 3.0146272390, 1.7138963510,
			0.0000000000, -3.0146272390, 1.7138963510,
		},
		{
			0.0000000000, 0.0000000000, -1.6929508795,
			0.0000000000, 1.2335354991, 0.6923003326,
			0.0000000000, -1.2335354991, 0.6923003326,
			0.0000000000, 2.9875928126, 1.7242445752,
			0.0000000000, -2.9875928126, 1.7242445752,
		},
		{
			0.0000000000, 0.0000000000, -1.6826636598,
			0.0000000000, 1.2382598141, 0.6926064201,
			0.0000000000, -1.2382598141, 0.6926064201,
			0.0000000000, 2.9956906742, 1.7187948778,
			0.0000000000, -2.9956906742, 1.7187948778,
		},
	}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestLoadEnergies(t *testing.T) {
	got := LoadEnergies("testfiles/min.dat")
	want := mat.NewDense(6, 1, []float64{
		0.000000000000,
		0.000453458157,
		0.000258906149,
		0.000268926371,
		0.000247426244,
		0.000251632099,
	})
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func compParams(a, b []Param, eps float64) bool {
	for i := range a {
		if a[i].Atom != b[i].Atom {
			fmt.Printf("atom: %q vs %q\n", a[i], b[i])
			return false
		}
		if a[i].Name != b[i].Name {
			fmt.Printf("name: %q vs %q\n", a[i], b[i])
			return false
		}
		if diff := math.Abs(a[i].Value - b[i].Value); diff > eps {
			fmt.Printf(
				"value: %q vs %q with difference %e\n",
				a[i], b[i], diff,
			)
			return false
		}
	}
	return true
}

func TestLoadParams(t *testing.T) {
	rc := LoadConfig("testfiles/test.in")
	got := rc.Params
	want := []Param{
		{"USS", "H", -11.246958000000},
		{"ZS", "H", 1.268641000000},
		{"BETAS", "H", -8.352984000000},
		{"GSS", "H", 14.448686000000},
		{"USS", "C", -51.089653000000},
		{"UPP", "C", -39.937920000000},
		{"ZS", "C", 2.047558000000},
		{"ZP", "C", 1.702841000000},
		{"BETAS", "C", -15.385236000000},
		{"BETAP", "C", -7.471929000000},
		{"GSS", "C", 13.335519000000},
		{"GPP", "C", 10.778326000000},
		{"GSP", "C", 11.528134000000},
		{"GP2", "C", 9.486212000000},
		{"HSP", "C", 0.717322000000},
	}
	if !compParams(got, want, 1e-12) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
