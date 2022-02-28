package main

import (
	"reflect"
	"testing"

	"gonum.org/v1/gonum/mat"
)

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

func TestLoadParams(t *testing.T) {
	tests := []struct {
		infile string
		want   []Param
	}{
		{
			infile: "testfiles/opt.out",
			want: []Param{
				{
					"H",
					[]string{
						"USS",
						"ZS",
						"BETAS",
						"GSS",
					},
					[]float64{
						-11.24695800,
						1.26864100,
						-8.35298400,
						14.44868600,
					},
				},
				{
					"C",
					[]string{
						"USS",
						"UPP",
						"ZS",
						"ZP",
						"BETAS",
						"BETAP",
						"GSS",
						"GPP",
						"GSP",
						"GP2",
						"HSP",
					},
					[]float64{
						-51.08965300,
						-39.93792000,
						2.04755800,
						1.70284100,
						-15.38523600,
						-7.47192900,
						13.33551900,
						10.77832600,
						11.52813400,
						9.48621200,
						0.71732200,
					},
				},
			},
		},
	}
	for _, test := range tests {
		got, _ := LoadParams(test.infile)
		if !reflect.DeepEqual(got, test.want) {
			t.Errorf("got %v, wanted %v\n", got, test.want)
		}
	}
}
