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
	delete(DERIVED_PARAMS, "GCore")
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
						"F0ss", "ZetaOverlap", "U", "Beta", "CoreKO",
						"GCore", "GCore", "GCore",
					},
					[]float64{
						0.5309794405, 1.268641, -0.4133181015,
						-0.3069665138, 0.9416560451, 0.0016794859,
						0.8557539975, 3.3750716455,
					},
				},
				{
					"C",
					[]string{"F0ss", "F0sp", "F0pp", "F2pp", "G1sp",
						"ZetaOverlap", "ZetaOverlap", "U", "U",
						"Beta", "Beta", "CoreKO", "GCore", "GCore",
						"GCore",
					},
					[]float64{
						0.490071306, 0.4236511293, 0.3644399818,
						0.1978513158, 0.0790832954, 2.047558,
						1.702841, -1.8775102017, -1.4676915546,
						-0.5653970197, -0.2745883383, 1.0202596926,
						0.003215496, 0.588117579, 2.5208171714,
					},
				},
			},
		},
		{
			infile: "testfiles/params.dat",
			want: []Param{
				{
					"H",
					[]string{
						"F0ss", "ZetaOverlap", "U", "Beta", "CoreKO",
						"GCore", "GCore", "GCore",
					},
					[]float64{
						0.5309794405, 1.268641, -0.4133181015,
						-0.3069665138, 0.9416560451, 0.0016794859,
						0.8557539975, 3.3750716455,
					},
				},
				{
					"C",
					[]string{"F0ss", "F0sp", "F0pp", "F2pp", "G1sp",
						"ZetaOverlap", "ZetaOverlap", "U", "U",
						"Beta", "Beta", "CoreKO", "GCore", "GCore",
						"GCore",
					},
					[]float64{
						0.490071306, 0.4236511293, 0.3644399818,
						0.1978513158, 0.0790832954, 2.047558,
						1.702841, -1.8775102017, -1.4676915546,
						-0.5653970197, -0.2745883383, 1.0202596926,
						0.003215496, 0.588117579, 2.5208171714,
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
