package main

import (
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestRelative(t *testing.T) {
	got := Relative(mat.NewDense(3, 1, []float64{1, 2, 3}))
	want := mat.NewDense(3, 1, []float64{0, 1, 2})
	n, _ := Norm(got, want)
	if  n != 0 {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestNorm(t *testing.T) {
	a := mat.NewDense(3, 1, []float64{1, 2, 3})
	b := mat.NewDense(3, 1, []float64{4, 5, 6})
	got, _ := Norm(a, b)
	want := 5.196152422706632 * htToCm
	if got != want {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestRMSD(t *testing.T) {
	a := mat.NewDense(3, 1, []float64{1, 2, 3})
	b := mat.NewDense(3, 1, []float64{4, 5, 6})
	got := RMSD(a, b)
	want := 3.0
	if got != want {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
