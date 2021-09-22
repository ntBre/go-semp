package main

import (
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestRMSD(t *testing.T) {
	a := mat.NewDense(3, 1, []float64{1, 2, 3})
	b := mat.NewDense(3, 1, []float64{4, 5, 6})
	got := Norm(a, b)
	want := 5.196152422706632
	if got != want {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
