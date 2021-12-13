package main

import (
	"bytes"
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

func TestWritePBS(t *testing.T) {
	var buf bytes.Buffer
	WritePBS(&buf, "the name", "input.com")
	got := buf.String()
	want := `#!/bin/bash
#SBATCH --job-name=the name
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o the name.pbs.out
#SBATCH --no-requeue
#SBATCH --mem=1gb

/home/qc/bin/g16b01.sh input.com
`
	if got != want {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
