package main

import (
	"bytes"
	"testing"
)

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

func TestSubmit(t *testing.T) {
	tmp := SUBMIT_CMD
	SUBMIT_CMD = "scripts/sbatch"
	defer func() {
		SUBMIT_CMD = tmp
	}()
	got := Submit("testfiles/file")
	want := "12345678"
	if got != want {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
