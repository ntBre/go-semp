package main

import (
	"bytes"
	"reflect"
	"testing"
)

func TestWritePBS(t *testing.T) {
	var buf bytes.Buffer
	WritePBS(&buf, "the name", []string{
		"input1",
		"input2",
	})
	got := buf.String()
	want := `#!/bin/bash
#SBATCH --job-name=semp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o the name.out
#SBATCH --no-requeue
#SBATCH --mem=1gb

/home/qc/bin/g16b01.sh input1.com
/home/qc/bin/g16b01.sh input2.com
`
	if got != want {
		t.Errorf("got\n%#+v, wanted\n%#+v\n", got, want)
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

func TestStat(t *testing.T) {
	tmp := STAT_CMD
	STAT_CMD = func() (string, []string) {
		return "cat", []string{
			"testfiles/squeue.dat",
		}
	}
	defer func() {
		STAT_CMD = tmp
	}()
	got := map[string]bool{
		"51009181": true,
		"51009182": true,
		"51009183": true,
		"51009184": true,
		"51009185": true,
		"51009191": true,
		"51009194": true,
		"51009203": true,
		"51009208": true,
		"51009210": true,
		"51009211": true,
	}
	Stat(&got)
	want := map[string]bool{
		"51009181": true,
		"51009182": true,
		"51009183": true,
		"51009184": true,
		"51009185": true,
		"51009191": true,
		"51009194": true,
		"51009203": true,
		"51009208": true,
		"51009210": true,
		"51009211": false,
	}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
