package main

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
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
#SBATCH -o the name.pbs.out
#SBATCH --no-requeue
#SBATCH --mem=1gb

export LD_LIBRARY_PATH=/home/qc/mopac2016/

hostname
date
/home/qc/mopac2016/MOPAC2016.exe input1.mop
/home/qc/mopac2016/MOPAC2016.exe input2.mop
date
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

func TestResubmit(t *testing.T) {
	tmp := SUBMIT_CMD
	SUBMIT_CMD, _ = filepath.Abs("testfiles/scripts/sbatch")
	setup()
	defer func() {
		SUBMIT_CMD = tmp
		test_takedown()
	}()
	src, err := os.Open("testfiles/job.mop")
	defer src.Close()
	if err != nil {
		panic(err)
	}
	inp := filepath.Join("inp/job.mop")
	dst, err := os.Create(inp)
	if err != nil {
		panic(err)
	}
	defer dst.Close()
	io.Copy(dst, src)
	got := Resubmit(Job{
		Filename: "job",
		Jobid:    "00000",
		I:        1,
		J:        2,
		Coeff:    12,
	})
	want := Job{
		Filename: "job_redo",
		I:        1,
		J:        2,
		Coeff:    12,
		Jobid:    "12345678",
	}
	if !reflect.DeepEqual(got, want) {
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
