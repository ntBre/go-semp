package main

import (
	"embed"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"text/template"
)

//go:embed pbs.tmpl
var Templates embed.FS

type PBS struct {
	Name   string
	Inputs []string
}

var PBS_TEMPLATE *template.Template

func init() {
	var err error
	PBS_TEMPLATE, err = template.ParseFS(Templates, "pbs.tmpl")
	if err != nil {
		panic(err)
	}
}

func WritePBS(w io.Writer, name string, infiles []string) {
	PBS_TEMPLATE.Execute(w, PBS{
		Name:   filepath.Base(name),
		Inputs: infiles,
	})
}

var SUBMIT_CMD string = "sbatch"

// Submit sends filename to the queue. The directory for the
// submission command is taken from the filename, so the full path
// needs to be present.
func Submit(filename string) string {
	dir := filepath.Dir(filename)
	base := filepath.Base(filename)
	cmd := exec.Command(SUBMIT_CMD, base)
	cmd.Dir = dir
	byts, err := cmd.Output()
	if err != nil {
		fmt.Printf("error on %q is: %v",
			cmd.String(), err,
		)
		panic(err)
	}
	// output like "Submitted batch job 49229449"
	fields := strings.Fields(string(byts))
	return fields[3]
}

func Resubmit(job Job) Job {
	src, err := os.Open(filepath.Join("inp", job.Filename+INFILE_SUFFIX))
	defer src.Close()
	if err != nil {
		panic(err)
	}
	inp := filepath.Join("inp", job.Filename+"_redo"+INFILE_SUFFIX)
	dst, err := os.Create(inp)
	if err != nil {
		panic(err)
	}
	defer dst.Close()
	io.Copy(dst, src)
	pbs := filepath.Join("inp", job.Filename+"_redo.pbs")
	f, err := os.Create(pbs)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	WritePBS(f, job.Filename, []string{job.Filename + "_redo"})
	return Job{
		Filename: job.Filename + "_redo",
		I:        job.I,
		J:        job.J,
		Jobid:    Submit(pbs),
		Coeff:    job.Coeff,
	}
}
