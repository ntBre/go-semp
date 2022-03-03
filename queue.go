package main

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
)

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
	fmt.Fprintf(os.Stderr, "resubmitting %s as %s\n",
		job.Filename, job.Filename+"_redo")
	return Job{
		Filename: job.Filename + "_redo",
		I:        job.I,
		J:        job.J,
		Jobid:    Submit(pbs),
		Coeff:    job.Coeff,
	}
}
