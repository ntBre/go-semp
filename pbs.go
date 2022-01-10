package main

import (
	"embed"
	"fmt"
	"io"
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

func Submit(filename string) string {
	dir := filepath.Dir(filename)
	base := filepath.Base(filename)
	cmd := exec.Command(SUBMIT_CMD, base)
	cmd.Dir = dir
	byts, err := cmd.Output()
	if err != nil {
		fmt.Println("the error is : ", cmd.String())
		panic(err)
	}
	// output like "Submitted batch job 49229449"
	fields := strings.Fields(string(byts))
	return fields[3]
}
