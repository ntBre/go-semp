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
	Name  string
	Input string
}

func WritePBS(w io.Writer, name, infile string) {
	// FIXME probably don't need to do this on every iteration
	f, err := template.ParseFS(Templates, "pbs.tmpl")
	if err != nil {
		panic(err)
	}
	f.Execute(w, PBS{
		Name:  name,
		Input: infile,
	})
}

func Submit(filename string) string {
	dir := filepath.Dir(filename)
	base := filepath.Base(filename)
	cmd := exec.Command("sbatch", base)
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
