package main

import (
	"embed"
	"io"
	"os/exec"
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
	cmd, err := exec.Command("sbatch", filename).Output()
	if err != nil {
		panic(err)
	}
	// output like "Submitted batch job 49229449"
	fields := strings.Fields(string(cmd))
	return fields[3]
}
