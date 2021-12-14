package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

func RunGaussian(i int, dir string, names []string, geom []float64,
	params []Param) (job Job) {
	// make file
	basename := filepath.Join(dir, fmt.Sprintf("job.%05d", i))
	comfile := basename + ".com"
	pbsfile := basename + ".pbs"
	input, err := os.Create(comfile)
	defer input.Close()
	if err != nil {
		panic(err)
	}
	WriteCom(input, names, geom, params)
	pbs, err := os.Create(pbsfile)
	defer pbs.Close()
	if err != nil {
		panic(err)
	}
	// submit
	WritePBS(pbs, filepath.Base(basename), filepath.Base(comfile))
	jobid := Submit(pbsfile)
	// return Job
	return Job{
		Filename: basename,
		Index:    i,
		Jobid:    jobid,
	}
}

func ParseGaussian(filename string) (ret float64, err error) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		return ret, ErrFileNotFound
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		if strings.Contains(scanner.Text(), "SCF Done") {
			fields := strings.Fields(scanner.Text())
			ret, err = strconv.ParseFloat(fields[4], 64)
			if err != nil {
				fmt.Printf("error parsing %q\n", fields)
				panic(err)
			}
			return
		}
	}
	return ret, ErrEnergyNotFound
}
