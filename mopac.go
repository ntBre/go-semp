package main

import (
	"bufio"
	"errors"
	"os"
	"strconv"
	"strings"
)

const (
	KCALHT = 627.5091809 // kcal/mol per hartree
)

func ReadOut(filename string) (energy float64, err error) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		err = ErrFileNotFound
		return
	}
	scanner := bufio.NewScanner(f)
	err = ErrEnergyNotFound
	var (
		line   string
		fields []string
		i      int
	)
	for i = 0; scanner.Scan(); i++ {
		line = scanner.Text()
		switch {
		case i == 0 && strings.Contains(strings.ToUpper(line), "PANIC"):
			panic("panic requested in output file")
		case i == 0 && strings.Contains(strings.ToUpper(line), "ERROR"):
			err = errors.New("file contains error")
			return
		case strings.Contains(line, "== MOPAC DONE =="):
			break
		}
	}
	if i == 0 {
		err = errors.New("blank output")
		return
	}
	// should I close old f first? what about deferring double
	// close?
	auxfile := TrimExt(filename) + ".aux"
	f, err = os.Open(auxfile)
	defer f.Close()
	if err != nil {
		err = ErrFileNotFound
		return
	}
	err = ErrEnergyNotFound
	scanner = bufio.NewScanner(f)
	for scanner.Scan() {
		line = scanner.Text()
		switch {
		case strings.Contains(line, "HEAT_OF_FORMATION"):
			fields = strings.Split(line, "=")
			strVal := fields[1]
			energy, err = strconv.ParseFloat(
				strings.Replace(strVal, "D", "E", -1),
				64,
			)
			energy /= KCALHT
		}
	}
	return
}
