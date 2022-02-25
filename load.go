package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"regexp"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/mat"
)

func LoadGeoms(filename string) (ret [][]float64) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	var line string
	row := make([]float64, 0, 3)
	scanner.Scan() // discard first "# GEOM" line
	for scanner.Scan() {
		line = scanner.Text()
		if strings.Contains(line, "#") {
			ret = append(ret, row)
			row = make([]float64, 0, len(row))
		} else {
			fields := strings.Fields(line)
			row = append(row, toFloat(fields)...)
		}
	}
	ret = append(ret, row)
	return
}

// LoadEnergies loads relative energies from filename and returns them
// as a column vector
func LoadEnergies(filename string) *mat.Dense {
	ret := make([]float64, 0)
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	var line string
	for scanner.Scan() {
		line = scanner.Text()
		fields := strings.Fields(line)
		if len(fields) == 1 {
			v, err := strconv.ParseFloat(fields[0], 64)
			if err != nil {
				panic(err)
			}
			ret = append(ret, v)
		}
	}
	return mat.NewDense(len(ret), 1, ret)
}

// LoadParams extracts semi-empirical parameters from a MOPAC output
// file
func LoadParams(filename string) (ret []Param, num int) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	var (
		line    string
		fields  []string
		first   bool
		inparam bool
		param   Param
		skip    int
	)
	blank := regexp.MustCompile(`^ *$`)
	gap := regexp.MustCompile(`(ALPB|XFAC)_ ?(\d+)`)
	// TODO might need to deduplicate symmetrical terms like 1
	// XFAC_6 vs 6 XFAC_1
	for scanner.Scan() {
		line = scanner.Text()
		fields = strings.Fields(line)
		if gap.MatchString(line) {
			atom_num := gap.FindStringSubmatch(line)[2]
			repl := fmt.Sprintf("${1}_%s", ATOMIC_NUMBERS[atom_num])
			line = gap.ReplaceAllString(line, repl)
			fields = strings.Fields(line)
		}
		switch {
		case skip > 0:
			skip--
		case strings.Contains(line, "MATRIX FROM HCORE"):
			return
		case strings.Contains(line, "PARAMETER VALUES USED"):
			skip = 4
			inparam = true
			first = true
			param = Param{}
		case inparam && blank.MatchString(line):
			if param.Atom != "" {
				ret = append(ret, param)
			}
			first = true
			param = Param{}
		case inparam && first:
			param.Atom = ATOMIC_NUMBERS[fields[0]]
			first = false
			fallthrough
		case inparam:
			if _, ok := DERIVED_PARAMS[fields[1]]; ok {
				continue
			}
			param.Names = append(param.Names, fields[1])
			v, err := strconv.ParseFloat(fields[2], 64)
			if err != nil {
				log.Fatalf("failed to parse %q as float",
					fields[2],
				)
			}
			param.Values = append(param.Values, v)
			num++
		}
	}
	return
}
