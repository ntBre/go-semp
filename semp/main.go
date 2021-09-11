package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
	"text/template"
)

var (
	// set representing derived semi-empirical parameters
	DERIVED_PARAMS = map[string]struct{}{
		"PQN":      {},
		"NValence": {},
		"DDN":      {},
		"KON":      {},
		"EISol":    {},
	}
	CHARGE = 0
	SPIN   = 1
)

func toFloat(strs []string) []float64 {
	ret := make([]float64, len(strs))
	var err error
	for i, s := range strs {
		ret[i], err = strconv.ParseFloat(s, 64)
		if err != nil {
			panic(err)
		}
	}
	return ret
}

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

func LoadEnergies(filename string) (ret []float64) {
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
	return
}

type Param struct {
	Atom   string
	Names  []string
	Values []float64
}

// LoadParams extracts semi-empirical parameters from a Gaussian
// output file
func LoadParams(filename string) (ret []Param) {
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
	)
	for scanner.Scan() {
		line = scanner.Text()
		fields = strings.Fields(line)
		switch {
		case line == " ****":
			inparam = true
			first = true
			if param.Atom != "" {
				ret = append(ret, param)
			}
			param = Param{}
		case inparam && line == " ":
			inparam = false
			return
		case inparam && first:
			param.Atom = fields[0]
			first = false
		case inparam:
			for _, f := range fields {
				split := strings.Split(f, "=")
				if _, ok := DERIVED_PARAMS[split[0]]; ok {
					continue
				}
				name := split[0]
				vals := strings.Split(split[1], ",")
				if name == "DCore" {
					name += "=" + vals[0] + "," + vals[1]
					vals = vals[2:]
				}
				for _, val := range vals {
					param.Names = append(param.Names, name)
					v, err := strconv.ParseFloat(val, 64)
					if err != nil {
						panic(err)
					}
					param.Values = append(param.Values, v)
				}
			}
		}
	}
	return
}

func DumpParams(params []Param, filename string) {
	f, err := os.Create(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	for _, param := range params {
		fmt.Fprintf(f, " ****\n%2s", param.Atom)
		var last, sep string
		for i, name := range param.Names {
			if strings.Contains(name, "DCore") {
				sep = ","
			} else {
				sep = "="
			}
			if name == last {
				fmt.Fprintf(f, ",%.10f", param.Values[i])
			} else {
				fmt.Fprintf(f, "\n%s%s%.10f",
					name, sep, param.Values[i])
				last = name
			}
		}
		fmt.Fprint(f, "\n")
	}
	fmt.Fprintln(f, " ****")
}

func WriteCom(names []string, coords []float64, paramfile, filename string) {
	f, err := os.Create(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	var geom strings.Builder
	for i := range names {
		fmt.Fprintf(&geom, "%s%20.12f%20.12f%20.12f\n",
			names[i], coords[i], coords[i+1], coords[i+2],
		)
	}

	t, err := template.New("com").Parse(`%mem=1000mb
%nprocs=1
#P PM6=(print,zero)

the title

{{.Charge}} {{.Spin}}
{{.Geom}}
@{{.Params}}

`)
	if err != nil {
		panic(err)
	}
	t.Execute(f, struct {
		Charge int
		Spin   int
		Geom   string
		Params string
	}{CHARGE, SPIN, geom.String(), paramfile})
}

func main() {
	labels := []string{"C", "C", "C", "H", "H"}
	fmt.Println(labels)
	geoms := LoadGeoms("file07")
	LoadEnergies("rel.dat")
	params := LoadParams("opt.out")
	DumpParams(params, "params.dat")
	WriteCom(labels, geoms[0], "params.dat", "out")
}
