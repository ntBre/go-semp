package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
	"text/template"

	"gonum.org/v1/gonum/mat"
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
	//  https://en.wikipedia.org/wiki/Numerical_differentiation
	//  recommends cube root of machine eps (~2.2e16) for step
	//  size
	DELTA = 6e-6
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

func Len(params []Param) (sum int) {
	for _, p := range params {
		sum += len(p.Names)
	}
	return sum
}

// SEnergy returns the relative semi-empirical energies in Ht
// corresponding to params
func SEnergy(names []string, geoms [][]float64, paramfile string) []float64 {
	// TODO
	return make([]float64, len(geoms))
}

func Sub(a, b []float64) []float64 {
	ret := make([]float64, len(a))
	for i := range a {
		ret[i] = a[i] - b[i]
	}
	return ret
}

func Scale(s float64, f []float64) []float64 {
	ret := make([]float64, len(f))
	for i := range f {
		ret[i] = s * f[i]
	}
	return ret
}

// NumJac computes the numerical Jacobian of energies vs params
func NumJac(names []string, geoms [][]float64, params []Param) *mat.Dense {
	rows := len(geoms)
	cols := Len(params)
	jac := mat.NewDense(rows, cols, nil)
	var col int
	for p := range params {
		for i := range params[p].Values {
			params[p].Values[i] += DELTA
			DumpParams(params, "params.dat")
			forward := SEnergy(names, geoms, "params.dat")

			params[p].Values[i] -= 2 * DELTA
			DumpParams(params, "params.dat")
			backward := SEnergy(names, geoms, "params.dat")

			// f'(x) ~ [ f(x+h) - f(x-h) ] / 2h ; where h
			// is DELTA here
			jac.SetCol(col, Scale(1/(2*DELTA), Sub(forward, backward)))

			// reset and move to next column
			params[p].Values[i] += DELTA
			col++
		}
	}
	return jac
}

func main() {
	labels := []string{"C", "C", "C", "H", "H"}
	fmt.Println(labels)
	geoms := LoadGeoms("file07")
	LoadEnergies("rel.dat")
	params := LoadParams("opt.out")
	DumpParams(params, "params.dat")
	WriteCom(labels, geoms[0], "params.dat", "out")
	NumJac(labels, geoms, params)
}
