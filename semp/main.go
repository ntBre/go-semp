package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"sync"
	"text/template"

	"os/exec"

	"gonum.org/v1/gonum/mat"
)

const (
	// from http://www.ilpi.com/msds/ref/energyunits.html
	htToCm = 219_474.5459784
	EPS    = 1e-14
)

var (
	// set representing derived semi-empirical parameters
	DERIVED_PARAMS = map[string]struct{}{
		"PQN":      {},
		"NValence": {},
		"DDN":      {},
		"KON":      {},
		"EISol":    {},
		// Adding these trying to make matrix not singular
		"DCore":  {},
		"EHeat":  {},
		"DipHyp": {},
		"GCore":  {},
	}
	CHARGE = 0
	SPIN   = 1
	//  https://en.wikipedia.org/wiki/Numerical_differentiation
	//  recommends cube root of machine eps (~2.2e16) for step
	//  size => 6e-6; adjusting down from there
	DELTA   = 1e-8
	LOGFILE io.Writer
)

// Flags
var (
	debug      = flag.Bool("debug", false, "toggle debugging information")
	cpuprofile = flag.String("cpu", "", "write a CPU profile")
	ncpus      = flag.Int("ncpus", 8, "number of cpus to use")
)

type Param struct {
	Atom   string
	Names  []string
	Values []float64
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

func Equal(a, b float64) bool {
	if math.Abs(a-b) > EPS {
		return false
	}
	return true
}

// LoadParams extracts semi-empirical parameters from a Gaussian
// output file
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
					v, err := strconv.ParseFloat(val, 64)
					if err != nil {
						panic(err)
					}
					// skip zero params for now
					if v != 0.0 {
						param.Names = append(param.Names, name)
						param.Values = append(param.Values, v)
						num++
					}
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
	fmt.Fprint(f, " ****\n\n")
}

var (
	counter int
)

func WriteCom(w io.Writer, names []string, coords []float64, paramfile string) {
	geom := ZipGeom(names, coords)
	t, err := template.New("com").Parse(`%mem=1000mb
%nprocs=1
#P PM6=(print,zero,input)

the title

{{.Charge}} {{.Spin}}
{{.Geom}}
@{{.Params}}

`)
	if err != nil {
		panic(err)
	}
	anon := struct {
		Charge int
		Spin   int
		Geom   string
		Params string
	}{CHARGE, SPIN, geom, paramfile}
	if *debug {
		f, _ := os.Create(fmt.Sprintf("debug/file%d.com", counter))
		counter++
		t.Execute(f, anon)
	}
	t.Execute(w, anon)
}

func RunGaussian(dir string, names []string,
	geom []float64, paramfile string) (ret float64) {
	var input bytes.Buffer
	WriteCom(&input, names, geom, paramfile)
	cmd := exec.Command("g16")
	cmd.Dir = dir
	cmd.Stdin = &input
	var r io.Reader
	r, _ = cmd.StdoutPipe()
	if *debug {
		outfile, _ := os.Create(fmt.Sprintf("debug/output%d.com", counter-1))
		r = io.TeeReader(r, outfile)
	}
	err := cmd.Start()
	if err != nil {
		fmt.Printf("error running %s:\n", cmd.String())
		panic(err)
	}
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		if strings.Contains(scanner.Text(), "SCF Done") {
			fields := strings.Fields(scanner.Text())
			// could return here, but need to shut down
			// command
			ret, err = strconv.ParseFloat(fields[4], 64)
			if err != nil {
				fmt.Printf("error parsing %q\n", fields)
				panic(err)
			}
		}
	}
	return
}

// SEnergy returns the relative semi-empirical energies in Ht
// corresponding to params
func PLSEnergy(dir string, names []string, geoms [][]float64, paramfile string) *mat.Dense {
	// parallel version, make sure the normal version works first
	ret := make([]float64, len(geoms))
	sema := make(chan struct{}, *ncpus)
	var wg sync.WaitGroup
	for i, geom := range geoms {
		sema <- struct{}{}
		wg.Add(1)
		go func(dir string, names []string,
			geom []float64, paramfile string, i int) {
			defer func() {
				<-sema
				wg.Done()
			}()
			ret[i] = RunGaussian(dir, names, geom, paramfile)
		}(dir, names, geom, paramfile, i)
	}
	wg.Wait()
	return mat.NewDense(len(ret), 1, ret)
}

// SEnergy returns the relative semi-empirical energies in Ht
// corresponding to params
func SEnergy(dir string, names []string, geoms [][]float64, paramfile string) []float64 {
	ret := make([]float64, len(geoms))
	for i, geom := range geoms {
		ret[i] = RunGaussian(dir, names, geom, paramfile)
		// hate to check this on every iteration
		if *debug {
			fmt.Printf("%5d%20.12f\n", i, ret[i])
		}
	}
	return ret
}

func CentralDiff(names []string, geoms [][]float64, params []Param, p, i int) []float64 {
	params[p].Values[i] += DELTA
	DumpParams(params, "params.dat")
	forward := SEnergy(".", names, geoms, "params.dat")

	params[p].Values[i] -= 2 * DELTA
	DumpParams(params, "params.dat")
	backward := SEnergy(".", names, geoms, "params.dat")

	// f'(x) ~ [ f(x+h) - f(x-h) ] / 2h ; where h is DELTA here
	return Scale(1/(2*DELTA), Sub(forward, backward))
}

func ForwardDiff(names []string, geoms [][]float64,
	params []Param, p, i int, energies *mat.Dense) []float64 {
	base := energies.RawMatrix().Data
	params[p].Values[i] += DELTA
	DumpParams(params, "params.dat")
	forward := SEnergy(".", names, geoms, "params.dat")

	// f'(x) ~ [ f(x+h) - f(x) ] / h ; where h is DELTA here
	return Scale(1/DELTA, Sub(forward, base))
}

// NumJac computes the numerical Jacobian of energies vs params
func NumJac(names []string, geoms [][]float64,
	params []Param, energies *mat.Dense) *mat.Dense {
	rows := len(geoms)
	cols := Len(params)
	jac := mat.NewDense(rows, cols, nil)
	var col int
	for p := range params {
		for i := range params[p].Values {
			data := CentralDiff(names, geoms, params, p, i)
			// data := ForwardDiff(names, geoms, params, p, i, energies)
			jac.SetCol(col, data)
			fmt.Fprintf(LOGFILE, "finished col %5d -> %s of %s\n", col,
				params[p].Names[i], params[p].Atom)
			col++
		}
	}
	return jac
}

func Relative(a *mat.Dense) *mat.Dense {
	min := mat.Min(a)
	r, c := a.Dims()
	if c != 1 {
		panic("too many columns")
	}
	ret := mat.NewDense(r, c, nil)
	for i := 0; i < r; i++ {
		ret.Set(i, 0, a.At(i, 0)-min)
	}
	return ret
}

func Norm(a, b *mat.Dense) float64 {
	var diff mat.Dense
	diff.Sub(a, b)
	return mat.Norm(&diff, 2)
}

func DumpVec(a *mat.Dense) {
	r, c := a.Dims()
	if c != 1 {
		panic("more than one column in expected vector")
	}
	for i := 0; i < r; i++ {
		fmt.Printf("%5d%20.12f\n", i, a.At(i, 0))
	}
}

func DumpMat(m mat.Matrix) {
	r, c := m.Dims()
	for i := 0; i < r; i++ {
		fmt.Printf("%5d", i)
		for j := 0; j < c; j++ {
			fmt.Printf("%12.8f", m.At(i, j))
		}
		fmt.Print("\n")
	}
	fmt.Print("\n")
}

func UpdateParams(params []Param, v *mat.Dense) []Param {
	ret := make([]Param, 0, len(params))
	var i int
	for _, p := range params {
		vals := make([]float64, 0, len(p.Values))
		for _, val := range p.Values {
			vals = append(vals, val+v.At(i, 0))
			i++
		}
		ret = append(ret, Param{
			Atom:   p.Atom,
			Names:  p.Names,
			Values: vals,
		})
	}
	return ret
}

func main() {
	flag.Parse()
	if *debug {
		os.Mkdir("debug", 0744)
	}
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			panic(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	LOGFILE, _ = os.Create("log")
	labels := []string{"C", "C", "C", "H", "H"}
	geoms := LoadGeoms("file07")
	ai := LoadEnergies("rel.dat")
	params, num := LoadParams("opt.out")
	fmt.Printf("loaded %d params\n", num)
	DumpParams(params, "params.dat")
	// BEGIN initial Norm
	baseEnergies := PLSEnergy(".", labels, geoms, "params.dat")
	energies := Relative(baseEnergies)
	var iter int
	fmt.Printf("Norm %5d: %20.4f cm-1\n", iter, Norm(ai, energies)*htToCm)
	iter++
	// END initial Norm

	// BEGIN first step
	for i := 0; i < 100; i++ {
		jac := NumJac(labels, geoms, params, baseEnergies)
		var prod mat.Dense
		prod.Mul(jac.T(), jac)
		var diff mat.Dense
		diff.Sub(ai, energies)
		var prod2 mat.Dense
		prod2.Mul(jac.T(), &diff)
		var step mat.Dense
		err := step.Solve(&prod, &prod2)
		if err != nil {
			fmt.Println("jacT")
			DumpMat(jac.T())
			fmt.Println("prod")
			DumpMat(&prod)
			fmt.Println("diff")
			DumpVec(&diff)
			fmt.Println("prod2")
			DumpMat(&prod2)
			fmt.Println("step")
			DumpVec(&step)
			panic(err)
		}
		newParams := UpdateParams(params, &step)
		DumpParams(newParams, "params.dat")
		energies = PLSEnergy(".", labels, geoms, "params.dat")
		energies = Relative(energies)
		fmt.Printf("Norm %5d: %20.4f cm-1\n", iter, Norm(ai, energies)*htToCm)
		iter++
		params = newParams
	}
}
