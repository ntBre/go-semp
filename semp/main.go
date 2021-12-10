package main

import (
	"bufio"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"runtime/pprof"
	"strconv"
	"strings"
	"text/template"
	"time"

	"gonum.org/v1/gonum/mat"
)

const (
	// from http://www.ilpi.com/msds/ref/energyunits.html
	htToCm = 219_474.5459784
	EPS    = 1e-14
	THRESH = 1.0
	NU     = 2.0
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
		// "GCore":  {},
	}
	CHARGE = 0
	SPIN   = 1
	//  https://en.wikipedia.org/wiki/Numerical_differentiation
	//  recommends cube root of machine eps (~2.2e16) for step
	//  size => 6e-6; adjusting down from there
	DELTA   = 1e-8
	LOGFILE io.Writer
)

/// Flags
var (
	debug      = flag.Bool("debug", false, "toggle debugging information")
	cpuprofile = flag.String("cpu", "", "write a CPU profile")
	ncpus      = flag.Int("ncpus", 8, "number of cpus to use")
	gauss      = flag.String("gauss", "g16", "command to run gaussian")
	lambda     = flag.Float64("lambda", 0.0, "initial lambda value for levmar")
	maxit      = flag.Int("maxit", 100, "maximum iterations")
	one        = flag.String("one", "", "run one iteration using the "+
		"params in params.dat and save the results in the argument")
)

/// Errors
var (
	ErrEnergyNotFound = errors.New("Energy not found in Gaussian output")
	ErrFileNotFound   = errors.New("Output file not found")
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

func WriteParams(params []Param, f io.Writer) {
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

func LogParams(w io.Writer, params []Param, iter int) {
	fmt.Fprintf(w, "Iter %5d\n", iter)
	WriteParams(params, w)
}

func DumpParams(params []Param, filename string) {
	f, err := os.Create(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	WriteParams(params, f)
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
	path, err := filepath.Abs(paramfile)
	if err != nil {
		panic(err)
	}
	anon := struct {
		Charge int
		Spin   int
		Geom   string
		Params string
	}{CHARGE, SPIN, geom, path}
	if *debug {
		f, _ := os.Create(fmt.Sprintf("debug/file%d.com", counter))
		counter++
		t.Execute(f, anon)
	}
	t.Execute(w, anon)
}

type Job struct {
	Filename string
	Index    int
	Jobid    string
}

func RunGaussian(i int, dir string, names []string, geom []float64,
	paramfile string) (job Job) {
	// make file
	basename := fmt.Sprintf("job.%05d", i)
	comfile := basename + ".com"
	pbsfile := basename + ".pbs"
	input, err := os.Create(comfile)
	defer input.Close()
	if err != nil {
		panic(err)
	}
	WriteCom(input, names, geom, paramfile)
	pbs, err := os.Create(pbsfile)
	defer pbs.Close()
	if err != nil {
		panic(err)
	}
	// submit
	WritePBS(pbs, basename, comfile)
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

// SEnergy returns the relative semi-empirical energies in Ht
// corresponding to params
func SEnergy(dir string, names []string, geoms [][]float64, paramfile string) *mat.Dense {
	ret := make([]float64, len(geoms))
	jobs := make([]Job, len(geoms))
	for i, geom := range geoms {
		jobs[i] = RunGaussian(i, dir, names, geom, paramfile)
	}
	for len(jobs) > 0 {
		for i := 0; i < len(jobs); i++ {
			job := jobs[i]
			energy, err := ParseGaussian(job.Filename + ".out")
			if err == nil {
				ret[job.Index] = energy
				l := len(jobs) - 1
				jobs[l], jobs[i] = jobs[i], jobs[l]
				jobs = jobs[:l]
			}
		}
		time.Sleep(1 * time.Second)
	}
	return mat.NewDense(len(ret), 1, ret)
}

func CentralDiff(names []string, geoms [][]float64, params []Param, p, i int) []float64 {
	params[p].Values[i] += DELTA
	DumpParams(params, "params.dat")
	forward := SEnergy(".", names, geoms, "params.dat")

	params[p].Values[i] -= 2 * DELTA
	DumpParams(params, "params.dat")
	backward := SEnergy(".", names, geoms, "params.dat")

	// f'(x) ~ [ f(x+h) - f(x-h) ] / 2h ; where h is DELTA here
	var diff mat.Dense
	diff.Sub(forward, backward)
	return Scale(1/(2*DELTA), diff.RawMatrix().Data)
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

func RMSD(a, b *mat.Dense) (ret float64) {
	as := a.RawMatrix().Data
	bs := b.RawMatrix().Data
	if len(as) != len(bs) {
		panic("dimension mismatch")
	}
	var count int
	for i := range as {
		// deviation
		diff := as[i] - bs[i]
		// square
		ret += diff * diff
		count++
	}
	// mean
	ret /= float64(count)
	// root
	return math.Sqrt(ret)
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

func Identity(n int) *mat.Dense {
	ret := mat.NewDense(n, n, nil)
	for i := 0; i < n; i++ {
		ret.Set(i, i, 1.0)
	}
	return ret
}

func UpdateParams(params []Param, v *mat.Dense, scale float64) []Param {
	ret := make([]Param, 0, len(params))
	var i int
	for _, p := range params {
		vals := make([]float64, 0, len(p.Values))
		for _, val := range p.Values {
			vals = append(vals, val+v.At(i, 0)*scale)
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

func LevMar(jac, ai, se *mat.Dense, params []Param, scale float64) []Param {
	// LHS
	var prod mat.Dense
	prod.Mul(jac.T(), jac)
	r, _ := prod.Dims()
	eye := Identity(r)
	var Leye mat.Dense
	Leye.Scale(*lambda, eye)
	var sum mat.Dense
	sum.Add(&prod, &Leye)
	// RHS
	var diff mat.Dense
	diff.Sub(ai, se)
	var prod2 mat.Dense
	prod2.Mul(jac.T(), &diff)
	var step mat.Dense
	err := step.Solve(&sum, &prod2)
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
	return UpdateParams(params, &step, scale)
}

func OneIter(labels []string, geoms [][]float64, paramfile, outfile string) {
	se := Relative(SEnergy(".", labels, geoms, paramfile))
	out, _ := os.Create(outfile)
	for _, s := range se.RawMatrix().Data {
		fmt.Fprintf(out, "%20.12f\n", s)
	}
}

func main() {
	host, _ := os.Hostname()
	flag.Parse()
	fmt.Printf("running with %d cpus on host: %s\n", *ncpus, host)
	fmt.Printf("initial lambda: %.14f\n", *lambda)
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
	labels := []string{"C", "C", "C", "H", "H"}
	geoms := LoadGeoms("file07")
	if *one != "" {
		OneIter(labels, geoms, "params.dat", *one)
		os.Exit(0)
	}
	LOGFILE, _ = os.Create("log")
	paramLog, _ := os.Create("params.log")
	ai := LoadEnergies("rel.dat")
	params, num := LoadParams("opt.out")
	fmt.Printf("loaded %d params\n", num)
	DumpParams(params, "params.dat")
	baseEnergies := SEnergy(".", labels, geoms, "params.dat")
	se := Relative(baseEnergies)
	norm := Norm(ai, se) * htToCm
	rmsd := RMSD(ai, se) * htToCm
	var (
		iter     int
		lastNorm float64
		lastRMSD float64
	)
	fmt.Printf("%17s%12s%12s%12s%12s\n",
		"cm-1", "cm-1", "cm-1", "cm-1", "s")
	fmt.Printf("%5s%12s%12s%12s%12s%12s\n",
		"Iter", "Norm", "ΔNorm", "RMSD", "ΔRMSD", "Time")
	fmt.Printf("%5d%12.4f%12.4f%12.4f%12.4f%12.1f\n",
		iter, norm, norm-lastNorm, rmsd, rmsd-lastRMSD, 0.0)
	LogParams(paramLog, params, iter)
	iter++
	lastNorm = norm
	lastRMSD = rmsd
	start := time.Now()
	for iter < *maxit && norm > THRESH {
		jac := NumJac(labels, geoms, params, baseEnergies)
		*lambda /= NU
		// BEGIN copy-paste
		newParams := LevMar(jac, ai, se, params, 1.0)
		DumpParams(newParams, "params.dat")
		se = SEnergy(".", labels, geoms, "params.dat")
		se = Relative(se)
		norm = Norm(ai, se) * htToCm
		rmsd = RMSD(ai, se) * htToCm
		// END copy-paste

		// case ii. and iii. of levmar, first case is ii. from
		// Marquardt63
		var bad bool
		for i := 0; norm > lastNorm; i++ {
			*lambda *= NU
			newParams := LevMar(jac, ai, se, params, 1.0)
			DumpParams(newParams, "params.dat")
			se = SEnergy(".", labels, geoms, "params.dat")
			se = Relative(se)
			norm = Norm(ai, se) * htToCm
			rmsd = RMSD(ai, se) * htToCm
			fmt.Fprintf(LOGFILE,
				"\tλ_%d to %g\n", i, *lambda)
			// give up after 5 increases
			if i > 4 {
				// case iii. failed, try footnote
				bad = true
				*lambda *= math.Pow(NU, float64(-(i + 1)))
				break
			}
		}
		var prev float64
		// just break if decreasing k doesnt change the norm
		for i := 2; bad && norm > lastNorm && norm-prev > 1e-6; i++ {
			k := 1.0 / float64(int(1)<<i)
			newParams := LevMar(jac, ai, se, params, k)
			DumpParams(newParams, "params.dat")
			se = SEnergy(".", labels, geoms, "params.dat")
			se = Relative(se)
			norm = Norm(ai, se) * htToCm
			rmsd = RMSD(ai, se) * htToCm
			fmt.Fprintf(LOGFILE, "\tk_%d to %g with Δ = %f\n",
				i, k, norm-lastNorm)
			prev = norm
		}
		fmt.Printf("%5d%12.4f%12.4f%12.4f%12.4f%12.1f\n",
			iter, norm, norm-lastNorm, rmsd, rmsd-lastRMSD,
			float64(time.Since(start))/1e9)
		start = time.Now()
		lastNorm = norm
		lastRMSD = rmsd
		params = newParams
		LogParams(paramLog, params, iter)
		iter++
	}
}
