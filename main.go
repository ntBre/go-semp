package main

import (
	"bufio"
	"bytes"
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

func WriteCom(w io.Writer, names []string, coords []float64, params []Param) {
	geom := ZipGeom(names, coords)
	t, err := template.New("com").Parse(`%mem=1000mb
%nprocs=1
#P PM6=(print,zero,input)

the title

{{.Charge}} {{.Spin}}
{{.Geom}}
{{.Params}}

`)
	if err != nil {
		panic(err)
	}
	var b bytes.Buffer
	WriteParams(params, &b)
	anon := struct {
		Charge int
		Spin   int
		Geom   string
		Params string
	}{CHARGE, SPIN, geom, b.String()}
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
	Target   *[]float64
}

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

// SEnergy returns the list of jobs needed to compute the desired
// semi-empirical energies after writing all of the input files
func SEnergy(dir string, names []string, geoms [][]float64, params []Param, dest *[]float64) []Job {
	os.Mkdir("inp", 0755)
	jobs := make([]Job, len(geoms))
	for i, geom := range geoms {
		jobs[i] = RunGaussian(i, dir, names, geom, params)
		jobs[i].Target = dest
	}
	return jobs
}

func CentralDiff(names []string, geoms [][]float64, params []Param, p, i int) []float64 {
	lgeom := len(geoms)
	params[p].Values[i] += DELTA
	fwd := make([]float64, lgeom)
	fwdJobs := SEnergy("inp", names, geoms, params, &fwd)

	params[p].Values[i] -= 2 * DELTA
	bwd := make([]float64, lgeom)
	bwdJobs := SEnergy("inp", names, geoms, params, &bwd)

	// have to restore the value
	params[p].Values[i] += DELTA

	jobs := append(fwdJobs, bwdJobs...)
	RunJobs(jobs)
	forward := mat.NewDense(lgeom, 1, fwd)
	backward := mat.NewDense(lgeom, 1, bwd)
	// f'(x) ~ [ f(x+h) - f(x-h) ] / 2h ; where h is DELTA here
	var diff mat.Dense
	diff.Sub(forward, backward)
	return Scale(1/(2*DELTA), diff.RawMatrix().Data)
}

// NumJac computes the numerical Jacobian of energies vs params
func NumJac(names []string, geoms [][]float64, params []Param) *mat.Dense {
	rows := len(geoms)
	cols := Len(params)
	jac := mat.NewDense(rows, cols, nil)
	var col int
	// these two loops are over params
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

// RunJobs runs jobs and stores the results in the jobs' Targets
func RunJobs(jobs []Job) {
	// SEnergy makes the directory, this deletes it
	defer os.RemoveAll("inp")
	qstat := make(map[string]bool)
	// initialize to true
	for _, j := range jobs {
		qstat[j.Jobid] = true
	}
	var shortened int
	for len(jobs) > 0 {
		for i := 0; i < len(jobs); i++ {
			job := jobs[i]
			energy, err := ParseGaussian(job.Filename + ".out")
			if err == nil {
				shortened++
				(*job.Target)[job.Index] = energy
				l := len(jobs) - 1
				jobs[l], jobs[i] = jobs[i], jobs[l]
				jobs = jobs[:l]
			} else if !qstat[job.Jobid] {
				jobs[i].Jobid = Submit(job.Filename + ".pbs")
			}
		}
		if shortened < 1 {
			// check the queue if no jobs finish
			Stat(&qstat)
		}
		time.Sleep(1 * time.Second)
		fmt.Fprintf(LOGFILE, "%d jobs remaining\n", len(jobs))
		shortened = 0
	}
}

func OneIter(labels []string, geoms [][]float64, params []Param, outfile string) {
	l := len(geoms)
	nrg := make([]float64, l)
	jobs := SEnergy("inp", labels, geoms, params, &nrg)
	RunJobs(jobs)
	se := Relative(mat.NewDense(l, 1, nrg))
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
		params, _ := LoadParams("params.dat")
		OneIter(labels, geoms, params, *one)
		os.Exit(0)
	}
	LOGFILE, _ = os.Create("log")
	paramLog, _ := os.Create("params.log")
	ai := LoadEnergies("rel.dat")
	params, num := LoadParams("opt.out")
	fmt.Printf("loaded %d params\n", num)
	nrg := make([]float64, len(geoms))
	jobs := SEnergy("inp", labels, geoms, params, &nrg)
	RunJobs(jobs)
	se := Relative(mat.NewDense(len(geoms), 1, nrg))
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
		jac := NumJac(labels, geoms, params)
		*lambda /= NU
		// BEGIN copy-paste
		newParams := LevMar(jac, ai, se, params, 1.0)
		jobs = SEnergy("inp", labels, geoms, newParams, &nrg)
		RunJobs(jobs)
		se = Relative(mat.NewDense(len(geoms), 1, nrg))
		norm = Norm(ai, se) * htToCm
		rmsd = RMSD(ai, se) * htToCm
		// END copy-paste

		// case ii. and iii. of levmar, first case is ii. from
		// Marquardt63
		var bad bool
		for i := 0; norm > lastNorm; i++ {
			*lambda *= NU
			newParams := LevMar(jac, ai, se, params, 1.0)
			jobs = SEnergy("inp", labels, geoms, newParams, &nrg)
			RunJobs(jobs)
			se = Relative(mat.NewDense(len(geoms), 1, nrg))
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
			DumpParams(newParams, "inp/params.dat")
			jobs = SEnergy("inp", labels, geoms, newParams, &nrg)
			RunJobs(jobs)
			se = Relative(mat.NewDense(len(geoms), 1, nrg))
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
