package main

import (
	"bytes"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"runtime/pprof"
	"strings"
	"text/template"
	"time"

	"gonum.org/v1/gonum/mat"
)

const (
	// from http://www.ilpi.com/msds/ref/energyunits.html
	htToCm        = 219_474.5459784
	EPS           = 1e-14
	THRESH        = 1.0
	NU            = 2.0
	CHUNK         = 128
	INFILE_SUFFIX = ".com"
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

type Params []Param

func (p Params) Values() (ret []float64) {
	for i := range p {
		ret = append(ret, p[i].Values...)
	}
	return
}

// WriteParams formats params for use in a Gaussian input file and
// writes them to w
func WriteParams(w io.Writer, params []Param) {
	for _, param := range params {
		fmt.Fprintf(w, " ****\n%2s", param.Atom)
		var last, sep string
		for i, name := range param.Names {
			if strings.Contains(name, "DCore") {
				sep = ","
			} else {
				sep = "="
			}
			if name == last {
				fmt.Fprintf(w, ",%.10f", param.Values[i])
			} else {
				fmt.Fprintf(w, "\n%s%s%.10f",
					name, sep, param.Values[i])
				last = name
			}
		}
		fmt.Fprint(w, "\n")
	}
	fmt.Fprint(w, " ****\n\n")
}

func LogParams(w io.Writer, params []Param, iter int) {
	fmt.Fprintf(w, "Iter %5d\n", iter)
	WriteParams(w, params)
}

func DumpParams(params []Param, filename string) {
	f, err := os.Create(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	WriteParams(f, params)
}

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
	WriteParams(&b, params)
	anon := struct {
		Charge int
		Spin   int
		Geom   string
		Params string
	}{CHARGE, SPIN, geom, b.String()}
	t.Execute(w, anon)
}

// monotonically increasing counter for job names
var counter int

type Type int

const (
	Fwd Type = iota
	Bwd
	None
)

// SEnergy returns the list of jobs needed to compute the desired
// semi-empirical energies after writing the Gaussian input files
func SEnergy(names []string, geoms [][]float64, params []Param, col int, calc Type) []Job {
	jobs := make([]Job, len(geoms))
	for i, geom := range geoms {
		name := fmt.Sprintf("inp/job.%010d", counter)
		input, err := os.Create(name + INFILE_SUFFIX)
		if err != nil {
			panic(err)
		}
		counter++
		WriteCom(input, names, geom, params)
		input.Close()
		jobs[i] = Job{
			Filename: filepath.Base(name),
			I:        i,
			J:        col,
		}
		switch calc {
		case Fwd:
			jobs[i].Coeff = 1 / (2 * DELTA)
		case Bwd:
			jobs[i].Coeff = -1 / (2 * DELTA)
		default:
			jobs[i].Coeff = 1
		}
	}
	return jobs
}

// Job is a type for containing all the information about a quantum
// chemical calculation. Filename should not include any directory
// information.
type Job struct {
	Filename string
	I, J     int
	Jobid    string
	Coeff    float64
}

// CentralDiff returns a list of Jobs needed to compute the ith column
// of the numerical Jacobian
func CentralDiff(names []string, geoms [][]float64, params []Param, p, i, col int) []Job {
	params[p].Values[i] += DELTA
	fwdJobs := SEnergy(names, geoms, params, col, Fwd)

	params[p].Values[i] -= 2 * DELTA
	bwdJobs := SEnergy(names, geoms, params, col, Bwd)

	// have to restore the value
	params[p].Values[i] += DELTA

	return append(fwdJobs, bwdJobs...)
}

// NumJac computes the numerical Jacobian of energies vs params
func NumJac(names []string, geoms [][]float64, params []Param) *mat.Dense {
	rows := len(geoms)
	cols := Len(params)
	jac := mat.NewDense(rows, cols, nil)
	var col int
	// Set up all the Jobs for the Jacobian
	jobs := make([]Job, 0)
	for p := range params {
		for i := range params[p].Values {
			if *debug {
				fmt.Printf("before:%12.9f\n", Params(params).Values())
			}
			jobs = append(jobs,
				CentralDiff(names, geoms, params, p, i, col)...)
			if *debug {
				fmt.Printf(" after:%12.9f\n", Params(params).Values())
			}
			col++
		}
	}
	RunJobs(jobs, jac)
	// fmt.Fprintf(LOGFILE, "finished col %5d -> %s of %s\n", col,
	// 	params[p].Names[i], params[p].Atom)
	cleanup()
	return jac
}

// Remove the `inp` directory and recreate it
func cleanup() {
	if err := os.RemoveAll("inp"); err != nil {
		panic(err)
	}
	if err := os.Mkdir("inp", 0744); err != nil {
		panic(err)
	}
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

func Filenames(jobs []Job) []string {
	ret := make([]string, len(jobs))
	for j := range jobs {
		ret[j] = jobs[j].Filename
	}
	return ret
}

// RunJobs runs jobs and stores the results in target
func RunJobs(jobs []Job, target *mat.Dense) {
	var check *mat.Dense
	if *debug {
		r, c := target.Dims()
		check = mat.NewDense(r, c, nil)
	}
	var pbs int
	chunkees := make([]Job, 0, CHUNK)
	runJobs := make([]Job, 0, len(jobs))
	// submit the jobs in groups of size CHUNK, then store the
	// updated jobs in runJobs
	for j, job := range jobs {
		chunkees = append(chunkees, job)
		if (j > 0 && j%CHUNK == 0) || j == len(jobs)-1 {
			name := fmt.Sprintf("inp/%d.pbs", pbs)
			pbs++
			f, err := os.Create(name)
			if err != nil {
				panic(err)
			}
			WritePBS(f, name, Filenames(chunkees))
			f.Close()
			jobid := Submit(name)
			for c := range chunkees {
				chunkees[c].Jobid = jobid
			}
			runJobs = append(runJobs, chunkees...)
			chunkees = make([]Job, 0, CHUNK)
		}
	}
	qstat := make(map[string]bool)
	// initialize to true
	for _, j := range runJobs {
		qstat[j.Jobid] = true
	}
	var shortened int
	for len(runJobs) > 0 {
		for i := 0; i < len(runJobs); i++ {
			job := runJobs[i]
			energy, err := ParseGaussian(
				filepath.Join("inp", job.Filename+".out"),
			)
			if err == nil {
				shortened++
				if *debug {
					check.Set(job.I, job.J, check.At(job.I, job.J)+1)
				}
				target.Set(job.I, job.J,
					target.At(job.I, job.J)+
						job.Coeff*energy)
				l := len(runJobs) - 1
				runJobs[l], runJobs[i] = runJobs[i], runJobs[l]
				runJobs = runJobs[:l]
			} else if !qstat[job.Jobid] {
				// delete the old job, resubmit, and
				// default qstat to true
				delete(qstat, job.Jobid)
				runJobs[i] = Resubmit(job)
				qstat[runJobs[i].Jobid] = true
			}
		}
		if shortened < 1 {
			// check the queue if no runJobs finish
			Stat(&qstat)
		}
		time.Sleep(1 * time.Second)
		fmt.Fprintf(LOGFILE, "%d jobs remaining\n", len(runJobs))
		shortened = 0
	}
	if *debug {
		DumpMat(check)
	}
}

func Resubmit(job Job) Job {
	src, err := os.Open(filepath.Join("inp", job.Filename+INFILE_SUFFIX))
	defer src.Close()
	if err != nil {
		panic(err)
	}
	inp := filepath.Join("inp", job.Filename+"_redo"+INFILE_SUFFIX)
	dst, err := os.Create(inp)
	if err != nil {
		panic(err)
	}
	defer dst.Close()
	io.Copy(dst, src)
	pbs := filepath.Join("inp", job.Filename+"_redo.pbs")
	f, err := os.Create(pbs)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	WritePBS(f, job.Filename, []string{job.Filename + "_redo"})
	fmt.Printf("resubmitting %s as %s\n", job.Filename, job.Filename+"_redo")
	return Job{
		Filename: job.Filename + "_redo",
		I:        job.I,
		J:        job.J,
		Jobid:    Submit(pbs),
		Coeff:    job.Coeff,
	}
}

func OneIter(names []string, geoms [][]float64, params []Param, outfile string) {
	l := len(geoms)
	nrg := mat.NewDense(l, 1, nil)
	jobs := SEnergy(names, geoms, params, 0, None)
	RunJobs(jobs, nrg)
	se := Relative(nrg)
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
	// TODO take these from input file
	labels := []string{"C", "C", "C", "H", "H"}
	geoms := LoadGeoms("file07")
	if *one != "" {
		params, _ := LoadParams("params.dat")
		OneIter(labels, geoms, params, *one)
		os.Exit(0)
	}
	os.RemoveAll("inp")
	os.Mkdir("inp", 0755)
	LOGFILE, _ = os.Create("log")
	paramLog, _ := os.Create("params.log")
	ai := LoadEnergies("rel.dat")
	params, num := LoadParams("opt.out")
	fmt.Printf("loaded %d params\n", num)
	nrg := mat.NewDense(len(geoms), 1, nil)
	jobs := SEnergy(labels, geoms, params, 0, None)
	RunJobs(jobs, nrg)
	se := Relative(nrg)
	norm, max := Norm(ai, se)
	rmsd := RMSD(ai, se) * htToCm
	var (
		iter     int
		lastNorm float64
		lastRMSD float64
	)
	fmt.Printf("%17s%12s%12s%12s%12s%12s\n",
		"cm-1", "cm-1", "cm-1", "cm-1", "cm-1", "s")
	fmt.Printf("%5s%12s%12s%12s%12s%12s%12s\n",
		"Iter", "Norm", "ΔNorm", "RMSD", "ΔRMSD", "Max", "Time")
	fmt.Printf("%5d%12.4f%12.4f%12.4f%12.4f%12.4f%12.1f\n",
		iter, norm, norm-lastNorm, rmsd, rmsd-lastRMSD, max, 0.0)
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
		jobs = SEnergy(labels, geoms, newParams, 0, None)
		nrg.Zero()
		RunJobs(jobs, nrg)
		se = Relative(nrg)
		norm, max = Norm(ai, se)
		rmsd = RMSD(ai, se) * htToCm
		// END copy-paste

		// case ii. and iii. of levmar, first case is ii. from
		// Marquardt63
		var bad bool
		for i := 0; norm > lastNorm; i++ {
			*lambda *= NU
			newParams := LevMar(jac, ai, se, params, 1.0)
			jobs = SEnergy(labels, geoms, newParams, 0, None)
			nrg.Zero()
			RunJobs(jobs, nrg)
			se = Relative(nrg)
			norm, max = Norm(ai, se)
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
			jobs = SEnergy(labels, geoms, newParams, 0, None)
			nrg.Zero()
			RunJobs(jobs, nrg)
			se = Relative(nrg)
			norm, max = Norm(ai, se)
			rmsd = RMSD(ai, se) * htToCm
			fmt.Fprintf(LOGFILE, "\tk_%d to %g with Δ = %f\n",
				i, k, norm-lastNorm)
			prev = norm
		}
		fmt.Printf("%5d%12.4f%12.4f%12.4f%12.4f%12.4f%12.1f\n",
			iter, norm, norm-lastNorm, rmsd, rmsd-lastRMSD, max,
			float64(time.Since(start))/1e9)
		start = time.Now()
		lastNorm = norm
		lastRMSD = rmsd
		params = newParams
		LogParams(paramLog, params, iter)
		iter++
	}
}
