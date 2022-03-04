package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
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
	INFILE_SUFFIX = ".mop"
	GAMMA0        = math.Pi / 4
	MAX_TRIES     = 5
)

var (
	CHARGE = 0
	SPIN   = 1
	//  https://en.wikipedia.org/wiki/Numerical_differentiation
	//  recommends cube root of machine eps (~2.2e16) for step
	//  size => 6e-6; adjusting down from there
	DELTA = 1e-8
)

// Flags
var (
	debug   = flag.String("debug", "", "enable debug info (params|jac)")
	cpuprof = flag.String("cpu", "", "write a CPU profile")
	one     = flag.Bool("one", false, "compute initial SE energies and exit")
	version = flag.Bool("v", false, "print version number and exit")
)

// WriteParams formats params for use in a MOPAC input file and
// writes them to w
func WriteParams(w io.Writer, params []Param) error {
	nw := bufio.NewWriter(w)
	for _, p := range params {
		fmt.Fprintf(nw, "%-8s%8s%20.12f\n",
			p.Name, p.Atom, p.Value,
		)
	}
	fmt.Fprint(nw, "\n")
	return nw.Flush()
}

func LogParams(w io.Writer, params []Param, iter int) {
	fmt.Fprintf(w, "Iter %5d\n", iter)
	WriteParams(w, params)
}

func WriteCom(w io.Writer, names []string, coords []float64, params []Param) {
	f, err := os.CreateTemp("tmparam", "")
	if err != nil {
		panic(err)
	}
	geom := ZipGeom(names, coords)
	/*
	   XYZ - Cartesian geometry
	   A0 - use atomic units (bohr)
	   precise - use higher precision
	   relscf - tighten precision by an additional multiplicative factor
	*/
	t, err := template.New("com").Parse(
		`XYZ A0 scfcrt=1.D-21 aux(precision=14) external={{.Name}} 1SCF charge={{.Charge}} PM6
blank line
blank line
{{.Geom}}

`)
	if err != nil {
		panic(err)
	}
	WriteParams(f, params)
	err = f.Close()
	if err != nil {
		panic(err)
	}
	abs, err := filepath.Abs(f.Name())
	if err != nil {
		panic(err)
	}
	t.Execute(w, struct {
		Geom   string
		Name   string
		Charge int
	}{
		Charge: CHARGE,
		Geom:   geom,
		Name:   abs,
	})
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
func SEnergy(names []string, geoms [][]float64, params []Param, col int,
	calc Type) []Job {
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
	Jobid    string
	I        int
	J        int
	Coeff    float64
}

// CentralDiff returns a list of Jobs needed to compute the ith column
// of the numerical Jacobian
func CentralDiff(names []string, geoms [][]float64, params []Param,
	i, col int) []Job {
	params[i].Value += DELTA
	fwdJobs := SEnergy(names, geoms, params, col, Fwd)

	params[i].Value -= 2 * DELTA
	bwdJobs := SEnergy(names, geoms, params, col, Bwd)

	// have to restore the value
	params[i].Value += DELTA

	return append(fwdJobs, bwdJobs...)
}

// NumJac computes the numerical Jacobian of energies vs params
func NumJac(names []string, geoms [][]float64, params []Param) *mat.Dense {
	rows := len(geoms)
	cols := len(params)
	jac := mat.NewDense(rows, cols, nil)
	var col int
	// Set up all the Jobs for the Jacobian
	jobs := make([]Job, 0)
	for i := range params {
		if *debug == "params" {
			fmt.Printf("before:%12.9f\n",
				Values(params))
		}
		jobs = append(jobs,
			CentralDiff(names, geoms, params, i, col)...)
		if *debug == "params" {
			fmt.Printf(" after:%12.9f\n",
				Values(params))
		}
		col++
	}
	RunJobs(jobs, jac)
	// fmt.Fprintf(LOGFILE, "finished col %5d -> %s of %s\n", col,
	// 	params[p].Names[i], params[p].Atom)
	cleanup()
	if *debug == "jac" {
		DumpMat(jac)
	}
	return jac
}

// Remove the `inp` directory and recreate it
func cleanup() {
	if err := os.RemoveAll("inp"); err != nil {
		panic(err)
	}
	if err := os.RemoveAll("tmparam"); err != nil {
		panic(err)
	}
	if err := os.Mkdir("inp", 0744); err != nil {
		panic(err)
	}
	if err := os.Mkdir("tmparam", 0744); err != nil {
		panic(err)
	}
}

func UpdateParams(params []Param, v *mat.Dense, scale float64) []Param {
	ret := make([]Param, 0, len(params))
	for i, p := range params {
		ret = append(ret, Param{
			Atom:  p.Atom,
			Name:  p.Name,
			Value: p.Value + v.At(i, 0)*scale,
		})
	}
	return ret
}

func LevMar(jac, ai, se *mat.Dense, params []Param, scale, lambda float64) (
	newParams []Param, gamma float64) {
	// LHS
	var a mat.Dense
	a.Mul(jac.T(), jac)
	// construct A* as a_rc = a_rc / [ sqrt(a_rr) * sqrt(a_cc) ], where a is
	// jacTjac
	r, c := a.Dims()
	Astar := mat.NewDense(r, c, nil)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			Astar.Set(i, j,
				a.At(i, j)/
					(math.Sqrt(a.At(i, i))*
						math.Sqrt(a.At(j, j))))
		}
	}
	eye := Identity(r)
	var Leye mat.Dense
	Leye.Scale(lambda, eye)
	var lhs mat.Dense
	lhs.Add(Astar, &Leye)
	// RHS
	var diff mat.Dense
	diff.Sub(ai, se)
	var g mat.Dense
	g.Mul(jac.T(), &diff)
	r, _ = g.Dims()
	// rhs = g*
	rhs := mat.NewDense(r, 1, nil)
	for i := 0; i < r; i++ {
		rhs.Set(i, 0, g.At(i, 0)/math.Sqrt(a.At(i, i)))
	}
	// d = δ
	var d mat.Dense
	err := d.Solve(&lhs, rhs)
	if err != nil {
		fmt.Fprintln(os.Stderr, "jac")
		DumpMat(jac)
		fmt.Fprintln(os.Stderr, "prod")
		DumpMat(&a)
		fmt.Fprintln(os.Stderr, "diff")
		DumpVec(&diff)
		fmt.Fprintln(os.Stderr, "prod2")
		DumpMat(&g)
		fmt.Fprintln(os.Stderr, "step")
		DumpVec(&d)
		if strings.Contains(err.Error(), "Inf") {
			panic(err)
		}
		fmt.Fprintf(os.Stderr, "Warning: %v\n", err)
	}
	// step is d* here, convert back to d with eqn 30
	r, _ = d.Dims()
	for i := 0; i < r; i++ {
		d.Set(i, 0, d.At(i, 0)/math.Sqrt(a.At(i, i)))
	}
	// gamma = acos( [delta(transpose) * g] / [mag(delta) * mag(g)] )
	gam := mat.NewDense(1, 1, nil)
	gam.Mul(d.T(), &g)
	magD := mat.Norm(&d, 2)
	magG := mat.Norm(&g, 2)
	gamma = math.Acos(gam.At(0, 0) / (magD * magG))
	newParams = UpdateParams(params, &d, scale)
	return
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
			energy, err := ReadOut(
				filepath.Join("inp", job.Filename+".out"),
			)
			if err == nil {
				shortened++
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
				log.Printf("resubmitting %s as %s for %v\n",
					job.Filename, job.Filename+"_redo",
					err,
				)
				runJobs[i] = Resubmit(job)
				qstat[runJobs[i].Jobid] = true
			}
		}
		if shortened < 1 {
			// check the queue if no runJobs finish
			Stat(&qstat)
			time.Sleep(1 * time.Second)
			fmt.Fprintf(os.Stderr,
				"%d jobs remaining\n", len(runJobs))
		}
		shortened = 0
	}
	fmt.Fprintln(os.Stderr, "jobs done")
}

func setup() {
	os.RemoveAll("tmparam")
	os.MkdirAll("tmparam", 0744)
	os.RemoveAll("inp")
	os.Mkdir("inp", 0744)
}

func takedown() {
	os.RemoveAll("tmparam")
}

func inner(atoms []string, geoms [][]float64, jac, ai, se *mat.Dense,
	params []Param, scale, lambda float64) (newParams []Param,
	newSe *mat.Dense, norm, max, gamma float64) {
	newParams, gamma = LevMar(jac, ai, se, params, scale, lambda)
	jobs := SEnergy(atoms, geoms, newParams, 0, None)
	nrg := mat.NewDense(len(geoms), 1, nil)
	RunJobs(jobs, nrg)
	newSe = Relative(nrg)
	norm, max = Norm(ai, newSe)
	return
}

func printHeader() {
	fmt.Printf("%17s%12s%12s%12s%12s%12s\n",
		"cm-1", "cm-1", "cm-1", "cm-1", "cm-1", "s")
	fmt.Printf("%5s%12s%12s%12s%12s%12s%12s\n",
		"Iter", "Norm", "ΔNorm", "RMSD", "ΔRMSD", "Max", "Time")
}
func printStep(iter int, norm, lastNorm, rmsd, lastRMSD, max, time float64) {
	fmt.Printf("%5d%12.4f%12.4f%12.4f%12.4f%12.4f%12.1f\n",
		iter, norm, norm-lastNorm, rmsd, rmsd-lastRMSD, max, time)
}

func work(infile string) (norms []float64) {
	conf := LoadConfig(infile)
	geoms := LoadGeoms(conf.GeomFile)
	paramLog, _ := os.Create("params.log")
	ai := LoadEnergies(conf.EnergyFile)
	params := conf.Params
	lambda := conf.Lambda
	fmt.Printf("loaded %d params\n", len(params))
	nrg := mat.NewDense(len(geoms), 1, nil)
	jobs := SEnergy(conf.Atoms, geoms, params, 0, None)
	RunJobs(jobs, nrg)
	if *one {
		DumpVec(nrg)
		os.Exit(0)
	}
	se := Relative(nrg)
	norm, max := Norm(ai, se)
	rmsd := RMSD(ai, se) * htToCm
	var (
		lastNorm float64
		lastRMSD float64
	)
	printHeader()
	printStep(0, norm, lastNorm, rmsd, lastRMSD, max, 0.0)
	LogParams(paramLog, params, 0)
	lastNorm = norm
	lastRMSD = rmsd
	start := time.Now()
	norms = []float64{norm}
	var (
		newParams []Param
		newSe     *mat.Dense
		gamma     float64
		delNorm   float64 = 1
	)
	for iter := 1; iter <= conf.MaxIt && norm > THRESH &&
		math.Abs(delNorm) > 1e-4; iter++ {
		jac := NumJac(conf.Atoms, geoms, params)
		lambda /= NU
		newParams, newSe, norm, max, gamma = inner(
			conf.Atoms, geoms, jac, ai, se, params, 1.0, lambda,
		)

		// case ii. and iii. of levmar, first case is ii. from
		// Marquardt63
		var (
			bad       bool
			lastGamma float64
		)
		for i := 0; norm > lastNorm; i++ {
			lambda *= NU
			dNorm := norm - lastNorm
			newParams, newSe, norm, max, gamma = inner(
				conf.Atoms, geoms, jac, ai,
				se, params, 1.0, lambda,
			)
			fmt.Fprintf(os.Stderr,
				"\tλ_%d to %g with ΔNorm = %f\n",
				i, lambda, norm-lastNorm,
			)
			switch {
			case i > 0 && gamma-lastGamma > 1e-6:
				log.Println("breaking gamma")
				bad = true
			case norm-lastNorm > dNorm:
				log.Println("breaking norm")
				bad = true
			case i > MAX_TRIES:
				log.Println("breaking i")
				bad = true
			}
			if bad {
				break
			}
		}
		var k float64 = 1
		for i := 2; bad && norm > lastNorm && k > 1e-14; i++ {
			k = 1.0 / math.Pow(10, float64(i))
			newParams, newSe, norm, max, gamma = inner(
				conf.Atoms, geoms, jac, ai,
				se, params, k, lambda,
			)
			fmt.Fprintf(os.Stderr,
				"\tk_%d to %g with ΔNorm = %f\n",
				i, k, norm-lastNorm)
		}
		rmsd = RMSD(ai, newSe) * htToCm
		printStep(iter, norm, lastNorm, rmsd, lastRMSD, max,
			float64(time.Since(start))/1e9)
		norms = append(norms, norm)
		start = time.Now()
		delNorm = norm - lastNorm
		lastNorm = norm
		lastRMSD = rmsd
		params = newParams
		se = newSe
		LogParams(paramLog, params, iter)
	}
	return
}

func main() {
	host, _ := os.Hostname()
	flag.Parse()
	if *version {
		fmt.Printf("semp version %s\n", VERSION)
		os.Exit(0)
	}
	args := flag.Args()
	infile := "semp"
	if len(args) >= 1 {
		infile = args[0]
	}
	DupOutErr(infile)
	fmt.Printf("semp version %s\n", VERSION)
	fmt.Printf("running on host: %s\n", host)
	setup()
	defer takedown()
	if *cpuprof != "" {
		f, err := os.Create(*cpuprof)
		if err != nil {
			panic(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	iters := work("semp.in")
	fmt.Printf("Convergence reached after %d iterations\n", len(iters))
}
