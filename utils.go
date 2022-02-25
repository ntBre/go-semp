package main

import (
	"fmt"
	"io"
	"math"
	"os"
	"path"
	"strconv"
	"strings"
	"syscall"

	"gonum.org/v1/gonum/mat"
)

const (
	// from https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
	toAng = 0.5291_772_109_03
)

func Len(params []Param) (sum int) {
	for _, p := range params {
		sum += len(p.Names)
	}
	return sum
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

// ZipGeom combines a list of names with a list of coordinates to
// yield a string geometry
func ZipGeom(names []string, coords []float64) string {
	var geom strings.Builder
	for i := range names {
		fmt.Fprintf(&geom, "%s%20.12f%20.12f%20.12f\n",
			names[i],
			coords[3*i],
			coords[3*i+1],
			coords[3*i+2],
		)
	}
	return geom.String()
}

// toFloat converts a list of strings to a float64 using
// strconv.ParseFloat
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

func Equal(a, b float64) bool {
	if math.Abs(a-b) > EPS {
		return false
	}
	return true
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
	WriteMat(os.Stdout, m)
}

func DumpJac(m mat.Matrix) {
	f, err := os.Create("jac")
	defer f.Close()
	if err != nil {
		panic(err)
	}
	WriteMat(f, m)
}

func WriteMat(w io.Writer, m mat.Matrix) {
	r, c := m.Dims()
	for i := 0; i < r; i++ {
		fmt.Fprintf(w, "%5d", i)
		for j := 0; j < c; j++ {
			fmt.Fprintf(w, "%12.8f", m.At(i, j))
		}
		fmt.Fprint(w, "\n")
	}
	fmt.Fprint(w, "\n")
}

func Identity(n int) *mat.Dense {
	ret := mat.NewDense(n, n, nil)
	for i := 0; i < n; i++ {
		ret.Set(i, i, 1.0)
	}
	return ret
}

// DupOutErr uses syscall.Dup2 to direct the stdout and stderr streams
// to files
func DupOutErr(infile string) {
	// set up output and err files and dup their fds to stdout and stderr
	// https://github.com/golang/go/issues/325
	base := infile[:len(infile)-len(path.Ext(infile))]
	outfile, _ := os.Create(base + ".out")
	errfile, _ := os.Create(base + ".log")
	syscall.Dup2(int(outfile.Fd()), 1)
	syscall.Dup2(int(errfile.Fd()), 2)
}
