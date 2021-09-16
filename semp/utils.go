package main

import (
	"fmt"
	"strconv"
	"strings"
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
			names[i], coords[3*i]*toAng,
			coords[3*i+1]*toAng, coords[3*i+2]*toAng,
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