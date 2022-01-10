package main

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// Relative makes the values in a relative to its minimum
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

// Norm computes the Euclidean norm between vectors a and b
func Norm(a, b *mat.Dense) float64 {
	var diff mat.Dense
	diff.Sub(a, b)
	return mat.Norm(&diff, 2)
}

// RMSD computes the root-mean-square deviation between vectors a and
// b
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
