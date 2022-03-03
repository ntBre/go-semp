package main

import (
	"fmt"
	"math"
	"os"

	"gonum.org/v1/gonum/mat"
)

func charComp(got, want string) {
	for c := range got {
		if len(got) <= c {
			fmt.Println("got too short")
			return
		} else if len(want) <= c {
			fmt.Println("want too short")
			return
		}
		if got[c] != want[c] {
			fmt.Printf("got\n%q, wanted\n%q\n",
				got[:c+1], want[:c+1])
			return
		}
	}
}

func test_takedown() {
	takedown()
	os.RemoveAll("inp")
}

// return the norm of the difference between got and want and whether or not it
// is greater than eps
func vecNorm(got, want *mat.Dense, eps float64) (float64, bool) {
	var diff mat.Dense
	diff.Sub(got, want)
	norm := mat.Norm(&diff, 2)
	return norm, norm > eps
}

// print the difference between column vectors got and want
func vecDiff(got, want *mat.Dense) {
	l, _ := got.Dims()
	fmt.Printf("\n%20s%20s%20s\n", "Got", "Want", "Diff")
	for i := 0; i < l; i++ {
		g := got.At(i, 0)
		w := want.At(i, 0)
		fmt.Printf("%20.12f%20.12f%20.12f\n",
			g, w, g-w,
		)
	}
}

func compMat(a, b *mat.Dense, eps float64) bool {
	var diff mat.Dense
	diff.Sub(a, b)
	return mat.Norm(&diff, 2) < eps
}

func compFloat(a, b []float64, eps float64) bool {
	for i := range a {
		if math.Abs(a[i]-b[i]) > eps {
			return false
		}
	}
	return true
}
