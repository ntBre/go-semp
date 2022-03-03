package main

import "fmt"

type Param struct {
	Name  string
	Atom  string
	Value float64
}

func (p Param) String() string {
	return fmt.Sprintf(
		"%-8s%4s%20.12f\n",
		p.Name, p.Atom, p.Value,
	)
}

func Values(params []Param) (ret []float64) {
	for _, p := range params {
		ret = append(ret, p.Value)
	}
	return
}
