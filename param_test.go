package main

import (
	"reflect"
	"testing"
)

func TestString(t *testing.T) {
	got := Param{
		Name:  "USS",
		Atom:  "H",
		Value: 1.0,
	}.String()
	want := "USS        H      1.000000000000\n"
	if got != want {
		t.Errorf("got\n%q, wanted\n%q\n", got, want)
	}
}

func TestValues(t *testing.T) {
	got := Values([]Param{
		{
			Name:  "USS",
			Atom:  "H",
			Value: 1.0,
		},
	})
	want := []float64{1.0}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
