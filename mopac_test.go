package main

import (
	"reflect"
	"testing"
)

func TestReadOut(t *testing.T) {
	got, _ := ReadOut("testfiles/job.out")
	want := 0.97127947459164715838e+02 / KCALHT
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
