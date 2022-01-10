package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

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
