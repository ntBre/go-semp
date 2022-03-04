package main

import (
	"bufio"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/BurntSushi/toml"
	"gonum.org/v1/gonum/mat"
)

type RawConf struct {
	// Intial semi-empirical parameters to be optimized, in MOPAC input
	// format. See MOPAC documentation for the external keyword for details
	Params string
	// Space-delimited list of atom labels like "C C C H H"
	Atoms string
	// INTDER-style file containg a list of geometries, defaults to file07
	GeomFile string
	// File containing a list of ab initio/training energies, defaults to
	// rel.dat
	EnergyFile string
	// The maximum number of iterations the optimization should run,
	// defaults to 250
	MaxIt int
	// The initial value for the lambda parameter in the Levenberg-Marquardt
	// algorithm, defaults to 1e-8
	Lambda float64
	// Whether or not the rank-one Broyden update should be used to
	// approximate the Jacobian
	Broyden bool
	// The interval at which to compute a full numerical Jacobian if Broyden
	// is used, defaults to 5
	BroydInt int
}

func (rc RawConf) ToConfig() (conf Config) {
	conf.Params = LoadParams(rc.Params)
	conf.Atoms = strings.Fields(rc.Atoms)
	conf.GeomFile = rc.GeomFile
	conf.EnergyFile = rc.EnergyFile
	conf.MaxIt = rc.MaxIt
	conf.Lambda = rc.Lambda
	conf.Broyden = rc.Broyden
	conf.BroydInt = rc.BroydInt
	return
}

type Config struct {
	GeomFile   string
	EnergyFile string
	Params     []Param
	Atoms      []string
	MaxIt      int
	Lambda     float64
	Broyden    bool
	BroydInt   int
}

func LoadConfig(filename string) Config {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	cont, err := io.ReadAll(f)
	if err != nil {
		panic(err)
	}
	// Defaults
	rc := RawConf{
		GeomFile:   "file07",
		EnergyFile: "rel.dat",
		MaxIt:      250,
		Lambda:     1e-8,
		Broyden:    false,
		BroydInt:   5,
	}
	err = toml.Unmarshal(cont, &rc)
	if err != nil {
		panic(err)
	}
	return rc.ToConfig()
}

func LoadGeoms(filename string) (ret [][]float64) {
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	var line string
	row := make([]float64, 0, 3)
	scanner.Scan() // discard first "# GEOM" line
	for scanner.Scan() {
		line = scanner.Text()
		if strings.Contains(line, "#") {
			ret = append(ret, row)
			row = make([]float64, 0, len(row))
		} else {
			fields := strings.Fields(line)
			row = append(row, toFloat(fields)...)
		}
	}
	ret = append(ret, row)
	return
}

// LoadEnergies loads relative energies from filename and returns them
// as a column vector
func LoadEnergies(filename string) *mat.Dense {
	ret := make([]float64, 0)
	f, err := os.Open(filename)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	scanner := bufio.NewScanner(f)
	var line string
	for scanner.Scan() {
		line = scanner.Text()
		fields := strings.Fields(line)
		if len(fields) == 1 {
			v, err := strconv.ParseFloat(fields[0], 64)
			if err != nil {
				panic(err)
			}
			ret = append(ret, v)
		}
	}
	return mat.NewDense(len(ret), 1, ret)
}

// LoadParams extracts semi-empirical parameters from an input parameter string.
// The expected format is that of a MOPAC input file, not output file
func LoadParams(params string) (ret []Param) {
	scanner := bufio.NewScanner(
		strings.NewReader(params),
	)
	var (
		line   string
		fields []string
	)
	// TODO might need to deduplicate symmetrical terms like 1
	// XFAC_6 vs 6 XFAC_1
	for scanner.Scan() {
		line = scanner.Text()
		fields = strings.Fields(line)
		switch {
		case len(fields) != 3:
			continue
		default:
			v, err := strconv.ParseFloat(fields[2], 64)
			if err != nil {
				log.Fatalf("failed to parse %q as float",
					fields[2],
				)
			}
			ret = append(ret,
				Param{
					Name:  fields[0],
					Atom:  fields[1],
					Value: v,
				})
		}
	}
	return
}
