package main

import (
	"bytes"
	"fmt"
	"os"
	"regexp"
	"testing"
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

func TestWriteParams(t *testing.T) {
	var b bytes.Buffer
	p, _ := LoadParams("testfiles/opt.out")
	WriteParams(&b, p)
	got := b.String()
	want := `USS            H    -11.246958000000
ZS             H      1.268641000000
BETAS          H     -8.352984000000
GSS            H     14.448686000000
FN11           H      0.024184000000
FN21           H      3.055953000000
FN31           H      1.786011000000
ALPB_1         H      3.540942000000
XFAC_1         H      2.243587000000
ALPB_6         H      1.027806000000
XFAC_6         H      0.216506000000
USS            C    -51.089653000000
UPP            C    -39.937920000000
ZS             C      2.047558000000
ZP             C      1.702841000000
BETAS          C    -15.385236000000
BETAP          C     -7.471929000000
GSS            C     13.335519000000
GPP            C     10.778326000000
GSP            C     11.528134000000
GP2            C      9.486212000000
HSP            C      0.717322000000
FN11           C      0.046302000000
FN21           C      2.100206000000
FN31           C      1.333959000000
ALPB_1         C      1.027806000000
XFAC_1         C      0.216506000000
ALPB_6         C      2.613713000000
XFAC_6         C      0.813510000000

`
	if got != want {
		charComp(got, want)
		t.Errorf("got\n%#v, wanted\n%#v\n", got, want)
	}
}

func TestWriteCom(t *testing.T) {
	var b bytes.Buffer
	p, _ := LoadParams("testfiles/opt.out")
	os.MkdirAll("tmparam", 0744)
	defer os.RemoveAll("tmparam")
	WriteCom(&b, []string{"C", "C", "C", "H", "H"},
		[]float64{
			0.0000000000, 0.0000000000, -1.6794733900,
			0.0000000000, 1.2524327590, 0.6959098120,
			0.0000000000, -1.2524327590, 0.6959098120,
			0.0000000000, 3.0146272390, 1.7138963510,
			0.0000000000, -3.0146272390, 1.7138963510,
		}, p)
	got := b.String()
	want := `XYZ A0 scfcrt=1.D-21 aux(precision=14) external=tmparam/1598109253 1SCF charge=0 PM6
blank line
blank line
C      0.000000000000      0.000000000000     -1.679473390000
C      0.000000000000      1.252432759000      0.695909812000
C      0.000000000000     -1.252432759000      0.695909812000
H      0.000000000000      3.014627239000      1.713896351000
H      0.000000000000     -3.014627239000      1.713896351000


`
	random := regexp.MustCompile(`(external=)/?.*tmparam/[^ ]+`)
	got = random.ReplaceAllString(got, "$1")
	want = random.ReplaceAllString(want, "$1")
	if got != want {
		charComp(got, want)
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}
