package main

import (
	"bufio"
	"os"
	"os/exec"
	"strings"
)

// Stat generates an updated map of job names to their queue
// status. The map value is true if the job is either pending (PD),
// queued (Q) or running (R) and false otherwise
func Stat(qstat *map[string]bool) {
	status, _ := exec.Command("squeue", "-u", os.Getenv("USER")).CombinedOutput()
	scanner := bufio.NewScanner(strings.NewReader(string(status)))
	var (
		line   string
		fields []string
		header = true
	)
	// initialize them all to false and set true if run
	for key := range *qstat {
		(*qstat)[key] = false
	}
	for scanner.Scan() {
		line = scanner.Text()
		if strings.Contains(line, "JOBID") {
			header = false
			continue
		} else if header {
			continue
		}
		fields = strings.Fields(line)
		if _, ok := (*qstat)[fields[0]]; ok {
			// jobs are initially put in PD = pending
			// state
			if strings.Contains("PDQR", fields[4]) {
				(*qstat)[fields[0]] = true
			}
		}
	}
}
