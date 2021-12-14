DEST = 'woods:semp/testing/.'

semp : main.go utils.go pbs.go load.go slurm.go stats.go
	go build .

files :
	scp -C file07 rel.dat opt.out $(DEST)

deploy: semp
	scp -C semp $(DEST)

eland: semp
	scp -C semp 'eland:programs/semp/.'

clean:
	rm -f params.dat out

test:
	go test .

cover:
	go test . -v -coverprofile=/tmp/pbqff.out; go tool cover -html /tmp/pbqff.out
