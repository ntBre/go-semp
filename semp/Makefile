DEST = 'woods:semp/testing/.'

semp : main.go utils.go
	go build .

files :
	scp -C file07 rel.dat opt.out $(DEST)

deploy: semp
	scp -C semp $(DEST)

eland: semp
	scp -C semp 'eland:semp/.'

clean:
	rm -f params.dat out
