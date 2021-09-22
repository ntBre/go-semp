DEST = 'woods:semp/testing/.'
semp : main.go utils.go
	go build .

files :
	scp -C file07 rel.dat $(DEST)

deploy: semp
	scp -C semp $(DEST)

clean:
	rm -f params.dat out
