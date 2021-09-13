semp : main.go utils.go
	go build .

deploy: semp
	scp -C semp 'woods:semp/testing/.'

clean:
	rm -f params.dat out
