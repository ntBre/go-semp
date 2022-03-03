DEST = 'woods:semp/testing/.'
SHORT = 0
TESTFLAGS = -v -failfast

ifeq ($(SHORT),1)
TESTFLAGS += -short
endif

semp : *.go *.tmpl
	go build .

files :
	scp -C file07 rel.dat opt.out $(DEST)

deploy: semp
	scp -C semp $(DEST)

eland: semp scripts/convert.py
	scp -C $? 'eland:programs/semp/.'
	date >> eland

clean:
	rm -f params.dat out

test:
	go test . ${TESTFLAGS}

prof:
	go test . -run=TestMain -cpuprofile=/tmp/semp.prof.out ${TESTFLAGS}

cover:
	go test . -v -coverprofile=/tmp/pbqff.out
	go tool cover -html /tmp/pbqff.out
