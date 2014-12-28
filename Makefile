all: code docs theory

code:
	make -C src analytics
	make -C src mc
	cp src/analytics ./
	cp src/monte-carlo ./

docs:
	make -C doc docs
	cp doc/документация.pdf ./

theory:
	make -C doc theory
	cp doc/теория.pdf ./

clean:
	rm -f документация.pdf теория.pdf spectrum *.dat *.gp
