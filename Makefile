all: code docs theory

code:
	make -C src release
	cp src/main ./

docs:
	make -C doc docs
	cp doc/документация.pdf ./

theory:
	make -C doc theory
	cp doc/теория.pdf ./

clean:
	rm -f документация.pdf теория.pdf main
