code:
	make -C src analytics
	make -C src mc
	cp src/analytics ./
	cp src/monte-carlo ./

docs:
	make -C doc docs

clean:
	rm -f *.dat *.gp monte-carlo analytics test
