code:
	make -C src analytics
	make -C src mc
	make -C src test
	cp src/analytics ./
	cp src/monte-carlo ./
	cp src/test ./

clean:
	rm -f *.dat *.gp monte-carlo analytics test
