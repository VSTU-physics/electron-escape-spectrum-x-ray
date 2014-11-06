DOC=документация.tex main.tex
THEORY=теория.tex

all: docs theory

docs: $(DOC)
	latexmk -pdf документация.tex

theory: $(THEORY)
	latexmk -pdf теория.tex
