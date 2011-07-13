.SUFFIXES: .pdf .tex

%.pdf: %.tex
	@echo "        PDF $<"
	pdflatex -shell-escape $<
	bibtex $(shell echo $< | sed -e 's/.tex//')
	pdflatex -shell-escape $<
	pdflatex -shell-escape $<

default: design.pdf
