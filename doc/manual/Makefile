.SUFFIXES: .pdf .tex

%.pdf: %.tex
	@echo "        PDF $<"
	pdflatex -shell-escape $<
	bibtex $(shell echo $< | sed -e 's/.tex//')
	pdflatex -shell-escape $<
	pdflatex -shell-escape $<

default: manual.pdf

manual.tex: version.tex
version.tex:
	echo $(shell bzr revision-info | awk '{print $$1}') > version.tex
