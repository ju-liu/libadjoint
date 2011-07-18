.SUFFIXES: .pdf .tex

%.pdf: %.tex
	@echo "        PDF $<"
	pdflatex -shell-escape $<
	bibtex $(shell echo $< | sed -e 's/.tex//')
	pdflatex -shell-escape $<
	pdflatex -shell-escape $<

default: manual.pdf

manual.tex: version.tex overview.tex about.tex introduction.tex
	touch manual.tex

version.tex:
	echo $(shell bzr revision-info | awk '{print $$1}') > version.tex

bib:
	cp /data/pfarrell/svn/literature/literature.bib .
