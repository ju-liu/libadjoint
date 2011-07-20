.SUFFIXES: .pdf .tex

%.pdf: %.tex
	@echo "        PDF $<"
	pdflatex -shell-escape $<
	bibtex $(shell echo $< | sed -e 's/.tex//')
	makeindex $(shell echo $< | sed -e 's/.tex//')
	pdflatex -shell-escape $<
	pdflatex -shell-escape $<

default: manual.pdf

manual.tex: version.tex $(shell grep include manual.tex | grep -v ^% | awk -F '{' '{print $$2}' | awk -F '}' '{print $$1}' | sed 's/$$/.tex/')
	touch manual.tex

version.tex:
	echo $(shell bzr revision-info | awk '{print $$1}') > version.tex

bib:
	cp /data/pfarrell/svn/literature/literature.bib .
