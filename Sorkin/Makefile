all: writeup.pdf

clean:
	rm *.aux *.log *.toc *.out

writeup.pdf: writeup.tex $(logbookentries) $(entries_inputs)
	pdflatex -draftmode $<
	bibtex $(basename $<).aux
	pdflatex -draftmode $<
	pdflatex $<
	rm $(basename $<).log $(basename $<).out $(basename $<).toc $(basename $<).aux
	rm $(basename $<).bbl $(basename $<).blg

input:
	mkdir $@
