FILE=CASOslo_Muff

all: CASOslo_Muff.Rnw makefile 
	echo "library('knitr'); knit2pdf('CASOslo_Muff.Rnw')" | R --no-save --no-restore

	pdflatex CASOslo_Muff
	bibtex CASOslo_Muff
	pdflatex CASOslo_Muff
	pdfjam --vanilla --quiet --nup 1x2 --no-landscape CASOslo_Muff.pdf

clean: 
	rm $(FILE).bbl $(FILE).aux $(FILE).log $(FILE).out $(FILE).blg 
