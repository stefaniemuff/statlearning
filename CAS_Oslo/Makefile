FILE=CASOslo_Muff

all:	
	R CMD Sweave $(FILE).Rnw;\
	pdflatex $(FILE);\
	bibtex $(FILE);\
	cp aa.bbl CASOslo_Muff.bbl ;\
	pdflatex $(FILE);\
	pdflatex $(FILE);\
	evince $(FILE).pdf &


# 
