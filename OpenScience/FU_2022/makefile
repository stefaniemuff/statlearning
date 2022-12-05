FILE=CBD_talk_steeringBoard

all: CBD_talk_steeringBoard.Rnw makefile 
	echo "library('knitr'); knit2pdf('CBD_talk_steeringBoard.Rnw')" | R --no-save --no-restore

	pdfjam --vanilla --quiet --nup 1x2 --no-landscape CBD_talk_steeringBoard.pdf

clean: 
	rm $(FILE).bbl $(FILE).aux $(FILE).log $(FILE).out $(FILE).blg 
