VIGS= intro.html changeADI.html

all: $(VIGS)

clean:
	rm -rf *.md figure/

%.html: %.Rmd
	$(R_HOME)/bin/Rscript -e "knitr::knit2html('$*.Rmd', options = c('toc', markdown::markdownHTMLOptions(TRUE))"
