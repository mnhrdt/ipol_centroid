DOC	= article
TEX	= $(wildcard *.tex)
BIB	= $(wildcard *.bib)
LATEX	= pdflatex -interaction nonstopmode -halt-on-error
BIBTEX	= bibtex
DPILR	= 96
PS2PDF	= ps2pdf -dCompatibilityLevel=1.4 -dEmbedAllFonts=true \
		-dAutoRotatePages=/None -dUseCIEColor=true
PS2PDFLR	= $(PS2PDF) -dPDFSETTINGS=/printer \
		-dColorImageResolution=$(DPILR) \
		-dGrayImageResolution=$(DPILR) \
		-dMonoImageResolution=$(DPILR)

default:	$(DOC).pdf

# LaTeX->PDF compilation
%.pdf	: %.tex $(TEX) $(BIB)
	$(LATEX) $*.tex
	-$(BIBTEX) $*.aux
	$(LATEX) $*.tex
	$(LATEX) $*.tex
	$(MAKE) $@~tag

# PDF filtering
%_lr.pdf	: %.pdf
	$(PS2PDFLR) $< $@
	pdfopt $@ $@.opt && mv $@.opt $@
	$(MAKE) $@~tag

# PDF tagging
# TODO: use XMP Metadata tagging instead
%.pdf~tag	: %.pdf meta.txt
	pdftk $*.pdf update_info_utf8 meta.txt output _$*.pdf \
		&& mv -f _$*.pdf $*.pdf

# cleanup
.PHONY	: clean distclean
clean	:
	$(RM) $(DOC).aux $(DOC).log $(DOC).out $(DOC).blg $(DOC).bbl
distclean	: clean
	$(RM) $(DOC).pdf $(DOC)_lr.pdf
