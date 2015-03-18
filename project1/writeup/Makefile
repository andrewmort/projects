M4_DIR := circuits
TEX_FILE := $(wildcard *.tex)
NAME := $(patsubst %.tex, %, $(TEX_FILE))
M4_FILES := $(wildcard $(M4_DIR)/*.m4)
MPS_FILES := $(patsubst %.m4, %.mps, $(M4_FILES))
LOG_FILES := $(patsubst %.m4, %.log, $(M4_FILES))

LIBDIR := /usr/local/texlive/texmf-local/circuit_macros
M4MACROS := mpost.m4
M4 := m4 -I $(LIBDIR)
DPIC := dpic
TR := | tr '=' '\134'

.PHONY: all clean

all: $(NAME).pdf

$(NAME).pdf: $(MPS_FILES) $(NAME).tex
	pdflatex $(NAME).tex
	pdflatex $(NAME).tex

%.mps: %.m4
	$(eval FILE_DIR = $(dir $*))
	$(eval FILE_NAME = $(notdir $*))
	@echo "Writing $*.mp header"
	@(printf 'verbatimtex\n' ;\
     printf '=documentclass[11pt]{article}\n' ;\
     printf '=usepackage{times,boxdims,graphicx}\n' ;\
     printf '=begin{document} etex\n' ) $(TR) > $*.mp
	cd $(FILE_DIR) && $(M4) $(M4MACROS) $(FILE_NAME).m4 | $(DPIC) -s >> $(FILE_NAME).mp
	@printf 'end\n' >> $*.mp
	cd $(FILE_DIR) && mpost --tex=latex $(FILE_NAME).mp
	mv $*.1 $*.mps
	-rm $*.mp $*.mpx

clean:
	@-rm $(MPS_FILES) $(LOG_FILES) $(NAME).out $(NAME).log $(NAME).aux $(NAME).toc &> /dev/null
