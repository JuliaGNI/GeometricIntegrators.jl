
all:
	xelatex -shell-escape commuting-diagram
	xelatex -shell-escape omega-orthogonal-projection
	xelatex -shell-escape post-projection
	xelatex -shell-escape post-projection-unstable
	xelatex -shell-escape standard-projection
	xelatex -shell-escape symmetric-projection-minus
	xelatex -shell-escape symmetric-projection-plus
	xelatex -shell-escape symplectic-projection-unstable
	xelatex -shell-escape variation-continuous
	xelatex -shell-escape variation-discrete

clean:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.pdf

empty: clean
	rm -f *.png
	rm -f *.svg
