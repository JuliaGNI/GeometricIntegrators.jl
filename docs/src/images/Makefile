pdf:
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

png:
	pdftoppm commuting-diagram.pdf              commuting-diagram.png              -png -r 150 -transp
	pdftoppm omega-orthogonal-projection.pdf    omega-orthogonal-projection.png    -png -r 150 -transp
	pdftoppm post-projection.pdf                post-projection.png                -png -r 150 -transp
	pdftoppm post-projection-unstable.pdf       post-projection-unstable.png       -png -r 150 -transp
	pdftoppm standard-projection.pdf            standard-projection.png            -png -r 150 -transp
	pdftoppm symmetric-projection-minus.pdf     symmetric-projection-minus.png     -png -r 150 -transp
	pdftoppm symmetric-projection-plus.pdf      symmetric-projection-plus.png      -png -r 150 -transp
	pdftoppm symplectic-projection-unstable.pdf symplectic-projection-unstable.png -png -r 150 -transp
	pdftoppm variation-continuous.pdf           variation-continuous.png           -png -r 150 -transp
	pdftoppm variation-discrete.pdf             variation-discrete.png             -png -r 150 -transp

clean:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.pdf

empty: clean
	rm -f *.png
	rm -f *.svg

all: pdf
