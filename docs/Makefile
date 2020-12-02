
.PHONY: documenter tutorial tutorial-html

all: tutorial documenter

documenter:
	julia --color=yes --project=. make.jl

tutorial:
	julia --color=yes --project=. tutorial.jl

tutorial-html:
	julia --color=yes --project=. tutorial-html.jl

clean:
	rm -Rf build
	rm -Rf src/tutorial
