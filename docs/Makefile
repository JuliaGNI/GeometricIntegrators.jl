
.PHONY: documenter images tutorial tutorial-html

all: tutorial images documenter

documenter:
	julia --color=yes --project=. make.jl

images:
	$(MAKE) all   -C src/images
	$(MAKE) clean -C src/images

tutorial:
	julia --color=yes --project=. tutorial.jl

tutorial-html:
	julia --color=yes --project=. tutorial-html.jl

clean:
	$(MAKE) empty -C src/images
	rm -Rf build
	rm -Rf src/tutorial
