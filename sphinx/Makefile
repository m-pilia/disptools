# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = disptools
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	rm -rf _build/*
	sphinx-apidoc -f -o source ../disptools
	python3 -m m2r ../README.md # convert GitHub README to rst
	mv ../README.rst readme.rst
	sed -i 's,https://.*sphinx/\([^?]\+\)?raw=true,sphinx/\1,g' readme.rst # convert image paths
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	for f in `ls _build/html | grep html`; do sed -i 's/_images/images/g;s/_modules/modules/g;s/_sources/sources/g;s/_static/static/g' _build/html/$$f; done
	for d in images modules sources static; do mv _build/html/_$$d _build/html/$$d; done
	rm -rf ../{images,modules,sources,static}/* || :
	cp -R _build/html/* ../.
