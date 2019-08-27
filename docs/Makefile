## ============================================================================
## @file    Makefile
## @brief   Makefile for creating docs
## @date    2019-08-27
## @license Please see the file named LICENSE in the project directory
## @website https://github.com/GillesPy2/GillesPy2
##
## The following are special commands recognized by this Makefile.
##
##   make help
##       Print this summary of the available commands.
##
##   make html
##   make docs
##       (Alternatives for the same thing.)  Regenerates the formatted HTML
##       output if necessary.  If none of the input files have been changed,
##       it does nothing.  Unlike "make autobuild" below, this only builds
##       the documentation once and then exits.
##
##   make autobuild
##
##   make OTHER
##       Pass OTHER directly to Sphinx.  For example, "make singlehtml"
##       runs the Sphinx command "singlehtml".
## ============================================================================

# Variables.  (These can also be set from the command line.)
# .............................................................................

SPHINXOPTS    = -j auto
SPHINXBUILD   = python3 -m sphinx
SOURCEDIR     = .
BUILDDIR      = build


# Main build targets.
# .............................................................................

help:
	@echo 'This is the Makefile for building the GillesPy2 documenation.'
	@echo 'Available commands:'
	@echo ''
	@echo '  make help'
	@echo '    Print this summary of available commands.'
	@echo ''
	@echo '  make html'
	@echo '  make docs'
	@echo '    (Both do the same action.) Regenerates the formatted HTML'
	@echo '    output if necessary. If none of the input files have been'
	@echo '    changed, it does nothing. Unlike "make autobuild", this'
	@echo '    only builds the documentation once and then exits.'
	@echo ''
	@echo '  make autobuild'
	@echo '    Runs sphinx-autobuild, which watches the file system and'
	@echo '    automatically regenerates the HTML documentation if any'
	@echo '    files are changed. Open your browser to localhost:8000'
	@echo '    to see the output refreshed automatically.'
	@echo ''
	@echo '  make OTHER'
	@echo '    Pass OTHER directly to Sphinx.  For example,'
	@echo '    "make singlehtml" runs the Sphinx command "singlehtml".'
	@echo ''

html docs:
	sphinx-apidoc -n -f -o "$(SOURCEDIR)" ../
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

autobuild:
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html"

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


# Miscellaneous directives.
# .............................................................................

.PHONY: help html docs build Makefile
