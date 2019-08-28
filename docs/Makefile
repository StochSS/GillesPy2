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
##       Runs sphinx-autobuild, which watches the file system and
##       automatically regenerates the HTML documentation if any
##       files are changed. Open your browser to localhost:8000
##       to see the output refreshed automatically.
##  
##  make clean
##       Delete the build directory.
##  
##  make publish
##       FOR MAIN REPO MAINTAINERS ONLY: do a "make clean", then
##       a "make html", then commit the changes to the gh-pages
##       branch, and finally, force-push that branch to GitHub.
##
## ============================================================================

# Preliminary setup.
# .............................................................................
# The following setup is needed to enable shell traps to work.
# However, the .ONESHELL directive requires GNU make version 3.82 or later.
# The clever test below is from https://stackoverflow.com/a/12231321/743730

ifeq ($(filter undefine,$(value .FEATURES)),)
    $(error Unsupported version of Make. Please use GNU Make 3.82 or above.)
endif

export SHELL:=/bin/bash
export SHELLOPTS:=$(if $(SHELLOPTS),$(SHELLOPTS):)pipefail:errexit
.ONESHELL:


# Variables.  (These can be overriden via command-line args.)
# .............................................................................

SPHINXOPTS  = -j auto
SPHINXBUILD = python3 -m sphinx
SOURCEDIR   = .
FORMAT      = html
BUILDDIR    = build/$(FORMAT)


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
	@echo '  make clean'
	@echo '    Delete the build directory.'
	@echo ''
	@echo '  make publish'
	@echo '    FOR MAIN REPO MAINTAINERS ONLY: do a "make clean", then'
	@echo '    a "make html", then commit the changes to the gh-pages'
	@echo '    branch, and finally, force-push that branch to GitHub.'
	@echo ''

$(FORMAT) docs: apidocs
	@$(SPHINXBUILD) -b $(FORMAT) "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

autobuild: $(FORMAT)
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/$(FORMAT)"

clean:
	rm -rf "$(BUILDDIR)"

publish: | clean $(FORMAT)
	function restoreGitBranch {
	    git checkout @{-1}
	}
	git checkout gh-pages
	-git add -f "$(BUILDDIR)"
	-git commit -m "Latest build." "$(BUILDDIR)"
	git push origin gh-pages -f
	trap restoreGitBranch EXIT

apidocs: build-dir
	sphinx-apidoc -f -T -o classes ../gillespy2

build-dir:
	mkdir -p build/$(FORMAT)


# Miscellaneous directives.
# .............................................................................

.PHONY: help html docs autobuild publish clean apidocs Makefile
