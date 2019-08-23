#!/bin/sh

sphinx-apidoc -f -o . ../
make html
