# Pylint coverage badge
# https://stackoverflow.com/a/47342642
#!/bin/bash
# pip3 install pylint anybadge

RATING=$(pylint ./gillespy2 | sed -n 's/^Your code has been rated at \([-0-9.]*\)\/.*/\1/p')
anybadge -o -l pylint -v $RATING -s '/10' -f .graphics/pylint.svg 2=red 4=orange 8=yellow 10=green
