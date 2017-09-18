export PYTHONPATH=`pwd`:$PYTHONPATH
if [ -d build_docs ]; then
    echo "build_docs directory exists, delete before running"
    exit 1
fi
sphinx-apidoc gillespy --full -o build_docs -H 'GillesPy' -A 'John H. Abel, Brian Drawert, Andreas Hellander' -V '1.1'
cd build_docs && make html
sed -i '' 's/gillespy.gillespy/gillespy/g' _build/html/*.html
cp -r _build/html/* ../docs/
