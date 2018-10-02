#!/bin/bash
for i in $(find gbasis tools scripts | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
(cd doc; make clean)
rm -v MANIFEST
rm -vr dist
rm -vr build
rm -vr doctrees
rm -v gbasis/cext.cpp
rm -v gbasis/*/cext.cpp
rm -v .coverage
rm -v doc/xml/*
rm -vr doc/cxxapi/*
rm -vr doc/pyapi/*
rm -vr doc/_build/*
exit 0
