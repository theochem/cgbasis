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
exit 0
