#!/bin/bash
for i in $(find gbasis | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
rm -v gbasis/cext*.cpp
rm -v gbasis/*/cext*.cpp
exit 0
