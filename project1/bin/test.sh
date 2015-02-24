#!/bin/bash

for i in `find $2 -name "*.cnf" -name "*sanity*"`; do
    echo
    echo "$i:"
    echo
    $1 $i
    echo
done

