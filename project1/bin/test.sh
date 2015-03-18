#!/bin/bash

#for i in `find $2 -name "*.cnf" -name "*sanity*"`; do
for i in `find $2 -name "*.cnf"`; do
    echo
    echo "$i:"
    echo
    time $1 $i
    echo
done

