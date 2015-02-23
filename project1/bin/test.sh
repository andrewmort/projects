#!/bin/bash

for i in `find $2 -name "*.cnf"`; do
    echo
    echo "$i:"
    echo
    $1 $i
    echo
done

