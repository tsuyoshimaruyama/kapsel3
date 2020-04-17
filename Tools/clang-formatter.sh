#!/bin/bash
for f in *.cxx; do
    echo $f
    clang-format -style=file $f > .tempfile
    mv .tempfile $f
done

for f in *.h; do
    echo $f
    clang-format -style=file $f > .tempfile
    mv .tempfile $f
done
