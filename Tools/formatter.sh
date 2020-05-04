#!/bin/bash
ls *.h | xargs dos2unix
ls *.cxx | xargs dos2unix
clang-format -style=file -i -verbose *.h
clang-format -style=file -i -verbose *.cxx

