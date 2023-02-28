#!/bin/sh
wget https://raw.githubusercontent.com/catchorg/Catch2/master/single_include/catch2/catch.hpp -O ../include/catch.hpp
make clean -C ../test
