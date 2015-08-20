#!/bin/sh
echo "library(tools); Sweave(\"$1\")" | Rterm --vanilla
