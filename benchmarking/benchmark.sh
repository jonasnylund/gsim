#!/bin/bash

common_args="-p 500 -n 300 --epsilon 0.05 -q"
name=$1
args=$2
loops=$3

echo "${name}: ${args}" ${common_args}
mkdir -p $name
pushd $name

# echo "${name}: ${args}" ${common_args} >> "stats.txt"
for ((i=1; i<=$loops; i++)); do
  echo -ne "\r${i}"
  ../../../gsim $common_args $args >> "stats.txt"
  sleep 5
done
echo -e "\r  "

popd