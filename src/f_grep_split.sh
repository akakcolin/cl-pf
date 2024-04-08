#!/bin/bash

new_dir=Pbed-N2done

for x in `grep Thank *out | awk '{print $1}' | perl -ne '{ @file=split(/\./,$_); print "$file[0]\n";}'`
  do
    echo $x
    echo "mv $x* $new_dir"
  done

