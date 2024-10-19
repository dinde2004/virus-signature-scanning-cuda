#!/bin/bash

for i in {0..9}
do
    diff tests/out_${i}.out tests/out_${i}.ans
done