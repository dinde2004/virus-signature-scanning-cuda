#!/bin/bash

for i in {1..30}
do
    diff tests/out_${i}.out tests/out_${i}.ans
done