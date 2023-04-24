#!/bin/bash

for chr in {1..19} 
do
    sbatch scripts/GMvSTITCH_genoprobs.sh chr
done