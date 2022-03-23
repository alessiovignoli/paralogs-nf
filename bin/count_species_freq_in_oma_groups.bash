#!/bin/bash

# Usage: count_species_freq_in_oma_groups.bash <OMA GROUPS>


sort -V $1 | uniq | awk '{$1=""; $2=""; $3=""; sub("  ", " "); print}' | tr "\n" " " | tr " " "\n" | awk '{print substr($0,1,5)}' | awk 'NF' | sort | uniq -c | sort -rnk 1