#!/bin/bash

grep ">" "$1" | tr -d ">" | tr "/" "\t" | cut -f1

