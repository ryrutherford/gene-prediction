#!/bin/bash
grep -w 'CDS' $1 | cut -f1-5 > parsed_gff3_1.txt