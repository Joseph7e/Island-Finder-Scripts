# Island-Finder-Scripts

This repository contains python and shell scripts that are used to locate pathogenecity islands from genome assemblies. The input to the program is a reference island FASTA and FASTA sequences for one or more genome assemblies. Specific settings can be changed within the scripts.

island_blast.sh - automates the blast with a specific format and runs the other two scripts.
island_finder.py - parses the BLAST and performs quality control
arrange_nodes.py - deals with things like recombination, and broken contigs. This script produces final island sequences and report file.

parse_island_blast.py -  a seperate script used to run BLAST and construct a xls table with precense/absence info about the pathogenetic islands.
