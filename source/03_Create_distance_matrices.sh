#!/usr/bin/env bash
FastTreeMP -makematrix data/archaea_msa.fasta > intermediate/archaea.dist
FastTreeMP -makematrix data/bacteria_msa.fasta > intermediate/bacteria.dist
