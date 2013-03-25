#!/bin/bash
samtools view -bS -o accepted_hits_round2.bam Aligned.out.sam
rm Aligned.out.sam
python ../../AssignReads2.py accepted_hits_round2.bam
samtools sort assigned_dmel.bam assigned_dmel_sorted
samtools reheader ../mel_only.header.sam assigned_dmel_sorted.bam > assigned_dmel_round2.bam
rm assigned_dmel_sorted.bam
