# Data files used

## LTI_mutations.xlsx
Contains five sheets of the SARS-CoV-2 defining mutations and one sheet with the references from which the mutations were selected.
- NSPs
- Spike
- N
- E
- M
- References

## NSPs.csv
Dataset consisting of the defining mutations affecting the following genes:  ORF1b, ORF1ab, ORF3a, ORF6, ORF7a and ORF8. This dataset was extracted from the first sheet of the file LTI_mutations.xlsx.

## Spike.csv
Dataset consisting of the defining mutations affecting the genes encoding for the S protein. This dataset was extracted from the second sheet of the file LTI_mutations.xlsx.

## N.csv
Dataset consisting of the defining mutations affecting the genes encoding for the N protein. This dataset was extracted from the third sheet of the file LTI_mutations.xlsx.

## E.csv
Dataset consisting of the defining mutations affecting the genes encoding for the E protein. This dataset was extracted from the fourth sheet of the file LTI_mutations.xlsx.

## M.csv
Dataset consisting of the defining mutations affecting the genes encoding for the M protein. This dataset was extracted from the fifth sheet of the file LTI_mutations.xlsx.

# Data frames containing mutation pairs and relevant statistics (e.g. pearson correlation) generated from the Athena query of SARS-CoV-2 variant calls
## mutation_pairs.spike.nodatefilter.csv
Subset of mutation pairs from the spike protein only. 

## mutation_pairs.spike.non_syn.nodatefilter.csv
Subset of mutation pairs from non-synonymous mutations from the spike protein only.

## mutation_pairs.wholegenome.nodatefilter.csv
Mutation pairs across the whole SARS-CoV-2 genome.

## mutation_pairs.wholegenome.non_syn.nodatefilter.csv
Non-synonymous mutations pairs across the whole SARS-CoV-2 genome.
