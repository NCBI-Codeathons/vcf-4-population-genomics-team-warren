# Team Linkage

## Team Members
List of participants and affiliations:
- Andrew Warren, University of Virginia (Team Leader)
- Cristina Delgado, Association of Public Health Laboratories
- David Stern, National Biodefense Analysis and Countermeasures Center
- Efe Sezgin, Izmir Institute of Technology
- Kaitlin Zajac, Rutgers University
- Katherine Ramos, Gibbons Lab
- Saber Tadros, National Institutes of Health

## Project Goals
Creation of a linkage landscape [dataset + visualization] for SARS-COV-2 using methods for identification of linked tuples of positions with an associated significance score.

## Approach
**Prerequisites**
- Connection to the Athena database.
- Python environment/python wrapper to make a SQL query against the Athena database.
- Pandas package.

**Methods**
Usage of the Pearson correlation coefficient value of ___ as a cutoff for the pairs that are found. Visualization with a graph of the number of samples involved per pair.
Conversion of the sets of positions into a GFF3 tracks.
Utilization of GFF3 tracks to visualize the linked positions per track on JBrowse, e.g., loading 10 clusters into a GFF3 to see them in the Jbrowser.

## Results

## Future Work
Calculation of additional subtractions sets from haplotypes to produce a better drop filter that would result in sets with higher quality. 
Development of a method to subtract/explain the relationship between lineage defining mutations and the clusters of paired mutations found in the query. 
Application of the Markov Clustering algorithm.

