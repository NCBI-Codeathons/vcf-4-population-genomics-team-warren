# Team Linkage
![Logo](logo.png)

## Table of Contents

1. [TeamMembers](#teammembers)
2. [Methods](#methods)
3. [ProjectGoals](#projectgoals)
4. [Methods](#methods)
5. [Results](#results)
6. [FutureWork](#futurework)


## Team Members
List of participants and affiliations:
- Andrew Warren, University of Virginia (Team Leader)
- Cristina Delgado, Association of Public Health Laboratories (Writer)
- David Stern, National Biodefense Analysis and Countermeasures Center (Tech Lead)
- Efe Sezgin, Izmir Institute of Technology,
- Kaitlin Zajac, Rutgers University,
- Katherine Ramos, Gibbons Lab at the Institute for Systems Biology,
- Saber Tadros, National Institutes of Health,

## Project Goals
Creation of a linkage landscape [dataset + visualization] for SARS-COV-2 using methods for identification of linked tuples of positions with an associated significance score.

## Methods
**Prerequisites**
- Connection to the Athena database.
- Python environment to make a SQL query against the Athena database.
- Pandas package.

**Procedure**
1. Usage of the Pearson correlation coefficient value of ___ as a cutoff for the pairs that are found. Visualization with a graph of the number of samples involved per pair.
2. Conversion of the sets of positions into a GFF3 tracks.
3. Utilization of GFF3 tracks to visualize the linked positions per track on JBrowse, e.g., loading 10 clusters into a GFF3 to see them in the Jbrowser.

**Running the code**

Installing the packages
```
# Install packages
%pip install boto3 pandas pyathena
import pandas as pd
import boto3
import pyathena
import time
```

## Results


## Future Work
Calculation of additional subtractions sets from haplotypes to produce a better drop filter that would result in sets with higher quality. 
Development of a method to subtract/explain the relationship between lineage defining mutations and the clusters of paired mutations found in the query. 
Application of the Markov Clustering algorithm.

