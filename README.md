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
- Installation of the following packges: pandas, boto3, pyathena and time.

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

Configurating AWS
```
# User must import their AWS credentials here
%export AWS_ACCESS_KEY_ID="your_access_key_id"
%export AWS_SECRET_ACCESS_KEY="your_secret_access_key"
%export AWS_DEFAULT_REGION="your_aws_region"

CLIENT = boto3.client("athena", aws_access_key_id=AWS_ACCESS_KEY,
    aws_secret_access_key=AWS_SECRET_KEY,
    region_name=AWS_REGION,
)

DATABASE_NAME = "ncbi-vcf-codeathon-rc-db1"
RESULT_OUTPUT_LOCATION = "s3://ncbi-vcf-codeathon-rc-athena/"
```

Creating and executing the query
```
code
```


## Results


## Future Work
Calculation of additional subtractions sets from haplotypes to produce a better drop filter that would result in sets with higher quality. 
Development of a method to subtract/explain the relationship between lineage defining mutations and the clusters of paired mutations found in the query. 
Application of the Markov Clustering algorithm.

