# Team Linkage
![Logo](logo.png)

## Table of Contents

1. [Team](#team)
2. [Goals](#goals)
3. [Methods](#methods)
4. [Results](#results)
5. [FutureWork](#futurework)


## Team
List of participants and affiliations:
- Andrew Warren, University of Virginia (Team Leader)
- Cristina Delgado, Association of Public Health Laboratories (Writer)
- David Stern, National Biodefense Analysis and Countermeasures Center (Tech Lead)
- Efe Sezgin, Izmir Institute of Technology,
- Kaitlin Zajac, Rutgers University,
- Katherine Ramos, Gibbons Lab at the Institute for Systems Biology,
- Saber Tadros, National Institutes of Health,

## Goals
Creation of a linkage landscape [dataset + visualization] for SARS-COV-2 using methods for identification of linked tuples of positions with an associated significance score.

## Methods
### Prerequisites
- Understanding of the concepts of: viral haplotypes, quasi-species, linkage disequilibrium, lineage defining mutations, epistatic mutations, etc.
- Connection to the Athena database.
- Python environment to make a SQL query against the Athena database.
- Installation of the following packges: pandas, boto3, pyathena, time, markov_clustering, xlrd, numpy, scikit-learn and scipy.

### Procedure
**Step 1. VCF DB query development**
1. Query to find meaningful correlation between mutations
2. Connect to athena.
3. Convert big query SQL to athena SQL.
4. Convert query to find generic position linkages.

**Step 2. Python connection to Athena (Boto)**
1. Demonstrated query to pandas data table.
2. Use converted big query to get pandas stored result.
3. Generate statistics on number of samples.
4. Graph of the number of samples involved per pair.

**Step 3. Cluster pairs of correlated mutations**
1. Markov clustering based on Pearson >0.1 < -0.1.
2. Deploy Silhouette score for understanding cluster significance.

**Step 4. Develop a method to subtract / explain the relationship between lineage defining mutations and the clusters of paired mutations (from query)**
1. Bipartite matching based on best overlap with lineage defining mutations.
2. Set threshold of overlap (heuristic based on distribution).

**Step 5. Check literature / group knowledge of “famous” mutations that are present in significant clusters**
1. Immunocompromised hallmarks.
2. Immune escape mutations.
3. SGTF (will deletion impact our variation calculation).

**Step 6. Visualization development**
1. Graph cluster visualization.
2. Convert python sets into GFF3 tracks for viewing on JBrowse.


### Query

**Running the code**

Installing the packages
```
# Installing necessary packages
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
# Creating the query
query="""WITH meta AS (
    SELECT acc AS run, biosample, bioproject
    FROM "metadata"
    WHERE cast(collection_date_sam[1] AS date) >= date_parse('2022-01-01', '%Y-%m-%d')
),
variations AS (
    SELECT run, POS, concat(ref, cast(pos AS varchar), alt) AS variation, concat(protein_name, cast(':' AS varchar), variation) AS aa_change
    FROM "annotated_variations"
    JOIN meta USING(run)
    WHERE G_AD_2/DP >= 0.5 AND G_AD_2 >= 50 AND DP >= 100 AND protein_name = concat('S','') AND alt_aa IS NOT NULL
),
records AS (
    SELECT count(distinct run) AS runs
    FROM variations
),
variation_totprobs AS (
    SELECT row_number() OVER(PARTITION BY NULL ORDER BY variation) AS rown, 
           variation,
           count(distinct run) AS var_tot,
           count(distinct run)/(SELECT runs FROM records) AS var_prob
    FROM variations
    GROUP BY variation
    ORDER BY variation    
),
var_pairs AS (
    SELECT v1.run, 
           v1.variation AS variation_1, 
           v2.variation AS variation_2, 
           v1.aa_change AS aa_change_1, 
           v2.aa_change AS aa_change_2, 
           v3.var_tot AS var_tot_1, 
           v4.var_tot AS var_tot_2, 
           v3.var_prob AS var_prob_1, 
           v4.var_prob AS var_prob_2
    FROM variations v1
    JOIN variations v2 ON v1.run = v2.run
    JOIN variation_totprobs v3 ON v1.variation = v3.variation
    JOIN variation_totprobs v4 ON v2.variation = v4.variation
    WHERE v1.variation != v2.variation AND v3.rown < v4.rown
),
results as (
    select v0.variation_1,
       v0.variation_2,
       aa_change_1,
       aa_change_2,
       var_tot_1 as var_1_count,
       var_tot_2 as var_2_count,
       least(var_tot_1, var_tot_2) as max_records_possible,
       count(distinct v0.run) as record_count,
       count(distinct biosample) as samples,
       count(distinct bioproject) as projects,
       count(distinct v0.run) / (select runs from records) as record_freq,
       (select runs from records) * count(distinct v0.run) / (var_tot_1 * var_tot_2) as rel_record_freq,
       count(distinct v0.run) / (var_tot_1) as con_prob_1,
       count(distinct v0.run) / (var_tot_2) as con_prob_2,
       count(distinct v0.run) / (var_tot_1 + var_tot_2 - count(distinct v0.run)) as jaccard,
       (((var_tot_1 - count(distinct v0.run))*((1 - var_prob_1)*(-1*var_prob_2)))+((var_tot_2 - count(distinct v0.run))*((-1*var_prob_1)*(1-var_prob_2)))+(count(distinct v0.run)*((1-var_prob_1)*(1-var_prob_2))))/sqrt((((var_tot_1 - count(distinct v0.run))*((1 - var_prob_1)*(1 - var_prob_2))+((var_tot_2 - count(distinct v0.run))*((-1*var_prob_1)*(-1*var_prob_2)))+(count(distinct v0.run)*((1-var_prob_1)*(1-var_prob_2)))))*((var_tot_1 - count(distinct v0.run))*((-1*var_prob_1)*(-1*var_prob_2))+((var_tot_2 - count(distinct v0.run))*((1-var_prob_1)*(1-var_prob_2)))+(count(distinct v0.run)*((1-var_prob_1)*(1-var_prob_2))))) as pearson
    from var_pairs v0
    join meta m on m.run = v0.run
    where var_prob_1 != 1
    group by variation_1, variation_2, var_tot_1, var_tot_2, var_prob_1, var_prob_2, aa_change_1, aa_change_2
)
select *
from results
where max_records_possible >100 and
        samples > 100 and
        projects >2 and 
        ((rel_record_freq > 1.5) or (rel_record_freq < 1/1.5) or 
        (jaccard > 1.5*con_prob_1*con_prob_2 or jaccard < 0.5*con_prob_1*con_prob_2 or con_prob_1>0.99 or con_prob_1<0.01 or con_prob_2>0.99 or con_prob_2<0.01 or jaccard>0.9 or jaccard<0.1 ) or
        (pearson > 1/10 or pearson < -1/10))
      """

# Getting the query execution ID
query_execution_id = response['QueryExecutionId']

# Function to check the status of the query execution
def is_query_running(query_execution_id):
    response = CLIENT.get_query_execution(QueryExecutionId=query_execution_id)
    state = response['QueryExecution']['Status']['State']
    return state in ['QUEUED', 'RUNNING']

# Waiting for the query to complete
while is_query_running(query_execution_id):
    time.sleep(5)  # Wait for 5 seconds before checking again
```

Obtaining the results
```
# Getting the results from Athena using get_query_results
response = CLIENT.get_query_results(QueryExecutionId=query_execution_id)

# Converting the results to a Pandas DataFrame
column_names = [col['Name'] for col in response['ResultSet']['ResultSetMetadata']['ColumnInfo']]
data_rows = [list(row['Data']) for row in response['ResultSet']['Rows'][1:]]
df = pd.DataFrame(data_rows, columns=column_names)

# Displaying the DataFrame
print(df)
```

### Markov Clustering

**Running the code**

Installing the packages
```
# Installing necessary packages
import pandas as pd
!pip install markov_clustering
!pip install xlrd
!pip install numpy
!pip install scikit-learn
import numpy as np
import markov_clustering as mc
import networkx as nx
import random
!pip install scipy
import scipy
```

## Results


## FutureWork
Calculation of additional subtractions sets from haplotypes to produce a better drop filter that would result in sets with higher quality. 
Development of a method to subtract/explain the relationship between lineage defining mutations and the clusters of paired mutations found in the query. 
Application of the Markov Clustering algorithm.

