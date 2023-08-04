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
),
variations0 AS (
    SELECT run, CONCAT(ref, CAST(pos as varchar), alt) AS ntvariation, CONCAT(protein_name, ': ', variation) AS aa_change, POS, REF, ALT, ref_aa, alt_aa, protein_name, variation
    FROM "annotated_variations"
    WHERE G_AD_2 / DP >= 0.5 AND G_AD_2 >= 50 AND DP >= 100
        AND run IN (SELECT run FROM meta)
        AND protein_name = 'S'
        AND alt_aa IS NOT NULL
),
variations AS (
    SELECT variations0.run, POS, ntvariations AS variation, aa_change
    FROM variations0
), records AS (
    SELECT COUNT(DISTINCT run) AS runs
    FROM variations
),
variation_totprobs AS (
    SELECT variation,
           COUNT(DISTINCT run) AS var_tot,
           COUNT(DISTINCT run) / (SELECT runs FROM records) AS var_prob
    FROM variations
    GROUP BY variation
)
      """

response = CLIENT.start_query_execution(
    QueryString=query,
    QueryExecutionContext={
        'Database': DATABASE_NAME
    },
    ResultConfiguration={
        'OutputLocation': RESULT_OUTPUT_LOCATION
    }
)

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

## Results


## FutureWork
Calculation of additional subtractions sets from haplotypes to produce a better drop filter that would result in sets with higher quality. 
Development of a method to subtract/explain the relationship between lineage defining mutations and the clusters of paired mutations found in the query. 
Application of the Markov Clustering algorithm.

