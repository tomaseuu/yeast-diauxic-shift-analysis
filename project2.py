import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np
import os
from scipy.stats import binom
np.random.seed(42)

os.environ["OMP_NUM_THREADS"] = "1"  # visual code told me to put this in case of memory leaks

#1. Data: in the github data/diauxic_raw_ratios.txt

diauxic_data = pd.read_csv('data/diauxic_raw_ratios.txt', sep='\t')

print(diauxic_data.head())

print("Shape of data:", diauxic_data.shape)



# 2. Plot the raw data as a heatmap

expression_vals = diauxic_data.drop(columns= ["ORF", "Name"])

plt.figure(figsize=(12, 8))
heatmap = sns.heatmap(expression_vals, cmap="viridis", vmin=0, vmax=5, yticklabels=False)

plt.title("Gene Expression Data Heatmap", fontsize= 20, fontweight= "bold")
plt.xlabel("Time Points", fontsize= 16, fontweight= "bold")
plt.ylabel("Genes", fontsize= 16, fontweight= "bold")
plt.xticks(fontsize=12, fontweight='bold')
colorbar = heatmap.collections[0].colorbar
colorbar.set_label('Expression Rate', fontsize=14, fontweight='bold')
plt.show()


# 3. Cluster the genes and plot the clustered matrix

# choose how many clusters
num_clusters = 4

kmeans = KMeans(n_clusters=num_clusters, random_state= 42)
gene_clusters = kmeans.fit_predict(expression_vals)

expression_vals["Cluster"] = gene_clusters

# calculate silhouette score for unnormalized data
score = silhouette_score(expression_vals.drop(columns=["Cluster"]), gene_clusters)
print(f"Silhouette Score: {score: .3f}")

# create clustered heatmap
clustered_data = expression_vals.sort_values(by= "Cluster")
clustered_data_no_cluster = clustered_data.drop(columns=["Cluster"])

plt.figure(figsize = (12, 8))

heatmap = sns.heatmap(clustered_data_no_cluster, cmap="viridis", vmin=0, vmax=5, yticklabels=False)

plt.title("Clustered Gene Expression Heatmap", fontsize= 20, fontweight= "bold")
plt.xlabel("Time Points", fontsize = 16, fontweight= "bold")
plt.ylabel("Genes (Grouped by Cluster)", fontsize = 16, fontweight= "bold")
plt.xticks(fontsize= 12, fontweight= "bold")

colorbar = heatmap.collections[0].colorbar
colorbar.set_label("Expression Rate", fontsize=14, fontweight= "bold")

plt.show()


# 4. Remove genes not of interest (non-varying)

std_devs = expression_vals.drop(columns=["Cluster"]).std(axis=1)

expression_vals["StdDev"] = std_devs
sorted_std = expression_vals.sort_values(by= "StdDev", ascending=False)

top_genes = sorted_std.head(230)    # keep 230 most variable genes
top_expressions = top_genes.drop(columns=["StdDev", "Cluster"])

print(top_expressions.shape)
print(top_expressions.head())

# 5. Find the authors' set of 230 most variable genes

authors_data =pd.read_csv('data/230genes_log_expression.txt', sep='\t')

print(authors_data.head())
print("Shape of authors' data:", authors_data.shape)

# obtain the sets of gene IDs
my_genes = set(diauxic_data.loc[top_expressions.index, "ORF"])
authors_genes = set(authors_data["ORF"])

intersection = my_genes.intersection(authors_genes)
union = my_genes.union(authors_genes)

# Jaccard coefficient
jaccard_coeff = len(intersection) / len(union)

print(f"Number of overlapping genes: {len(intersection)}")
print(f"Jaccard Coefficient: {jaccard_coeff:.3f}")


# 6. Redo clustering analysis using the authors' 230 genes

authors_expressions = authors_data.drop(columns=["Name"])

plt.figure(figsize=(12, 8))
heatmap = sns.heatmap(authors_expressions.drop(columns=['ORF']), cmap="flare", yticklabels=False)

plt.title("Authors' 230 Genes - Raw Heatmap", fontsize=20, fontweight="bold")
plt.xlabel("Time Points", fontsize=16, fontweight="bold")
plt.ylabel("Genes", fontsize=16, fontweight="bold")
plt.xticks(fontsize=12, fontweight='bold')
colorbar = heatmap.collections[0].colorbar
colorbar.set_label('Expression Rate', fontsize=14, fontweight='bold')
plt.show()


#REDO clustering

num_clusters_authors = 4

kmeans_authors = KMeans(n_clusters=num_clusters_authors, random_state=42)
gene_clusters_authors = kmeans_authors.fit_predict(authors_expressions.drop(columns=['ORF']))

authors_expressions['Cluster'] = gene_clusters_authors
clustered_authors_data = authors_expressions.sort_values(by='Cluster')

clustered_authors_no_cluster = clustered_authors_data.drop(columns=['Cluster', 'ORF'])

plt.figure(figsize=(12, 8))
heatmap = sns.heatmap(clustered_authors_no_cluster, cmap="flare", yticklabels=False)

plt.title("Authors' 230 Genes - Clustered Heatmap", fontsize=20, fontweight="bold")
plt.xlabel("Time Points", fontsize=16, fontweight="bold")
plt.ylabel("Genes (Grouped by Cluster)", fontsize=16, fontweight="bold")
plt.xticks(fontsize=12, fontweight="bold")
colorbar = heatmap.collections[0].colorbar
colorbar.set_label('Expression Rate', fontsize=14, fontweight='bold')
plt.show()




# 2.7 Re-clustering with my 230 genes!

# using MY top 230 genes
num_clusters_mine = 4

# cluster my 230-gene expression matrix
kmeans_mine = KMeans(n_clusters=num_clusters_mine, random_state=42)
gene_clusters_mine = kmeans_mine.fit_predict(top_expressions)

top_expressions["Cluster"] = gene_clusters_mine
clustered_mine_data = top_expressions.sort_values(by='Cluster')

clustered_mine_no_cluster = clustered_mine_data.drop(columns=["Cluster"])

plt.figure(figsize=(12, 8))
heatmap = sns.heatmap(clustered_mine_no_cluster, cmap="vlag", yticklabels=False, robust=True)

plt.title("My 230 Genes - Clustered Heatmap", fontsize=20, fontweight="bold")
plt.xlabel("Time Points", fontsize=16, fontweight="bold")
plt.ylabel("Genes (Grouped by Cluster)", fontsize=16, fontweight="bold")
plt.xticks(fontsize=12, fontweight="bold")
colorbar = heatmap.collections[0].colorbar
colorbar.set_label('Expression Rate', fontsize=14, fontweight='bold')
plt.show()











# PART 2: Evaluation and statistical intuition

# 1. Test if your list of genes contains GO terms that are overrepresented using the established GO tool: https://geneontology.org/
# obtain gene names of top 230 genes
gene_names = diauxic_data.loc[top_expressions.index, 'Name']
# drop missing names
gene_names = gene_names.dropna()
for name in gene_names:
    print(name)

# 2. Calculate the probability that your list of items (genes) has more Green genes than expected by chance
n = 50    
p = 0.71

p_val = binom.sf(37, n, p)   # sf = survival function = probability of getting more than 37 successes
print("P-value:", p_val)

# 3.1 Run a simulation to estimate p-value

n_genes_total = 6000  # posslbe genes
green_prob = 0.71
sample_size = 50      # randomly draw 50 genes
threshold = 38
X = 10000   # lets repeat it 10000 times

# create fake genome
# 71% chance of being Green (1 means Green)
# 29% chance of being not green (0 means not green)
# but within 6000 entries

n_green_total = int(n_genes_total * green_prob)
genome = np.array([1] * int(n_genes_total * green_prob) + [0] * int(n_genes_total * (1 - green_prob)))

successes = 0 # counter : how many times did I get 38 or more green genes?
rng = np.random.default_rng(42)  # reproducible random generator so its consistent

for _ in range(X):
    sample = np.random.choice(genome, size=sample_size, replace=False)
    green_count = np.sum(sample)
    if green_count >= threshold:
        successes +=1

simulated_p_val = successes / X
print(f"Simulated P-value after {X} trials: {simulated_p_val:.5f}")


# 4.2 How many simulations (X) do you need to run to be certain of the probability

simulation_steps = [100, 500, 1000, 2000, 5000, 10000, 50000, 80000]
estimated_p_vals = []

for num_trials in simulation_steps:
    successes = 0
    for _ in range (num_trials):
        sample = np.random.choice(genome, size=sample_size, replace=False)
        green_count = np.sum(sample)
        if green_count >= threshold:
            successes += 1
    estimated_p = successes / num_trials
    estimated_p_vals.append(estimated_p)
    print(f"Trials: {num_trials}, Estimated p-value: {estimated_p:.5f}")

plt.figure(figsize=(10,6))
plt.plot(simulation_steps, estimated_p_vals, marker='o') # marker places a circle at each data point
plt.xlabel('Number of Simulations (X)', fontsize=14)
plt.ylabel('Estimated P-value', fontsize=14)
plt.title('P-value Stability with Increasing Simulations', fontsize=16, fontweight="bold")
plt.grid(True)
plt.show()


# Bonus:

# 2.3 Is Green appearing more frequently than expected? 
# What if there were 35 or 40 Green items in the list rather than 38?

n = 50
p = 0.71
p_vals_38 = binom.sf(37, n, p)
print("P-value for 38 Greens:", p_vals_38)

p_vals_35 = binom.sf(34, n, p)
print("P-value for 35 Greens:", p_vals_35)

p_vals_40 = binom.sf(39, n, p)
print("P-value for 40 Greens:", p_vals_40)

# 2.4 Find a pair of numbers (n,k) such that if n 
# (number of items in the background) and 
# k (number of items in the research list) are annotated with Purple, 
# a student working on this problem will calculate a p-value 
# between (0.003 and 0.01) for Purple

p = 0.1 
n_range = range(30, 101)  
k_range = range(1, 30)    

for n in n_range:
    for k in k_range:
        p_value = binom.sf(k-1, n, p)  # P(X >= k)
        if 0.003 <= p_value <= 0.01:
            print(f"n = {n}, k = {k}, p-value = {p_value:.5f}")


