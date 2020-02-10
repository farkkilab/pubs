import pandas as pd
import numpy as np
import os
import time
import glob
from multiprocessing import Process, freeze_support
## TRIBUS - spatial analysis of class-class attraction and avoidance scores
# Adapted from Jake's TMA analysis https://github.com/labsyspharm/hms-lincs-tma/blob/master/neighborhood_analysis_of_clusters.py
# /julia.casado
#

def get_neighbor_fractions(neighbor_idx, data):
    neighbors = data.reindex(neighbor_idx).dropna(
        how='all').cluster.value_counts()
    total_neighbors = neighbors.sum()
    neighbor_fractions = neighbors / total_neighbors
    return neighbor_fractions.sort_index()

def get_neighbors(data, exclude_self=True):
    """Adapted from Jake's TMA analysis code
    """
    neighbor_cols = [x for x in data.columns if 'neighbour' in x]

    cell_ids = data[neighbor_cols].values.flatten()
    cell_ids = cell_ids[(~np.isnan(cell_ids)) & (cell_ids != 0)]
    cell_ids = np.unique(cell_ids).astype(int)

    neighbors_idx = cell_ids
    if exclude_self:
        neighbor_idx = [
            x for x in neighbors_idx if x not in data.index]
    return neighbors_idx

def permutation_neighborhood(subset, data, num_permutations=1000, verbose=True):
    """Permutate cluster labels in the data table to get a null distribution of
       observing the neighbor by cluster profile by chance.
    """
    name = subset.name
    neighbors = get_neighbors(subset, exclude_self=True)
    true_neighbor_fractions = get_neighbor_fractions(neighbors, data)
    pval_higher = pd.Series(0, index=true_neighbor_fractions.index)
    pval_lower = pd.Series(0, index=true_neighbor_fractions.index)
    permutated_fractions = pd.Series(0, index=true_neighbor_fractions.index)

    t = time.time()
    for i in range(num_permutations):
        permutated_cluster_labels = np.random.permutation(data.cluster.values)
        permutated_data = data.copy()
        permutated_data.loc[:, 'cluster'] = permutated_cluster_labels
        permutated_neighbor_fractions = get_neighbor_fractions(
            neighbors, permutated_data)
        # this line might not be necessary any more
        permutated_neighbor_fractions = permutated_neighbor_fractions.reindex(
            true_neighbor_fractions.index, fill_value=0)
        # what is the likelihood that true data is higher than random?
        pval_higher += (true_neighbor_fractions <= permutated_neighbor_fractions)
        # what is the likelihood that true data is lower than random?
        pval_lower += (true_neighbor_fractions >= permutated_neighbor_fractions)
        permutated_fractions +=  permutated_neighbor_fractions
        if (i % 100 == 0) & (verbose):
            print('Iteration: {}'.format(str(i)))
            print(permutated_neighbor_fractions)
    print('Time', name, ':', time.time() - t,'seconds')

    fold_change = np.log2( true_neighbor_fractions / ( permutated_fractions / num_permutations ) )
    pvals = pval_lower / num_permutations
    tmp_pvals_higher = pval_higher / num_permutations
    pvals[fold_change > 0] = tmp_pvals_higher[fold_change > 0]

    return pvals, fold_change, true_neighbor_fractions

def neighborhood_analysis_single_sample(data, cluster_col='cluster', **kwargs):
    """Check out that apply 8-)
    """
    t = time.time()
    permutation_result = data.groupby('cluster').apply(lambda group: permutation_neighborhood(group,data))
    print('Total time',time.time() - t,'seconds')
    # split groupby result into different dataframes
    pval_report = pd.concat( [s[0] for s in permutation_result], axis=1, keys=permutation_result.index, sort=False)
    fc_report = pd.concat( [s[1] for s in permutation_result], axis=1, keys=permutation_result.index, sort=False)
    fraction_report = pd.concat( [s[2] for s in permutation_result], axis=1, keys=permutation_result.index, sort=False)

    return pval_report, fc_report, fraction_report

def processFile(fname, path_out):
    t = time.time()
    df = pd.read_csv(fname)
    df.index += 1
    df.index =  df.index.astype('int64')
    pvals, fc, fractions = neighborhood_analysis_single_sample(df)
    key = fname.split('/')[1]
    pvals.to_csv(path_out + "/pvals_" + key)
    fc.to_csv(path_out + "/foldchanges_" + key)
    fractions.to_csv(path_out + "/fractions_" + key)

def main():
    t0 = time.time()
    path_in = 'neighbor_data'
    path_out = 'neighbor_results'
    files = [filename for filename in glob.iglob(path_in+'/*.csv', recursive=True)]
    n_files = len(files)
    Pros = []
    for fname in files:
        p = Process(target=processFile, args=(fname, path_out))
        Pros.append(p)
        freeze_support()
        p.start()

    for t in Pros:
        t.join()

    t1 = time.time() - t0
    print("Finished all ", n_files, " files in ", t1, " seconds.\n")

if __name__ == "__main__":
    main()
