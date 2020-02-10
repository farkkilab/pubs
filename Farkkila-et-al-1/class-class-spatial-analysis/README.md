# Permutation testing

 * R: One single column with "cluster labels" aka cell types
 * Steps:
    1. getNeighbors.m to create extra columns with neighbor IDs.
    For all cells or just for a given cell type.
    2. permutation_test.py to test of interaction significance and produce neighbor class fractions
