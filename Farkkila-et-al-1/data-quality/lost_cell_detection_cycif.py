import os
import pandas as pd
import cycifsuite.detect_lost_cells as dlc
import numpy as np
import time
import glob
from multiprocessing import Process, freeze_support

## First clone Jake Wang's repository from https://github.com/yunguan-wang/cycif_analysis_suite

def processFile(fname, path_out, manual_threshold=None):
    t = time.time()
    df = pd.read_csv(fname)
    qc_cols = [x for x in df.columns if 'DNA' in x]
    df_qc = df[qc_cols]
    n_cycles = len(qc_cols)
    df_qc = np.power(np.e,df_qc)
    fig = path_out + "/thresholding_" + fname.split('\\')[1] + '.png'
    if not manual_threshold:
        _,_,threshold = dlc.ROC_lostcells(df_qc,0,1,steps=50,n_cycles=n_cycles, filtering_method = 'cycle_diff',fld_stat_method='overall', automatic=True, figname=fig)
    else:
        threshold = manual_threshold
    lc,_ = dlc.get_lost_cells(df_qc, threshold,n_cycles=n_cycles, filtering_method='cycle_diff')
    df.loc[lc.index,'lost'] = True
    df.lost.fillna(False,inplace=True)
    df.to_csv(path_out + "/annotated_" + fname.split('\\')[1])
    timeInterval = time.time() - t
    print("[Done]", fname, "\tElapsed time: ", timeInterval, " seconds.\n")

def main():
    t0 = time.time()
    path_in = 'Z:/ANNIINA/Julia/2019/histoCAT/output'
    path_out = 'D:/julia/dataFromcellRingMask'
    os.chdir(path_in)
    files = [filename for filename in glob.iglob('**/*.csv', recursive=True)]
    n_files = len(files)
    print(n_files)
    manual_thresholds ={"HG0098":0.25, "C0046":0.285, "C0050":0.19, "C0085":0.18}
    Pros = []
    for fname in files:
        key = fname.split('\\')[0]
        if key in list(manual_thresholds.keys()):
            p = Process(target=processFile, args=(fname, path_out, manual_thresholds[key]))
            Pros.append(p)
            freeze_support()
            p.start()
    
    for t in Pros:
        t.join()
    
    t1 = time.time() - t0
    print("Finished all ", n_files, " files in ", t1, " seconds.\n")

if __name__ == "__main__":
    main()


""" for fname in files:
    print(fname)
    t = time.time()
    print("\tLoading data...")
    df = pd.read_csv(fname)
    qc_cols = [x for x in df.columns if 'DNA' in x]
    df_qc = df[qc_cols]
    n_cycles = len(qc_cols)
    df_qc = np.power(np.e,df_qc)
    print("\tCalculating threshold...")
    fig = path_out + "thresholding_" + fname + '.png'
    _,_,threshold = dlc.ROC_lostcells(df_qc,0,1,steps=50,n_cycles=n_cycles, filtering_method = 'cycle_diff',fld_stat_method='overall', automatic=True, figname=fig)
    print("\t",threshold)
    lc,_ = dlc.get_lost_cells(df_qc, threshold,n_cycles=n_cycles, filtering_method='cycle_diff')
    df.loc[lc.index,'lost'] = True
    df.lost.fillna(False,inplace=True)
    print("\tWriting annotation column...")
    df.to_csv(path_out + "/annot_" + fname)
    timeInterval = time.time() - t
    print("\tElapsed time: ",timeInterval, " seconds.\n") """