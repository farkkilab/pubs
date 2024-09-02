# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 17:27:58 2024

@author: ingamari
"""
#umaps for t-CycIF cells

import matplotlib.pyplot as plt
import scimap as sm
import anndata as ad
import pandas as pd
import sys
import os
import numpy as np

import scanpy as sc
import seaborn as sns
sns.set(color_codes=True)

#the data is already log transformed
ddata = sm.pp.mcmicro_to_scimap("E:/Sciset/data_14062022.csv", unique_CellId=True, split="X.Position", log=False, CellId='ID')

#subset and select markers
bdata = sc.pp.subsample(ddata, n_obs=500000, random_state=0, copy=True)
bdata = bdata[:,bdata.var_names.isin(['Ecadherin', 'Vimentin','SMA', 'MHCII', 'CD3d', 'CD4', 'CD8a', 'IBA1', 'CD11c', 'CD163', 'CK7', 'FOXP3', 'CD31'])]


#compute
sc.pp.neighbors(bdata, n_neighbors=30, n_pcs=10)
sc.tl.umap(bdata)
sc.pl.umap(bdata, color=['Ecadherin', 'Vimentin','SMA', 'MHCII', 'CD3d', 'CD4', 'CD8a', 'IBA1', 'CD11c', 'CD163', 'CK7', 'FOXP3', 'CD31'], size=2, use_raw=False)

#plot one at a time
sc.pl.umap(bdata, color=['Ecadherin'], size=2, use_raw=False)
sc.pl.umap(bdata, color=[ 'Vimentin'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['SMA'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['MHCII'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['CD3d'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['CD4'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['CD8a'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['IBA1'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['CD11c'], size=2, use_raw=False)
sc.pl.umap(bdata, color=[ 'CD163'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['CK7'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['FOXP3'], size=2, use_raw=False)
sc.pl.umap(bdata, color=['CD31'], size=2, use_raw=False)


#sc.pl.umap(bdata, color=['Global_celltype_tumor_stroma_merged2'], size=1, palette=palette)

#rename cells to merge categories
rename = {'Tumor': ['EMT', 'Proliferating.EMT', 'Epithelial', 'Proliferating.epithelial'], 'Myeloid': ["IBA1.CD11c.Macrophages", "IBA1.CD163.Macrophages", "CD163.Macrophages", "CD11c.myeloid"], 'Stroma': ["Desmin.positive.cell", "Myofibroblast", "Fibroblast", "SMA.Desmin.positive.cell"], 'Endothelial.cell' : ["SMA.CD31.positive.cell"]}
bdata = sm.hl.rename(bdata, rename, from_column='GlobalCellType2', to_column='GlobalCellType4')
sns.set(rc={'figure.figsize': (11, 11), 'axes.facecolor':'white', 'figure.facecolor':'white'})

plt.rcParams.update({'axes.labelsize' : 'large'})

palette = {"CD8.T.cells": "#ff5a5f","CD4.T.cells": "#57cc99", "FOXP3.CD4.Tregs": "#bfd7ea", "CD163.Macrophages": "#d0c9ea",'IBA1.Macrophages':"#a786db", 'IBA1.CD11c.Macrophages':"#b59ce0", 'IBA1.CD163.Macrophages':"#c2b3e5",  "Bcell":"#0b3954", "CD11c.myeloid":"#996fd6", "Stroma":"#d5bdaf", "Tumor":"#747483", 'Endothelial.cell':"#9E8E6E"}

sc.pl.umap(bdata, color=['GlobalCellType4'], size=2, palette={"Tumor":"#f6114a","Stroma":"#0aa0bf","Endothelial.cell":"#fedb39","Myeloid":"#78b177","CD8.T.cells":"#9862a2","CD4.T.cells":"#f05006","FOXP3.CD4.Tregs":"#f36e98","Bcell":"#fca00c"})

#plot markers also with this size

sc.pl.umap(bdata, color=['Ecadherin', 'Vimentin','SMA', 'MHCII', 'CD3d', 'CD4', 'CD8a', 'IBA1', 'CD11c', 'CD163', 'CK7', 'FOXP3', 'CD31'], size=2, use_raw=False)
