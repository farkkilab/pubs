# -*- coding: utf-8 -*-
"""
Created on Thu May  9 18:39:08 2024

@author: ingamari
"""
import anndata as ad
import pandas as pd
import sys
import os
import numpy as np

import scanpy as sc
import seaborn as sns
sns.set(color_codes=True)

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.cluster import hierarchy as sch
from matplotlib import pyplot as plt
from copy import deepcopy
%matplotlib inline

import scimap as sm

adata = sm.pp.mcmicro_to_scimap("E:\sciset/data_20231125.csv", CellId = "ID", split="X.Position", log=False)


s01 = adata[adata.obs["imageid"] == "s01"]
s02 = adata[adata.obs["imageid"] == "s02"]
s03 = adata[adata.obs["imageid"] == "s03"]
s04 = adata[adata.obs["imageid"] == "s04"]
s05 = adata[adata.obs["imageid"] == "s05"]
s06 = adata[adata.obs["imageid"] == "s06"]
s07 = adata[adata.obs["imageid"] == "s07"]
s08 = adata[adata.obs["imageid"] == "s08"]
s09 = adata[adata.obs["imageid"] == "s09"]
s11 = adata[adata.obs["imageid"] == "s11"]
s12 = adata[adata.obs["imageid"] == "s12"]
s13 = adata[adata.obs["imageid"] == "s13"]
s14 = adata[adata.obs["imageid"] == "s14"]
s15 = adata[adata.obs["imageid"] == "s15"]
s16 = adata[adata.obs["imageid"] == "s16"]
s17 = adata[adata.obs["imageid"] == "s17"]
s18 = adata[adata.obs["imageid"] == "s18"]
s19 = adata[adata.obs["imageid"] == "s19"]
s20 = adata[adata.obs["imageid"] == "s20"]
s21 = adata[adata.obs["imageid"] == "s21"]
s22 = adata[adata.obs["imageid"] == "s22"]
s23 = adata[adata.obs["imageid"] == "s23"]

sq.gr.spatial_neighbors(s01, spatial_key='spatial', library_key='imageid')
sq.gr.centrality_scores(s01, 'GlobalCellType2')
sq.gr.centrality_scores(s01, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s01, 'neighbordood_cluster2')

sq.gr.spatial_neighbors(s02, spatial_key='spatial', coord_type="generic")
sq.gr.centrality_scores(s02, 'GlobalCellType2')
sq.gr.centrality_scores(s02, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s02, 'neighbordood_cluster2')
sq.gr.centrality_scores(s02, 'classify')


sq.gr.centrality_scores(s03, 'GlobalCellType2')
sq.gr.centrality_scores(s03, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s03, 'neighbordood_cluster2')

sq.gr.centrality_scores(s04, 'GlobalCellType2')
sq.gr.centrality_scores(s04, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s04 'neighbordood_cluster2')

sq.gr.centrality_scores(s05, 'GlobalCellType2')
sq.gr.centrality_scores(s05, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s05, 'neighbordood_cluster2')

sq.gr.centrality_scores(s06, 'GlobalCellType2')
sq.gr.centrality_scores(s06, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s06, 'neighbordood_cluster2')

sq.gr.centrality_scores(s07, 'GlobalCellType2')
sq.gr.centrality_scores(s07, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s07, 'neighbordood_cluster2')

sq.gr.centrality_scores(s08, 'GlobalCellType2')
sq.gr.centrality_scores(s08, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s08, 'neighbordood_cluster2')

sq.gr.centrality_scores(s09, 'GlobalCellType2')
sq.gr.centrality_scores(s09, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s09, 'neighbordood_cluster2')

sq.gr.spatial_neighbors(s11, spatial_key='spatial', coord_type="generic")
sq.gr.centrality_scores(s11, 'GlobalCellType2')
sq.gr.centrality_scores(s11, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s11, 'neighbordood_cluster2')
sq.gr.centrality_scores(s11, 'classify')

sq.gr.centrality_scores(s12, 'GlobalCellType2')
sq.gr.centrality_scores(s12, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s12, 'neighbordood_cluster2')

sq.gr.centrality_scores(s13, 'GlobalCellType2')
sq.gr.centrality_scores(s13, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s13, 'neighbordood_cluster2')

sq.gr.spatial_neighbors(s14, spatial_key='spatial', coord_type="generic")
sq.gr.centrality_scores(s14, 'GlobalCellType2')
sq.gr.centrality_scores(s14, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s14, 'neighbordood_cluster2')
sq.gr.centrality_scores(s14, 'classify')


sq.gr.centrality_scores(s15, 'GlobalCellType2')
sq.gr.centrality_scores(s15, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s15, 'neighbordood_cluster2')

sq.gr.centrality_scores(s16, 'GlobalCellType2')
sq.gr.centrality_scores(s16, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s16, 'neighbordood_cluster2')

sq.gr.centrality_scores(s17, 'GlobalCellType2')
sq.gr.centrality_scores(s17, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s17, 'neighbordood_cluster2')

sq.gr.centrality_scores(s18, 'GlobalCellType2')
sq.gr.centrality_scores(s18, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s18, 'neighbordood_cluster2')

sq.gr.centrality_scores(s19, 'GlobalCellType2')
sq.gr.centrality_scores(s19, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s19, 'neighbordood_cluster2')

sq.gr.spatial_neighbors(s20, spatial_key='spatial', coord_type="generic")
sq.gr.centrality_scores(s20, 'GlobalCellType2')
sq.gr.centrality_scores(s20, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s20, 'neighbordood_cluster2')
sq.gr.centrality_scores(s20, 'classify')


sq.gr.centrality_scores(s21, 'GlobalCellType2')
sq.gr.centrality_scores(s21, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s21, 'neighbordood_cluster2')

sq.gr.centrality_scores(s22, 'GlobalCellType2')
sq.gr.centrality_scores(s22, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s22, 'neighbordood_cluster2')

sq.gr.centrality_scores(s23, 'GlobalCellType2')
sq.gr.centrality_scores(s23, 'GlobalCellType_cd8')
sq.gr.centrality_scores(s23, 'neighbordood_cluster2')


s01.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s01_globalcelltype2.csv")
s01.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s01_globalcelltypecd8_v2.csv")
s01.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s01_neighbordood_cluster2.csv")

s02.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s02_globalcelltype2.csv")
s02.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s02_globalcelltypecd8.csv")
s02.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s02_neighbordood_cluster2.csv")
s02.uns['classify_centrality_scores'].to_csv("D:/Sciset/closeness_scores/pSTAT1/s02_classify.csv")

s03.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s03_globalcelltype2.csv")
s03.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s03_globalcelltypecd8.csv")
s03.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s03_neighbordood_cluster2.csv")

s04.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s04_globalcelltype2.csv")
s04.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s04_globalcelltypecd8.csv")
s04.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s04_neighbordood_cluster2.csv")

s05.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s05_globalcelltype2.csv")
s05.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s05_globalcelltype_cd8.csv")
s05.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s05_neighbordood_cluster2.csv")

s06.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s06_globalcelltype2.csv")
s06.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s06_globalcelltypecd8.csv")
s06.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s06_neighbordood_cluster2.csv")

s07.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s07_globalcelltype2.csv")
s07.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s07_globalcelltypecd8.csv")
s07.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s07_neighbordood_cluster2.csv")

s08.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s08_globalcelltype2.csv")
s08.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s08_globalcelltypecd8.csv")
s08.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s08_neighbordood_cluster2.csv")

s09.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s09_globalcelltype2.csv")
s09.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s09_globalcelltypecd8.csv")
s09.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s09_neighbordood_cluster2.csv")

s11.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s11_globalcelltype2.csv")
s11.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s11_globalcelltypecd8.csv")
s11.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s11_neighbordood_cluster2.csv")
s11.uns['classify_centrality_scores'].to_csv("D:/Sciset/closeness_scores/pSTAT1/s11_classify.csv")

s12.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s12_globalcelltype2.csv")
s12.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s12_globalcelltypecd8.csv")
s12.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s12_neighbordood_cluster2.csv")

s13.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s13_globalcelltype2.csv")
s13.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s13_globalcelltypecd8.csv")
s13.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s13_neighbordood_cluster2.csv")

s14.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s14_globalcelltype2.csv")
s14.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s14_globalcelltypecd8.csv")
s14.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s14_neighbordood_cluster2.csv")
s14.uns['classify_centrality_scores'].to_csv("D:/Sciset/closeness_scores/pSTAT1/s14_classify.csv")

s15.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s15_globalcelltype2.csv")
s15.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s15_globalcelltypecd8.csv")
s15.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s15_neighbordood_cluster2.csv")

s16.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s16_globalcelltype2.csv")
s16.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s16_globalcelltypecd8.csv")
s16.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s16_neighbordood_cluster2.csv")

s17.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s17_globalcelltype2.csv")
s17.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s17_globalcelltypecd8.csv")
s17.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s17_neighbordood_cluster2.csv")

s18.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s18_globalcelltype2.csv")
s18.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s18_globalcelltypecd8.csv")
s18.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s18_neighbordood_cluster2.csv")

s19.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s19_globalcelltype2.csv")
s19.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s19_globalcelltypecd8.csv")
s19.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s19_neighbordood_cluster2.csv")

s20.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s20_globalcelltype2.csv")
s20.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s20_globalcelltypecd8.csv")
s20.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s20_neighbordood_cluster2.csv")
s20.uns['classify_centrality_scores'].to_csv("D:/Sciset/closeness_scores/pSTAT1/s20_classify.csv")


s21.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s21_globalcelltype2.csv")
s21.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s21_globalcelltypecd8.csv")
s21.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s21_neighbordood_cluster2.csv")

s22.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s22_globalcelltype2.csv")
s22.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s22_globalcelltypecd8.csv")
s22.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s22_neighbordood_cluster2.csv")

s23.uns['GlobalCellType2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s23_globalcelltype2.csv")
s23.uns['GlobalCellType_cd8_centrality_scores'].to_csv("E:/Sciset/closeness_scores/GlobalCellType2/s23_globalcelltypecd8.csv")
s23.uns['neighbordood_cluster2_centrality_scores'].to_csv("E:/Sciset/closeness_scores/neighbordood_cluster2/s23_neighbordood_cluster2.csv")


