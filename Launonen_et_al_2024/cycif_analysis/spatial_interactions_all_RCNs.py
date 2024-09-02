# -*- coding: utf-8 -*-
"""
Created on Sun May 12 19:51:01 2024

@author: ingamari
"""
import scimap as sm
import anndata as ad
import pandas as pd
import sys
import os
import numpy as np
import scanpy as sc
import seaborn as sns
sns.set(color_codes=True)

#spatial interactions per RCN

cdata_RCN = sm.pp.mcmicro_to_scimap("P:/h345/afarkkilab/Data/Sciset/Epithelial_EMT_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Epithelial_EMT_interactions.csv')

#this is not run yet
cdata_RCN = sm.pp.mcmicro_to_scimap("P:/h345/afarkkilab/Data/Sciset/EMT_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/EMT_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/Proliferating.EMT_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Proliferating.EMT_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/Epithelial_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Epithelial_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/Proliferating.epithelial_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Proliferating.epithelial_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/epithelial_and_proliferating_epithelial_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/epithelial_and_proliferating_epithelial_interactions.csv')

#
cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/Fibroblast_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
cdata_RCN = sm.hl.dropFeatures(cdata_RCN, drop_markers=None, drop_cells=None, drop_meta_columns=None, drop_groups='s23', groups_column='imageid', subset_raw=True, verbose=True)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Fibroblast_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/Myofibroblast_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Myofibroblast_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/SMA.CD31.positive_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/SMA.CD31.positive_interactions.csv')

#here
cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/Immune_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Immune_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/CD8_CD4.T.cells_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/CD8_CD4.T.cells_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/Macrophages_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Macrophages_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/Desmin.positive_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/Desmin.positive_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/CD11c.myeloid_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/CD11c.myeloid_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/stroma_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/stroma_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/IBA1.CD163.Macrophages_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/IBA1.CD163.Macrophages_interactions.csv')


cdata_RCN = sm.pp.mcmicro_to_scimap("E:/sciset/SMA.Desmin.myofibroblast_neighborhood.csv", unique_CellId=True, split="X.Position", log=False)
sm.tl.spatial_interaction(cdata_RCN, x_coordinate="X.Position", y_coordinate="Y.Position", phenotype="GlobalCellType2", method="radius", radius=45)
interactions = cdata_RCN.uns['spatial_interaction']
interactions.to_csv('D:/Sciset/SMA.Desmin.myofibroblast_interactions.csv')
















