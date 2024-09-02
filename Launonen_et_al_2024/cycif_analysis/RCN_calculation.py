# -*- coding: utf-8 -*-
"""
Created on Thu May  9 19:46:52 2024

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

#sciset RCN calculation

cdata = sm.pp.mcmicro_to_scimap("D:/Sciset/data_030123.csv", unique_CellId=True, split="X.Position", log=False)

cdata = sm.tl.spatial_count(cdata, x_coordinate='X.Position', y_coordinate='Y.Position', phenotype='GlobalCellType2', method='radius', radius=100, knn=10, imageid='imageid', subset=None, label='spatial_count')

cdata = sm.tl.spatial_cluster(cdata, df_name='spatial_count', method='kmeans', k=18, n_pcs=None, resolution=1, phenograph_clustering_metric='euclidean', nearest_neighbors=30, random_state=0, label=None, output_dir=None)

sm.pl.stacked_barplot(cdata, x_axis='spatial_kmeans', y_axis='GlobalCellType2', subset_xaxis=None, subset_yaxis=None, order_xaxis=None, order_yaxis=None, method='percent', plot_tool='matplotlib', matplotlib_cmap=None, matplotlib_bbox_to_anchor=(1, 1.02), matplotlib_legend_loc=2, return_data=False)

#view neighborhoods on top of images
cdata.obs['spatial_kmeans'].astype('category')

image_path="E:/sciset/ometiffs/Sample_14.ome.tif"
sm.pl.image_viewer(image_path, cdata, subset="s14", imageid='imageid', channel_names=['DNA1', 'Rabbit', 'Goat', 'Mouse', 'DNA2', 'TAZ', 'CD207', 'SNAT1', 'DNA3', 'CD163', 'CD57', 'CD20', 'DNA4', 'Annexin', 'pSTAT1', 'KRAS', 'DNA5', 'CD4', 'pERK', 'CD8a', 'DNA6', 'CD45RO', 'FOXP3', 'CD3d', 'DNA7', 'TIM3', 'oldCD68', 'Desmin', 'DNA8', 'CD15','oldCD11b', 'FOXOA3', 'DNA9', 'HE4', 'CD11c', 'yH2AX', 'DNA10', 'Ki67', 'Vimentin', 'MHCII', 'DNA11', 'LaminB1', 'CK7', 'MHCI', 'DNA12', 'Ecadherin', 'SMA', 'CD31', 'DNA13', 'IBA1', 'CD68', 'CD11b'], point_color='white', point_size=30, flip_y=False, overlay='spatial_kmeans', x_coordinate="X.Position", y_coordinate="Y.Position")


