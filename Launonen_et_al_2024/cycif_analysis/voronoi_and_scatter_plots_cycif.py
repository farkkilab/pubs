# -*- coding: utf-8 -*-
"""
Created on Thu May  9 18:35:26 2024

@author: ingamari
"""
#voronoi images of neighborhoods and scatter plots of cell types
#for figure 3 and 4

import scimap as sm
import anndata as ad
import pandas as pd
import sys
import os
import numpy as np

import scanpy as sc
import seaborn as sns
sns.set(color_codes=True)

cdata = ad.read("E:/sciset/data_20231125.h5ad")

#for viewing images in Napari with cell type calls
image_path="E:/sciset/ometiffs/Sample_12.ome.tif"
sm.pl.image_viewer(image_path, cdata, subset="s12", imageid='imageid', channel_names=['DNA1', 'Rabbit', 'Goat', 'Mouse', 'DNA2', 'TAZ', 'CD207', 'SNAT1', 'DNA3', 'CD163', 'CD57', 'CD20', 'DNA4', 'Annexin', 'pSTAT1', 'KRAS', 'DNA5', 'CD4', 'pERK', 'CD8a', 'DNA6', 'CD45RO', 'FOXP3', 'CD3d', 'DNA7', 'TIM3', 'oldCD68', 'Desmin', 'DNA8', 'CD15','oldCD11b', 'FOXOA3', 'DNA9', 'HE4', 'CD11c', 'yH2AX', 'DNA10', 'Ki67', 'Vimentin', 'MHCII', 'DNA11', 'LaminB1', 'CK7', 'MHCI', 'DNA12', 'Ecadherin', 'SMA', 'CD31', 'DNA13', 'IBA1', 'CD68', 'CD11b'], point_color='white', point_size=30, flip_y=False, overlay='GlobalCellType2', x_coordinate="X.Position", y_coordinate="Y.Position")


#Voronoi plots of whole slides images
colors = {"Proliferating.epithelial": "#100EAF","epithelial_and_proliferating_epithelial":"#0A087F", "Epithelial":"#29279F","Epithelial_EMT" :"#3F3596","EMT":"#261B83" ,"Proliferating.EMT" :"#110575", "tumor-stroma-interface":"#A020F0", "Desmin.positive":"#383f51","SMA.Desmin.myofibroblast":"#dddbf1","stroma":"#929487", "SMA.CD31.positive" :"#D3D4D9", "Myofibroblast":"#d1beb0",  "Fibroblast":"#ab9f9d", "Macrophages":"#247ba0","IBA1.CD163.Macrophages": "#70c1b3","CD11c.myeloid":"#b2dbbf", "CD8_CD4_T.cells":"#ff1654","Immune":"#f3ffbd"}

sm.pl.voronoi(cdata, imageid = "imageid", subset="s04",  colors = colors, figsize=(30, 30), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)
sm.pl.voronoi(cdata, imageid = "imageid", subset="s14",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)
sm.pl.voronoi(cdata, imageid = "imageid", subset="s07",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)
sm.pl.voronoi(cdata, imageid = "imageid", subset="s03",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)
sm.pl.voronoi(cdata, imageid = "imageid", subset="s13",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)


colors = {'Epithelial' : "#627D87", 'EMT' : "#467181", 'Proliferating.epithelial' :"#747483", 'Proliferating.EMT' : "#235A6F", 'CD8.T.cells' : "#ff5a5f", 'CD4.T.cells':"#57cc99", 'Bcell':"#0b3954", 'FOXP3.CD4.Tregs':"#bfd7ea", 'CD11c.myeloid':"#996fd6", 'IBA1.Macrophages':"#a786db", 'IBA1.CD11c.Macrophages':"#b59ce0", 'IBA1.CD163.Macrophages':"#c2b3e5", 'CD163.Macrophages':"#d0c9ea",'SMA.Desmin.positive.cell':"#edede9", 'Desmin.positive.cell':"#d6ccc2", 'SMA.CD31.positive.cell':"#f5ebe0",'Fibroblast':"#e3d5ca", 'Myofibroblast':"#d5bdaf", 'Endothelial.cell':"#9E8E6E"}

sm.pl.voronoi(cdata, imageid = "imageid", subset="s07",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='GlobalCellType2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)

#voronoi plots of spesific regions

#subsetting s14 - for fig 2
bdata = cdata[cdata.obs['imageid'] == "s14"]
bdata = bdata[bdata.obs['X.Position'] > 8293]
bdata = bdata[bdata.obs['X.Position'] < 10624]
bdata = bdata[bdata.obs['Y.Position'] > 13524]
bdata = bdata[bdata.obs['Y.Position'] < 16686]

sm.pl.voronoi(bdata, imageid = "imageid", subset="s14",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)

#subsetting s15 figure 2
bdata = cdata[cdata.obs['imageid'] == "s15"]
bdata = bdata[bdata.obs['X.Position'] < 24087]
bdata = bdata[bdata.obs['X.Position'] > 23436]
bdata = bdata[bdata.obs['Y.Position'] > 8979]
bdata = bdata[bdata.obs['Y.Position'] < 9870]

sm.pl.voronoi(bdata, imageid = "imageid", subset="s15",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)


bdata = cdata[cdata.obs['imageid'] == "s20"]
bdata = bdata[bdata.obs['X.Position'] < 9363]
bdata = bdata[bdata.obs['X.Position'] > 8806]
bdata = bdata[bdata.obs['Y.Position'] > 12591]
bdata = bdata[bdata.obs['Y.Position'] < 13221]

sm.pl.voronoi(bdata, imageid = "imageid", subset="s20",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)

bdata = cdata[cdata.obs['imageid'] == "s13"]
bdata = bdata[bdata.obs['X.Position'] < 21832]
bdata = bdata[bdata.obs['X.Position'] > 21292]
bdata = bdata[bdata.obs['Y.Position'] > 6876]
bdata = bdata[bdata.obs['Y.Position'] < 7481]

sm.pl.voronoi(bdata, imageid = "imageid", subset="s13",  colors = colors, figsize=(15, 15), x_coordinate='X.Position', y_coordinate='Y.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=False)

#for fig 2

bdata = cdata[cdata.obs['imageid'] == "s04"]
bdata = bdata[bdata.obs['Y.Position'] > 7009]
bdata = bdata[bdata.obs['Y.Position'] < 10644]
bdata = bdata[bdata.obs['X.Position'] > 5682]
bdata = bdata[bdata.obs['X.Position'] < 7997]



colors = {"Proliferating.epithelial": "#100EAF","epithelial_and_proliferating_epithelial":"#0A087F", "Epithelial":"#29279F","Epithelial_EMT" :"#3F3596","EMT":"#261B83" ,"Proliferating.EMT" :"#110575", "tumor-stroma-interface" :"#3d5a80", "Desmin.positive":"#383f51","SMA.Desmin.myofibroblast":"#dddbf1",
           "stroma":"#929487",  "SMA.CD31.positive" :"#D3D4D9", "Myofibroblast":"#d1beb0",  "Fibroblast":"#ab9f9d",
           "Macrophage":"#247ba0","IBA1.CD163.Macrophages": "#70c1b3","CD11c.myeloid":"#b2dbbf", 
           "CD8_CD4_T.cells":"#ff1654","Immune":"#f3ffbd" }

sm.pl.voronoi(bdata, imageid = "imageid",  colors = colors, x_coordinate='Y.Position', y_coordinate='X.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=True, dpi=600, outputDir="E:/Sciset/scimap/plots/scatter/", outputFileName='voronoi_plot_s04_subset.png')


#scatter plots

colors = {'Epithelial' : "#627D87", 'EMT' : "#467181", 'Proliferating.epithelial' :"#747483", 'Proliferating.EMT' : "#235A6F", 'CD8.T.cells' : "#ff5a5f", 'CD4.T.cells':"#57cc99", 'Bcell':"#0b3954", 'FOXP3.CD4.Tregs':"#bfd7ea", 'CD11c.myeloid':"#996fd6", 'IBA1.Macrophages':"#a786db", 'IBA1.CD11c.Macrophages':"#b59ce0", 'IBA1.CD163.Macrophages':"#c2b3e5", 'CD163.Macrophages':"#d0c9ea",'SMA.Desmin.positive.cell':"#edede9", 'Desmin.positive.cell':"#d6ccc2", 'SMA.CD31.positive.cell':"#f5ebe0",'Fibroblast':"#e3d5ca", 'Myofibroblast':"#d5bdaf", 'Endothelial.cell':"#9E8E6E"}

sm.pl.spatial_scatterPlot(bdata, colorBy = "GlobalCellType2", x_coordinate='Y.Position', y_coordinate='X.Position', imageid='imageid', layer=None, subset=None, s=10, ncols=None, alpha=1, dpi=600, fontsize=None, plotLegend=True, cmap='RdBu_r',figsize=(10, 5), vmin=None, vmax=None, customColors=colors, invert_yaxis=True)


#next: CD8CD4 rich area
bdata = cdata[cdata.obs['imageid'] == "s04"]
bdata = bdata[bdata.obs['Y.Position'] > 7419]
bdata = bdata[bdata.obs['Y.Position'] < 7797]
bdata = bdata[bdata.obs['X.Position'] > 7392]
bdata = bdata[bdata.obs['X.Position'] < 7749]



colors = {'Epithelial' : "#627D87", 'EMT' : "#467181", 'Proliferating.epithelial' :"#747483", 'Proliferating.EMT' : "#235A6F", 'CD8.T.cells' : "#ff5a5f", 'CD4.T.cells':"#57cc99", 'Bcell':"#0b3954", 'FOXP3.CD4.Tregs':"#bfd7ea", 'CD11c.myeloid':"#996fd6", 'IBA1.Macrophages':"#a786db", 'IBA1.CD11c.Macrophages':"#b59ce0", 'IBA1.CD163.Macrophages':"#c2b3e5", 'CD163.Macrophages':"#d0c9ea",'SMA.Desmin.positive.cell':"#edede9", 'Desmin.positive.cell':"#d6ccc2", 'SMA.CD31.positive.cell':"#f5ebe0",'Fibroblast':"#e3d5ca", 'Myofibroblast':"#d5bdaf", 'Endothelial.cell':"#9E8E6E"}

sm.pl.spatial_scatterPlot(bdata, colorBy = "GlobalCellType2", x_coordinate='Y.Position', y_coordinate='X.Position', imageid='imageid', layer=None, subset=None, s=80, ncols=None, alpha=1, dpi=600, fontsize=None, plotLegend=False, cmap='RdBu_r',figsize=(5, 5), vmin=None, vmax=None, customColors=colors, invert_yaxis=True)



#next: tsinterface rich area


bdata = cdata[cdata.obs['imageid'] == "s04"]
bdata = bdata[bdata.obs['Y.Position'] > 9547]
bdata = bdata[bdata.obs['Y.Position'] < 9877]
bdata = bdata[bdata.obs['X.Position'] > 7331]
bdata = bdata[bdata.obs['X.Position'] < 7561]



colors = {'Epithelial' : "#627D87", 'EMT' : "#467181", 'Proliferating.epithelial' :"#747483", 'Proliferating.EMT' : "#235A6F", 'CD8.T.cells' : "#ff5a5f", 'CD4.T.cells':"#57cc99", 'Bcell':"#0b3954", 'FOXP3.CD4.Tregs':"#bfd7ea", 'CD11c.myeloid':"#996fd6", 'IBA1.Macrophages':"#a786db", 'IBA1.CD11c.Macrophages':"#b59ce0", 'IBA1.CD163.Macrophages':"#c2b3e5", 'CD163.Macrophages':"#d0c9ea",'SMA.Desmin.positive.cell':"#edede9", 'Desmin.positive.cell':"#d6ccc2", 'SMA.CD31.positive.cell':"#f5ebe0",'Fibroblast':"#e3d5ca", 'Myofibroblast':"#d5bdaf", 'Endothelial.cell':"#9E8E6E"}

sm.pl.spatial_scatterPlot(bdata, colorBy = "GlobalCellType2", x_coordinate='Y.Position', y_coordinate='X.Position', imageid='imageid', layer=None, subset=None, s=80, ncols=None, alpha=1, dpi=600, fontsize=None, plotLegend=False, cmap='RdBu_r',figsize=(5, 5), vmin=None, vmax=None, customColors=colors, invert_yaxis=True)


#next: macrophage rich area


bdata = cdata[cdata.obs['imageid'] == "s04"]
bdata = bdata[bdata.obs['Y.Position'] > 7706]
bdata = bdata[bdata.obs['Y.Position'] < 8025]
bdata = bdata[bdata.obs['X.Position'] > 5682]
bdata = bdata[bdata.obs['X.Position'] < 5959]



colors = {'Epithelial' : "#627D87", 'EMT' : "#467181", 'Proliferating.epithelial' :"#747483", 'Proliferating.EMT' : "#235A6F", 'CD8.T.cells' : "#ff5a5f", 'CD4.T.cells':"#57cc99", 'Bcell':"#0b3954", 'FOXP3.CD4.Tregs':"#bfd7ea", 'CD11c.myeloid':"#996fd6", 'IBA1.Macrophages':"#a786db", 'IBA1.CD11c.Macrophages':"#b59ce0", 'IBA1.CD163.Macrophages':"#c2b3e5", 'CD163.Macrophages':"#d0c9ea",'SMA.Desmin.positive.cell':"#edede9", 'Desmin.positive.cell':"#d6ccc2", 'SMA.CD31.positive.cell':"#f5ebe0",'Fibroblast':"#e3d5ca", 'Myofibroblast':"#d5bdaf", 'Endothelial.cell':"#9E8E6E"}

sm.pl.spatial_scatterPlot(bdata, colorBy = "GlobalCellType2", x_coordinate='Y.Position', y_coordinate='X.Position', imageid='imageid', layer=None, subset=None, s=80, ncols=None, alpha=1, dpi=600, fontsize=None, plotLegend=False, cmap='RdBu_r',figsize=(5, 5), vmin=None, vmax=None, customColors=colors, invert_yaxis=True)


#s07 full


colors = {'Epithelial' : "#627D87", 'EMT' : "#467181", 'Proliferating.epithelial' :"#747483", 'Proliferating.EMT' : "#235A6F", 'CD8.T.cells' : "#ff5a5f", 'CD4.T.cells':"#57cc99", 'Bcell':"#0b3954", 'FOXP3.CD4.Tregs':"#bfd7ea", 'CD11c.myeloid':"#996fd6", 'IBA1.Macrophages':"#a786db", 'IBA1.CD11c.Macrophages':"#b59ce0", 'IBA1.CD163.Macrophages':"#c2b3e5", 'CD163.Macrophages':"#d0c9ea",'SMA.Desmin.positive.cell':"#edede9", 'Desmin.positive.cell':"#d6ccc2", 'SMA.CD31.positive.cell':"#f5ebe0",'Fibroblast':"#e3d5ca", 'Myofibroblast':"#d5bdaf", 'Endothelial.cell':"#9E8E6E"}

sm.pl.spatial_scatterPlot(cdata, colorBy = "GlobalCellType2",subset="s07", x_coordinate='Y.Position', y_coordinate='X.Position', imageid='imageid', layer=None, s=10, ncols=None, alpha=1, dpi=600, fontsize=None, plotLegend=False, cmap='RdBu_r',figsize=(6, 5), vmin=None, vmax=None, customColors=colors, invert_yaxis=True)



#next EMT in s13


bdata = cdata[cdata.obs['imageid'] == "s13"]
bdata = bdata[bdata.obs['Y.Position'] > 8517]
bdata = bdata[bdata.obs['Y.Position'] < 10197]
bdata = bdata[bdata.obs['X.Position'] > 15016]
bdata = bdata[bdata.obs['X.Position'] < 16277]


colors = {"Proliferating.epithelial": "#100EAF","epithelial_and_proliferating_epithelial":"#0A087F", "Epithelial":"#29279F","Epithelial_EMT" :"#3F3596","EMT":"#261B83" ,"Proliferating.EMT" :"#110575", "tumor-stroma-interface" :"#3d5a80", "Desmin.positive":"#383f51","SMA.Desmin.myofibroblast":"#dddbf1",
           "stroma":"#929487",  "SMA.CD31.positive" :"#D3D4D9", "Myofibroblast":"#d1beb0",  "Fibroblast":"#ab9f9d",
           "Macrophage":"#247ba0","IBA1.CD163.Macrophages": "#70c1b3","CD11c.myeloid":"#b2dbbf", 
           "CD8_CD4_T.cells":"#ff1654","Immune":"#f3ffbd" }

sm.pl.voronoi(bdata, imageid = "imageid",  colors = colors, x_coordinate='Y.Position', y_coordinate='X.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=True, dpi=600, outputDir="E:/Sciset/scimap/plots/scatter/", outputFileName='voronoi_plot_s04_subset.png')


#

image_path = "D:/sciset/ometiffs/Sample_18.ome.tif"
sm.pl.image_viewer(image_path, cdata, subset="s18", imageid='imageid', channel_names=['DNA1', 'Rabbit', 'Goat', 'Mouse', 'DNA2', 'TAZ', 'CD207', 'SNAT1', 'DNA3', 'CD163', 'CD57', 'CD20', 'DNA4', 'Annexin', 'pSTAT1', 'KRAS', 'DNA5', 'CD4', 'pERK', 'CD8a', 'DNA6', 'CD45RO', 'FOXP3', 'CD3d', 'DNA7', 'TIM3', 'oldCD68', 'Desmin', 'DNA8', 'CD15','oldCD11b', 'FOXOA3', 'DNA9', 'HE4', 'CD11c', 'yH2AX', 'DNA10', 'Ki67', 'Vimentin', 'MHCII', 'DNA11', 'LaminB1', 'CK7', 'MHCI', 'DNA12', 'Ecadherin', 'SMA', 'CD31', 'DNA13', 'IBA1', 'CD68', 'CD11b'], point_color='white', point_size=30, flip_y=False, overlay='neighbordood_cluster2', x_coordinate="X.Position", y_coordinate="Y.Position")


bdata = cdata[cdata.obs['imageid'] == "s18"]
bdata = bdata[bdata.obs['Y.Position'] > 9493]
bdata = bdata[bdata.obs['Y.Position'] < 15622]
bdata = bdata[bdata.obs['X.Position'] > 9322]
bdata = bdata[bdata.obs['X.Position'] < 13242]


#colors = {"Proliferating.epithelial": "#6b1650","epithelial_and_proliferating_epithelial":"#9c1057", "Epithelial":"#cc095d","Epithelial_EMT" :"#fd0363","EMT":"#ffc8dd" ,"Proliferating.EMT" :"#ffafcc", "tumor-stroma-interface" :"#3d5a80", "Desmin.positive":"#0077b6","SMA.Desmin.myofibroblast":"#90e0ef",  "stroma":"#caf0f8",  "SMA.CD31.positive" :"#118ab2", "Myofibroblast":"#3a86ff",  "Fibroblast":"#90e0ef", "Macrophage":"#c7f9cc","IBA1.CD163.Macrophages": "#a9fdac","CD11c.myeloid":"#44cf6c", "CD8_CD4_T.cells":"#32a287","Immune":"#38a3a5" }


colors = {"Proliferating.epithelial": "#100EAF","epithelial_and_proliferating_epithelial":"#0A087F", "Epithelial":"#29279F","Epithelial_EMT" :"#3F3596","EMT":"#261B83" ,"Proliferating.EMT" :"#110575", "tumor-stroma-interface" :"#3d5a80", "Desmin.positive":"#383f51","SMA.Desmin.myofibroblast":"#dddbf1",
           "stroma":"#929487",  "SMA.CD31.positive" :"#D3D4D9", "Myofibroblast":"#d1beb0",  "Fibroblast":"#ab9f9d",
           "Macrophage":"#247ba0","IBA1.CD163.Macrophages": "#70c1b3","CD11c.myeloid":"#b2dbbf", 
           "CD8_CD4_T.cells":"#ff1654","Immune":"#f3ffbd" }

sm.pl.voronoi(bdata, imageid = "imageid",  colors = colors, x_coordinate='Y.Position', y_coordinate='X.Position', color_by='neighbordood_cluster2',voronoi_edge_color = 'black',voronoi_line_width = 0.01, voronoi_alpha = 0.9, overlay_points=None,plot_legend=True, legend_size=6, flip_y=True, dpi=600, outputDir="E:/Sciset/scimap/plots/scatter/", outputFileName='voronoi_plot_s04_subset.png', figsize=(10,6))




#joku muu myös

#scatter plots of neighborhoods for reviewers

colors = {"Proliferating.epithelial": "#100EAF","epithelial_and_proliferating_epithelial":"#0A087F", "Epithelial":"#29279F","Epithelial_EMT" :"#3F3596","EMT":"#261B83" ,"Proliferating.EMT" :"#110575", "tumor-stroma-interface":"#A020F0", "Desmin.positive":"#383f51","SMA.Desmin.myofibroblast":"#dddbf1","stroma":"#929487", "SMA.CD31.positive" :"#D3D4D9", "Myofibroblast":"#d1beb0",  "Fibroblast":"#ab9f9d", "Macrophages":"#247ba0","IBA1.CD163.Macrophages": "#70c1b3","CD11c.myeloid":"#b2dbbf", "CD8_CD4_T.cells":"#ff1654","Immune":"#f3ffbd"}
sm.pl.spatial_scatterPlot(cdata, subset="s04",colorBy = "neighbordood_cluster2", x_coordinate='Y.Position', y_coordinate='X.Position', imageid='imageid', layer=None, s=80, ncols=None, alpha=1, dpi=600, fontsize=None, plotLegend=False, cmap='RdBu_r',figsize=(5, 5), vmin=None, vmax=None, customColors=colors, invert_yaxis=True)












