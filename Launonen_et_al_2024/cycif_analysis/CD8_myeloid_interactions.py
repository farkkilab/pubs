#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 11:51:36 2024

@author: oncosys
"""

#myelonets

#calculate which CD8+T-cells interact with different myeloid cells

paths = '/media/oncosys/T7/sciset/data_with_myelonets.csv'
bdata = sm.pp.mcmicro_to_scimap(paths, unique_CellId=True, split="X.Position", log=False, CellId = 'ID')

sm.tl.spatial_pscore(bdata, proximity=['IBA1.CD163.Macrophages','CD8.T.cells'], score_by='imageid', x_coordinate='X.Position', y_coordinate='Y.Position', z_coordinate=None, phenotype='GlobalCellType2', method='radius', radius=45, knn=3, imageid='imageid', subset=None, label='spatial_pscore_IBA1.CD163.Macrophages')
sm.tl.spatial_pscore(bdata, proximity=['IBA1.CD11c.Macrophages','CD8.T.cells'], score_by='imageid', x_coordinate='X.Position', y_coordinate='Y.Position', z_coordinate=None, phenotype='GlobalCellType2', method='radius', radius=45, knn=3, imageid='imageid', subset=None, label='spatial_pscore_IBA1.CD11c.Macrophages')
sm.tl.spatial_pscore(bdata, proximity=['CD163.Macrophages','CD8.T.cells'], score_by='imageid', x_coordinate='X.Position', y_coordinate='Y.Position', z_coordinate=None, phenotype='GlobalCellType2', method='radius', radius=45, knn=3, imageid='imageid', subset=None, label='spatial_pscore_CD163.Macrophages')
sm.tl.spatial_pscore(bdata, proximity=['CD11c.myeloid','CD8.T.cells'], score_by='imageid', x_coordinate='X.Position', y_coordinate='Y.Position', z_coordinate=None, phenotype='GlobalCellType2', method='radius', radius=45, knn=3, imageid='imageid', subset=None, label='spatial_pscore_CD11c.myeloid')


#save the interacting cells as csv


proximities = bdata.obs['spatial_pscore_IBA1.CD163.Macrophages']
proximities.to_csv('/media/oncosys/T7/sciset/IBA1.CD163_interactions.csv')


proximities = bdata.obs['spatial_pscore_IBA1.CD11c.Macrophages']
proximities.to_csv('/media/oncosys/T7/sciset/IBA1.CD11c_interactions.csv')


proximities = bdata.obs['spatial_pscore_CD163.Macrophages']
proximities.to_csv('/media/oncosys/T7/sciset/CD163_interactions.csv')


proximities = bdata.obs['spatial_pscore_CD11c.myeloid']
proximities.to_csv('/media/oncosys/T7/sciset/CD11c.myeloid_interactions.csv')








