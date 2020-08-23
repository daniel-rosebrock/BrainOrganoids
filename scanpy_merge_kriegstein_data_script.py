import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Load Kriegstein datasets
adata = sc.read_text("/project/elkabetz_lab_data/public_data/kriegstein_2020/oragnoid_scrna/exprMatrix.tsv.gz")
adata = adata.transpose()
meta = pd.read_csv("/project/elkabetz_lab_data/public_data/kriegstein_2020/oragnoid_scrna/meta.tsv", sep="\t")
adata.obs = meta

genes = list(adata.var_names)
sc.pp.highly_variable_genes(adata,n_bins=50)
min_mean = 0.0125
top_2000 = sorted(adata.var['dispersions_norm'][adata.var['means'] >= min_mean],reverse=True)[3000]
adata.var['highly_variable'] = (adata.var['dispersions_norm'] >= top_2000) & (adata.var['means'] >= min_mean)
adata.var[adata.var['highly_variable'] == True].sort_values(['dispersions_norm'],ascending=False)[0:50]
adata = adata[:, adata.var[adata.var['highly_variable'] == True].index]

# Load and Normalize Mutukula et al. dataset
adata_org = sc.read('/project/elkabetz_lab_data/sequencing/scrna/d50_organoids/cortical_cells_count_matrix.csv')
sc.pp.normalize_total(adata_org)
sc.pp.log1p(adata_org)

genes_org = list(adata_org.var_names)
genes_overlap = list(set(adata.var_names).intersection(adata_org.var_names))
adata_org = adata_org[:,genes_overlap]
sc.pp.scale(adata_org)

# Normalize Kriegstein organoid data
adata = adata[:,genes_overlap]
sc.pp.regress_out(adata,['Batch'])
sc.pp.scale(adata, max_value=10)

# Scale Mutukula et al. dataset accordingly
sc.pp.scale(adata_org, max_value=10)

# add meta data
meta_dat = pd.read_csv('/project/elkabetz_lab_data/sequencing/scrna/d50_organoids/cortical_cells_d50_cell_ident.txt',sep='\t')
cols_name = [x if x != 'Unnamed: 0' else 'Cell' for x in meta_dat.columns]
adata_org.obs = meta_dat

## merge datasets
adata_merge = adata.concatenate(adata_org,batch_key='sample')

## combat
sc.pp.combat(adata_merge,key='sample')


sc.tl.pca(adata_merge)
sc.pp.neighbors(adata_merge)
sc.tl.umap(adata_merge)

tx_dict = {'nt':'Inhibitor-free','sbnx':'Triple-i','sbn':'Dual SMAD-i'}

adata_merge.obs['Protocol'] = [x if x in ['Sasai','Xiang','Pasca'] else tx_dict[adata_merge.obs['Cell'][j].split(".")[0]] for j,x in enumerate(adata_merge.obs['Protocol'])]
adata_merge.obs['Age'] = [x if x in [6.0,10.0,14.0,18.0,24.0] else tx_dict[adata_merge.obs['Cell'][j].split(".")[0]] for j,x in enumerate(adata_merge.obs['Age'])]
adata_merge.obs.astype({'Age':'category','Cluster':'category'}).dtypes



adata_merge.obs['Index'] = range(len(adata_merge.obs))

color_d = {'Sasai':'lightcoral','Xiang':'paleturquoise','Pasca':'lightgreen','Inhibitor-free':'green','Dual SMAD-i':'purple','Triple-i':'yellow'}
size_d = {'Sasai':0.5,'Xiang':0.5,'Pasca':0.5,'Inhibitor-free':1,'Dual SMAD-i':1,'Triple-i':1}
alpha_d = {'Sasai':0.1,'Xiang':0.1,'Pasca':0.1,'Inhibitor-free':0.5,'Dual SMAD-i':0.5,'Triple-i':0.5}


## kriegstein fig
fig = plt.figure(figsize=(10,10))

for tx in ['Sasai','Xiang','Pasca']:
    plt.scatter(x=[x[0] for x in adata_merge.obsm['X_umap'][adata_merge.obs[adata_merge.obs['Protocol']==tx]['Index']]],y=[x[1] for x in adata_merge.obsm['X_umap'][adata_merge.obs[adata_merge.obs['Protocol']==tx]['Index']]],color=color_d[tx],s=size_d[tx],alpha=alpha_d[tx])

plt.xlim(plt.xlim())
plt.ylim(plt.ylim())

for tx in ['Sasai','Xiang','Pasca']:
    plt.plot(-10000,-10000,'o',label=tx,color=color_d[tx])

plt.legend()
plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/umap.kriegstein.pdf')
plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/umap.kriegstein.png',dpi=500)


## overall fig
fig = plt.figure(figsize=(10,10))

for tx in ['Sasai','Xiang','Pasca','Inhibitor-free','Dual SMAD-i','Triple-i']:
    plt.scatter(x=[x[0] for x in adata_merge.obsm['X_umap'][adata_merge.obs[adata_merge.obs['Protocol']==tx]['Index']]],y=[x[1] for x in adata_merge.obsm['X_umap'][adata_merge.obs[adata_merge.obs['Protocol']==tx]['Index']]],color=color_d[tx],s=size_d[tx],alpha=alpha_d[tx])

plt.xlim(plt.xlim())
plt.ylim(plt.ylim())

for tx in ['Sasai','Xiang','Pasca']:
    plt.plot(-10000,-10000,'o',label=tx,color=color_d[tx])

plt.legend()
plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/merged_umap.pdf')
plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/merged_umap.png',dpi=500)

## indiv clust figs
clust_dict = {0:4,1:1,2:9,3:3,4:7,5:5,6:10,7:2,8:8,9:6,10:12,11:11}
for clust in clust_dict.keys():
    fig = plt.figure(figsize=(10,10))
    print('here')
    for tx in ['Sasai','Xiang','Pasca']:
        plt.scatter(x=[x[0] for x in adata_merge.obsm['X_umap'][adata_merge.obs[adata_merge.obs['Protocol']==tx]['Index']]],y=[x[1] for x in adata_merge.obsm['X_umap'][adata_merge.obs[adata_merge.obs['Protocol']==tx]['Index']]],color=color_d[tx],s=size_d[tx],alpha=alpha_d[tx])
    for tx in ['Inhibitor-free','Dual SMAD-i','Triple-i']:
        plt.scatter(x=[x[0] for x in adata_merge.obsm['X_umap'][adata_merge.obs[(adata_merge.obs['Protocol']==tx) & (adata_merge.obs['res.0.6']==clust)]['Index']]],y=[x[1] for x in adata_merge.obsm['X_umap'][adata_merge.obs[(adata_merge.obs['Protocol']==tx) & (adata_merge.obs['res.0.6']==clust)]['Index']]],color=color_d[tx],s=size_d[tx],alpha=alpha_d[tx])
    plt.xlim(plt.xlim())
    plt.ylim(plt.ylim())
    for tx in ['Sasai','Xiang','Pasca','Inhibitor-free','Dual SMAD-i','Triple-i']:
        plt.plot(-10000,-10000,'o',label=tx,color=color_d[tx])
    plt.legend()
    plt.title('Cluster '+str(clust_dict[clust]))
    plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/merged_umap.C'+str(clust_dict[clust])+'.pdf')
    plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/merged_umap.C'+str(clust_dict[clust])+'.png',dpi=500)

## age
age_color_dict = {3.:'paleturquoise', 5.:'dodgerblue', 8.:'navy', 10.:'maroon', 15.:'tomato', 24.:'lightsalmon'}
fig = plt.figure(figsize=(10,10))

for age in [3.,5.,8.,10.,15.,24.]:
    plt.scatter(x=[x[0] for x in adata_merge.obsm['X_umap'][adata_merge.obs[adata_merge.obs['Age']==age]['Index']]],y=[x[1] for x in adata_merge.obsm['X_umap'][adata_merge.obs[adata_merge.obs['Age']==age]['Index']]],color=age_color_dict[age],s=0.5,alpha=0.1)

plt.xlim(plt.xlim())
plt.ylim(plt.ylim())

for age in [3.,5.,8.,10.,15.,24.]:
    plt.plot(-10000,-10000,'o',label='W'+str(int(age)),color=age_color_dict[age])

plt.legend()
plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/umap.kriegstein.age.pdf')
plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/umap.kriegstein.age.png',dpi=500)

gene = 'FABP7'
fig = sc.pl.umap(adata_merge,color=[gene],color_map='Purples',return_fig=True)
plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/umap_merged.'+gene+'.png',dpi=500)
plt.savefig('/project/elkabetz_lab/Daniel/Manuscript/Figures/kriegstein_merge/umap_merged.'+gene+'.pdf')
