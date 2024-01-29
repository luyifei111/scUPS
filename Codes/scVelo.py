#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.chdir("/cluster/home/yflu/UPS/scvelo")
os.getcwd()


# In[2]:


import scanpy as sc
import anndata 
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd


# In[3]:


import scvelo as scv
import scanpy as sc
#import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad


# In[4]:


import matplotlib.pyplot as plt


# In[3]:


X = io.mmread(f"counts.mtx")


# In[4]:


adata = anndata.AnnData(
    X=X.transpose().tocsr()
)


# In[5]:


adata


# In[6]:


with open(f"gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()


# In[9]:


cell_meta = pd.read_csv(f"metadata.csv")


# In[10]:


adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']


# In[11]:


adata.var.index = gene_names


# In[12]:


pca = pd.read_csv(f"pca.csv")
pca.index = adata.obs.index


# In[13]:


adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T


# In[14]:


sc.pl.umap(adata, color='celltype', frameon=False, save=True)


# In[5]:


import scvelo as scv
import scanpy as sc
#import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad


# In[6]:


import matplotlib.pyplot as plt


# In[17]:


scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, 
                               fontsize=6, color_map='viridis',
                               frameon=False)


# In[19]:


sample_loom_1="/cluster/home/yflu/UPS/scvelo/loom/ZL-1.loom"
ldata1 = scv.read(sample_loom_1, cache=True)
ldata1
ldata1.var_names_make_unique()


# In[20]:


ldata1.obs.index[0:2]


# In[21]:


adata.obs.index[0:2]


# In[22]:


barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_1' for bc in barcodes]
ldata1.obs.index = barcodes
ldata1.obs.index[0:5]


# In[23]:


sample_loom_2="/cluster/home/yflu/UPS/scvelo/loom/ZL-2.loom"
ldata2 = scv.read( sample_loom_2, cache=True)
ldata2.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_2' for bc in barcodes]
ldata2.obs.index = barcodes
ldata2.obs.index[0:5]


# In[25]:


sample_loom_3="/cluster/home/yflu/UPS/scvelo/loom/ZL-9.loom"
ldata3 = scv.read( sample_loom_3, cache=True)
ldata3.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_3' for bc in barcodes]
ldata3.obs.index = barcodes
ldata3.obs.index[0:5]


# In[26]:


sample_loom_4="/cluster/home/yflu/UPS/scvelo/loom/ZL-10.loom"
ldata4 = scv.read( sample_loom_4, cache=True)
ldata4.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in ldata4.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_4' for bc in barcodes]
ldata4.obs.index = barcodes
ldata4.obs.index[0:5]


# In[27]:


sample_loom_5="/cluster/home/yflu/UPS/scvelo/loom/ZL-13.loom"
ldata5 = scv.read( sample_loom_5, cache=True)
ldata5.var_names_make_unique()
barcodes = [bc.split(':')[1] for bc in ldata5.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_5' for bc in barcodes]
ldata5.obs.index = barcodes
ldata5.obs.index[0:5]


# In[28]:


ldata = ldata1.concatenate([ldata2,ldata3,ldata4,ldata5])
ldata


# In[29]:


barcodes = [bc[0:len(bc)-2] for bc in ldata.obs.index.tolist()]
barcodes[0:5]


# In[30]:


ldata.obs.index = barcodes
ldata.obs.index[0:5]


# In[31]:


adata.obs.index[0:5]


# In[32]:


ldata.obs.index[0:5]


# In[33]:


ldata = ldata[np.isin(ldata.obs.index,adata.obs.index)]


# In[34]:


adata.obs.index[0:5]


# In[35]:


ldata.obs.index[0:5]


# In[36]:


adata = scv.utils.merge(adata, ldata)
adata


# In[38]:


adata.obs['orig.ident'].value_counts()


# In[39]:


adata.obs.index.name = "cells"
adata.write_h5ad(f'UPS.adata_ldata.h5ad')


# In[5]:


adata = sc.read(f'UPS.adata_ldata.h5ad')
adata


# In[6]:


scv.pp.filter_and_normalize(adata)


# In[7]:


scv.pp.moments(adata, n_neighbors=30, n_pcs=30)


# In[8]:


scv.pl.proportions(adata, groupby='celltype')


# In[5]:


type="UPS"


# In[6]:


import gc
gc.collect()
#
temp_pre= f"{type}_nue.in_process2"
if False==os.path.exists(f"{temp_pre}.velo.gz.h5ad"):
    scv.tl.recover_dynamics(adata, var_names='all', n_jobs=20)
    scv.tl.velocity(adata, mode='dynamical')
    adata.write(f"{temp_pre}.velo.gz.h5ad", compression='gzip')
    print(">>Write to file")
else:
    adata = sc.read(f"{temp_pre}.velo.gz.h5ad", compression='gzip', ext="h5ad")
    print(">>read from file")


# In[7]:


scv.tl.velocity_graph(adata)
scv.tl.velocity_embedding(adata, basis="umap")


# In[10]:


scv.pl.velocity_embedding(adata, basis = 'umap', title="UPS",arrow_length=3, arrow_size=2,
                          save=f"{type}.embedding.pdf",
                          color="celltype")


# In[8]:


scv.settings.set_figure_params('scvelo', dpi=900, dpi_save=900)
fig, ax = plt.subplots()
ax.set_aspect(1)

scv.pl.velocity_embedding_grid(adata, basis='umap',color='celltype', title='UPS',
                               arrow_size=1, arrow_length=2, arrow_color="#D2691E",
                               alpha=0.1,
                               #density=0.9,
                               legend_loc='right margin',legend_fontsize=5,
                               show=True,
                               save=f"{type}.grid.pdf", #figsize=(10,10),
                               xlim=[-10,10],ylim=[-10,10], ax=ax)


# In[11]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
scv.settings.set_figure_params('scvelo', dpi=900, dpi_save=900)
#fig, ax = plt.subplots()
#ax.set_aspect(1)
scv.pl.velocity_embedding_stream(adata, basis='umap',color='celltype', title='UPS',
                               #arrow_size=1, ##arrow_length=2, 
                               #arrow_color="#D2691E",
                               #alpha=0.01, density=0.9,
                               legend_loc='right margin',legend_fontsize=5,
                               show=True,
                               save=f"{type}.stream.pdf")
                                #, #figsize=(10,10),
                               #xlim=[-10,10],ylim=[-10,10], 
                               #  ax=ax)


# In[132]:


scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', 
               save=f"{type}.latent_time.pdf",# figsize=(10,9),
               color_map='gnuplot', size=80)


# In[128]:


scv.tl.velocity_pseudotime(adata)


# In[130]:


help(scv.tl.latent_time)


# In[129]:


scv.pl.scatter(adata, color='velocity_pseudotime',save=f"{type}.pseudotime.pdf",# figsize=(10,9),
               color_map='gnuplot', size=80)


# In[33]:


celltype_MSC.index


# In[16]:


adata


# In[22]:


adata.obs.celltype.isin(['MSC'])


# In[63]:


celltype_MSC = adata.obs.celltype[adata.obs.celltype.isin(["MSC"])].copy()


# In[57]:


celltype_DLK1 = adata.obs.celltype[adata.obs.celltype.isin(["DLK1 UPS"])].copy()


# In[93]:


celltype_DLK1


# In[96]:


celltype_MSC.index[2]


# In[114]:


MSC = ""


# In[115]:


for i in range(1,len(celltype_MSC.index)):
    MSC = MSC + celltype_MSC.index[i] + "|"


# In[116]:


MSC_1 = MSC[:len(MSC)-1]


# In[117]:


MSC_1


# In[118]:


DLK1 = ""


# In[119]:


for i in range(1,len(celltype_DLK1.index)):
    DLK1 = DLK1 + celltype_DLK1.index[i] + "|"


# In[120]:


DLK1_1 = DLK1[:len(DLK1)-1]


# In[121]:


DLK1_1


# In[64]:


celltype_MSC.index


# In[65]:


np.savetxt('MSC_index.csv', celltype_MSC.index, fmt='%s')


# In[66]:


np.savetxt('DLK1_index.csv', celltype_DLK1.index, fmt='%s')


# In[89]:


MSC_str = pd.read_csv('MSC_STR.csv',dtype = "str")


# In[91]:


DLK1_str = pd.read_csv('DLK1_STR.csv',dtype = "str")


# In[92]:


DLK1_str


# In[84]:


help(pd.read_csv)


# In[90]:


MSC_str


# In[50]:


test = celltype_MSC.index.tolist().to_csv(MSC.csv)


# In[45]:


test


# In[37]:


f"index_{str(celltype_MSC.index)}


# In[13]:


help(scv.tl.latent_time)

