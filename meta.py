

import pandas as pd
import numpy as np

def ensembl2entrez(ensembl):
    import mygene
    import math
    mg = mygene.MyGeneInfo()
    ens = ensembl
    ginfo = mg.querymany(ens, scopes='ensembl.gene', returnall=False,fields='entrezgene',as_dataframe=True)
    ginfo['index'] = ginfo.index
    ginfo.drop_duplicates(subset='index', keep='first', inplace=True)
    entrez1 = []
    for item in ensembl:
        entrez1.append(ginfo.loc[item,'entrezgene'])
    entrez2 = []
    for item in entrez1:
        try:
            math.isnan(item)
            entrez2.append('***')
        except TypeError:  # item is not a real number, means it is legit
            entrez2.append(item)
    return entrez2

def check(y):
    type_ = [type(i) for i in y ]
    truth = [ False if i == str else True for i in type_]
    y = y[truth]
    y = np.where(y<=0,0.005,y)
    return y

def hedges_g(ye,yc):  # ye is a ndarray of experimental group, yc is a ndarray of control group
    import math
    ye_bar = np.mean(ye)
    yc_bar = np.mean(yc)
    ne = len(ye)
    nc = len(yc)
    m = ne + nc -2
    cm = 1 - 3/(4*m-1)
    ye_sd = np.std(ye)
    yc_sd = np.std(yc)
    S_pool = math.sqrt( ((ne-1)*ye_sd**2 + (nc-1)*yc_sd**2) / m )
    T = cm * (ye_bar - yc_bar) / S_pool
    V = (ne+nc)/(ne*nc) + (T**2)/(2*(ne+nc))
    return ye_bar,yc_bar,S_pool,ne,nc,cm,T,V


# gtex skin TPM matrix
skin = set(list(pd.read_csv('./skin.txt',sep='\t',header=None)[0]))
print('loading gtex')
gtex = pd.read_csv('./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct',sep='\t',
                   skiprows=[0,1])
print('common sample')
coln = set(gtex.columns.tolist())
common = ['Name','Description'] + list(coln.intersection(skin)) 
gtex[common].to_csv('./gtex_skin.txt',sep='\t',index=None)

# gtex transfer to entrez
gtex = pd.read_csv('./gtex_skin.txt',sep='\t')
print('finished loading')
id =  [item[0:15] for item in gtex['Name']]
entrez = ensembl2entrez(id)
print(len(entrez))
gtex['Name'] = entrez
gtex.to_csv('./gtex_skin_entrez.txt',sep='\t',index=None)

# find common genes
prefix = './melanoma1/'
study1 = pd.read_csv(prefix + 'skcm_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt',sep='\t')
study2 = pd.read_csv(prefix + 'mel_tsam_liang_2017/data_RNA_Seq_expression_median.txt',sep='\t')
study3 = pd.read_csv(prefix + 'skcm_dfci_2015/data_RNA_Seq_expression_median.txt',sep='\t')
study4 = pd.read_csv(prefix + 'skcm_mskcc_2014/data_RNA_Seq_expression_median.txt',sep='\t')
study5 = pd.read_csv(prefix + 'skcm_tcga/data_RNA_Seq_v2_expression_median.txt',sep='\t')
gtex = pd.read_csv('./gtex_skin_entrez.txt',sep='\t')

study1_set = set(list(study1['Entrez_Gene_Id']))
study2_set = set(list(study2['Entrez_Gene_Id']))
study3_set = set(list(study3['Entrez_Gene_Id']))
study4_set = set(list(study4['Entrez_Gene_Id']))
study5_set = set(list(study5['Entrez_Gene_Id']))
gtex_set = set(list(gtex['Name']))
gtex_set.remove('***')
gtex_set = [int(item) for item in gtex_set]


common = study1_set.intersection(study2_set,study3_set,study4_set,study5_set,gtex_set)

# retain common genes for each study
study1_fine = study1.loc[[True if item in common else False for item in study1['Entrez_Gene_Id']]]
study1_fine.drop_duplicates(subset='Entrez_Gene_Id',inplace=True)
study1_fine = study1_fine.set_index(pd.Index(np.arange(study1_fine.shape[0])))
study2_fine = study2.loc[[True if item in common else False for item in study2['Entrez_Gene_Id']]]
study2_fine.drop_duplicates(subset='Entrez_Gene_Id',inplace=True)
study2_fine = study2_fine.set_index(pd.Index(np.arange(study2_fine.shape[0])))
study3_fine = study3.loc[[True if item in common else False for item in study3['Entrez_Gene_Id']]]
study3_fine.drop_duplicates(subset='Entrez_Gene_Id',inplace=True)
study3_fine = study3_fine.set_index(pd.Index(np.arange(study3_fine.shape[0])))
study4_fine = study4.loc[[True if item in common else False for item in study4['Entrez_Gene_Id']]]
study4_fine.drop_duplicates(subset='Entrez_Gene_Id',inplace=True)
study4_fine = study4_fine.set_index(pd.Index(np.arange(study4_fine.shape[0])))
study5_fine = study5.loc[[True if item in common else False for item in study5['Entrez_Gene_Id']]]
study5_fine.drop_duplicates(subset='Entrez_Gene_Id',inplace=True)
study5_fine = study5_fine.set_index(pd.Index(np.arange(study5_fine.shape[0])))
common_str = [str(item) for item in common]
gtex_fine = gtex.loc[[True if item in common_str else False for item in gtex['Name']]]
gtex_fine.drop_duplicates(subset='Name',inplace=True)
gtex_fine = gtex_fine.set_index(pd.Index(np.arange(gtex_fine.shape[0])))

# save
study1_fine.to_csv('./melanoma1/study1_fine.txt',sep='\t',index=None)
study2_fine.to_csv('./melanoma1/study2_fine.txt',sep='\t',index=None)
study3_fine.to_csv('./melanoma1/study3_fine.txt',sep='\t',index=None)
study4_fine.to_csv('./melanoma1/study4_fine.txt',sep='\t',index=None)
study5_fine.to_csv('./melanoma1/study5_fine.txt',sep='\t',index=None)
gtex_fine.to_csv('./gtex_skin_entrez_fine.txt',sep='\t',index=None)

import pickle
with open('./common_entrez.p','wb') as f:
    pickle.dump(common,f)

# hedges g estimator per study per gene
import pickle
with open('./common_entrez.p','rb') as f:
    common = pickle.load(f)
common = list(common)

study1 = pd.read_csv('./melanoma1/study1_fine.txt',sep='\t')
study2 = pd.read_csv('./melanoma1/study2_fine.txt',sep='\t')
study3 = pd.read_csv('./melanoma1/study3_fine.txt',sep='\t')
study4 = pd.read_csv('./melanoma1/study4_fine.txt',sep='\t')
study5 = pd.read_csv('./melanoma1/study5_fine.txt',sep='\t')
gtex = pd.read_csv('./gtex_skin_entrez_fine.txt',sep='\t')

result = np.empty([len(common),5,3])
study_all = [study1,study2,study3,study4,study5]

for i in range(len(common)):
    gene = common[i]
    gtex.index = gtex['Name']
    yc = gtex.loc[gene,:].values[2:]
    for j in range(len(study_all)):
        study = study_all[j]
        study.index = study['Entrez_Gene_Id']
        ye = study.loc[gene,:].values[2:]

        ye,yc = check(ye),check(yc)
        ye = ye.astype(np.float64)
        yc = yc.astype(np.float64)
        if j != 1:    # liang genome research paper seems to have log transform the FPKM
            ye = np.log2(ye)
        yc = np.log2(yc)
        try:
            ye_bar,yc_bar,S_pool,ne,nc,cm, t,v = hedges_g(ye,yc)
        except:
            print(gene)
            raise Exception
        result[i,j,:] = [int(j+1),t,v]

with open('./skin_result.p','wb') as f:
    pickle.dump(result,f)


# Let's move to glioblastoma
# gtex skin TPM matrix
brain = set(list(pd.read_csv('./glioblastoma/brain.txt',sep='\t',header=None)[0]))
print('loading gtex')
gtex = pd.read_csv('./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct',sep='\t',
                   skiprows=[0,1])
print('common sample')
coln = set(gtex.columns.tolist())
common = ['Name','Description'] + list(coln.intersection(brain)) 
gtex[common].to_csv('./glioblastoma/gtex_brain.txt',sep='\t',index=None)

# gtex transfer to entrez
gtex = pd.read_csv('./glioblastoma/gtex_brain.txt',sep='\t')
print('finished loading')
id =  [item[0:15] for item in gtex['Name']]
entrez = ensembl2entrez(id)
print(len(entrez))
gtex['Name'] = entrez
gtex.to_csv('./glioblastoma/gtex_brain_entrez.txt',sep='\t',index=None)

# find common genes
prefix = './glioblastoma/'
study1 = pd.read_csv(prefix + 'gbm_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt',sep='\t')
study2 = pd.read_csv(prefix + 'gbm_tcga_pub2013/data_RNA_Seq_v2_expression_median.txt',sep='\t')
study3 = pd.read_csv(prefix + 'gbm_tcga/data_RNA_Seq_v2_expression_median.txt',sep='\t')
gtex = pd.read_csv('./glioblastoma/gtex_brain_entrez.txt',sep='\t')

study1_set = set(list(study1['Entrez_Gene_Id']))
study2_set = set(list(study2['Entrez_Gene_Id']))
study3_set = set(list(study3['Entrez_Gene_Id']))
gtex_set = set(list(gtex['Name']))
gtex_set.remove('***')
gtex_set = [int(item) for item in gtex_set]


common = study1_set.intersection(study2_set,study3_set,gtex_set)

# retain common genes for each study
study1_fine = study1.loc[[True if item in common else False for item in study1['Entrez_Gene_Id']]]
study1_fine.drop_duplicates(subset='Entrez_Gene_Id',inplace=True)
study1_fine = study1_fine.set_index(pd.Index(np.arange(study1_fine.shape[0])))
study2_fine = study2.loc[[True if item in common else False for item in study2['Entrez_Gene_Id']]]
study2_fine.drop_duplicates(subset='Entrez_Gene_Id',inplace=True)
study2_fine = study2_fine.set_index(pd.Index(np.arange(study2_fine.shape[0])))
study3_fine = study3.loc[[True if item in common else False for item in study3['Entrez_Gene_Id']]]
study3_fine.drop_duplicates(subset='Entrez_Gene_Id',inplace=True)
study3_fine = study3_fine.set_index(pd.Index(np.arange(study3_fine.shape[0])))
common_str = [str(item) for item in common]
gtex_fine = gtex.loc[[True if item in common_str else False for item in gtex['Name']]]
gtex_fine.drop_duplicates(subset='Name',inplace=True)
gtex_fine = gtex_fine.set_index(pd.Index(np.arange(gtex_fine.shape[0])))

# save
study1_fine.to_csv('./glioblastoma/study1_fine.txt',sep='\t',index=None)
study2_fine.to_csv('./glioblastoma/study2_fine.txt',sep='\t',index=None)
study3_fine.to_csv('./glioblastoma/study3_fine.txt',sep='\t',index=None)
gtex_fine.to_csv('./glioblastoma/gtex_brain_entrez_fine.txt',sep='\t',index=None)

import pickle
with open('./glioblastoma/common_entrez.p','wb') as f:
    pickle.dump(common,f)


# hedges g estimator per study per gene
import pickle
with open('./glioblastoma/common_entrez.p','rb') as f:
    common = pickle.load(f)
common = list(common)

study1 = pd.read_csv('./glioblastoma/study1_fine.txt',sep='\t')
study2 = pd.read_csv('./glioblastoma/study2_fine.txt',sep='\t')
study3 = pd.read_csv('./glioblastoma/study3_fine.txt',sep='\t')

gtex = pd.read_csv('./glioblastoma/gtex_brain_entrez_fine.txt',sep='\t')

result = np.empty([len(common),3,3])
study_all = [study1,study2,study3]

for i in range(len(common)):
    gene = common[i]
    gtex.index = gtex['Name']
    yc = gtex.loc[gene,:].values[2:]
    for j in range(len(study_all)):
        study = study_all[j]
        study.index = study['Entrez_Gene_Id']
        ye = study.loc[gene,:].values[2:]

        ye,yc = check(ye),check(yc)
        ye = ye.astype(np.float64)
        yc = yc.astype(np.float64)

        ye = np.log2(ye)
        yc = np.log2(yc)
        try:
            ye_bar,yc_bar,S_pool,ne,nc,cm, t,v = hedges_g(ye,yc)
        except:
            print(gene)
            raise Exception
        result[i,j,:] = [int(j+1),t,v]

with open('./glioblastoma/brain_result.p','wb') as f:
    pickle.dump(result,f)



# back from cluster
import pickle
with open('/Users/ligk2e/Desktop/meta/new/skin_result.p','rb') as f:
    skin_result = pickle.load(f)
with open('/Users/ligk2e/Desktop/meta/new/common_entrez.p','rb') as f:
    common = pickle.load(f)

common = list(common)
skin_result_inter = skin_result.reshape(skin_result.shape[0],-1)
np.savetxt('/Users/ligk2e/Desktop/meta/new/skin_result_inter.txt',skin_result_inter)
np.savetxt('/Users/ligk2e/Desktop/meta/new/skin_common.txt',np.array(list(common)).reshape(-1,1))

# back from R
crystal = pd.read_csv('/Users/ligk2e/Desktop/meta/new/skin_result_R.txt',sep=' ')
crystal_diff_df_m = crystal.loc[crystal['V3'] < 0.05]

# draw two representative forest plot
# let's pick entrez ID
common.index(16)  # 9
r = skin_result[58]


# let's draw forest plot
import matplotlib.pyplot as plt
import seaborn as sns

store = {}
for i in range(r.shape[0]):
    data = np.random.normal(r[i,1],r[i,2],1000)
    store['study{}'.format(i+1)] = data
store['combine'] = np.random.normal(4.13,0.88,1000)
df = pd.DataFrame(store)
df = df.stack().to_frame()
df.reset_index(level=1,inplace=True)
fig,ax = plt.subplots()
sns.pointplot(x=0,y='level_1',data=df,ax=ax,join=False)
ax.set_xlabel('Gene expression(Hedges\' g estimator)')
ax.set_ylabel('Studies')


crystal_diff_m = crystal.loc[crystal['V3'] < 0.05]['V5'].astype(np.int)
import mygene
mg = mygene.MyGeneInfo()
ginfo = mg.querymany(list(crystal_diff), scopes='entrezgene', returnall=False, fields='symbol', as_dataframe=True)
crystal_diff_symbol = ginfo['symbol']
crystal_diff_symbol.to_csv('/Users/ligk2e/Desktop/meta/new/skin_diff.txt',sep='\n',index=None)







### how about glioblastoma
import pickle
with open('/Users/ligk2e/Desktop/meta/new/glioblastoma/brain_result.p','rb') as f:
    brain_result = pickle.load(f)
with open('/Users/ligk2e/Desktop/meta/new/glioblastoma/common_entrez.p','rb') as f:
    common = pickle.load(f)
brain_result_inter = brain_result.reshape(brain_result.shape[0],-1)
np.savetxt('/Users/ligk2e/Desktop/meta/new/glioblastoma/brain_result_inter.txt',brain_result_inter)
np.savetxt('/Users/ligk2e/Desktop/meta/new/glioblastoma/brain_common.txt',np.array(list(common)).reshape(-1,1))

# back from R
crystal = pd.read_csv('/Users/ligk2e/Desktop/meta/new/brain_result_R.txt',sep=' ')
crystal_diff_df_g = crystal.loc[crystal['V3'] < 0.05]
common = list(common)
len(set(crystal_diff_df['V5'].values).intersection(set(crystal_diff_df_m['V5'].values)))

# draw two representative forest plot
# let's pick entrez ID
common.index(24)  # 9
r = brain_result[17]


# let's draw forest plot
import matplotlib.pyplot as plt
import seaborn as sns

store = {}
for i in range(r.shape[0]):
    data = np.random.normal(r[i,1],r[i,2],1000)
    store['study{}'.format(i+1)] = data
store['combine'] = np.random.normal(15.31,5.5,1000)
df = pd.DataFrame(store)
df = df.stack().to_frame()
df.reset_index(level=1,inplace=True)
fig,ax = plt.subplots()
sns.pointplot(x=0,y='level_1',data=df,ax=ax,join=False)
ax.set_xlabel('Gene expression(Hedges\' g estimator)')
ax.set_ylabel('Studies')


crystal_diff_g = crystal.loc[crystal['V3'] < 0.05]['V5'].astype(np.int)
import mygene
mg = mygene.MyGeneInfo()
ginfo = mg.querymany(list(crystal_diff), scopes='entrezgene', returnall=False, fields='symbol', as_dataframe=True)
crystal_diff_symbol = ginfo['symbol']
crystal_diff_symbol.to_csv('/Users/ligk2e/Desktop/meta/new/skin_diff.txt',sep='\n',index=None)


# shared in melanoma and glioblastoma
share = list(set(crystal_diff_g.values).intersection(set(crystal_diff_m)))
mela = crystal_diff_df_m.loc[ [True if int(item) in share else False for item in crystal_diff_df_m['V5']] ]

glio = crystal_diff_df_g.loc[ [True if int(item) in share else False for item in crystal_diff_df_g['V5']] ]

mela.sort_values(by='V5',inplace=True)
glio.sort_values(by='V5',inplace=True)

mela.set_index(pd.Index(np.arange(mela.shape[0])),inplace=True)
glio.set_index(pd.Index(np.arange(glio.shape[0])),inplace=True)

fig,ax = plt.subplots()
ax.scatter(mela['V1'],glio['V1'],color='blue',alpha=0.2)
ax.plot([0,60],[0,60],color='orange')
ax.set_xlabel('effect size in melanoma')
ax.set_ylabel('effect size in glioblastoma')

np.where(glio['V1'].values > 60)













