#sh run_all_simple_benchmark.sh
#source activate fastai
import pickle as pk
import xgboost as xgb
import numpy as np
import pandas as pd
import matplotlib as plt
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import Binarizer
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import classification_report
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
import time
import random

import os
from scipy.io import mmread

def load_mtx(genome_dir):
        barcodes_tsv = os.path.join(genome_dir, "barcodes.tsv")
        genes_tsv = os.path.join(genome_dir, "genes.tsv")
        matrix_mtx = os.path.join(genome_dir, "matrix.mtx")
        for filepath in [barcodes_tsv, genes_tsv, matrix_mtx]:
            if not os.path.exists(filepath):
                raise IOError("Required file not found: %s" % filepath)
        barcodes = pd.read_csv(barcodes_tsv, delimiter='\t', header=None, usecols=[0]).values.squeeze()
        genes = pd.read_csv(genes_tsv, delimiter='\t', header=None, usecols=[1]).values.squeeze()
        #genes = [cr_constants.Gene(gene_id, None, None, None, None) for gene_id in genes]
        matrix = mmread(matrix_mtx)
        
        B = matrix.todense()
        ge = pd.DataFrame(B, range(1, B.shape[0] + 1), range(1, B.shape[1] + 1))
        
        ge.index = genes
        return ge

def ge_transform(df_GE, genes, doScale=True):
    scaler = MinMaxScaler()
    print(len(df_GE))
    binarizer = Binarizer(threshold=0.01)
    df_features = df_GE.transpose()
    print(len(df_features))
    #df_features = df_features.groupby(df_features.columns, axis=1).agg(max)
    
    #missing_genes = set(genes).difference( set( df_GE.index )  )
    #print(missing_genes)
    
    #df_features.set_index([list(missing_genes)], append=True)
    #for i in missing_genes:
    #    df_features[i] = 0
    
    df_features = df_features[genes]
    
    if doScale == True:
        scaler.fit(df_features)
        df_features = scaler.transform(df_features)
    
    binarizer.fit(df_features)
    df_features = binarizer.transform(df_features)
    print(len(df_features))
    df_features = pd.DataFrame(df_features)
    
    df_features.columns = genes
    return df_features

def simplifyPredTsub(y_pred): 
    y_pred_simple = [1 if x==0 else x for x in y_pred]
    y_pred_simple = [8 if x==3 else x for x in y_pred_simple]
    y_pred_simple = [8 if x==9 else x for x in y_pred_simple]
    y_pred_simple = [8 if x==10 else x for x in y_pred_simple] #DPT should be considered as cytotoxic
    y_pred_simple = [6 if x==7 else x for x in y_pred_simple] #cytotoxic
    y_pred_simple = [6 if x==11 else x for x in y_pred_simple]
    #y_pred_simple = [6 if x==13 else x for x in y_pred_simple]
    
    y_pred_simple = [2 if x==12 else x for x in y_pred_simple]
    y_pred_simple = [2 if x==14 else x for x in y_pred_simple]
    y_pred_simple = [2 if x==15 else x for x in y_pred_simple]
    y_pred_simple = [1 if x==16 else x for x in y_pred_simple]
    y_pred_simple = [1 if x==17 else x for x in y_pred_simple]
    return y_pred_simple

def simplifyPredCoarse(y_pred):
    y_pred_simple = [1 if x==0 else x for x in y_pred]
    y_pred_simple = [6 if x==3 else x for x in y_pred_simple]
    y_pred_simple = [6 if x==9 else x for x in y_pred_simple]
    y_pred_simple = [6 if x==10 else x for x in y_pred_simple]
    y_pred_simple = [6 if x==7 else x for x in y_pred_simple]
    y_pred_simple = [6 if x==11 else x for x in y_pred_simple]
    #y_pred_simple = [6 if x==13 else x for x in y_pred_simple]
    y_pred_simple = [6 if x==8 else x for x in y_pred_simple]
    
    y_pred_simple = [2 if x==12 else x for x in y_pred_simple]
    y_pred_simple = [2 if x==14 else x for x in y_pred_simple]
    y_pred_simple = [2 if x==15 else x for x in y_pred_simple]    
    
    y_pred_simple = [1 if x==16 else x for x in y_pred_simple]
    y_pred_simple = [1 if x==17 else x for x in y_pred_simple]

    return y_pred_simple

param = {'booster':'gbtree', 'objective':'multi:softmax','num_class':10, 
            'eval_metric':'mlogloss',
            #'learning_rate':0.05,
            #'max_depth':3,
            #'n_estimators':50,
            #'min_child_weight':1,
            #'colsample_bylevel':1,
            #'subsample':0.5,
            #'reg_alpha':0,
            #'reg_lambda':0,
            #'gamma':0,
            #'random_state':10,
            #'tree_method':'gpu_hist'
            }


#
print("Loading selected genes")

genes_raw = pd.read_csv('mountArvados/by_id/fa4af28a57cad7c66d23b9fd5d702a28+1405/dea_selected_genes_0.05_5_gex.tsv') #, header=None
set_new = genes_raw.x.values.flatten().tolist()


availablegenes = pd.read_csv('mountArvados/by_id/fa4af28a57cad7c66d23b9fd5d702a28+1405/availablegenes.tsv', header=None)
set_new = set(set_new).intersection( set(availablegenes.values.flatten().tolist() ) )


## Load data sets

print("Loading data sets")


#df_GE_CBMC = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection/cbmc_filtered_scran.csv')
#set_new = set(set_new).intersection( set(df_GE_CBMC.index.values.flatten().tolist() ) )
#del df_GE_CBMC

#df_GE_SORT = pd.read_csv('/pstore/scratch/u/ouyangt/data/benchmark/pbmcsort_multigenes_dea_scran.csv')
#set_new = set(set_new).intersection( set(df_GE_SORT.columns.values.flatten().tolist() ) )

#print(len(set_new))

#df_LABEL9_SORT = pd.read_csv('/pstore/scratch/u/ouyangt/data/benchmark/pbmcsort_label.csv')
#df_GE_SORT_TRANSFORMED = ge_transform(np.log1p(df_GE_SORT.fillna(0).iloc[0:85423,:].transpose()), list(set_new))
#del df_GE_SORT
#d_SORT = xgb.DMatrix(df_GE_SORT_TRANSFORMED , label=df_LABEL9_SORT.truth[0:85423])

df_GE_TENX1K = pd.read_csv('mountArvados/by_id/59876256d61cf0596bd625a9f195660a+2629/pbmc_1k_v3_filtered_scran.csv') #_scran norm by size later
df_LABEL9_TENX1K = pd.read_csv('mountArvados/by_id/fa4af28a57cad7c66d23b9fd5d702a28+1405/pbmc1k_celltype_clean.csv')
df_GE_TENX1K_TRANSFORMED = ge_transform(df_GE_TENX1K, list(set_new))
del df_GE_TENX1K
d_TENX1K = xgb.DMatrix(df_GE_TENX1K_TRANSFORMED, label=df_LABEL9_TENX1K)

#df_GE_TENX10gex = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection/pbmc10k_gex_filtered_scran.csv') #5k_pbmc_v3_filtered.csv
#df_LABEL9_TENX10gex = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection2/pbmc_10kgex_celltype_clean.csv')
#df_GE_TENX10gex_TRANSFORMED = ge_transform(df_GE_TENX10gex, list(set_new))
#del df_GE_TENX10gex
#d_TENX10gex = xgb.DMatrix(df_GE_TENX10gex_TRANSFORMED, label=df_LABEL9_TENX10gex)

df_GE_TENX10v3 = pd.read_csv('mountArvados/by_id/59876256d61cf0596bd625a9f195660a+2629/pbmc10k_v3_filtered_scran.csv')
df_LABEL9_TENX10v3 = pd.read_csv('mountArvados/by_id/fa4af28a57cad7c66d23b9fd5d702a28+1405/pbmc_10kv3_celltype_clean.csv')
df_GE_TENX10v3_TRANSFORMED = ge_transform(df_GE_TENX10v3, list(set_new))
del df_GE_TENX10v3
d_TENX10v3 = xgb.DMatrix(df_GE_TENX10v3_TRANSFORMED, label=df_LABEL9_TENX10v3)

#df_GE_CBMC = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection/cbmc_filtered_scran.csv')
#df_LABEL_CBMC = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection2/cbmc_celltype_clean_asnk.csv')
#df_GE_CBMC_TRANSFORMED = ge_transform(df_GE_CBMC, list(set_new)) 
#del df_GE_CBMC
#d_CBMC = xgb.DMatrix(df_GE_CBMC_TRANSFORMED, label=df_LABEL_CBMC)

#df_GE_MALT = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection/malt_filtered_scran.csv')
#df_LABEL_MALT = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection2/malt_celltype_clean.csv')
#df_GE_MALT_TRANSFORMED = ge_transform(df_GE_MALT, list(set_new))
#del df_GE_MALT
#d_MALT = xgb.DMatrix(df_GE_MALT_TRANSFORMED , label=df_LABEL_MALT)


#df_GE_TENX5v3 = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection/5k_pbmc_v3_filtered_scran.csv') #5k_pbmc_v3_filtered.csv
#df_GE_TENX5v3_TRANSFORMED = ge_transform(df_GE_TENX5v3, list(set_new))
#del df_GE_TENX5v3
#df_LABEL9_TENX5v3 = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection2/pbmc5k_celltype_clean.csv')
#d_TENX5v3 = xgb.DMatrix(df_GE_TENX5v3_TRANSFORMED , label=df_LABEL9_TENX5v3)


#df_GE_UNKNOWN = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection/unknown_50.csv')
#df_GE_UNKNOWN_TRANSFORMED = ge_transform(df_GE_UNKNOWN, list(set_new), doScale = False )
#del df_GE_UNKNOWN
#df_LABEL9_UNKNOWN = pd.read_csv('/pstore/scratch/u/ouyangt/data/feature_selection/unknown_50_celltype.csv')


#df_MULTI_TRANSFORMED = pd.concat([df_GE_TENX10v3_TRANSFORMED, #Must have all three -> PBMCsort
#                                  df_GE_TENX10gex_TRANSFORMED,
#                                  df_GE_CBMC_TRANSFORMED
#                                  ], ignore_index=False, axis=0)  #, df_GE_UNKNOWN_TRANSFORMED

df_MULTI_TRANSFORMED = df_GE_TENX1K_TRANSFORMED

#labels_MULTI_TRANSFORMED = pd.DataFrame(np.concatenate(( df_LABEL9_TENX10v3,
#                                                        df_LABEL9_TENX10gex,
#                                                        df_LABEL_CBMC)), columns=range(1)) #, df_LABEL9_UNKNOWN

labels_MULTI_TRANSFORMED = df_LABEL9_TENX1K
#labels_MULTI_TRANSFORMED = simplifyPredTsub(labels_MULTI_TRANSFORMED)

d_MULTI = xgb.DMatrix(df_MULTI_TRANSFORMED, label=labels_MULTI_TRANSFORMED)

## Model training 

print("Training models")

print("Training XGBoost")

start = time.time()
clfXgb = xgb.train(param, d_MULTI)
end = time.time()
print(end - start)

file = open('xgb_model.pkl', 'wb')
pk.dump(clfXgb, file)
file.close()

print("Training Logistic")

clfLog = LogisticRegression(random_state=822, 
                         solver='lbfgs', 
                         multi_class='multinomial')

#clfLog = CalibratedClassifierCV(ClassifierLog)

start = time.time()
clfLog.fit(df_MULTI_TRANSFORMED, labels_MULTI_TRANSFORMED)
end = time.time()
print(end - start)

file = open('logistic_model_clean.pkl', 'wb')
pk.dump(clfLog, file)
file.close()

#AttributeError: 'CalibratedClassifierCV' object has no attribute 'coef_'
pd.DataFrame(clfLog.coef_).to_csv("logistic_model_coef_scran_clean_dea.csv")

print("Training SVC")

ClassifierSVC = LinearSVC() #C=0.8
clfSVC = CalibratedClassifierCV(ClassifierSVC)
start = time.time()
clfSVC.fit(df_MULTI_TRANSFORMED, labels_MULTI_TRANSFORMED)
end = time.time()
print(end - start)

file = open('svm_model_clean.pkl', 'wb')
pk.dump(clfSVC, file)
file.close()

## Prediction 

score_list_dataset = []

score_list_xgb = []
score_list_log = []
score_list_svm = []

score_list_xgb_sub = []
score_list_log_sub = []
score_list_svm_sub = []

score_list_xgb_coarse = []
score_list_log_coarse = []
score_list_svm_coarse = []

#10kv3

print("10kv3")
label_this = df_LABEL9_TENX10v3

start = time.time()
pred_this = clfXgb.predict(d_TENX10v3)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("xgb_pred_10kv3_scran_clean_dea.csv")


score_list_dataset.insert(len(score_list_dataset), 'pbmc10kv3')

score_list_xgb.insert(len(score_list_xgb), accuracy_score(label_this, pred_this))
score_list_xgb_sub.insert(len(score_list_xgb_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_xgb_coarse.insert(len(score_list_xgb_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfLog.predict(df_GE_TENX10v3_TRANSFORMED)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("log_pred_10kv3_scran_clean_dea.csv")

score_list_log.insert(len(score_list_log), accuracy_score(label_this, pred_this))
score_list_log_sub.insert(len(score_list_log_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_log_coarse.insert(len(score_list_log_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfSVC.predict(df_GE_TENX10v3_TRANSFORMED)
end = time.time()
print(end - start)

score_list_svm.insert(len(score_list_svm), accuracy_score(label_this, pred_this))
score_list_svm_sub.insert(len(score_list_svm_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_svm_coarse.insert(len(score_list_svm_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )


pd.DataFrame(pred_this).to_csv("svm_pred_10kv3_scran_clean_dea.csv")
'''
#10gex

print("10kng")
label_this = df_LABEL9_TENX10gex

start = time.time()
pred_this = clfXgb.predict(d_TENX10gex)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("xgb_pred_10gex_scran_clean_dea.csv")


score_list_dataset.insert(len(score_list_dataset), 'pbmc10kng')

score_list_xgb.insert(len(score_list_xgb), accuracy_score(label_this, pred_this))
score_list_xgb_sub.insert(len(score_list_xgb_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_xgb_coarse.insert(len(score_list_xgb_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfLog.predict(df_GE_TENX10gex_TRANSFORMED)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("log_pred_10gex_scran_clean_dea.csv")

score_list_log.insert(len(score_list_log), accuracy_score(label_this, pred_this))
score_list_log_sub.insert(len(score_list_log_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_log_coarse.insert(len(score_list_log_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfSVC.predict(df_GE_TENX10gex_TRANSFORMED)
end = time.time()
print(end - start)


score_list_svm.insert(len(score_list_svm), accuracy_score(label_this, pred_this))
score_list_svm_sub.insert(len(score_list_svm_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_svm_coarse.insert(len(score_list_svm_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

pd.DataFrame(pred_this).to_csv("svm_pred_10gex_scran_clean_dea.csv")

#CBMC
print("CBMC")
label_this = df_LABEL_CBMC

start = time.time()
pred_this = clfXgb.predict(d_CBMC)
end = time.time()
print(end - start)

score_list_dataset.insert(len(score_list_dataset), 'cbmc')

score_list_xgb.insert(len(score_list_xgb), accuracy_score(label_this, pred_this))
score_list_xgb_sub.insert(len(score_list_xgb_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_xgb_coarse.insert(len(score_list_xgb_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

pd.DataFrame(pred_this).to_csv("xgb_pred_cbmc_scran_clean_dea.csv")

start = time.time()
pred_this = clfLog.predict(df_GE_CBMC_TRANSFORMED)
end = time.time()
print(end - start)

score_list_log.insert(len(score_list_log), accuracy_score(label_this, pred_this))
score_list_log_sub.insert(len(score_list_log_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_log_coarse.insert(len(score_list_log_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

pd.DataFrame(pred_this).to_csv("log_pred_cbmc_scran_clean_dea.csv")

start = time.time()
pred_this = clfSVC.predict(df_GE_CBMC_TRANSFORMED)
end = time.time()
print(end - start)

score_list_svm.insert(len(score_list_svm), accuracy_score(label_this, pred_this))
score_list_svm_sub.insert(len(score_list_svm_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_svm_coarse.insert(len(score_list_svm_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

pd.DataFrame(pred_this).to_csv("svm_pred_cbmc_scran_clean_dea.csv")

#1K 
print("1k")

label_this = df_LABEL9_TENX1K

start = time.time()
pred_this = clfXgb.predict(d_TENX1K)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("xgb_pred_1k_scran_clean_dea.csv")

score_list_dataset.insert(len(score_list_dataset), 'pbmc1k')

score_list_xgb.insert(len(score_list_xgb), accuracy_score(label_this, pred_this))
score_list_xgb_sub.insert(len(score_list_xgb_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_xgb_coarse.insert(len(score_list_xgb_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfLog.predict(df_GE_TENX1K_TRANSFORMED)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("log_pred_1k_scran_clean_dea.csv")

score_list_log.insert(len(score_list_log), accuracy_score(label_this, pred_this))
score_list_log_sub.insert(len(score_list_log_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_log_coarse.insert(len(score_list_log_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfSVC.predict(df_GE_TENX1K_TRANSFORMED)
end = time.time()
print(end - start)

score_list_svm.insert(len(score_list_svm), accuracy_score(label_this, pred_this))
score_list_svm_sub.insert(len(score_list_svm_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_svm_coarse.insert(len(score_list_svm_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

pd.DataFrame(pred_this).to_csv("svm_pred_1k_scran_clean_dea.csv")


# pbmc5kv3

print("5kv3")

label_this = df_LABEL9_TENX5v3

start = time.time()
pred_this = clfXgb.predict(d_TENX5v3)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("xgb_pred_5k_scran_clean_dea.csv")

score_list_dataset.insert(len(score_list_dataset), 'pbmc5k')

score_list_xgb.insert(len(score_list_xgb), accuracy_score(label_this, pred_this))
score_list_xgb_sub.insert(len(score_list_xgb_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_xgb_coarse.insert(len(score_list_xgb_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfLog.predict(df_GE_TENX5v3_TRANSFORMED)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("log_pred_5k_scran_clean_dea.csv")

score_list_log.insert(len(score_list_log), accuracy_score(label_this, pred_this))
score_list_log_sub.insert(len(score_list_log_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_log_coarse.insert(len(score_list_log_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfSVC.predict(df_GE_TENX5v3_TRANSFORMED)
end = time.time()
print(end - start)

score_list_svm.insert(len(score_list_svm), accuracy_score(label_this, pred_this))
score_list_svm_sub.insert(len(score_list_svm_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_svm_coarse.insert(len(score_list_svm_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

pd.DataFrame(pred_this).to_csv("svm_pred_5k_scran_clean_dea.csv")


#MALT
print("MALT")

label_this = df_LABEL_MALT

start = time.time()
pred_this = clfXgb.predict(d_MALT)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("xgb_pred_malt_scran_clean_dea.csv")


score_list_dataset.insert(len(score_list_dataset), 'malt')

score_list_xgb.insert(len(score_list_xgb), accuracy_score(label_this, pred_this))
score_list_xgb_sub.insert(len(score_list_xgb_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_xgb_coarse.insert(len(score_list_xgb_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfLog.predict(df_GE_MALT_TRANSFORMED)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("log_pred_malt_scran_clean_dea.csv")


score_list_log.insert(len(score_list_log), accuracy_score(label_this, pred_this))
score_list_log_sub.insert(len(score_list_log_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_log_coarse.insert(len(score_list_log_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfSVC.predict(df_GE_MALT_TRANSFORMED)
end = time.time()
print(end - start)

score_list_svm.insert(len(score_list_svm), accuracy_score(label_this, pred_this))
score_list_svm_sub.insert(len(score_list_svm_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_svm_coarse.insert(len(score_list_svm_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

pd.DataFrame(pred_this).to_csv("svm_pred_malt_scran_clean_dea.csv")


# SORT

print("Sort")

label_this = df_LABEL9_SORT[0:85423]

start = time.time()
pred_this = clfXgb.predict(d_SORT)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("xgb_pred_sort_scran_clean_dea.csv")

score_list_dataset.insert(len(score_list_dataset), 'pbmcsort')


score_list_xgb.insert(len(score_list_xgb), accuracy_score(label_this.truth, pred_this))
score_list_xgb_sub.insert(len(score_list_xgb_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_xgb_coarse.insert(len(score_list_xgb_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfLog.predict(df_GE_SORT_TRANSFORMED)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("log_pred_sort_scran_clean_dea.csv")

score_list_log.insert(len(score_list_log), accuracy_score(label_this.truth, pred_this))
score_list_log_sub.insert(len(score_list_log_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_log_coarse.insert(len(score_list_log_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )

start = time.time()
pred_this = clfSVC.predict(df_GE_SORT_TRANSFORMED)
end = time.time()
print(end - start)

pd.DataFrame(pred_this).to_csv("svm_pred_sort_scran_clean_dea.csv")


score_list_svm.insert(len(score_list_svm), accuracy_score(label_this.truth, pred_this))
score_list_svm_sub.insert(len(score_list_svm_sub), accuracy_score(simplifyPredTsub(label_this.truth), simplifyPredTsub(pred_this) ) )
score_list_svm_coarse.insert(len(score_list_svm_coarse), accuracy_score(simplifyPredCoarse(label_this.truth), simplifyPredCoarse(pred_this) ) )
'''

#Save f-1 scores

pd.DataFrame(score_list_dataset).to_csv("dataset_id_scran_clean_dea.csv")

pd.DataFrame(score_list_xgb).to_csv("xgb_f1scores_t3_scran_clean_dea.csv")
pd.DataFrame(score_list_log).to_csv("log_f1scores_t3_scran_clean_dea.csv")
pd.DataFrame(score_list_svm).to_csv("svm_f1scores_t3_scran_clean_dea.csv")

#pd.DataFrame(score_list_xgb_sub).to_csv("xgb_f1scores_t2_scran_clean_dea.csv")
#pd.DataFrame(score_list_log_sub).to_csv("log_f1scores_t2_scran_clean_dea.csv")
#pd.DataFrame(score_list_svm_sub).to_csv("svm_f1scores_t2_scran_clean_dea.csv")

#pd.DataFrame(score_list_xgb_coarse).to_csv("xgb_f1scores_t1_scran_clean_dea.csv")
#pd.DataFrame(score_list_log_coarse).to_csv("log_f1scores_t1_scran_clean_dea.csv")
#pd.DataFrame(score_list_svm_coarse).to_csv("svm_f1scores_t1_scran_clean_dea.csv")

print("Done")