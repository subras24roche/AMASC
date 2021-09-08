#AMASC Marker Selection

import xgboost as xgb
import numpy as np
import pandas as pd
import matplotlib as plt
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import Binarizer

import time
import random
import sys

## Parameters ## 

RNA_DATA_PATH = sys.argv[1]
CELLTYPE_DATA_PATH = sys.argv[2]
OUTPUT_GENESET_PATH = sys.argv[3]
OUTPUT_ACCURACY_PATH = sys.argv[4]
num_exp = int(sys.argv[5])
data_split = 0.8
threshold_binarize = 0.01


print("\nLoading data")

df_LABEL = pd.read_csv(CELLTYPE_DATA_PATH, index_col = 0)
df_LABEL['truth'] = df_LABEL['truth'].replace(8,6)
#evens = df_LABEL[df_LABEL.truth.isin([2,4,6])].index
#df_LABEL= df_LABEL[df_LABEL.index.isin(evens)]
df_GE = pd.read_csv(RNA_DATA_PATH, index_col = 0)
df_GE.columns = df_GE.columns.str.replace('.','-')
#df_GE = df_GE[[col for col in df_GE.columns if col in evens]]

lamprecht_ge_raw = pd.read_csv("Lamprecht_1a_normalized_filt.csv", index_col = 0)
lamprecht_label = pd.read_csv("Lamprecht_1a_celltype_clean.csv", index_col = 0)
lamprecht_ge_raw = lamprecht_ge_raw[~lamprecht_ge_raw.index.duplicated(keep='first')]
lamprecht_ge_raw = lamprecht_ge_raw[lamprecht_ge_raw.index.isin(df_GE.index)]
df_GE = df_GE[df_GE.index.isin(lamprecht_ge_raw.index)]
lamprecht_ge_raw = lamprecht_ge_raw.reindex(df_GE.index)

print("\nData loaded")

def ge_transform(df_GE, genes):
    scaler = MinMaxScaler()
    print(len(df_GE))
    binarizer = Binarizer(threshold=threshold_binarize) 
    df_features = df_GE.transpose()
    print(len(df_features))
    df_features = df_features.groupby(df_features.columns, axis=1).agg(max)
    
    df_features = df_features[genes]
    scaler.fit(df_features)
    df_features = scaler.transform(df_features)
    binarizer.fit(df_features)
    df_features = binarizer.transform(df_features)
    print(len(df_features))
    df_features = pd.DataFrame(df_features)
    
    df_features.columns = genes
    return df_features

list_GEs_raw = df_GE

geneset = list(set(list_GEs_raw.index))

geneset_clean = [x for x in geneset if not x.startswith('MT-') and not x.startswith('#')]
geneset_clean = [x for x in geneset_clean if not x.startswith('RPL') and not x.startswith('#')]
geneset_clean = [x for x in geneset_clean if not x.startswith('RPS') and not x.startswith('#')]

list_GEs = ge_transform(list_GEs_raw, genes=geneset_clean)

del list_GEs_raw
del df_GE

list_labels = df_LABEL

#random.seed(102)
#split_sets = np.random.rand(len(list_labels)) < data_split
labels_train_global = list_labels
GEs_train_global = list_GEs

geneset = list(set(lamprecht_ge_raw.index))
geneset_clean = [x for x in geneset if not x.startswith('MT-') and not x.startswith('#')]
geneset_clean = [x for x in geneset_clean if not x.startswith('RPL') and not x.startswith('#')]
geneset_clean = [x for x in geneset_clean if not x.startswith('RPS') and not x.startswith('#')]
lamprecht_ge = ge_transform(lamprecht_ge_raw, genes=geneset_clean)
labels_test_global = lamprecht_label
GEs_test_global = lamprecht_ge

param = {'booster':'gbtree', 'objective':'multi:softmax','num_class':17, 
            'eval_metric':'mlogloss',
            'learning_rate':0.05,
            'max_depth':3,
            'n_estimators':100,
            'min_child_weight':1,
            'colsample_bylevel':1,
            'subsample':0.5,
            'reg_alpha':0,
            'reg_lambda':0,
            'gamma':0,
            }

random.seed(26)

print("Number of sampling:"+str(num_exp)+"\n")
accuracy_matrix_train = np.zeros(num_exp)


t0 = time.time()
for sub_itr in range(num_exp):

    split_sets_train = np.random.rand(len(labels_train_global)) < data_split
    split_sets_test = np.random.rand(len(labels_test_global)) < data_split

    labels_train = labels_train_global[split_sets_train]
    labels_test = labels_test_global[split_sets_test]
    GEs_train = GEs_train_global[split_sets_train]
    GEs_test = GEs_test_global[split_sets_test]

    print(sub_itr, GEs_train.shape)
    print(sub_itr, GEs_test.shape)

    dtrain_this = xgb.DMatrix(GEs_train, label=labels_train)
    dvalid_this = xgb.DMatrix(GEs_test, label=labels_test)
    bst_this = xgb.train(param, dtrain_this)

    y_pred_this = bst_this.predict(dvalid_this)
    accuracy_matrix_train[sub_itr] = accuracy_score(labels_test, y_pred_this)
    print(sub_itr, accuracy_matrix_train[sub_itr] )
    set_this = set(bst_this.get_score(importance_type='weight').keys())

    set_score = bst_this.get_score(importance_type='weight')
    set_score_combined = {'itr':[sub_itr]*len(set_score), 'genes':list(set_score.keys()), 'weight':list(set_score.values())}
    pd.DataFrame(set_score_combined).to_csv('lamp_genelist_weights.csv', mode='a', header=False)

    output_set_all_combined = {'itr': [sub_itr]*len(set_this), 'genes': list(set_this)}
    pd.DataFrame(output_set_all_combined).to_csv(OUTPUT_GENESET_PATH, mode='a', header=False)
t1 = time.time()

total_run_time = t1-t0

#output_set = pd.read_csv(OUTPUT_GENESET_PATH, header = None)
#genes = pd.DataFrame(output_set.groupby(2)[1].agg(lambda x: x.count()))
#genes = genes[genes[1] == num_exp]
#genelist = list(genes.index)
#pd.DataFrame(genelist).to_csv('pbmc_genelist.csv', header = False)

output_set = pd.read_csv('lamp_genelist_weights.csv', header = None)
genes = pd.DataFrame(output_set.groupby(2).agg({1:'count',3:'mean'}))
genes = genes[genes[1] >= num_exp*0.90]
pd.DataFrame(genes).to_csv('lamp_genelist.csv', header = False)

pd.DataFrame(accuracy_matrix_train).to_csv(OUTPUT_ACCURACY_PATH, mode='a', header=False)

print(total_run_time)
