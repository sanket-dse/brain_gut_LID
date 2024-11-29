# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 14:56:26 2024

@author: sanke
"""

## Importing libraries
import pandas as pd
from sklearn.linear_model import HuberRegressor,LinearRegression,Ridge,LogisticRegression
import numpy as np
from sklearn.impute import KNNImputer,SimpleImputer
from sklearn.metrics import r2_score
import statsmodels.api as sm
from statistics import mean,stdev
from sklearn import preprocessing
from math import log
from sklearn.preprocessing import PolynomialFeatures,OneHotEncoder,KBinsDiscretizer
from patsy import dmatrices,dmatrix
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import scipy.stats as stats
import random
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score, accuracy_score, f1_score, \
    roc_auc_score,precision_score,recall_score,balanced_accuracy_score, \
    confusion_matrix, matthews_corrcoef, roc_curve, auc
from sklearn.model_selection import StratifiedKFold, KFold, RepeatedKFold, RepeatedStratifiedKFold
# Importing custom module
import sys
sys.path.append("C:/Users/sanke/Downloads/PD_models")
import important_modules as ip
 
## Reading the csv file
#df1 = pd.read_csv("C:/Users/sanke/Downloads/PD_models/Complete_Features.csv")
#df2 = pd.read_csv('C:/Users/sanke/Downloads/PD_models/post-pre_Complete_Features.csv')
df_hrv = pd.read_csv('C:/Users/sanke/Downloads/PD_models/heart_features.csv')
df1 = pd.read_csv("C:/Users/sanke/Downloads/PD_models/Complete_Features_v3.csv")

## Getting the pre and post-pre features
idx = [0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27,32,33,
       34,35,40,41,42,43,48,49,50,51,56,57,58,59,64,65,66,67]
idx2 = list(range(72,108))
idx = idx + idx2

df = pd.concat([df1.iloc[:,idx],df1[['PD_Healthy','DyskinesiaScore','Age','Gender']]],axis=1) # Combining the dataframes

## Normalizing the dyskinesia score
dys_high = 8 
age_high = 80
post_hrv_high = 3
df.DyskinesiaScore = df.DyskinesiaScore / dys_high
df.Age = df.Age / age_high

# Getting the features and the target variable
#X_dys = df.loc[df['PD_Healthy']==1,df.columns[0:-4]] #Only PD patients & PAC features
#X_nondys = df.loc[df['PD_Healthy']==1,df.columns[0:-4]] 
X = df.loc[df['PD_Healthy'] == 1,df.columns[0:-4]] 
y = np.array(df.loc[df['PD_Healthy'] == 1,['DyskinesiaScore']])
y = np.nan_to_num(y)
y = y.flatten()


#Imputing with the mean values
imputer = SimpleImputer(strategy='mean')
#X_dys_imp = imputer.fit_transform(X_dys)
#X_nondys_imp = imputer.fit_transform(X_nondys)
X_imputed = imputer.fit_transform(X)

## Log transforming and scaling the PAC features
pac_high = 4
X_imputed[:,0:36] = np.log(X_imputed[:,0:36])
X_imputed[:,0:36] = X_imputed[:,0:36] / pac_high

# Creating a dataframe to crossverify the feature names that are selected
imputed_df = pd.DataFrame.from_records(X_imputed)
imputed_df.columns = df.columns[0:-4]

############################# Sequential feature loop ####################################

X = X_imputed
'''
## Discretize the target variable
binner = KBinsDiscretizer(n_bins = 3, encode = 'ordinal', strategy = 'quantile')
y_labels = binner.fit_transform(y.reshape(-1,1))

features = []
splits = 8
repeats = 4
rkf = RepeatedStratifiedKFold(n_splits = splits, n_repeats = repeats, random_state= 1)
#rkf = RepeatedKFold(n_splits = 8, n_repeats = 3, random_state = 1)

for i, (train_index, test_index) in enumerate(rkf.split(X,y_labels)):
    ## Splitting into train dataset
    X_train = X[train_index,:]    
    y_train = y[train_index]
    imputed_df_train = imputed_df.iloc[train_index,:]
    # Sequential feature selection
    features_to_consider,features_list,R2_list = ip.SFS_alt(X_train,y_train,imputed_df_train,train_index)
    
    features.append(features_list)
'''

## Randomly choosing indices 
features = []
repeats = 100
perc = 0.9 

for i in range(repeats) : 
    train_index = random.sample(range(len(y)), round(perc*len(y)))
    train_index.sort()
    X_train = X[train_index,:]    
    y_train = y[train_index]
    imputed_df_train = imputed_df.iloc[train_index,:]
    # Sequential feature selection
    features_to_consider,features_list,R2_list = ip.SFS_alt(X_train,y_train,imputed_df_train,train_index)

    features.append(features_list)
    
## Checking the feature frequency
features_list = []
for i in range(len(features)) :
    
    features_temp = features[i]
    features_temp = np.array(features_temp).tolist()

    for j in range(len(features_temp)):
        name = features_temp[j]
        features_temp[j] = name[0]
        
    for k in range(len(features_temp)) :
        features_list.append(features_temp[k])


uniq, counts = np.unique(features_list, return_counts=True) 

## Selecting the prominent features
thresh_count = round((1/3)*(repeats))
prom_idx = [i for i in range(len(counts)) if counts[i] > thresh_count]
prom_features = uniq[prom_idx].tolist()
features_to_consider = np.array(imputed_df.loc[:,prom_features])

#################### Creating the model with the prominent features ########################

prom_features = ['pac_normo_pre_F_T2','pac_normo_F_T3','pac_tachy_TO_T4','pac_brady_pre_CP_T4',\
            'pac_brady_pre_F_T4']

features_to_consider = np.array(imputed_df.loc[:,prom_features])
    
#Using statsmodels
model = sm.RLM(y,features_to_consider,M=sm.robust.norms.HuberT())
reg = model.fit(scale_est = sm.robust.scale.HuberScale())
reg_coefs = reg.params
y_pred = reg.predict(features_to_consider)

#print(reg.summary(yname="Dyskinesia Score", xname=["var_%d" % i for i in range(len(reg.params))] ))
print(reg.summary(yname="Dyskinesia Score", xname= prom_features ))
print('Adjusted R2 = ', ip.adjusted_R2(reg,features_to_consider,y))
print('RMSE = ',mean_squared_error(y, y_pred, squared=False))
print('p-value for entire model = ', ip.f_test_model(reg,features_to_consider,y))
print('R2 value = ', ip.R2(reg,features_to_consider,y))

## Saving the features in a dataframe and writing them to a csv file
features_df = pd.DataFrame.from_records(np.hstack((features_to_consider,y.reshape(-1,1))))
features_df.columns = [str(prom_features[i]) for i in range(len(features_to_consider[0])) ] + ['y']

features_df.to_csv('Selected_features_v3.csv', index=False, header=True)

################ Calculating cohen's F2 for the significant features ###############

# PS - the significant features are pac_normo_pre_F_T2, pac_normo_F_T3

features_df = imputed_df.loc[:,prom_features]
print('Cohen\'s F2 for pac_normo_pre_F_T2 : ', ip.cohensf2(features_df, y, \
                                                           'pac_normo_pre_F_T2'))
print('Cohen\'s F2 for pac_normo_F_T3 : ', ip.cohensf2(features_df, y, \
                                                           'pac_normo_F_T3'))
    
################ Plotting the feature importance ############################

#Sorting the coefs
data = {'Coefs' : reg_coefs, 'Features' : prom_features }
coef_df = pd.DataFrame(data)
coef_df = coef_df.sort_values(by=['Coefs'], key=abs, ascending=False)

fig, ax = plt.subplots()
y_pos = np.arange(len(coef_df.Coefs))
ax.barh(y_pos,coef_df.Coefs, align='center')
ax.set_yticks(y_pos, labels=np.array(coef_df.Features))
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Regression coefficients')
ax.set_title('Dyskinesia score')

plt.show()

################ Permutation feature importance ############################################

model = HuberRegressor()
reg = model.fit(features_to_consider,y)
r = permutation_importance(reg, features_to_consider, y,
                           n_repeats=30, scoring = 'r2',
                           random_state=0)

ax = sns.boxplot(np.transpose(r.importances),orient='h',linecolor='skyblue',color='white',
            medianprops = dict(color='green'), showfliers=False)
ax.axvline(x = 0,color='black',linestyle='dashed')
ax.set_yticklabels(features_list)
plt.xlabel('Permutation feature importance')

##################### Feature frequency ###########################################

prom_freq = counts[prom_idx]
#Sorting the coefs
data = {'Coefs' : prom_freq, 'Features' : prom_features }
coef_df = pd.DataFrame(data)
coef_df = coef_df.sort_values(by=['Coefs'], key=abs, ascending=False)

fig, ax = plt.subplots()
y_pos = np.arange(len(coef_df.Coefs))
ax.barh(y_pos,coef_df.Coefs, align='center')
ax.set_yticks(y_pos, labels=coef_df.Features)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Frequency of selection')
ax.set_title('Features Importance')

plt.show()

################### Feature interaction with heart features and demographics #############

sig_features_list = ['pac_normo_pre_F_T2','pac_normo_F_T3']

features_to_consider = np.array(imputed_df.loc[:,sig_features_list])

#Imputing the HRV measures
hrv_imputed = imputer.fit_transform(df_hrv)

# Sqrt transforming and scaling the variables and accessing only the pre med variables
hrv_imputed = np.sqrt(hrv_imputed) 
hrv_imputed = hrv_imputed / post_hrv_high
hrv_imputed = hrv_imputed[:,[0,1,2,3,8,9,10,11]]

 
#One hot encoding the Gender variable
encoder = OneHotEncoder(sparse_output=False)
gender = encoder.fit_transform(np.array(df['Gender']).reshape(-1,1))
gender = np.delete(gender,0,1)


#Scaling the variables
#hrv_scaled = standard_scaler.fit_transform(hrv_imputed)
#age_scaled = standard_scaler.fit_transform(np.array(df.Age).reshape(-1,1))
hrv_scaled = hrv_imputed
age_scaled = np.array(df.Age).reshape(-1,1)

# Finding the task averaged HRV and LFHF ratio
pre_hrv = []
pre_hflf = []

for i in range(len(hrv_scaled)) : 
    pre_hrv.append(mean(hrv_scaled[:,[0,1,2,3]][i]))
    pre_hflf.append(mean(hrv_scaled[:,[4,5,6,7]][i]))

#Combining all the demographic variables
demographics = np.hstack((np.array(pre_hrv).reshape(-1,1),np.array(pre_hflf).reshape(-1,1),\
                          age_scaled,gender, np.array(df.PD_Healthy).reshape(-1,1)))

#Creating a dataframe with the required features
demographics_df = pd.DataFrame.from_records(demographics)
demographics_df.columns = np.hstack(('pre_HRV','pre_HFLF','Age','Gender_M','PD_Healthy'))
new_demo_df = demographics_df.loc[demographics_df['PD_Healthy'] == 1].reset_index() #For PD Patients
#new_demo_df.to_csv('Demographics_heart_features.csv', index=False, header=True)

#new_demo_df = demographics_df.drop('PD_Healthy',axis=1) #For all subjects
features_df = pd.DataFrame.from_records(np.hstack((features_to_consider,y.reshape(-1,1))))
#features_df.columns = ['x' + str(i) for i in range(len(features_to_consider[0]))] + ['y']
features_df.columns = [str(sig_features_list[i]) for i in range(len(features_to_consider[0])) ] + ['y']

merged_df = pd.concat([features_df,new_demo_df],axis=1)

#creating interaction terms
#des_form = ['y ~ ('] + ['x'+str(i)+'+' for i in range(len(features_to_consider[0]))] + [')*pre_HRV_T1']
des_form = ['y ~ ('] + [str(sig_features_list[i]) + '+' for i in range(len(features_to_consider[0]))] + [')*Age']
des_form[-2] = des_form[-2][0:-1]
des_form = ''.join(des_form)
interaction_terms = dmatrices(des_form,data=merged_df, return_type = 'dataframe')
#features = np.delete(np.array(interaction_terms[1]),np.s_[0:len(features_to_consider[0])+1],1)
features = np.array(interaction_terms[1])
y = np.array(interaction_terms[0]).ravel()
feat_list = interaction_terms[1].columns.tolist() #Names of the features


## Adding demographics as features
# PS : Run this for adding all demographic features as covariates
features = pd.concat([merged_df.Age,merged_df.Gender_M,merged_df.pre_HRV,merged_df.pre_HFLF,\
                      merged_df.iloc[:,0:len(features_to_consider[0])]],
                     axis=1)
feat_list = features.columns.tolist()
features = np.array(features)

#Using statsmodels
model = sm.RLM(y,features,M=sm.robust.norms.HuberT())
reg = model.fit(scale_est = sm.robust.scale.HuberScale())
reg_coefs = reg.params
y_pred = reg.predict(features)

#print(reg.summary(yname="Dyskinesia Score", xname=["var_%d" % i for i in range(len(reg.params))] ))
print(reg.summary(yname="Dyskinesia Score", xname= feat_list ))
print('Adjusted R2 = ', ip.adjusted_R2(reg,features,y))
print('RMSE = ',mean_squared_error(y, y_pred, squared=False))
print('p-value for entire model = ', ip.f_test_model(reg,features,y))
print('R2 value = ', ip.R2(reg,features,y))

#################### Using these features to do logistic regression #################

X = df.loc[:,df.columns[0:-4]] 
#y = np.array(df.PD_Healthy)
y = np.array(df.DyskinesiaScore) #Target variable
y = np.nan_to_num(y)
y = y.flatten()
y = np.array([1 if x>0 else 0 for x in y])#

#Imputing with the mean values
imputer = SimpleImputer(strategy='mean')
#X_dys_imp = imputer.fit_transform(X_dys)
#X_nondys_imp = imputer.fit_transform(X_nondys)
X_imputed = imputer.fit_transform(X)

## Log transforming and scaling the PAC features
pac_high = 4
X_imputed[:,0:36] = np.log(X_imputed[:,0:36])
X_imputed[:,0:36] = X_imputed[:,0:36] / pac_high

# Creating a dataframe to crossverify the feature names that are selected
imputed_df = pd.DataFrame.from_records(X_imputed)
imputed_df.columns = df.columns[0:-4]

features_to_consider = np.array(imputed_df.loc[:,prom_features])

## Giving cross validation score 
skf = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 20, random_state= 1)
cls_coefs = []
f1_train = []
roc_train = []
acc_train = []
prec_train = []
recall_train = []
spec_train = []
mcc_train = []

f1_test = []
roc_test = []
acc_test = []
prec_test = []
recall_test = []
spec_test = []
mcc_test = []
r_val = []
r_train = []

for i, (train_index, test_index) in enumerate(skf.split(features_to_consider,y)):
    ## Splitting into train and test
    X_train = features_to_consider[train_index,:]
    X_test = features_to_consider[test_index,:]
    y_train = y[train_index]
    y_test = y[test_index]
    
    #Fitting the model on training data 
    try :
        #model = sm.Logit(y_train,X_train).fit(method='bfgs',maxiter=100)
        model = LogisticRegression(C=1e30).fit(X_train,y_train)
        #model = DecisionTreeClassifier().fit(X_train,y_train)
        #model = RandomForestClassifier(max_depth = 4, n_jobs = -1,random_state =1).fit(X_train,y_train)
    except :
        continue
    

    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)
    
    cls_coefs.append(model.coef_.transpose().flatten())
    
    '''
    # Permutation importance
    r = permutation_importance(model, X_test, y_test,
                               n_repeats=30, scoring = 'accuracy',
                               random_state=0)
    if i==0 :
        r_val = r.importances
    else :
        r_val = np.hstack((r_val,r.importances))
        
    # Permutation importance
    r = permutation_importance(model, X_train, y_train,
                               n_repeats=30, scoring = 'accuracy',
                               random_state=0)
    if i==0 :
        r_train = r.importances
    else :
        r_train = np.hstack((r_train,r.importances))
    '''
    
    #Calculating metrics
    acc_test.append(accuracy_score(y_test,y_pred_test))
    f1_test.append(f1_score(y_test,y_pred_test))
    roc_test.append(roc_auc_score(y_test, model.predict_proba(X_test)[:,-1]))
    prec_test.append(precision_score(y_test,y_pred_test))
    recall_test.append(recall_score(y_test,y_pred_test))
    tn, fp, fn, tp = confusion_matrix(y_test,y_pred_test).ravel()
    spec = tn / (tn+fp)
    spec_test.append(spec)
    mcc_test.append(matthews_corrcoef(y_test,y_pred_test))
    
    acc_train.append(accuracy_score(y_train,y_pred_train))
    f1_train.append(f1_score(y_train,y_pred_train))
    roc_train.append(roc_auc_score(y_train, model.predict_proba(X_train)[:,-1]))
    prec_train.append(precision_score(y_train,y_pred_train))
    recall_train.append(recall_score(y_train,y_pred_train))
    tn, fp, fn, tp = confusion_matrix(y_train,y_pred_train).ravel()
    spec = tn / (tn+fp)
    spec_train.append(spec)
    mcc_train.append(matthews_corrcoef(y_train, y_pred_train))
    
print('Train Accuracy : ',mean(acc_train))
print('Test Accuracy : ',mean(acc_test))
print('Train Recall : ', mean(recall_train))
print('Test Recall : ', mean(recall_test))
print('Train Precision : ',mean(prec_train))
print('Test Precision : ',mean(prec_test))
print('Train F1 score : ',mean(f1_train))
print('Test F1 score : ',mean(f1_test))
print('Train ROC AUC score : ',mean(roc_train))
print('Test ROC AUC score : ',mean(roc_test))
print('Train specifictiy score : ',mean(spec_train))
print('Test specificity score : ',mean(spec_test))
print('Train MCC score : ',mean(mcc_train))
print('Test MCC score : ',mean(mcc_test))


################ Plotting the classification coefficients ############################
'''
model = LogisticRegression(C=1e30).fit(features_to_consider,y)
cls_coefs = model.coef_.transpose().flatten()

fig, ax = plt.subplots()
y_pos = np.arange(len(cls_coefs))
ax.barh(y_pos,cls_coefs, align='center')
ax.set_yticks(y_pos, labels=np.array(prom_features).tolist())
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Feature importance')
ax.set_title('Dyskinesia score')

plt.show()


#Sorting the coefs
features_list = np.array(prom_features).tolist()
data = {'Coefs' : cls_coefs, 'Features' : features_list }
coef_df = pd.DataFrame(data)
coef_df = coef_df.sort_values(by=['Coefs'], key=abs, ascending=False)
'''
mean_coefs = []
std_coefs = []
for i in range(len(cls_coefs[0])) :
    temp= []
    for j in range(len(cls_coefs)):
        temp.append(cls_coefs[j][i])
        
    mean_coefs.append(mean(temp))
    std_coefs.append(stdev(temp))

#Sorting the coefs
features_list = np.array(prom_features).tolist()
data = {'Coefs' : mean_coefs, 'Features' : features_list , 'Error' : std_coefs}
coef_df = pd.DataFrame(data)
coef_df = coef_df.sort_values(by=['Coefs'], key=abs, ascending=False)

# Plotting the coefficients
fig, ax = plt.subplots()
y_pos = np.arange(len(coef_df.Coefs))
ax.barh(y_pos,coef_df.Coefs, align='center')
ax.set_yticks(y_pos, labels=np.array(coef_df.Features))
ax.errorbar(coef_df.Coefs,y_pos,xerr=coef_df.Error, fmt="o", color="r")
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Classification coefficients')
ax.set_title('Dyskinesia score')

plt.show()
