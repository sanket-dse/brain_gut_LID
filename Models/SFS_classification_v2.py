# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 10:02:14 2024

@author: sanke
"""

## Importing libraries
import pandas as pd
from sklearn.feature_selection import SequentialFeatureSelector, VarianceThreshold 
from sklearn.linear_model import HuberRegressor,LinearRegression,Ridge, LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.metrics import r2_score, accuracy_score, f1_score, \
    roc_auc_score,precision_score,recall_score,balanced_accuracy_score, \
    confusion_matrix, matthews_corrcoef, roc_curve, auc
import statsmodels.api as sm
from statistics import mean, stdev, median
from sklearn import preprocessing
from math import log
from sklearn.preprocessing import PolynomialFeatures,OneHotEncoder
from patsy import dmatrices,dmatrix
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import scipy.stats as stats
from scipy.stats import spearmanr, median_abs_deviation
from statistics import median
import random
from sklearn.model_selection import StratifiedKFold, KFold, RepeatedKFold, RepeatedStratifiedKFold
import pickle
from mpl_toolkits.mplot3d import Axes3D 
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser'
# Importing custom module
import sys
sys.path.append("C:/Users/sanke/Downloads/PD_models")
import important_modules as ip

import scipy.cluster.hierarchy as shc
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from sklearn.metrics import mutual_info_score, adjusted_rand_score
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.manifold import TSNE, MDS
import umap


## Reading the csv file
df1 = pd.read_csv("C:/Users/sanke/Downloads/PD_models/Complete_Features_v3.csv")
df_hrv = pd.read_csv('C:/Users/sanke/Downloads/PD_models/heart_features.csv')

## Getting the pre and post-pre features
idx = [0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27,32,33,
       34,35,40,41,42,43,48,49,50,51,56,57,58,59,64,65,66,67]
idx2 = list(range(72,108))
idx = idx + idx2

#df = pd.concat([df1.iloc[:,0:-13],df1[['PD_Healthy','DyskinesiaScore','Age','Gender']]],axis=1) # Combining the dataframes
df = pd.concat([df1.iloc[:,idx],df1[['PD_Healthy','DyskinesiaScore','Age','Gender']]],axis=1) # Combining the dataframes

# Defining high values for all features
age_high = 80
dys_high = 8 
post_hrv_high = 3
pre_hflf_high = 3
pac_high = 4


# Scaling the variables
df.Age = df.Age / age_high
df.DyskinesiaScore = df.DyskinesiaScore / dys_high

# Getting features and the scores
X = df.iloc[:,0:-4] #For all subjects
#y = np.array(df.PD_Healthy) #1-> PD 0-> Healthy
y = np.array(df.DyskinesiaScore) #Target variable
y = np.nan_to_num(y)
y = y.flatten()
y = np.array([1 if x>0 else 0 for x in y])

# Imputing missing values
#imputer = KNNImputer(n_neighbors=2, weights="uniform")
#X_imputed = imputer.fit_transform(X)
imputer = SimpleImputer(strategy='mean')
X_imputed = imputer.fit_transform(X)


## Log transforming and scaling the PAC features
X_imputed[:,0:36] = np.log(X_imputed[:,0:36])
X_imputed[:,0:36] = X_imputed[:,0:36] / pac_high


# Creating a dataframe to crossverify the feature names that are selected
imputed_df = pd.DataFrame.from_records(X_imputed)
imputed_df.columns = df.columns[0:-4]

#################################### Feature Selection #####################################

X = X_imputed
features = []

'''
splits = 10
repeats = 4
rkf = RepeatedStratifiedKFold(n_splits = splits, n_repeats = repeats, random_state= 1)
for i, (train_index, test_index) in enumerate(rkf.split(X,y)):  
    
    ## Splitting into train and test
    X_train = X[train_index,:]
    X_test = X[test_index,:]
    y_train = y[train_index]
    y_test = y[test_index]
    
    #Creating a dataframe with the imputed features
    imputed_df_train = imputed_df.iloc[train_index,:]
    
    #Sequential feature selection on the training dataset
    features_to_consider,features_list, _ = ip.SFS_logist(X_train, y_train, 
                                                             imputed_df_train, 20)
    
    features.append(features_list)  
    '''
    
## Randomly choosing indices 
repeats = 100
perc = 0.9 

for i in range(repeats) : 
    train_index = random.sample(range(len(y)), round(perc*len(y)))
    train_index.sort()
    
    ## Splitting into train and test
    X_train = X[train_index,:]   
    y_train = y[train_index]  
    
    #Creating a dataframe with the imputed features
    imputed_df_train = imputed_df.iloc[train_index,:]
    
    #Sequential feature selection on the training dataset
    features_to_consider,features_list, R2_list = ip.SFS_logist(X_train, y_train, 
                                                             imputed_df_train, 20)
    
    features.append(features_list)  

'''
    # Getting the features
    features_list = np.array(features_list).tolist()
    for j in range(len(features_list)):
        name = features_list[j]
        features_list[j] = name[0]
        
    #Train the model on these selected features on the training dataset
    model = LogisticRegression().fit(features_to_consider,y_train)
    
    # Predicting on the test dataset
    y_pred_test = model.predict(imputed_df.loc[test_index,features_list])
    y_pred_train = model.predict(imputed_df.loc[train_index,features_list])
    
    #Calculating metrics
    acc_test.append(accuracy_score(y_test,y_pred_test))
    f1_test.append(f1_score(y_test,y_pred_test))
    roc_test.append(roc_auc_score(y_test, model.predict_proba(
        imputed_df.loc[test_index,features_list])[:,-1]))
    prec_test.append(precision_score(y_test,y_pred_test))
    recall_test.append(recall_score(y_test,y_pred_test))
    tn, fp, fn, tp = confusion_matrix(y_test,y_pred_test).ravel()
    spec = tn / (tn+fp)
    spec_test.append(spec)
    mcc_test.append(matthews_corrcoef(y_test,y_pred_test))
    
    acc_train.append(accuracy_score(y_train,y_pred_train))
    f1_train.append(f1_score(y_train,y_pred_train))
    roc_train.append(roc_auc_score(y_train, model.predict_proba(
        imputed_df.loc[train_index,features_list])[:,-1]))
    prec_train.append(precision_score(y_train,y_pred_train))
    recall_train.append(recall_score(y_train,y_pred_train))
    tn, fp, fn, tp = confusion_matrix(y_train,y_pred_train).ravel()
    spec = tn / (tn+fp)
    spec_train.append(spec)
    mcc_train.append(matthews_corrcoef(y_train, y_pred_train))
    '''
######################### Selecting prominent features #####################################

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

prom_features = ['pac_normo_pre_F_T4','pac_normo_TO_T3','pac_tachy_pre_CP_T3',\
    'pac_normo_pre_TO_T2','pac_normo_pre_F_T2']

'''
prom_features = ['pac_brady_pre_F_T1', 'pac_brady_pre_F_T2', 'pac_brady_pre_F_T3', 'pac_brady_pre_F_T4',\
                 'pac_brady_pre_CP_T1', 'pac_brady_pre_CP_T2', 'pac_brady_pre_CP_T3', 'pac_brady_pre_CP_T4',\
                 'pac_brady_pre_TO_T1', 'pac_brady_pre_TO_T2', 'pac_brady_pre_TO_T3', 'pac_brady_pre_TO_T4',\
                 'pac_normo_pre_F_T1', 'pac_normo_pre_F_T2', 'pac_normo_pre_F_T4',]
'''

features_to_consider = np.array(imputed_df.loc[:,prom_features])

########################## Running the model with validation ###############################
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

skf = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 2, random_state= 1)

for i, (train_index, test_index) in enumerate(skf.split(features_to_consider,y)): 
    
    ## Splitting into train and test
    X_train = features_to_consider[train_index,:]
    X_test = features_to_consider[test_index,:]
    y_train = y[train_index]
    y_test = y[test_index]
    
    #Train the model on these selected features on the training dataset
    model = LogisticRegression(C=1e30, max_iter=1000).fit(X_train,y_train)
    cls_coefs.append(model.coef_.transpose().flatten())
   
    # Predicting on the test dataset
    y_pred_test = model.predict(X_test)
    y_pred_train = model.predict(X_train)
   
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

##################### Feature frequency ###########################################

prom_freq = counts[prom_idx]

fig, ax = plt.subplots()
y_pos = np.arange(len(prom_freq))
ax.barh(y_pos,prom_freq, align='center')
ax.set_yticks(y_pos, labels=prom_features)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Frequency of selection')
ax.set_title('Features Importance')

plt.show()

#Sorting the coefs
features_list = np.array(prom_features).tolist()
data = {'Coefs' : prom_freq, 'Features' : features_list }
coef_df = pd.DataFrame(data)
coef_df = coef_df.sort_values(by=['Coefs'], key=abs, ascending=False)

fig, ax = plt.subplots()
y_pos = np.arange(len(coef_df.Coefs))
ax.barh(y_pos,coef_df.Coefs, align='center')
ax.set_yticks(y_pos, labels=np.array(coef_df.Features))
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Frequency of selection')
ax.set_title('Feature importance')

plt.show()

###################### Classification coefficients #####################################

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

fig, ax = plt.subplots()
y_pos = np.arange(len(coef_df.Coefs))
ax.barh(y_pos,coef_df.Coefs, align='center')
ax.errorbar(coef_df.Coefs,y_pos,xerr=coef_df.Error, fmt="o", color="r")
ax.set_yticks(y_pos, labels=np.array(coef_df.Features))
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Logistic Regression coefficients')
ax.set_title('PD v/s Healthy')

plt.show()

######################### Feature coefficients ############################################

# These were the ones that got significant. Comment these out when running a new instance
# of feature selection
#features_list = ['pac_brady_pre_CP_T1','pac_normo_TO_T3','pac_tachy_pre_CP_T3',\
#    'pac_normo_pre_TO_T2','pac_normo_pre_F_T2']
#features_to_consider = np.array(imputed_df.loc[:,features_list])
#new_features_list = ['Intercept','pac_brady_pre_CP_T1','pac_normo_TO_T3','pac_tachy_pre_CP_T3',\
#    'pac_normo_pre_TO_T2','pac_normo_pre_F_T2']

new_features_list = ['Intercept'] + prom_features
features_to_consider = imputed_df.loc[:,prom_features]
reg = sm.Logit(y,sm.add_constant(features_to_consider)).fit(method='bfgs',maxiter = 1000)

print(reg.summary(yname="Dyskinesia Score", xname= new_features_list ))
print(reg.get_margeff(at='overall',method='dydx').summary())  #Marginal effects

#Cohen's F2
for i in features_to_consider.columns : 
    print('Feature :',i)
    ip.cohensf2_logist(features_to_consider, y, i)
    
##################### Adding interactions with significant features ##################

sig_features_list = ['pac_normo_TO_T3','pac_tachy_pre_CP_T3',\
    'pac_normo_pre_F_T2']

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
#new_demo_df = demographics_df.loc[demographics_df['PD_Healthy'] == 1].reset_index() #For PD Patients
#new_demo_df.to_csv('Demographics_heart_features.csv', index=False, header=True)

#new_demo_df = demographics_df.drop('PD_Healthy',axis=1) #For all subjects
features_df = pd.DataFrame.from_records(np.hstack((features_to_consider,y.reshape(-1,1))))
#features_df.columns = ['x' + str(i) for i in range(len(features_to_consider[0]))] + ['y']
features_df.columns = [str(sig_features_list[i]) for i in range(len(features_to_consider[0])) ] + ['y']

merged_df = pd.concat([features_df,demographics_df],axis=1)

#creating interaction terms
#des_form = ['y ~ ('] + ['x'+str(i)+'+' for i in range(len(features_to_consider[0]))] + [')*pre_HRV_T1']
des_form = ['y ~ ('] + [str(sig_features_list[i]) + '+' for i in range(len(features_to_consider[0]))] + [')*pre_HRV']
des_form[-2] = des_form[-2][0:-1]
des_form = ''.join(des_form)
interaction_terms = dmatrices(des_form,data=merged_df, return_type = 'dataframe')
#features = np.delete(np.array(interaction_terms[1]),np.s_[0:len(features_to_consider[0])+1],1)
features = np.array(interaction_terms[1])
y = np.array(interaction_terms[0]).ravel()
feat_list = interaction_terms[1].columns.tolist() #Names of the features

'''
## Adding demographics as features
# PS : Run this for adding all demographic features as covariates
features = pd.concat([merged_df.Age,merged_df.Gender_M,merged_df.pre_HRV,merged_df.pre_HFLF,\
                      merged_df.iloc[:,0:len(features_to_consider[0])]],
                     axis=1)
feat_list = features.columns.tolist()
feat_list = ['Intercept'] + feat_list
'''

## Run this for just adding Age and Gender as covariates to the model
#features = pd.concat([merged_df.Age,merged_df.Gender_M,merged_df.iloc[:,0:len(features_to_consider[0])]],axis=1)
#feat_list = features.columns.tolist()


## Fitting the model

reg = sm.Logit(y,features).fit(method='bfgs',maxiter = 1000)
features_list = np.array(feat_list).tolist()
print(reg.summary(yname="Dyskinesia Score", xname= features_list))

#Cohen's F2
temp_df = pd.DataFrame(features, columns=features_list) #Creating a temporary dataframe
for i in temp_df.columns : 
    print('Feature :',i)
    ip.cohensf2_logist(temp_df, y, i)

####################### Using marginal effects command #################################

heart_idx = 4 #The index of heart measure   
pvals = reg.pvalues
mean_arr = np.array([interaction_terms[1].mean()])
reg_params = reg.params

### dydx w.r.t HRV for various values of PAC 

for i in range(heart_idx+1,len(pvals)) : 
    
    if pvals[i] < 0.05 :
        
        #Getting indices of the pac features to be analyzed

        idx = i #Index of the interaction of heart measure with pac value
        colname = interaction_terms[1].columns[idx].split(':')[0] #name of the PAC value
        pac_idx = interaction_terms[1].columns.values.tolist().index(colname)        
         
        
        #Creating evenly spaced out values for the PAC feature and heart meas
        pac_values = np.linspace(min(features[:,pac_idx]), max(features[:,pac_idx]),10) 
        heart_values = np.linspace(min(features[:,heart_idx]), max(features[:,heart_idx]),10)
        
        #Creating exog values matrix and calculating marginal effects at different values of 
        # PAC feature
        exog_values = np.tile(mean_arr,[10,1])
        marginal_effects = np.tile(0.0,[10,])
        standard_error = np.tile(0.0,[10,])        
        
        for j in range(10) : 
            for k in range(10) : 
                exog_values[k,pac_idx] = pac_values[j] #Replacing the PAC features
                exog_values[k,i] = heart_values[k] * pac_values[j] #Replacing the interaction of PAC and heart meas
                exog_values[k,heart_idx] = heart_values[k]        
        
            
            # Calculate marginal effects 
            marg_eff = reg.get_margeff(at='overall',atexog = exog_values)

            # Get the marginal effect and standard error for pre_HRV
            marginal_effects[j] = marg_eff.margeff[heart_idx-1]
            standard_error[j] = marg_eff.margeff_se[heart_idx-1]
        
        ## Plotting the marginal effects 
       
        plt.figure(figsize=(15, 15))
        
        plt.plot(pac_values, marginal_effects)
        plt.errorbar(pac_values, marginal_effects, yerr = standard_error, \
                     fmt='o', linestyle='-', capsize=4)
        plt.axhline(y = 0, color = 'b', linestyle = '--') 
        plt.title(''.join(['Marginal Effects of pre-medication HRV measure on Dyskinesia at different levels of ',colname]), \
                  fontsize = 10)
        plt.xlabel(colname, fontsize=10)
        plt.ylabel('Marginal Effect of pre-medication HRV measure on Probability of Dyskinesia', \
                   fontsize=10)
        plt.tick_params(axis='both', which='major', labelsize=30)
        plt.grid(True)

        # Show plot
        plt.show()
        

### dydx w.r.t PAC for various values of HRV 

for i in range(heart_idx+1,len(pvals)) : 
    
    if pvals[i] < 0.05 :
        
        #Getting indices of the pac features to be analyzed

        idx = i #Index of the interaction of heart measure with pac value
        colname = interaction_terms[1].columns[idx].split(':')[0] #name of the PAC value
        pac_idx = interaction_terms[1].columns.values.tolist().index(colname)         
         
        
        #Creating evenly spaced out values for the PAC feature and heart meas
        pac_values = np.linspace(min(features[:,pac_idx]), max(features[:,pac_idx]),10) 
        heart_values = np.linspace(min(features[:,heart_idx]), max(features[:,heart_idx]),10)
        
        #Creating exog values matrix and calculating marginal effects at different values of 
        #pre hrv
        exog_values = np.tile(mean_arr,[10,1])
        marginal_effects = np.tile(0.0,[10,])
        standard_error = np.tile(0.0,[10,])        
        
        for j in range(10) : 
            for k in range(10) : 
                exog_values[k,pac_idx] = pac_values[k] #Replacing the PAC features
                exog_values[k,i] = heart_values[j] * pac_values[k] #Replacing the interaction of PAC and heart meas
                exog_values[k,heart_idx] = heart_values[j]        
        
            
            # Calculate marginal effects 
            marg_eff = reg.get_margeff(at='overall',atexog = exog_values)

            # Get the marginal effect and standard error for pre_HRV
            marginal_effects[j] = marg_eff.margeff[pac_idx-1]
            standard_error[j] = marg_eff.margeff_se[pac_idx-1]
        
        ## Plotting the marginal effects 
       
        plt.figure(figsize=(15, 15))
        
        plt.plot(heart_values, marginal_effects)
        plt.errorbar(heart_values, marginal_effects, yerr = standard_error, \
                     fmt='o', linestyle='-', capsize=4)
        plt.axhline(y = 0, color = 'b', linestyle = '--') 
        # Plot settings
        plt.title(''.join(['Marginal Effects of ', colname,' on Dyskinesia at different levels of pre-medication HRV measures']), \
                  fontsize=10)
        plt.xlabel('Pre-HRV', fontsize=10)
        plt.ylabel(''.join(['Marginal Effect of ',colname,' on Probability of Dyskinesia']), \
                   fontsize=10)
        plt.tick_params(axis='both', which='major', labelsize=30)
        plt.grid(True)

        # Show plot
        plt.show()    
        

'''
###################### Plotting the marginal effects #############################

heart_idx = 4 #The index of heart measure   
pvals = reg.pvalues
mean_arr = np.array(interaction_terms[1].mean())
reg_params = reg.params
mean_heart_meas = mean(interaction_terms[1].loc[:,'pre_HRV'])
sd_heart_meas = stdev(interaction_terms[1].loc[:,'pre_HRV'])

# Define logistic function
def logistic(z):
     return 1 / (1 + np.exp(-z))

# Define the marginal effect function
def marginal_effect(X1, X2, beta_1, beta_3, probabilities):
     return probabilities * (1 - probabilities) * (beta_1 + beta_3 * X2)

### Varying levels of pre HRV measure

for i in range(heart_idx+1,len(pvals)) : 
    
    if pvals[i] < 0.05 :
        
        #Getting indices of the pac features to be analyzed

        idx = i #Index of the interaction of heart measure with pac value
        colname = interaction_terms[1].columns[idx].split(':')[0] #name of the PAC value
        pac_idx = interaction_terms[1].columns.values.tolist().index(colname) 
        indices = np.array(range(0,len(interaction_terms[1].columns.values)))
        const_idx = np.delete(indices, [idx,pac_idx])   
        
        #Setting up the equation to probe interaction
        
        intercept = 0
        for i in const_idx :
            intercept += reg_params[i]*mean_arr[i]
    
        x = np.linspace(min(interaction_terms[1].iloc[:,pac_idx]), 
                        max(interaction_terms[1].iloc[:,pac_idx]), 100)
        
    
        # We will calculate the marginal effects at different calorie intake levels for a fixed value of exercise (X1)
        fixed_heart_meas = [(mean_heart_meas-sd_heart_meas), mean_heart_meas,\
                                (mean_heart_meas+sd_heart_meas)]  
        
        # Create a plot
        plt.figure(figsize=(10, 6))

        # Loop over exercise hours and calculate marginal effects for different calorie intake levels
        for heart_meas in fixed_heart_meas:
            log_odds = intercept + reg_params[heart_idx]*heart_meas + (reg_params[idx]*heart_meas + reg_params[pac_idx])*x
            probabilities = logistic(log_odds)
            marginal_effects = marginal_effect(heart_meas, x, reg_params[heart_idx],\
                                               reg_params[i], probabilities)
            plt.plot(x, marginal_effects, label=f'HRV : {heart_meas}')
            
        

        # Plot settings
        plt.title(''.join(['Marginal Effects of pre-medication HRV measure on Dyskinesia at different levels of ',colname]))
        plt.xlabel(colname)
        plt.ylabel('Marginal Effect of pre-medication HRV measure on Probability of Dyskinesia')
        plt.legend(title='pre HRV')
        plt.grid(True)

        # Show plot
        plt.show()   

### Varying levels of PAC measures

for i in range(heart_idx+1,len(pvals)) : 
    
    if pvals[i] < 0.05 :
        
        #Getting indices of the pac features to be analyzed

        idx = i #Index of the interaction of heart measure with pac value
        colname = interaction_terms[1].columns[idx].split(':')[0] #name of the PAC value
        pac_idx = interaction_terms[1].columns.values.tolist().index(colname) 
        indices = np.array(range(0,len(interaction_terms[1].columns.values)))
        const_idx = np.delete(indices, [idx,heart_idx])   
        
        #Setting up the equation to probe interaction
        
        intercept = 0
        for i in const_idx :
            intercept += reg_params[i]*mean_arr[i]
    
        x = np.linspace(min(interaction_terms[1].iloc[:,heart_idx]), 
                        max(interaction_terms[1].iloc[:,heart_idx]), 100)

        mean_pac = mean(interaction_terms[1].iloc[:,pac_idx])
        sd_pac = stdev(interaction_terms[1].iloc[:,pac_idx])
    
        # We will calculate the marginal effects at different calorie intake levels for a fixed value of exercise (X1)
        fixed_pac = [(mean_pac-sd_pac), mean_pac,\
                                (mean_pac+sd_pac)]  
        
        # Create a plot
        plt.figure(figsize=(10, 6))

        # Loop over exercise hours and calculate marginal effects for different calorie intake levels
        for pac in fixed_pac:
            log_odds = intercept + reg_params[pac_idx]*pac + (reg_params[idx]*pac + reg_params[pac_idx])*x
            probabilities = logistic(log_odds)
            marginal_effects = marginal_effect(pac, x, reg_params[heart_idx],\
                                               reg_params[i], probabilities)
            plt.plot(x, marginal_effects, label=f'PAC : {pac}')
        

        # Plot settings
        plt.title(''.join(['Marginal Effects of ', colname,' on Dyskinesia at different levels of pre-medication HRV measures']))
        plt.xlabel('Pre-HRV')
        plt.ylabel(''.join(['Marginal Effect of ',colname,' on Probability of Dyskinesia']))
        plt.legend(title=colname)
        plt.grid(True)

        # Show plot
        plt.show()   

        
        
################# Visualizing interactions ######################################

# find indices with less than median
idx = imputed_df.index[imputed_df['pac_normo_TO_T3'] < median(imputed_df['pac_normo_TO_T3'])].tolist()
idx2 = imputed_df.index[imputed_df['pac_normo_TO_T3'] > median(imputed_df['pac_normo_TO_T3'])].tolist()

HRV = merged_df.loc[idx,'pre_HRV'].tolist()
probs = reg.predict(features)
selected_probs = probs[idx]
HRV2 = merged_df.loc[idx2,'pre_HRV'].tolist()
selected_probs2 = probs[idx2]

plt.scatter(HRV,selected_probs)
plt.scatter(HRV2, selected_probs2)
plt.show()

#####################

probs = reg.predict(features)
idx = [index for index,value in enumerate(probs) if value > 0.5]
idx2 = [index for index,value in enumerate(probs) if value < 0.5]

HRV = merged_df.loc[idx2,'pre_HRV'].tolist()
pac = imputed_df.loc[idx2,'pac_normo_TO_T3']

plt.scatter(HRV,pac)
plt.show()

#############################

# find indices with less than me
idx = merged_df.index[merged_df['pre_HRV'] < median(merged_df['pre_HRV'])].tolist()
idx2 = merged_df.index[merged_df['pre_HRV'] > median(merged_df['pre_HRV'])].tolist()

pac = imputed_df.loc[idx2,'pac_normo_TO_T3']
probs = reg.predict(features)
selected_probs = probs[idx2]

plt.scatter(pac,selected_probs)
plt.show()

##########################

HRV = merged_df['pre_HRV'].tolist()
probs = reg.predict(features)
#pac = merged_df['pac_normo_TO_T3'].tolist()
pac = merged_df['pac_tachy_pre_CP_T3'].tolist()


plt.scatter(HRV,probs,c=pac,cmap='Greens')
plt.colorbar()
plt.show()

plt.scatter(pac,probs,c=HRV,cmap='Greens')
plt.colorbar()
plt.show()
'''

############################ Clustering ##############################################

########### Hierarchial clustering
plt.figure()  
plt.title("Dendrograms")  
dend = shc.dendrogram(shc.linkage(features_to_consider, method='ward'))

########## Doing PCA 
pca = PCA()
x_new = pca.fit_transform(features_to_consider)

def myplot_2d(score,coeff,labels=None):
    xs = score[:,0]
    ys = score[:,1]
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    plt.scatter(xs * scalex,ys * scaley, c = y)
    for i in range(n):
        plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
        if labels is None:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+str(i+1), color = 'g', ha = 'center', va = 'center')
        else:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center')
    
    
def myplot_3d(score,coeff,labels=None):
    xs = score[:,0]
    ys = score[:,1]
    zs = score[:,2]
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    scalez = 1.0/(zs.max() - zs.min())
    
    ax = go.Figure(data=[go.Scatter3d(
    x=xs * scalex, y=ys * scaley, z=zs * scalez,
    mode='markers',
    marker=dict(size=8, color=y, colorscale='Viridis',   
        opacity=0.8))])
    
        
    for i in range(n):
        ax.add_trace(go.Scatter3d(x = [0,coeff[i,0]], y = [0, coeff[i,1]],z = [0, coeff[i,2]],
                                  mode='lines',line_color='black', 
                                  name = "Var"+str(i+1)))
        
    ax.update_layout(scene = dict(
                    xaxis_title='PC 1',
                    yaxis_title='PC 2',
                    zaxis_title='PC 3'))
    
    ax.show()
    
#Call the function. Use only the 2 PCs.
myplot_2d(x_new[:,0:2],np.transpose(pca.components_[0:2, :]))
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.xlabel("PC{}".format(1))
plt.ylabel("PC{}".format(2))
plt.grid()
plt.show()

## Using embeddings 
embedding = []
for i in range(len(mean_coefs)):
    if i==0:
        embedding = features_to_consider[:,i]*mean_coefs[i]
    else:     
        embedding =  np.c_[embedding,(features_to_consider[:,i]*mean_coefs[i])]


########### USing t-SNE for visualization
tsne = TSNE(n_components=2, verbose=1, perplexity=15, n_iter=5000)
X_tsne = tsne.fit_transform(embedding)

# Plotting the clustering
fig = plt.figure(figsize=(15, 15))
#ax = fig.add_subplot(projection='3d')
  
for i in range(2):
    plt.scatter(X_tsne[y == i , 0] , X_tsne[y == i , 1] ,label=i)

plt.xlabel('Dim 1')
plt.ylabel('Dim 2')
plt.legend()
plt.show()

##################### Using umap for visualization
X_umap = umap.UMAP(n_neighbors=4, min_dist=0.4,
                   metric = 'correlation').fit_transform(embedding)

# Plotting the clustering
fig = plt.figure(figsize=(15, 15))
#ax = fig.add_subplot(projection='3d')
  
for i in range(2):
    plt.scatter(X_umap[y == i , 0] , X_umap[y == i , 1] ,label=i)

plt.xlabel('Dim 1')
plt.ylabel('Dim 2')
plt.legend()
plt.show()

###################### USing MDS for visualization
mds = MDS(n_components=2,eps = 0.001)
X_mds = mds.fit_transform(embedding)

# Plotting the clustering
fig = plt.figure(figsize=(15, 15))
#ax = fig.add_subplot(projection='3d')
  
for i in range(2):
    plt.scatter(X_mds[y == i , 0] , X_mds[y == i , 1] ,label=i)

plt.xlabel('Dim 1')
plt.ylabel('Dim 2')
plt.legend()
plt.show()



