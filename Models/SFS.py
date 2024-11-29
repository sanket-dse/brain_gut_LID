# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:40:41 2024

@author: sanke
"""
## Importing libraries
import pandas as pd
from sklearn.feature_selection import SequentialFeatureSelector 
from sklearn.linear_model import HuberRegressor,LinearRegression,Ridge
import numpy as np
from sklearn.impute import KNNImputer,SimpleImputer
from sklearn.metrics import r2_score
import statsmodels.api as sm
from statistics import mean,stdev
from sklearn import preprocessing
from math import log
from sklearn.preprocessing import PolynomialFeatures,OneHotEncoder
from patsy import dmatrices,dmatrix
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import scipy.stats as stats
from scipy.stats import spearmanr
from sklearn.inspection import permutation_importance
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

#Slicing the dataframes
#df1_copy = df1.iloc[:,0:72]
#df1_copy2 = df1.iloc[:,96:-13]
#df2_copy = df2.iloc[:,0:36]
#df2_copy2 = df2.iloc[:,48:-13]
#df1_new = pd.concat([df1_copy,df1_copy2],axis=1)
#f2_new = pd.concat([df2_copy,df2_copy2],axis = 1)
#df = pd.concat([df1_new,df2_new,df1[['PD_Healthy','DyskinesiaScore','Age','Gender']]],axis=1) # Combining the dataframes

#df = pd.concat([df1.iloc[:,0:72],df2],axis=1) # Combining the dataframes\
#df = pd.concat([df1.iloc[:,0:108],df1[['PD_Healthy','DyskinesiaScore','Age','Gender']]],axis=1) # Combining the dataframes
idx = [0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27,32,33,
       34,35,40,41,42,43,48,49,50,51,56,57,58,59,64,65,66,67]
idx2 = list(range(72,108))
idx = idx +idx2

df = pd.concat([df1.iloc[:,idx],df1[['PD_Healthy','DyskinesiaScore','Age','Gender']]],axis=1) # Combining the dataframes

# Defining high values for all features
age_high = 80
dys_high = 8 
post_hrv_high = 3
#hf_high = 70
#lf_high = 70
pre_hflf_high = 3
pac_high = 4
#eeg_power_high = 0.2
#egg_power_high = 0.05

# Scaling the variables
df.Age = df.Age / age_high
df.DyskinesiaScore = df.DyskinesiaScore / dys_high
#df_hrv.iloc[:,0:4] = df_hrv.iloc[:,0:4] / hrv_high
#df_hrv.iloc[:,12:16] = df_hrv.iloc[:,12:16] / hflf_high


# Getting features and the scores
#X = df.iloc[:,0:72] #For all subjects
#y = np.array(df.DyskinesiaScore) #Target variable
#y = np.nan_to_num(y)
#y = y.flatten()
#X = df.loc[df['PD_Healthy']==1,df.columns[0:36]] 
X = df.loc[df['PD_Healthy']==1,df.columns[0:-4]] #Only PD patients & PAC features
y = np.array(df.loc[df['PD_Healthy']==1,['DyskinesiaScore']])
y = np.nan_to_num(y)
y = y.flatten()



# Imputing missing values
#imputer = KNNImputer(n_neighbors=2, weights="uniform")
#X_imputed = imputer.fit_transform(X)
imputer = SimpleImputer(strategy='mean')
X_imputed = imputer.fit_transform(X)


## Log transforming and scaling the PAC features
X_imputed[:,0:36] = np.log(X_imputed[:,0:36])
X_imputed[:,0:36] = X_imputed[:,0:36] / pac_high


'''
# Scaling the features and the target variable
standard_scaler = preprocessing.StandardScaler()
X_imputed = standard_scaler.fit_transform(X_imputed)
y = standard_scaler.fit_transform(y.reshape(-1,1))
y = np.ravel(y)



# Scaling the variables between 0 and 1
min_max_scaler = preprocessing.MinMaxScaler()
X_imputed = min_max_scaler.fit_transform(X_imputed)
'''

# Creating a dataframe to crossverify the feature names that are selected
imputed_df = pd.DataFrame.from_records(X_imputed)
imputed_df.columns = df.columns[0:-4]

'''
# Creating a feature correlation heatmap 
corr_mat = X.corr()
plt.figure(figsize=(25, 16))
sns.heatmap(corr_mat, vmin=-1, vmax=1,annot=False,cmap='BrBG')

# Checking if PAC features are correlated to EEG and EGG power
test = imputed_df[['pac_tachy_CP_T4','tachy_power_T4','eeg_power_CP_T4']]
corr_mat = test.corr(method='spearman')
plt.figure(figsize=(25, 16))
sns.heatmap(corr_mat, vmin=-1, vmax=1,annot=True,cmap='BrBG')
plt.scatter(imputed_df['pac_normo_F_T1'], imputed_df['eeg_power_F_T1'])
plt.xlabel('pac')
plt.ylabel('eeg power')
plt.show()

'''

## Defining adjusted R^2 and AIC as the scoring function
def adjusted_R2(model,X,y_true):
    y_pred = model.predict(X)
    #R2 = r2_score(y_true, y_pred)
    #adjusted_R2 = 1-(1 - R2)*((len(X)-1)/(len(X)-len(X[0])-1))
    SST = sum((y_true-mean(y_true))**2)
    SSE = sum((y_true - y_pred)**2)
    n = len(X) # No. of observations
    p = len(X[0]) # No. of features
    adjusted_R2 = 1 - ((n-1)/(n-p))*(SSE/SST)
    return adjusted_R2

def aic(model,X,y_true) :
    y_pred = model.predict(X)
    SSE = sum((y_true - y_pred)**2)
    n = len(X) # No. of observations
    p = len(X[0]) # No. of features
    MSE = SSE/n
    aic =  n*log(MSE) + 2*(p+1) + ((2*(p+1)*(p+2))/(n-p))
    return aic

# Defining f test for the significance of the model

def f_test_model(model,X,y_true) : 
    y_pred = model.predict(X)
    TSS = sum((y_true-mean(y_true))**2)
    ESS = sum((y_pred-mean(y_pred))**2)
    RSS = sum((y_true-y_pred)**2)
    df_model = len(X[0])
    df_residuals = len(X) - df_model
    
    MSR = ESS/df_model
    MSE = RSS/df_residuals
    F_stat = MSR/MSE
    p_value = 1 - (stats.f.cdf(F_stat, df_model, df_residuals))
    return p_value



######################################################################################

## Initializing some variables
features_to_hold = [] #Selected features
R2_list = [] #List of R2 values to calculate R2_diff
R2_diff = 0
tol = 5 #Play around with this value (0.03 for adj_R2) (5 for aic)
features_to_permute = [] #Features to choose from

################# Sequential Feature Selection loop ##########################

X = X_imputed

## For the first round of selection

count = len(X[0])
adj_R2 = []

for i in range(0,count):
    #Use these two lines for this model
    #model = HuberRegressor(max_iter=100)
    #reg = model.fit(X[:,i].reshape(-1,1),y)
    #Use these two lines for this model
    model = sm.RLM(y,X[:,i].reshape(-1,1),M=sm.robust.norms.HuberT())
    reg = model.fit(scale_est = sm.robust.scale.HuberScale())
    adj_R2.append(adjusted_R2(reg,X[:,i].reshape(-1,1),y))
    
           
idx = adj_R2.index(max(adj_R2)) # Get the one with highest R2 score
#idx = adj_R2.index(min(adj_R2))
features_to_hold = X[:,idx] #Store it in an array 
features_to_permute = np.delete(X,idx,1)
R2_list.append(max(adj_R2)) # Store all the R2 values
#R2_list.append(min(adj_R2))    
    
        
## For subsequent rounds
                    
while R2_diff < tol: 


    count = len(features_to_permute[0])
    adj_R2 = []

    
    for i in range(0,count):
        features = np.column_stack((features_to_permute[:,i],features_to_hold))
        #Use these two lines for this model
        #model = HuberRegressor(max_iter=100)
        #reg = model.fit(features,y)
        #Use these two lines for this model
        model = sm.RLM(y,features,M=sm.robust.norms.HuberT())
        reg = model.fit(scale_est = sm.robust.scale.HuberScale())
        
        adj_R2.append(aic(reg,features,y))
       
        
    #idx = adj_R2.index(max(adj_R2))
    idx = adj_R2.index(min(adj_R2))
    #adding to the features to be held
    features_to_hold = np.column_stack((features_to_hold,features_to_permute[:,idx])) 
    features_to_permute = np.delete(features_to_permute,idx,1)
    #R2_list.append(max(adj_R2))
    R2_list.append(min(adj_R2))
    
    try:
        #R2_diff = R2_list[-2] - R2_list[-1] #For adjusted R2
        R2_diff = R2_list[-1] - R2_list[-2] # For AIC
    except:
        R2_diff = 0
        
    
 
 ## Getting the indices of the selected features
 
#idx = R2_list.index(max(R2_list))
idx = R2_list.index(min(R2_list)) 
features_to_consider = features_to_hold[:,0:idx+1]
count = len(features_to_consider[0])
features_list = []
for i in range(0,count):
    features_list.append(imputed_df.columns[imputed_df.eq(features_to_consider[:,i], axis=0).all(0)])
        
features_list = np.array(features_list).tolist()

for j in range(len(features_list)):
    name = features_list[j]
    features_list[j] = name[0]

############################ Fitting the best model #############################

#Using sklearn
#model = HuberRegressor()
#reg = model.fit(features_to_consider,y)
#reg_coefs = reg.coef_
#Using statsmodels
model = sm.RLM(y,features_to_consider,M=sm.robust.norms.HuberT())
reg = model.fit(scale_est = sm.robust.scale.HuberScale())
reg_coefs = reg.params
y_pred = reg.predict(features_to_consider)

#print(reg.summary(yname="Dyskinesia Score", xname=["var_%d" % i for i in range(len(reg.params))] ))
print(reg.summary(yname="Dyskinesia Score", xname= features_list ))
print('Adjusted R2 = ', adjusted_R2(reg,features_to_consider,y))
print('RMSE = ',mean_squared_error(y, y_pred, squared=False))
print('p-value for entire model = ', f_test_model(reg,features_to_consider,y))

################ Plotting the feature importance ############################

#Sorting the coefs
data = {'Coefs' : reg_coefs, 'Features' : features_list }
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

################ Handling no. of features > 12 #########################

max_features = 12
if len(features_to_consider[0]) > max_features:
    feat_idx = np.array(coef_df.index[0:max_features]).tolist()
    reg_coefs = np.array(coef_df.Coefs[0:max_features]).tolist()
    features_to_consider = features_to_consider[:,feat_idx]
    features_list = np.array(features_list)[feat_idx].tolist()
    
    #features_to_consider = features_to_consider[:,0:max_features]
    #features_list = features_list[0:max_features]
    
    ## Building another model
    model = sm.RLM(y,features_to_consider,M=sm.robust.norms.HuberT())
    reg = model.fit(scale_est = sm.robust.scale.HuberScale())
    reg_coefs = reg.params
    y_pred = reg.predict(features_to_consider)

    #print(reg.summary(yname="Dyskinesia Score", xname=["var_%d" % i for i in range(len(reg.params))] ))
    print(reg.summary(yname="Dyskinesia Score", xname= features_list ))
    print('Adjusted R2 = ', adjusted_R2(reg,features_to_consider,y))
    print('RMSE = ',mean_squared_error(y, y_pred, squared=False))
    print('p-value for entire model = ', f_test_model(reg,features_to_consider,y))
    
    ## Plotting their regression coefs 
    
    #Sorting the coefs
    data = {'Coefs' : reg_coefs, 'Features' : features_list }
    coef_df = pd.DataFrame(data)
    coef_df = coef_df.sort_values(by=['Coefs'], key=abs, ascending=False)
    
    fig, ax = plt.subplots()
    y_pos = np.arange(len(coef_df.Coefs))
    ax.barh(y_pos,coef_df.Coefs, align='center')
    ax.set_yticks(y_pos, labels=np.array(coef_df.Features))
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('Feature importance')
    ax.set_title('Dyskinesia score')

    plt.show()
    
## Saving the features in a dataframe and writing them to a csv file
features_df = pd.DataFrame.from_records(np.hstack((features_to_consider,y.reshape(-1,1))))
features_df.columns = [str(features_list[i]) for i in range(len(features_to_consider[0])) ] + ['y']

features_df.to_csv('Selected_features_v3.csv', index=False, header=True)

####################### Adding the next set of EEG & EGG features #########################

# Getting features and scores
#X = df.iloc[:,72:120] #All subjects & EEG,EGG features
X = df.loc[df['PD_Healthy']==1,df.columns[72:120]] #Only PD patients & EEG,EGG features
#y = np.array(df.loc[df['PD_Healthy']==1,['PartIVScore']])

#Imputing the missing values
X_imputed = imputer.fit_transform(X)

## Sqrt transformation (No need of scaling)
X_imputed = np.sqrt(X_imputed)

'''
# Scaling the features 
min_max_scaler = preprocessing.MinMaxScaler()
X_imputed = min_max_scaler.fit_transform(X_imputed)


# Scaling the features and the target variable
standard_scaler = preprocessing.StandardScaler()
X_imputed = standard_scaler.fit_transform(X_imputed)
y = standard_scaler.fit_transform(y.reshape(-1,1))
y = np.ravel(y)
'''

# Creating a dataframe to crossverify the feature names that are selected
imputed_df2 = pd.DataFrame.from_records(X_imputed)
#imputed_df.columns = df.columns[0:-4]
imputed_df2.columns = df.columns[72:120]
final_imputed_df = pd.concat([imputed_df,imputed_df2],axis=1)

## Initializing some variables
features_to_hold = [] #Selected features
R2_list = [] #List of R2 values to calculate R2_diff
R2_diff = 0
tol = 5 #Play around with this value (0.03 for adj_R2)
features_to_permute = [] #Features to choose from

## SFS loop ##

X = X_imputed
    
## For the first round of selection

count = len(X[0])
adj_R2 = []

for i in range(0,count):
    features = np.column_stack((features_to_consider,X[:,i].reshape(-1,1)))
    model = sm.RLM(y,features,M=sm.robust.norms.HuberT())
    reg = model.fit(scale_est = sm.robust.scale.HuberScale())
    adj_R2.append(aic(reg,features,y))
    
           
#idx = adj_R2.index(max(adj_R2)) # Get the one with highest R2 score
idx = adj_R2.index(min(adj_R2))
features_to_hold = np.column_stack((features_to_consider,X[:,idx])) #Store it in an array 
features_to_permute = np.delete(X,idx,1)
#R2_list.append(max(adj_R2)) # Store all the R2 values
R2_list.append(min(adj_R2))    
    
        
## For subsequent rounds
                    
while R2_diff < tol: 


    count = len(features_to_permute[0])
    adj_R2 = []

    
    for i in range(0,count):
        features = np.column_stack((features_to_permute[:,i],features_to_hold))
        #Use these two lines for this model
        #model = HuberRegressor(max_iter=100)
        #reg = model.fit(features,y)
        #Use these two lines for this model
        model = sm.RLM(y,features,M=sm.robust.norms.HuberT())
        reg = model.fit(scale_est = sm.robust.scale.HuberScale())
        
        adj_R2.append(aic(reg,features,y))
       
        
    #idx = adj_R2.index(max(adj_R2))
    idx = adj_R2.index(min(adj_R2))
    #adding to the features to be held
    features_to_hold = np.column_stack((features_to_hold,features_to_permute[:,idx])) 
    features_to_permute = np.delete(features_to_permute,idx,1)
    #R2_list.append(max(adj_R2))
    R2_list.append(min(adj_R2))
    
    try:
        #R2_diff = R2_list[-2] - R2_list[-1] #For adjusted R2
        R2_diff = R2_list[-1] - R2_list[-2] # For AIC
    except:
        R2_diff = 0
        
    
 
 ## Getting the indices of the selected features
 
#idx = R2_list.index(max(R2_list))
idx = R2_list.index(min(R2_list))
features_to_consider_final = features_to_hold[:,0:len(features_to_consider[0])+idx+1]
count = len(features_to_consider_final[0])
features_list = []
for i in range(0,count):
    features_list.append(final_imputed_df.columns[final_imputed_df.eq(features_to_consider_final[:,i], axis=0).all(0)])
        

########################## Fitting the final best model ###############################

model = sm.RLM(y,features_to_consider_final,M=sm.robust.norms.HuberT())
reg = model.fit(scale_est = sm.robust.scale.HuberScale())
reg_coefs = reg.params
y_pred = reg.predict(features_to_consider_final)
print(reg.summary(yname="y", xname=["var_%d" % i for i in range(len(reg.params))] ))
print('Adjusted R2 = ', adjusted_R2(reg,features_to_consider_final,y))
print('RMSE = ',mean_squared_error(y, y_pred, squared=False))
print('p-value for entire model = ', f_test_model(reg,features_to_consider_final,y))


################ Plotting the feature importance ############################


fig, ax = plt.subplots()
y_pos = np.arange(len(reg_coefs))
ax.barh(y_pos,reg_coefs, align='center')
ax.set_yticks(y_pos, labels=np.array(features_list))
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Feature importance')
ax.set_title('Dyskinesia score')

plt.show()
       
############ Checking the interaction of these features with demographics ##########

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


#Combining all the demographic variables
demographics = np.hstack((hrv_scaled,age_scaled,gender,np.array(df.PD_Healthy).reshape(-1,1)))

#Creating a dataframe with the required features
demographics_df = pd.DataFrame.from_records(demographics)
demographics_df.columns = np.hstack((df_hrv.columns[[0,1,2,3,8,9,10,11]],'Age','Gender_M','PD_Healthy'))
new_demo_df = demographics_df.loc[demographics_df['PD_Healthy'] == 1].reset_index() #For PD Patients
demographics_df.to_csv('Demographics_heart_features_v2.csv', index=False, header=True)

#new_demo_df = demographics_df.drop('PD_Healthy',axis=1) #For all subjects
features_df = pd.DataFrame.from_records(np.hstack((features_to_consider,y.reshape(-1,1))))
#features_df.columns = ['x' + str(i) for i in range(len(features_to_consider[0]))] + ['y']
features_df.columns = [str(features_list[i]) for i in range(len(features_to_consider[0])) ] + ['y']

merged_df = pd.concat([features_df,new_demo_df],axis=1)

#creating interaction terms
#des_form = ['y ~ ('] + ['x'+str(i)+'+' for i in range(len(features_to_consider[0]))] + [')*pre_HRV_T1']
des_form = ['y ~ ('] + [str(features_list[i]) + '+' for i in range(len(features_to_consider[0]))] + [')*pre_HRV_T1']
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
features = pd.concat([merged_df.Age,merged_df.Gender_M,merged_df.pre_HRV_T1,merged_df.pre_HRV_T2,
                      merged_df.pre_HRV_T3,merged_df.pre_HRV_T4,merged_df.pre_HFLF_T1,merged_df.pre_HFLF_T2,
                      merged_df.pre_HFLF_T3,merged_df.pre_HFLF_T4,merged_df.iloc[:,0:len(features_to_consider[0])]],
                     axis=1)
feat_list = features.columns.tolist()
'''

## Run this for just adding Age and Gender as covariates to the model
#features = pd.concat([merged_df.Age,merged_df.Gender_M,merged_df.iloc[:,0:len(features_to_consider[0])]],axis=1)
#feat_list = features.columns.tolist()


## Fitting the model

huber_t = sm.RLM(y, features, M=sm.robust.norms.HuberT())
hub_results = huber_t.fit(scale_est = sm.robust.scale.HuberScale())
#print(hub_results.summary(yname="y", xname=["var_%d" % i for i in range(len(hub_results.params))] ))
print(hub_results.summary(yname="y", xname= feat_list ))
y_pred = hub_results.predict(features)
print('Adjusted R2 = ', adjusted_R2(hub_results,np.array(features),y))
print('RMSE = ',mean_squared_error(y, y_pred, squared=False))
print('p-value for entire model = ', f_test_model(hub_results,np.array(features),y))

############## Visualizing the interaction terms ####################################

######## At different levels of cardiac features ##############
pvals = hub_results.pvalues
mean_arr = np.array(interaction_terms[1].mean())
reg_params = hub_results.params
mean_heart_meas = mean(interaction_terms[1].iloc[:,10])
sd_heart_meas = stdev(interaction_terms[1].iloc[:,10])

for i in range(10,len(pvals)) : 
    
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

        y_mean = intercept + reg_params[10]*mean_heart_meas + (reg_params[idx]*mean_heart_meas + reg_params[pac_idx])*x
        y_mean_p_sd = intercept + reg_params[10]*(mean_heart_meas +  sd_heart_meas) + \
            (reg_params[idx]*(mean_heart_meas +  sd_heart_meas) + reg_params[pac_idx])*x
        y_mean_m_sd = intercept + reg_params[10]*(mean_heart_meas -  sd_heart_meas) + \
            (reg_params[idx]*(mean_heart_meas -  sd_heart_meas) + reg_params[pac_idx])*x
    
        # Plotting the lines 

        plt.plot(x, y_mean, label = ''.join(["Mean ",interaction_terms[1].columns[10]])) 
        plt.plot(x, y_mean_p_sd, label = ''.join(["(Mean + sd)  ",interaction_terms[1].columns[10]])) 
        plt.plot(x, y_mean_m_sd, label = ''.join(["(Mean - sd)  ",interaction_terms[1].columns[10]])) 
        plt.legend()
        plt.title(''.join(['Interaction of ', colname, ' with ', interaction_terms[1].columns[10]])) 
        plt.xlabel(colname)
        plt.ylabel('Dyskinesia score')
        plt.show()

