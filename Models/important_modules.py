# -*- coding: utf-8 -*-
"""
Created on Sat May  4 12:30:33 2024

@author: sanke
"""
from sklearn.linear_model import LogisticRegression, HuberRegressor, LinearRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from statistics import mean
from math import log
import scipy.stats as stats
import numpy as np
import statsmodels.api as sm
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error, accuracy_score, matthews_corrcoef,\
    confusion_matrix
from sklearn.model_selection import RepeatedKFold, RepeatedStratifiedKFold
from sklearn.preprocessing import KBinsDiscretizer

# Defining R2 score

def R2(model,X,y_true):
    y_pred = model.predict(X)
    r2 = r2_score(y_true,y_pred) 
    return r2

def rmse(model,X,y_true):
    y_pred = model.predict(X)
    rmse = mean_squared_error(y_true, y_pred, squared=False)
    return rmse

def rmse_alt(y_true,y_pred):
    #y_pred = model.predict(X)
    rmse = mean_squared_error(y_true, y_pred, squared=False)
    return rmse

## Defining cohen's f2 for linear regression model

def cohensf2(feat_df, y, feature_name):
    
    # Defining the full model
    model = sm.RLM(y,feat_df,M=sm.robust.norms.HuberT())
    reg = model.fit(scale_est = sm.robust.scale.HuberScale())
    y_pred = reg.predict(feat_df)
    R2ab = r2_score(y,y_pred)
    
    # Defining the model without the desired features
    red_feat_df = feat_df.drop(feature_name, axis=1)
    model = sm.RLM(y,red_feat_df,M=sm.robust.norms.HuberT())
    reg = model.fit(scale_est = sm.robust.scale.HuberScale())
    y_pred = reg.predict(red_feat_df)
    R2a = r2_score(y,y_pred)
    
    # Calculating the cohensf2
    f2 = (R2ab - R2a) / (1 - R2ab)
    return f2

## Defining Cohen's F2 for logistic regression model
# feat_df = Features dataframe with features
# y = The target variable
# feature_name = The feature name for which cohen's f2 has to be calculated

def cohensf2_logist(feat_df, y, feature_name) : 
    
    temp_df = feat_df.copy()
    
    # Full model: All predictors
    logit_full = sm.Logit(y, sm.add_constant(temp_df)).fit(method='bfgs',maxiter = 1000, disp=0)

    # Reduced model: Exclude one predictor (e.g., excluding the last one)
    temp_df = temp_df.drop(feature_name, axis=1)
    logit_reduced = sm.Logit(y, sm.add_constant(temp_df)).fit(method='bfgs',maxiter = 1000, disp=0)

    # McFadden's pseudo-R2
    R2_full = 1 - logit_full.llf / logit_full.llnull
    R2_reduced = 1 - logit_reduced.llf / logit_reduced.llnull

    # Compute f^2
    f2 = (R2_full - R2_reduced) / (1 - R2_full)
    print(f"Cohen's f^2: {f2:.4f}")
    return f2


# Defining the model evaluation metric adjusted R2
def adjusted_R2(model,X,y_true):
    y_pred = model.predict(X)
    #R2 = r2_score(y_true, y_pred)
    #adjusted_R2 = 1-(1 - R2)*((len(X)-1)/(len(X)-len(X[0])-1))
    SST = sum((y_true-mean(y_true))**2)
    SSE = sum((y_true - y_pred)**2)
    n = len(X)
    p = len(X[0])+1
    adjusted_R2 = 1 - ((n-1)/(n-p))*(SSE/SST)
    return adjusted_R2

# Defining the model evaluation metric adjusted R2(for the SFS_regress module)
def adjusted_R2_alt(X,y_true,y_pred):
    #y_pred = model.predict(X)
    #R2 = r2_score(y_true, y_pred)
    #adjusted_R2 = 1-(1 - R2)*((len(X)-1)/(len(X)-len(X[0])-1))
    SST = sum((y_true-mean(y_true))**2)
    SSE = sum((y_true - y_pred)**2)
    n = len(X)
    p = len(X[0])+1
    adjusted_R2 = 1 - ((n-1)/(n-p))*(SSE/SST)
    return adjusted_R2


# Defining the scoring metric AIC
def aic(model,X,y_true) :
    y_pred = model.predict(X)
    SSE = sum((y_true - y_pred)**2)
    n = len(X) # No. of observations
    p = len(X[0]) # No. of features
    MSE = SSE/n
    aic =  n*log(MSE) + 2*(p+1) + ((2*(p+1)*(p+2))/(n-p))
    return aic

# Defining the scoring metric AIC (for the SFS_regress module)
def aic_alt(X,y_true,y_pred) :
    #y_pred = model.predict(X)
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

def acc_cv(X,y_true) : 
    acc = []
    rkf = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 2, random_state= 1)
    for i, (train_index, test_index) in enumerate(rkf.split(X,y_true)):
        ## Splitting into train and test
        X_train = X[train_index,:]
        X_test = X[test_index,:]
        y_train = y_true[train_index]
        y_test = y_true[test_index]
        #Fitting the model on training data 
        try :
            #model = sm.Logit(y_train,X_train)
            #reg = model.fit()
            model = LogisticRegression()
            #model = DecisionTreeClassifier()
            #model = RandomForestClassifier(max_depth=4, n_jobs = -1,random_state =1)
            reg = model.fit(X_train,y_train) 
        except :
            continue
        #Predicting on test dataset
        #yhat = reg.predict(X_test) 
        #y_pred = list(map(round, yhat))
        
        y_pred = reg.predict(X_test)
        #Calculating accuracy score
        acc.append(accuracy_score(y_test,y_pred))
    
    try :
        acc_score = mean(acc)
    except:
        acc_score = 0
        
    return acc_score


def mcc_cv(X,y_true) : 
    mcc = []
    rkf = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 2, random_state= 1)
    for i, (train_index, test_index) in enumerate(rkf.split(X,y_true)):
        ## Splitting into train and test
        X_train = X[train_index,:]
        X_test = X[test_index,:]
        y_train = y_true[train_index]
        y_test = y_true[test_index]
        #Fitting the model on training data 
        try :
            #model = sm.Logit(y_train,X_train)
            #reg = model.fit()
            model = LogisticRegression()
            #model = DecisionTreeClassifier()
            #model = RandomForestClassifier(max_depth=4, n_jobs = -1,random_state =1)
            reg = model.fit(X_train,y_train) 
        except :
            continue
        #Predicting on test dataset
        #yhat = reg.predict(X_test) 
        #y_pred = list(map(round, yhat))
        y_pred = reg.predict(X_test)
        #Calculating accuracy score
        mcc.append(matthews_corrcoef(y_test,y_pred))
    
    try :
        mcc_score = mean(mcc)
    except:
        mcc_score = 0
        
    return mcc_score

def spec_cv(X,y_true) : 
    spec = []
    rkf = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 2, random_state= 1)
    for i, (train_index, test_index) in enumerate(rkf.split(X,y_true)):
        ## Splitting into train and test
        X_train = X[train_index,:]
        X_test = X[test_index,:]
        y_train = y_true[train_index]
        y_test = y_true[test_index]
        #Fitting the model on training data 
        try :
            #model = sm.Logit(y_train,X_train)
            #reg = model.fit()
            model = LogisticRegression()
            #model = DecisionTreeClassifier()
            #model = RandomForestClassifier(max_depth=4, n_jobs = -1,random_state =1)
            reg = model.fit(X_train,y_train) 
        except :
            continue
        #Predicting on test dataset
        #yhat = reg.predict(X_test) 
        #y_pred = list(map(round, yhat))
        y_pred = reg.predict(X_test)
        #Calculating accuracy score
        tn, fp, fn, tp = confusion_matrix(y_test,y_pred).ravel()
        spec_temp = tn / (tn+fp)
        spec.append(spec_temp)
    
    try :
        spec_score = mean(spec)
    except:
        spec_score = 0
        
    return spec_score

def rmse_cv(X,y_true) : 
    ## Discretize the target variable
    binner = KBinsDiscretizer(n_bins = 3, encode = 'ordinal', strategy = 'quantile')
    y_labels = binner.fit_transform(y_true.reshape(-1,1))
    
    aic_score = []
    rkf = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 1, random_state= 1)
    for i, (train_index, test_index) in enumerate(rkf.split(X,y_labels)):
        ## Splitting into train and test
        X_train = X[train_index,:]
        X_test = X[test_index,:]
        y_train = y_true[train_index]
        y_test = y_true[test_index]
        #Fitting the model on training data 
        try :
            model = sm.RLM(y_train,X_train,M=sm.robust.norms.HuberT())
            reg = model.fit(scale_est = sm.robust.scale.HuberScale())
            #model = HuberRegressor()
            #model = LinearRegression(fit_intercept=False)
            #reg = model.fit(X_train,y_train)
        except :
            continue
        #Predicting on test dataset
        y_pred = reg.predict(X_test)
        #Calculating accuracy score
        aic_score.append(rmse_alt(y_test,y_pred))
    
    try :
        mean_aic = mean(aic_score)
    except:
        mean_aic = 0
        
    return mean_aic

def adj_R2_cv(X,y_true) : 
    ## Discretize the target variable
    binner = KBinsDiscretizer(n_bins = 3, encode = 'ordinal', strategy = 'quantile')
    y_labels = binner.fit_transform(y_true.reshape(-1,1))
    
    R2_score = []
    rkf = RepeatedStratifiedKFold(n_splits = 6, n_repeats = 3, random_state= 2)
    for i, (train_index, test_index) in enumerate(rkf.split(X,y_labels)):
        ## Splitting into train and test
        X_train = X[train_index,:]
        X_test = X[test_index,:]
        y_train = y_true[train_index]
        y_test = y_true[test_index]
        #Fitting the model on training data 
        try :
            #model = sm.RLM(y_train,X_train,M=sm.robust.norms.HuberT())
            #reg = model.fit(scale_est = sm.robust.scale.HuberScale())
            #model = HuberRegressor()
            model = LinearRegression(fit_intercept=False)
            reg = model.fit(X_train,y_train)
        except :
            continue
        #Predicting on test dataset
        y_pred = reg.predict(X_test)
        #Calculating accuracy score
        R2_score.append(adjusted_R2_alt(X_test,y_test,y_pred))
    
    try :
        mean_R2 = mean(R2_score)
    except:
        mean_R2 = 0
        
    return mean_R2
        
## Pvalues for the logistic regression 
def logit_pvalue(model, x):
    """ Calculate z-scores for scikit-learn LogisticRegression.
    parameters:
        model: fitted sklearn.linear_model.LogisticRegression with intercept and large C
        x:     matrix on which the model was fit
    This function uses asymtptics for maximum likelihood estimates.
    """
    p = model.predict_proba(x)
    n = len(p)
    m = len(model.coef_[0]) + 1
    coefs = np.concatenate([model.intercept_, model.coef_[0]])
    x_full = np.matrix(np.insert(np.array(x), 0, 1, axis = 1))
    ans = np.zeros((m, m))
    for i in range(n):
        ans = ans + np.dot(np.transpose(x_full[i, :]), x_full[i, :]) * p[i,1] * p[i, 0]
    vcov = np.linalg.inv(np.matrix(ans))
    se = np.sqrt(np.diag(vcov))
    t =  coefs/se  
    p = (1 - norm.cdf(abs(t))) * 2
    return p

## Calculating average marginal effects for logistic regression

def AME(df,y,feature_name) :
    
    # 1. Get the mean of feature and calculate the small change
    calc_df = df.copy() # This is to save all the extra calculation
    predict_df = df.copy() # This is to save the increment of the desired features
    mean_feature = predict_df[feature_name].mean()
    std_feature = predict_df[feature_name].std()
    h = (abs(mean_feature) + 0.0001) * 0.0001
    print("Small change (h):", h)

    # 2. Fit the logistic regression model
    #model = sm.Logit(y,sm.add_constant(predict_df)).fit(method='bfgs',maxiter = 1000, disp=0) 
    model = LogisticRegression(C=1e30).fit(predict_df,y)
    
    # 3. Predict log-odds for the original data
    #calc_df['lw_0'] = model.predict(sm.add_constant(predict_df))
    calc_df['lw_0'] = model.predict_proba(predict_df)[:,1]
    
    # 4. Change feature by a small amount
    predict_df[feature_name] += h  # Increment cigs by the small change
    #calc_df['lw_1'] = model.predict(sm.add_constant(predict_df)) # Predict log-odds after the change
    calc_df['lw_1'] = model.predict_proba(predict_df)[:,1]
    
    # 5. Calculate the change in predicted probabilities
    calc_df['dydx'] = (calc_df['lw_1'] - calc_df['lw_0']) / h  # Marginal effect for each observation

    # 6. Average the marginal effects
    average_dydx = calc_df['dydx'].mean()
    print("Average marginal effect (dy/dx):", average_dydx)
    return average_dydx
    
    
# Sequential feature selection loop (normal)
## Input : Features, target variable labels, the dataframe with imputed features
## Output : reduced feature matrix, names of the features and the list of AIC values

def SFS(X_train,y_train,X_test,y_test,imputed_df,train_indices,test_indices) : 
    
    ## Initializing some variables
    features_to_hold = [] #Selected features
    R2_list = [] #List of R2 values to calculate R2_diff
    R2_diff = 0
    tol = 5 #Play around with this value (0.03 for adj_R2) (5 for aic)
    features_to_permute = [] #Features to choose from

    ################# Sequential Feature Selection loop ##########################
        
    ## For the first round of selection

    count = len(X_train[0])
    adj_R2 = []

    for i in range(0,count):
        #Use these two lines for this model
        #model = HuberRegressor(max_iter=100)
        #reg = model.fit(X[:,i].reshape(-1,1),y)
        #Use these two lines for this model
        model = sm.RLM(y_train,X_train[:,i].reshape(-1,1),M=sm.robust.norms.HuberT())
        reg = model.fit(scale_est = sm.robust.scale.HuberScale())
        adj_R2.append(rmse(reg,X_test[:,i].reshape(-1,1),y_test))
        
               
    #idx = adj_R2.index(max(adj_R2)) # Get the one with highest R2 score
    idx = adj_R2.index(min(adj_R2)) #Get the one with lowest AIC score
    
    features_to_hold_train = X_train[:,idx] #Store the selected features in an array
    features_to_hold_test = X_test[:,idx]
    
    features_to_permute_train = np.delete(X_train,idx,1) #Store the features that need to be permuted
    features_to_permute_test = np.delete(X_test,idx,1) #Store the features that need to be permuted
    
    #R2_list.append(max(adj_R2)) # Store all the R2 values
    R2_list.append(min(adj_R2))   #Store the AIC values 
        
            
    ## For subsequent rounds
                        
    while R2_diff < tol: 


        count = len(features_to_permute_train[0])
        adj_R2 = []

        
        for i in range(0,count):
            features_train = np.column_stack((features_to_permute_train[:,i],features_to_hold_train))
            features_test = np.column_stack((features_to_permute_test[:,i],features_to_hold_test))
            #Use these two lines for this model
            #model = HuberRegressor(max_iter=100)
            #reg = model.fit(features,y)
            #Use these two lines for this model
            model = sm.RLM(y_train,features_train,M=sm.robust.norms.HuberT())
            reg = model.fit(scale_est = sm.robust.scale.HuberScale())
            
            adj_R2.append(rmse(reg,features_test,y_test))
           
            
        #idx = adj_R2.index(max(adj_R2))
        idx = adj_R2.index(min(adj_R2))
        #adding to the features to be held
        features_to_hold_train = np.column_stack((features_to_hold_train,features_to_permute_train[:,idx]))
        features_to_hold_test = np.column_stack((features_to_hold_test,features_to_permute_test[:,idx]))
        
        features_to_permute_train = np.delete(features_to_permute_train,idx,1)
        features_to_permute_test = np.delete(features_to_permute_test,idx,1)
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
    features_to_consider = np.concatenate([features_to_hold_train[:,0:idx+1],
                                           features_to_hold_test[:,0:idx+1]],axis=0)
    
    # Create a dataframe with the indices and then sort the rows
    indices = np.concatenate([train_indices,test_indices],axis=0)
    temp_df = pd.DataFrame(data=features_to_consider,index=indices)
    temp_df = temp_df.sort_index()
    features_to_consider = temp_df.to_numpy()
    
    count = len(features_to_consider[0])
    features_list = []
    for i in range(0,count):
        features_list.append(imputed_df.columns[imputed_df.eq(features_to_consider[:,i], axis=0).all(0)])
        
    return features_to_consider,features_list,R2_list

# Sequential feature selection loop (normal)
## Input : Features, target variable labels, the dataframe with imputed features
## Output : reduced feature matrix, names of the features and the list of AIC values

def SFS_alt(X,y,imputed_df,train_indices) : 
    
    ## Initializing some variables
    features_to_hold = [] #Selected features
    R2_list = [] #List of R2 values to calculate R2_diff
    R2_diff = 0
    tol = 1 #Play around with this value (0.03 for adj_R2) (5 for aic)
    features_to_permute = [] #Features to choose from

    ################# Sequential Feature Selection loop ##########################
        
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
        adj_R2.append(aic(reg,X[:,i].reshape(-1,1),y))
        
               
    #idx = adj_R2.index(max(adj_R2)) # Get the one with highest R2 score
    idx = adj_R2.index(min(adj_R2)) #Get the one with lowest AIC score
    
    features_to_hold = X[:,idx] #Store the selected features in an array
    
    features_to_permute = np.delete(X,idx,1) #Store the features that need to be permuted
    
    #R2_list.append(max(adj_R2)) # Store all the R2 values
    R2_list.append(min(adj_R2))   #Store the AIC values 
        
            
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
    
    # Create a dataframe with the indices and then sort the rows
    temp_df = pd.DataFrame(data=features_to_consider,index=train_indices)
    temp_df = temp_df.sort_index()
    features_to_consider = temp_df.to_numpy()
    
    count = len(features_to_consider[0])
    features_list = []
    for i in range(0,count):
        features_list.append(imputed_df.columns[imputed_df.eq(features_to_consider[:,i], axis=0).all(0)])
        
    return features_to_consider,features_list,R2_list


# Sequential feature selection loop for classification (Logistic regression)
## Input : Features, target variable labels, the dataframe with imputed features
## Output : reduced feature matrix, names of the features and the list of accuracy values
## Note : Change the model in the acc_cv function 

def SFS_logist(X_train,y_train,imputed_df,max_num_features) : 
    
    ## Initializing some variables
    features_to_hold = [] #Selected features
    R2_list = [] #List of R2 values to calculate R2_diff
    R2_diff = 0
    tol = 0.03 #Play around with this value (0.03 for adj_R2) (5 for aic)
    features_to_permute = [] #Features to choose from

    ################# Sequential Feature Selection loop ##########################
        
    ## For the first round of selection

    count = len(X_train[0])
    adj_R2 = []

    for i in range(0,count):
        adj_R2.append(acc_cv(X_train[:,i].reshape(-1,1),y_train))
        
               
    idx = adj_R2.index(max(adj_R2)) # Get the one with highest accuracy score
    #idx = adj_R2.index(min(adj_R2)) #Get the one with lowest AIC score
    
    features_to_hold = X_train[:,idx] #Store the selected features in an array    
    features_to_permute = np.delete(X_train,idx,1) #Store the features that need to be permuted
   
    
    R2_list.append(max(adj_R2)) # Store all the acc values
    #R2_list.append(min(adj_R2))   #Store the AIC values 
        
            
    ## For subsequent rounds
                        
    while R2_diff < tol:       

        count = len(features_to_permute[0])
        adj_R2 = []

        
        for i in range(0,count):
            features = np.column_stack((features_to_permute[:,i],features_to_hold))
            adj_R2.append(acc_cv(features,y_train))
           
            
        idx = adj_R2.index(max(adj_R2))
        #idx = adj_R2.index(min(adj_R2))
        
        #adding to the features to be held
        features_to_hold = np.column_stack((features_to_hold,features_to_permute[:,idx]))     
        features_to_permute = np.delete(features_to_permute,idx,1)
        
        R2_list.append(max(adj_R2))
        #R2_list.append(min(adj_R2))
        
        try:
            R2_diff = R2_list[-2] - R2_list[-1] #For maximizing accuracy score
            #R2_diff = R2_list[-1] - R2_list[-2] # For AIC
        except:
            R2_diff = 0
            
        if len(features_to_hold[0]) > max_num_features:
            break
            
        
     
     ## Getting the indices of the selected features
     
    idx = R2_list.index(max(R2_list))
    #idx = R2_list.index(min(R2_list))
    features_to_consider = features_to_hold[:,0:idx+1]    
        
    count = len(features_to_consider[0])
    features_list = []
    for i in range(0,count):
        features_list.append(imputed_df.columns[imputed_df.eq(features_to_consider[:,i], axis=0).all(0)])
        
    return features_to_consider,features_list,R2_list

# Sequential feature selection loop for classification (Logistic regression)
## Input : Features, target variable labels, the dataframe with imputed features
## Output : reduced feature matrix, names of the features and the list of aic values
## Note : Change the model in the rmse_cv function 

def SFS_regress(X_train,y_train,imputed_df,max_num_features) : 
    
    ## Initializing some variables
    features_to_hold = [] #Selected features
    R2_list = [] #List of R2 values to calculate R2_diff
    R2_diff = 0
    tol = 0.03 #Play around with this value (0.03 for adj_R2) (5 for aic)
    features_to_permute = [] #Features to choose from

    ################# Sequential Feature Selection loop ##########################
        
    ## For the first round of selection

    count = len(X_train[0])
    adj_R2 = []

    for i in range(0,count):
        adj_R2.append(rmse_cv(X_train[:,i].reshape(-1,1),y_train))
        #adj_R2.append(adj_R2_cv(X_train[:,i].reshape(-1,1),y_train))
        
               
    #idx = adj_R2.index(max(adj_R2)) # Get the one with highest accuracy score
    idx = adj_R2.index(min(adj_R2)) #Get the one with lowest AIC score
    
    features_to_hold = X_train[:,idx] #Store the selected features in an array    
    features_to_permute = np.delete(X_train,idx,1) #Store the features that need to be permuted
   
    
    #R2_list.append(max(adj_R2)) # Store all the acc values
    R2_list.append(min(adj_R2))   #Store the AIC values 
        
            
    ## For subsequent rounds
                        
    while R2_diff < tol:       

        count = len(features_to_permute[0])
        adj_R2 = []

        
        for i in range(0,count):
            features = np.column_stack((features_to_permute[:,i],features_to_hold))
            adj_R2.append(rmse_cv(features,y_train))
            #adj_R2.append(adj_R2_cv(features,y_train))
           
            
        #idx = adj_R2.index(max(adj_R2))
        idx = adj_R2.index(min(adj_R2))
        
        #adding to the features to be held
        features_to_hold = np.column_stack((features_to_hold,features_to_permute[:,idx]))     
        features_to_permute = np.delete(features_to_permute,idx,1)
        
        #R2_list.append(max(adj_R2))
        R2_list.append(min(adj_R2))
        
        try:
            #R2_diff = R2_list[-2] - R2_list[-1] #For maximizing accuracy score
            R2_diff = R2_list[-1] - R2_list[-2] # For AIC
        except:
            R2_diff = 0
            
        if len(features_to_hold[0]) > max_num_features:
            break
            
        
     
     ## Getting the indices of the selected features
     
    #idx = R2_list.index(max(R2_list))
    idx = R2_list.index(min(R2_list))
    features_to_consider = features_to_hold[:,0:idx+1]    
        
    count = len(features_to_consider[0])
    features_list = []
    for i in range(0,count):
        features_list.append(imputed_df.columns[imputed_df.eq(features_to_consider[:,i], axis=0).all(0)])
        
    return features_to_consider,features_list,R2_list

