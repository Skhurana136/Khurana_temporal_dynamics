# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 17:21:42 2020

@author: khurana
"""

#temporal_variation_factor_analysis
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
import statsmodels.formula.api as smf

head_path = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/headatinlet.csv"
head = pd.read_csv(head_path, sep = ",")
cov1 = np.round(np.cov(head["H1"]),2)
cov2 = np.round(np.cov(head["H2"]),2)
cov3 = np.round(np.cov(head["H3"]),2)

def corrfunc (data):
    normdata = data - np.mean(data)
    autocorr_f = (np.correlate(normdata, normdata, mode='full'))/(np.std(data)**2 * np.shape(data)[0])
    return autocorr_f

cor1 = corrfunc(head["H1"])[5474:]
cor2 = corrfunc(head["H2"])[5474:]
cor3 = corrfunc(head["H3"])[5474:]

cor1tim = np.where(cor1 < 0.75)[0][0]
cor2tim = np.where(cor2 < 0.75)[0][0]
cor3tim = np.where(cor3 < 0.75)[0][0]
mincor = cor2.min()

cleanup_tim = {"cov" : {1 : cov1, 2 : cov2, 5 : cov3},
               "timtrace" : {1 : cor1tim, 2: cor2tim, 5 : cor3tim},
               "chem_factors" : {'DOC':0,'DO':1, 'Ammonium':2, 'Nitrate':3, 'Sulphate':4, 'Particulate organic matter':5,
                                 'Mobile active aerobic degraders':6, 'Mobile active nitrate reducers':7,
                                 'Mobile active sulphate reducers':8, 'Mobile active ammonia oxidizers':9,
                                 'Mobile inactive aerobic degraders':10, 'Mobile inactive nitrate reducers':11,
                                 'Mobile inactive sulphate reducers':12, 'Mobile inactive ammonia oxidizers':13,
                                 'Nitrogen':14, 'TOC':15},
                   "Pe" : {'Slow':2,'Medium':11, 'Fast':22}}

directory = r"Y:\Home\khurana\4. Publications\Restructuring\Paper2\Figurecodes\/"
filename = "Normalized_RMSamplitude_chem.csv"
data = pd.read_csv(directory + filename, sep="\t")
print(data.columns)
print(data.shape)
print(data.Chem.unique())
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")
chemstoplot = ["DOC", "DO", "TOC", "Nitrogen"]
data["cov"] = data ["Time_series"]

data  = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/mass_flux_sensitivity_generalized.csv", sep="\t")

data["timtrace"] = data ["Time_series"]
data["chem_factors"] = data ["Chem"]
data["Pe"] = data ["Regime"]
data.replace(cleanup_tim, inplace = True)

dadata = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/Conc_da_ss.csv", sep = "\t")

mdata = pd.merge(data, dadata[["Regime", "Trial", "Chem", "%reldelmassflux", "Breakthroughtime", "avgconc_in"]], on = ["Regime", "Trial", "Chem"])
mdata.loc[mdata["PeDa"] > 40, "PeDamark"] = 3
mdata.loc[(mdata["PeDa"] > 20) & (mdata["PeDa"] < 40), "PeDamark"] = 2
mdata.loc[(mdata["PeDa"] > 1) & (mdata["PeDa"] < 20), "PeDamark"] = 1
mdata.loc[mdata["PeDa"] < 1, "PeDamark"] = 0
labels = {0 : "Da/Pe < 1",
          1 : "1 < Da/Pe < 15",
          2 : "15 < Da/Pe < 40",
          3 : "Da/Pe > 40"}

#Transient
allfeatures = dadata.columns.tolist()
features = ['Variance', 'Anisotropy', 'Regime', 'Chem', 'Breakthroughtime', 'fraction', 'cov', 'PeDa', 'PeDamark']
yfeature = ['Sensitivity', 'Sensitivitybase']

mydata = mdata[features+yfeature]
chemstoplot = ["DOC", "DO", "Nitrogen", "TOC"]
mydata = mydata[mydata['Chem'].isin (chemstoplot)]
mydata["Variance"]= mydata["Variance"].astype(str)
mydata["Anisotropy"]= mydata["Anisotropy"].astype(str)

colsymdict = {"Slow" : "indianred", "Medium" : "g", "Fast" : "steelblue",
              "DOC": 'D', "DO" : '^', "Nitrogen" : 's', "TOC" : 'o'}

#Plot yfeature with breakthrough time divided by regime, variance and anisotropy
import seaborn as sns
sns.pairplot(mydata, kind="scatter", hue="Regime", palette=colsymdict)

#fixed effect - fraction
#Random effect - Regime, covariance, 
xplot = "fraction"
yplot = 'Sensitivitybase'
y = mydata[yplot]

row = []
f1 = yplot + ' ~ ' + xplot
for g in ["Regime", "Variance", "Anisotropy", "cov"]:
    print(g)
    mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g])
    mdf = mymd.fit()
    y_predict = mdf.fittedvalues
    RMSE = np.sqrt(((y-y_predict)**2).mean())
    row.append([yplot, xplot,g,RMSE, "Fixed"])

g = "Regime"
mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g], re_formula = '~' + xplot)
mdf = mymd.fit()
y_predict = mdf.fittedvalues
RMSE = np.sqrt(((y-y_predict)**2).mean())
row.append([yplot, xplot,g,RMSE, "RandomFixed"])

mydata["R_V"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str)
mydata["R_A"] = mydata["Regime"] + "_" + mydata["Anisotropy"].astype(str)
mydata["V_A"] = mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str)
mydata["R_V_A"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str)
mydata["R_V_C"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Chem"]
mydata["R_A_C"] = mydata["Regime"] + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata["V_A_C"] = mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata["R_V_A_C"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"]
mydata ["R_C"] = mydata["Regime"] + "_" + mydata["Chem"]
mydata["R_V_Co"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["cov"].astype(str)
mydata["R_A_Co"] = mydata["Regime"] + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["cov"].astype(str)
mydata["V_A_Co"] = mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["cov"].astype(str)
mydata["R_V_A_Co"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["cov"].astype(str)
mydata["R_V_C_Co"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Chem"] + "_" + mydata["cov"].astype(str)
mydata["R_A_C_Co"] = mydata["Regime"] + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"] + "_" + mydata["cov"].astype(str)
mydata["V_A_C_Co"] = mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"] + "_" + mydata["cov"].astype(str)
mydata["R_V_A_C_Co"] = mydata["Regime"] + "_" + mydata["Variance"].astype(str) + "_" + mydata["Anisotropy"].astype(str) + "_" + mydata["Chem"] + "_" + mydata["cov"].astype(str)
mydata ["R_C_Co"] = mydata["Regime"] + "_" + mydata["Chem"] + "_" + mydata["cov"].astype(str)
mydata ["R_Co"] = mydata["Regime"] + "_" + mydata["cov"].astype(str)

for g in ["R_Co", "R_V_A_C_Co"]:
    print(g)
    mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g], re_formula = '~' + xplot)
    mdf = mymd.fit()
    y_predict = mdf.fittedvalues
    RMSE = np.sqrt(((y-y_predict)**2).mean())
    row.append([yplot, xplot,g,RMSE, "RandomFixed"])

for g in ["PeDamark", "cov"]:
    print(g)
    mymd = smf.mixedlm(formula = f1, data = mydata, groups = mydata[g], re_formula = '~' + xplot)
    mdf = mymd.fit()
    y_predict = mdf.fittedvalues
    RMSE = np.sqrt(((y-y_predict)**2).mean())
    row.append([yplot, xplot,g,RMSE, "RandomFixed"])
    
results = pd.DataFrame.from_records(row, columns = ["response", "independent", "Details", "RMSE", "Model_Type"])
results2 = pd.DataFrame.from_records(row, columns = ["response", "independent", "Details", "RMSE", "Model_Type"])

#LMEM is not improving the results considering PeDamark alone
mydata.columns
sns.scatterplot(mydata["Breakthroughtime"]/mydata["cov"], "Sensitivity", hue = "PeDamark", data = mydata)

#Trying sci-kit learn for feature selection
#Transient
allfeatures = mdata.columns.tolist()
features = ['Variance', 'Anisotropy', 'Pe', 'Chem','conc_Da', 'fraction', 'cov', 'PeDa', 'PeDamark']
yfeature = ['Sensitivity', 'Sensitivitybase']

mydata = mdata[features+yfeature]
chemstoplot = ["DOC", "DO", "Nitrogen", "TOC"]
mydata = mydata[mydata['Chem'].isin (chemstoplot)]
mydata["Variance"]= mydata["Variance"].astype(str)
mydata["Anisotropy"]= mydata["Anisotropy"].astype(str)

from sklearn.ensemble import RandomForestRegressor as rf
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.svm import SVR
from sklearn import linear_model
from sklearn.feature_selection import RFE, SelectFromModel

feattofilt = ['Variance', 'Anisotropy', 'Pe', 'conc_Da', 'fraction', 'cov']
X = mydata[feattofilt]
y = mydata[yfeature[0]]
X_train, X_test, y_train, y_test = train_test_split (X, y, test_size = 0.2, random_state = 0)
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

regressor = rf(n_estimators = 20, random_state = 0)
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))

selector = RFE (regressor, n_features_to_select = 3, step = 1)
selector = selector.fit(X, y)
ranking = selector.ranking_

estimator = SVR(kernel = "linear")
estimator.fit(X_train, y_train)
y_pred = regressor.predict(X_test)
print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))
selector = RFE (estimator, n_features_to_select = 2, step = 1)
selector = selector.fit(X, y)
ranking = selector.ranking_

feattofilt = ['Variance', 'Anisotropy', 'Pe', 'conc_Da', 'fraction', 'cov']
X = mydata[feattofilt]
y = mydata[yfeature[0]]
X_train, X_test, y_train, y_test = train_test_split (X, y, test_size = 0.2, random_state = 0)
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

lasso = linear_model.Lasso(alpha=0.1)
ols = linear_model.LinearRegression()
svr = SVR(kernel = "linear")

estilist = [svr, ols, lasso]

for e in estilist:
    selector = RFE(e, n_features_to_select = 2, step = 1).fit(X_train, y_train)
    ranking = selector.ranking_
    maxidx = list(np.where(ranking==1))
    impfeatures = np.array(feattofilt)[maxidx].tolist()
    print(e, "Important features: ")
    print(impfeatures)

for e in estilist:
    selector = SelectFromModel(e).fit(X, y)
    ranking = selector.get_support()
    impfeatures = np.array(feattofilt)[ranking].tolist()
    print(e, " Important features: ")
    print(impfeatures)