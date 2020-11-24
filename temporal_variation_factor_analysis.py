# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 17:21:42 2020

@author: khurana
"""

#temporal_variation_factor_analysis
#!pip install factor_analyzer
from factor_analyzer import FactorAnalyzer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
import seaborn as sns


#Figure S1: Autocorrelation discussion for imposed temporal dynamics

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

mdata = pd.merge(data, dadata[["Regime", "Trial", "Chem", "%reldelmassflux", "%fraction_rel_delmf"]], on = ["Regime", "Trial", "Chem"])
mdata.loc[mdata["PeDa"] > 40, "PeDamark"] = 3
mdata.loc[(mdata["PeDa"] > 20) & (mdata["PeDa"] < 40), "PeDamark"] = 2
mdata.loc[(mdata["PeDa"] > 1) & (mdata["PeDa"] < 20), "PeDamark"] = 1
mdata.loc[mdata["PeDa"] < 1, "PeDamark"] = 0
labels = {0 : "Da/Pe < 1",
          1 : "1 < Da/Pe < 15",
          2 : "15 < Da/Pe < 40",
          3 : "Da/Pe > 40"}

corrMatrix = mdata.drop(['Unnamed: 0', 'Unnamed: 0.1', 'Sensitivity%', 'Sensitivitybase%', 'PeDamark', 'sensbase', 'timtrace', 'chem_factors', 'Time_series'], axis = 1).corr()
sns.heatmap(corrMatrix, annot = True)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    p2 = sns.heatmap(corrMatrix, mask=mask, square=True, cmap = "vlag")

for k in ['pearson', 'kendall', 'spearman']:
    corrMatrix = mdata.drop(['Unnamed: 0', 'Unnamed: 0.1', 'Sensitivity%', 'Sensitivitybase%', 'PeDamark', 'sensbase', 'timtrace', 'chem_factors', 'Time_series'], axis = 1).corr(method = k)
    plt.figure()
    sns.heatmap(corrMatrix, annot = True)


#Subset of the data, the 14 columns containing the survey answers

features = ['conc_Da', 'Pe', '%reldelmassflux', 'cov', 'timtrace','Time', 'PeDa', '%fraction_rel_delmf']
yfeature = 'chem_factors'
x = mdata[features+ [yfeature]] 
fa = FactorAnalyzer(2, rotation = 'varimax')
fa.fit(x[features] , y = x[yfeature])
#Get Eigen values and plot them
ev, v = fa.get_eigenvalues()
ev
plt.plot(range(1,x[features].shape[1]+1),ev)

from sklearn.decomposition import FactorAnalysis, PCA
X = mdata[[features]] 
transformer = FactorAnalysis(n_components=2, random_state=0)
X_fit = transformer.fit(X[features.remove(yfeature)], y = X[[yfeature]])
X_transformed = transformer.fit_transform(X[features.remove(yfeature)], y = X[[yfeature]])
X_transformed.shape
X_fit.get_covariance()
X_fit.score(X[features], y = X[[yfeature]])

x = mdata[features + [yfeature]]
fa = FactorAnalyzer(2, rotation = 'varimax')
fa.fit(x[features] , y = x[yfeature])
#Get Eigen values and plot them
ev, v = fa.get_eigenvalues()
ev
plt.plot(range(1,x[features].shape[1]+1),ev)

for nc in [4, 5]:
    pca = PCA(n_components = nc)
    pcadata = mdata[features]
    #X = pcadata[['conc_Da','cov', 'fraction', 'Time']]
    X = mdata[features+[yfeature]] 
    y = X[[yfeature]]
    scaler.fit(X[features])
    X_scaled = scaler.transform(X[features])
    pca.fit(X_scaled)
    X_pca = pca.transform(X_scaled)
    
    # Create a new dataset from principal components
    y1 = pd.concat([y, mdata[['PeDamark', 'Senssquared']]], axis = 1)
    result_df = pd.concat([pd.DataFrame(data = X_pca), y1], axis=1)
    result_df.head(5)
    
    # Visualize Principal Components with a scatter plot
    
    fig = plt.figure(figsize = (12,10))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('First Principal Component ', fontsize = 15)
    ax.set_ylabel('Second Principal Component ', fontsize = 15)
    ax.set_title('Principal Component Analysis with components:' + str(nc), fontsize = 20)
    ax.scatter(result_df[0], result_df[1], s = 50*result_df['Senssquared'], c = result_df['chem_factors'], alpha = 0.7)
    ax.legend()
    
    print("PCA component shape: {}".format(pca.components_.shape))
    print("PCA components:\n{}".format(pca.components_))
    print("PCA components:\n{}".format(pca.explained_variance_ratio_))
    
    plt.matshow(pca.components_, cmap='viridis')
    plt.yticks(list(range(nc)), ["First component", "Second component", "Third component", "Fourth component"])
    plt.colorbar()
    plt.xticks(range(len(X[features].columns.tolist())),
               X[features].columns.tolist(), rotation=60, ha='left')
    plt.xlabel("Feature")
    plt.ylabel("Principal components")
    plt.title('Principal Component Analysis with components:' + str(nc), fontsize = 20)



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(result_df[0], result_df[2], result_df[3], c= result_df['PeDamark'])

ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')

plt.show()

sns.pairplot(X)
    
from sklearn.decomposition import PCA, KernelPCA

np.random.seed(0)
X = mdata.drop(['Unnamed: 0', 'Unnamed: 0.1', 'Regime', 'Chem','Domain','Sensitivity', 'Sensitivitybase', 'sensbase', 
                'Variance', 'Anisotropy', 'Time_series', 'timtrace', 'PeDamark','Trial', 'chem_factors'], axis = 1)
X = mdata[features]
cols = X.columns.tolist()
#X['Trial']= X['Trial'].replace('H',0)
#X.Trial = X.Trial.astype(int)
#X.drop(['Trial'], axis = 1, inplace = True)
X = scaler.fit_transform(X)
kpca = KernelPCA(n_components = 4, kernel="rbf", fit_inverse_transform=True, gamma=10)
X_kpca = kpca.fit_transform(X)
X_back = kpca.inverse_transform(X_kpca)
pca = PCA(n_components = 4)
X_pca = pca.fit_transform(X)


plt.figure()
plt.title("Original space")
plt.scatter(mdata['fraction'], mdata['Sensitivity%'], c = mdata['PeDamark'], alpha = 0.5)
plt.xlabel ("x1")
plt.ylabel ("x2")
    
for arrays, labels in zip([X_kpca, X_pca, X_back], ["KPCA", "PCA", "inverse_transformed"]):
#    print(np.shape(arrays))
#    y1 = pd.concat([y, mdata[['PeDamark', 'Senssquared']]], axis = 1)
#    result_df = pd.concat([pd.DataFrame(data = X_pca), y1], axis=1)
#    result_df.head(5)
    compdf = pd.concat([pd.DataFrame(data = arrays), mdata[["PeDamark"]]], axis = 1)
    plt.figure()
    plt.title(labels)
    plt.scatter(compdf[0], compdf[1], s= 40, c = compdf['PeDamark'], alpha = 0.5)
    plt.xlabel ("x1")
    plt.ylabel ("x2")


#SVD
from sklearn.decomposition import TruncatedSVD
svd = TruncatedSVD(n_components = 7)
X_svd = svd.fit_transform(X)

nc = 7
plt.matshow(svd.components_, cmap='viridis')
plt.yticks(list(range(nc)), ["First component", "Second component", "Third component", "Fourth component"])
plt.colorbar()
plt.xticks(range(len(cols)),
              cols, rotation=60, ha='left')
plt.xlabel("Feature")
plt.ylabel("Principal components")
#plt.title('Principal Component Analysis with components:' + str(nc), fontsize = 20)

