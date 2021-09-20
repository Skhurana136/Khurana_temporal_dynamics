# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 14:00:44 2021

@author: khurana
"""

import os
import pandas as pd
import data_reader.data_processing as proc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()

from sklearn.cross_decomposition import CCA
ca = CCA()

raw_dir = "E:/Zenodo_temporal_dynamics"
out_dir = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes"

Regimes = ["Slow", "Equal", "Fast"]
vardict = proc.speciesdict("Saturated")
States = ["Active", "Inactive"]
biomass_vars = list(t for t in vardict.keys() if (vardict[t]["State"] in (States)) and ('sulphate' not in t))
chem_vars = ["DOC", "DO", "TOC", "Ammonium","Nitrogen"]

mf_path = os.path.join(raw_dir, "pca_chem.csv")
biomass_path = os.path.join(raw_dir, "pca.csv")

mf_data= pd.read_csv(mf_path, low_memory = False)
biomass_data = pd.read_csv(biomass_path)

for data in [mf_data, biomass_data]:
    print(data.shape)

X = mf_data[chem_vars+['Variance','Anisotropy','Time_series',"Normalised_vel"]].reset_index()
y = biomass_data[biomass_vars]

X_scaled = scaler.fit_transform(X)
y_scaled = scaler.fit_transform(y)

ca.fit(X_scaled, y_scaled)
X1, y1 = ca.transform(X_scaled, y_scaled)
print(X1.shape)
print(y1.shape)
X1_y1 = pd.DataFrame({"CCX_1": X1[:,0],
                       "CCY_1":y1[:, 0],
                       "CCX_2":X1[:, 1],
                       "CCY_2":y1[:, 1],
                       })

all_res = pd.concat([mf_data[["Regime","Trial","Variance","Anisotropy","Time_series", "Normalised_vel"]].reset_index(),
                              X1_y1], axis = 1)
all_res.head()
cc_res_m = pd.concat([mf_data.reset_index(),X1_y1[["CCX_1", "CCX_2"]]], axis = 1)
cc_res_b = pd.concat([y,X1_y1[["CCY_1", "CCY_2"]]], axis = 1)

np.corrcoef(X1[:, 0], y1[:, 0])

plt.figure(figsize=(10,8))
sns.scatterplot(x="CCX_1",
                y="CCY_1", 
                data=all_res, hue ="Normalised_vel", edgecolor="None")
plt.title('First Pair of Canonical Covariate, corr = %.2f' %
         np.corrcoef(X1[:, 0], y1[:, 0])[0, 1])

plt.figure(figsize=(10,8))
sns.boxplot(x="Regime",
                y="CCX_1", 
               data=all_res)
sns.stripplot(x="Regime",
                y="CCX_1", 
                 data=all_res)

plt.figure(figsize=(10,8))
sns.boxplot(x="Regime",
                y="CCY_1", 
               data=all_res)
sns.stripplot(x="Regime",
                y="CCY_1", 
                 data=all_res)

corr_X_df= cc_res_m.corr(method='pearson') 
corr_X_df.head()
X_df_lt = corr_X_df.where(np.tril(np.ones(corr_X_df.shape)).astype(np.bool))

plt.figure(figsize=(8,6))
sns.heatmap(X_df_lt,cmap="coolwarm",annot=True,fmt='.1f')
plt.tight_layout()
plt.savefig("Heatmap_Canonical_Correlates_from_X_and_data.jpg",
                    format='jpeg',
                    dpi=300)

corr_y_df= cc_res_b.corr(method='pearson') 
corr_y_df.head()
y_df_lt = corr_y_df.where(np.tril(np.ones(corr_y_df.shape)).astype(np.bool))

plt.figure(figsize=(10,8))
sns.heatmap(y_df_lt,cmap="coolwarm",annot=True,fmt='.1f')
plt.tight_layout()
plt.savefig("Heatmap_Canonical_Correlates_from_y_and_data.jpg",
                    format='jpeg',
                    dpi=300)

all_res_het = all_res[all_res.Variance != 0]
g = sns.FacetGrid(all_res_het, col = "Variance", row = "Regime")
g.map_dataframe(sns.scatterplot, x = "CCX_1",
                y = "CCY_1", hue = "Normalised_vel",
                style = "Anisotropy",
                cmap = "flare",
                alpha = 0.6,
                edgecolor = "None")
g.add_legend()

#Ridge Regression
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y_scaled, random_state=42
)

lmr = Ridge()
lmr.fit(X_train, y_train)
y_pred = lmr.predict(X_train)

mae = np.sum(y_train - y_pred)
string_score = f'MAE on training set: {mae:.2f}'
plt.scatter(y_train, y_pred)

y_pred = lmr.predict(X_test)
mae = np.sum(y_test - y_pred)
string_score += f'\nMAE on testing set: {mae:.2f}'
fig, ax = plt.subplots(figsize=(5, 5))
plt.scatter(y_test, y_pred)
ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls="--", c="red")
plt.text(string_score)
plt.title('Ridge model, small regularization')
plt.ylabel('Model predictions')
plt.xlabel('Truths')