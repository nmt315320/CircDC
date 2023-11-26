import shap
import xgboost
import matplotlib.pyplot as plt
import pandas as pd

#df=pd.read_csv('rcm_features.csv')
y=df['label']
X=df.drop('label',axis=1)
X,y=shap.datasets.boston()
model=xgboost.train({"learning_rate":0.01},xgboost.DMatrix(X,label=y),100)
explainer=shap.TreeExplainer(model)
shap_values=explainer.shap_values(X)
shap.force_plot(explainer.expected_value,shap_values[0,:],X.iloc[0,:],show=False)
plt.savefig('rcm/cun.png')
#plt.clear()
shap.force_plot(explainer.expected_value,shap_values,X,show=False)
#shap.force_plot(explainer.expected_value[0], shap_values_1[0], data_arry, matplotlib=True,show=False)
plt.savefig('rcm/R2.png')
#plt.clear()
shap.dependence_plot("RM",shap_values,X,show=False)
plt.savefig('rcm/R1.png')
shap.summary_plot(shap_values,X,show=False)
plt.savefig('rcm/R3.png')
#plt.clear()
shap.summary_plot(shap_values,X,plot_type="bar",show=False)
plt.savefig('rcm/R4.png')
#plt.clear()