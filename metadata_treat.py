import pandas as pd
import math
import numpy as np

df = pd.read_csv(r'C:/Users/jeremie/Desktop/SraRunTable.txt')
df=df[['Run','Age','disease','AvgSpotLen','Bases','sex','disease_stage','Library Name']]
df = df[df['Library Name'].str.endswith('PC')]


df['reads']=((df['Bases']/(df['AvgSpotLen']*2)))
df['reads']=df['reads'].apply(np.floor).astype(int)

df_cancer =df[df['disease']=='COLORECTAL CANCER']
df_cancer.to_csv(r'c:\users\jeremie\desktop\df_cancer_IC.txt',index=None)
print(df_cancer.head())

df_control=df[df['disease']=='CONTROL']
df_control.to_csv(r'c:\users\jeremie\desktop\df_control_IC.txt',index=None)
print(df_control.head())
