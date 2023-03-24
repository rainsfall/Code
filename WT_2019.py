import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from pmdarima import auto_arima
from joblib import Parallel, delayed
from statsmodels.tsa.seasonal import STL
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.stattools import arma_order_select_ic
from statsmodels.tsa.arima_model import ARIMA
from pmdarima.model_selection import train_test_split
df = pd.read_csv('C:/Users/Onlyone/OneDrive/桌面/数据/风电数据/2019_lianyungang.csv', parse_dates=['Time'], index_col='Time')


x = df.index
WT1_d = df['WT1'].resample('D',closed = 'right',label = 'right').mean()
WT1_h = df['WT1'].resample('H',closed = 'right',label = 'right').mean()
WT1_q = df['WT1'].resample('15T',closed = 'right',label = 'right').mean()
WT1_normal = df['WT1']
'''绘制数据图像'''
'''
fig, axs = plt.subplots(nrows=3, figsize=(24, 15))
axs[0].plot(x, df['WT1'],color = (46/255,89/255,167/255))
axs[1].plot(x, df['WT2'],color = (242/255,200/255,103/255))
axs[2].plot(x, df['WT3'],color = (237/255,109/255,61/255))
axs[0].set_title('WT1')
axs[1].set_title('WT2')
axs[2].set_title('WT3')
fig.suptitle('风机功率')
plt.show()
'''
'''进行平稳性检验
# ADF平稳性检验
def check_stationarity(column):
    result = sm.tsa.stattools.adfuller(column)
    p_value = result[1]
    return p_value
# 并行计算多个列的p值
p_values = Parallel(n_jobs=-1)(delayed(check_stationarity)(df[column]) for column in df.columns)
# 将结果整理成一个DataFrame
results = pd.DataFrame({'Column': df.columns, 'P-Value': p_values})
# 打印结果
print(results)
'''
'''
# 绘制ACF与PACF图
fig, ax = plt.subplots(8,figsize=(20,15))
fig.subplots_adjust(hspace = 2)
plot_acf(WT1_d, ax=ax[0],color = (46/255,89/255,167/255))
ax[0].set_title('ACF(Sampled a day)')
plot_pacf(WT1_d, ax=ax[1],color = (46/255,89/255,167/255))
ax[1].set_title('PACF(Sampled a day)')
plot_acf(WT1_h, ax=ax[2],color = (242/255,200/255,103/255))
ax[2].set_title('ACF(Sampled an hour)')
plot_pacf(WT1_h, ax=ax[3],color = (242/255,200/255,103/255))
ax[3].set_title('PACF(Sampled an day)')
plot_acf(WT1_q, ax=ax[4],color = (237/255,109/255,61/255))
ax[4].set_title('ACF(Sampled a 15 min)')
plot_pacf(WT1_q, ax=ax[5],color = (237/255,109/255,61/255))
ax[5].set_title('PACF(Sampled a 15 min)')
plot_acf(WT1_normal, ax=ax[6],color = (108/255,168/255,175/255))
ax[6].set_title('ACF(Sampled a 5 min)')
plot_pacf(WT1_normal, ax=ax[7],color = (108/255,168/255,175/255))
ax[7].set_title('PACF(Sampled a 5 min)')
'''

'''白噪声检验'''
'''
res_d = acorr_ljungbox(WT1_d,lags = [6,12,24],return_df = True)
res = acorr_ljungbox(WT1_normal,lags = [6,12,24],return_df = True)
'''

'''STL分解'''
'''
res_d = STL(WT1_d).fit()
res_normal = STL(WT1_normal,period = 288).fit()
fig, ax = plt.subplots(8,figsize=(20,15))
fig.subplots_adjust(hspace = 2)
res_d.plot()
res_normal.plot()
plt.show()
'''

'''数据分割'''
train, test = train_test_split(WT1_d, train_size=0.8)
'''模型定阶'''
model = auto_arima(WT1_d, seasonal=True, m=30)
print(model.summary())

'''模型拟合'''
fitted = model.predict_in_sample()

'''效果展示'''
fig, axs = plt.subplots(2, 2,figsize = (20,16))

model.resid.plot(ax=axs[0][0])
axs[0][0].set_title('residual')

model.resid.plot(kind='hist', ax=axs[0][1])
axs[0][1].set_title('histogram')

sm.qqplot(model.resid, line='45', fit=True, ax=axs[1][0])
axs[1][0].set_title('Normal Q-Q')

plot_acf(model.resid, ax=axs[1][1])

plt.show()