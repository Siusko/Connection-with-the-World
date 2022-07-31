import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

###Простой вызов файла
data= pd.read_csv("DefaultDataset.csv", delimiter=",")
data.columns = ["AAPL_x","AAPL_y"]
fig = px.line(data, x = 'AAPL_x', y = 'AAPL_y', log_x=True) #по x log шкала
fig.show()

###Вызов указывающий путь к файлу
data = pd.read_csv('/home/yesius/Documents/ФЛЕШКА/КОНФЕРЕНЦИИ И СТАТЬИ/2022/Ковтун/a.csv')
data.columns = ["AAPL_x","AAPL_y"]
fig = px.line(data, x = 'AAPL_x', y = 'AAPL_y') #по x log шкала
fig.show()

###открыть csv как массив
import numpy as np
with open("/home/yesius/Documents/ФЛЕШКА/КОНФЕРЕНЦИИ И СТАТЬИ/2022/Ковтун/a.csv") as file_name:
    data = np.loadtxt(file_name, delimiter=",")
#назначаем x и y
x=data[:,0]
y=data[:,1]
len(x)
from scipy.interpolate import interp1d
f= interp1d(x, y)
f2 = interp1d(x, y, kind='cubic')
xnew = np.linspace(data[0,0], data[46,0], num=10000, endpoint=True)

###Применение фильтра Савицкого-Голая к уже интерполированным данным
import matplotlib.pyplot as plt

plt.plot(xnew, f2(xnew))

from scipy.signal import savgol_filter
yhat = savgol_filter(f2(xnew),2000, 3) # window size 2000, polynomial order 3

plt.plot(x,y)
plt.plot(xnew,yhat, color='red')
plt.show()
