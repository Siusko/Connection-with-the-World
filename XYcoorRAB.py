# %% markdown
# #### Importing Python libraries
#
# The next cell activates the Matplotlib plotting library and loads several useful libraries
# %% codecell
#  displays plots directly in this notebook as images
%matplotlib inline
import numpy as np       # the NumPy library for fast matrix and vector data types and math operations
from numpy import sqrt, sin, cos, pi,e, arctan
import matplotlib.pyplot as plt   # functions for plotting, loaded under a convenient name alias 'plt'
import scipy
from scipy import constants
import math
import pandas as pd

# %% markdown
# #### CONSTANTS
# %% codecell
f=26*10**9 # частота зондирования
f1=37*10**9
f2=26*10**9
R=0.34 # радиус зондируемой плазмы
b=R+0.1 # к радиусу плазмы прибавляем размер антенны
pi=constants.pi
me=constants.m_e # electron mass
c=constants.speed_of_light # speed of light in vacuum
e0=constants.epsilon_0 # the electric constant (vacuum permittivity)
q=constants.e  # elementary charge
γ=3 # Степень функции
λ=constants.nu2lambda(f) # перевод частоты в длинну волны
ω=f*2*pi # угловая частота
Nc=me*(pi**2)*4*(f**2)*e0*(10**-6)/(q**2)#Плотность в зависимомти от частоты
Nc1=me*(pi**2)*4*(f1**2)*e0*(10**-6)/(q**2)
Nc2=me*(pi**2)*4*(f2**2)*e0*(10**-6)/(q**2)
Ncconst=me*(pi**2)*4*((37*10**9)**2)*e0*(10**-6)/(q**2)#меняем плотность расчитав для двух частот
Nm=0.8*Ncconst#меняем плотность расчитав для двух частот
ωm=math.sqrt(Nm*(q**2)/(me*e0))

numbers = [Nc1,Nc2,Ncconst,Nc1/Nc2] # распечатать все константы в столбец (in scientific notation)
for x1 in numbers:
    print("{:e}".format(x1))
Nc
# %% codecell
Nc
# %% markdown
# #### VARIABALES
# %% codecell
k=[0,1,2,3,4,5,6,7,8,9,10,11]
# %% codecell

# %% codecell
α1=[]
for n in range(len(k)) :
    a1=(n+1)*pi/180
    α1.append(a1)
α1
# %% codecell
α2=[]
for n in range(len(k)) :
    if b*(math.sin(a1))<R:
        a2=pi-math.asin(b*math.sin((n+1)*pi/180)/R)
        α2.append(a2)
    elif R>=b:
        a2=math.asin(b*math.sin((n+1)*pi/180)/R)
        α2.append(a2)
    else:
        α2.append(0)
α2

# %% codecell
k
# %% codecell
βk=[]
for n in range(len(k)) :
    β=pi-α1[n]-α2[n]
    βk.append(β)
βk
# %% codecell
θk=[]
for n in range(len(k)) :
    θ=pi-βk[n]
    θk.append(θ)
θk
# %% codecell
Φk=[]
for n in range(len(k)) :
    Φ=α1[n]+βk[n]
    Φk.append(Φ)
Φk
# %% codecell

# %% markdown
# import sympy
# from sympy import symbols, solve
# ak=[]
# for n in range(len(k)) :
#     r0 = symbols('r0')
#     expr = 1-(Nm*(1-(r0/R)**2)/Nc)-((R**2)/(r0**2))*(math.sin(Φk[n]))**2
#     sol = solve(expr)
#     ak.append(sol)
# ak1=[]
# for i in range(len(k)) :
#     ak[i] = [n for n in ak[i] if n.is_real]
#     ak[i] = [n for n in ak[i] if n.is_positive]
#     ak1.append(ak[i])
# ak1
# %% codecell
from scipy.optimize import fsolve
kak=[]
r00=0.0001 #The starting estimate for the roots of func(x) = 0.
for n in range(len(k)) :
    kak.append(fsolve(lambda r0: 1-(Nm*(1/(1+9*(r0/R)**γ))/Nc)-((R**2)/(r0**2))*(math.sin(Φk[n]))**2, r00))
ak = np.array(kak, dtype=np.float32)
# %% codecell
ak
# %% codecell

# %% codecell

# %% codecell
#result = integrate.quad(lambda x: x**7, 2, 0), где lambda x по сути dx, x**7 - это подинтегральная функция,
#2 - нижний придел, 0 - верхний придел интегрирования, result = integrate.quad - по сути знак интеграла (метод)


import scipy.integrate as integrate #вызываем метод интегрирования
import scipy.special as special
Θk=[]
for n in range(len(k)) :
    Θ = integrate.quad(lambda r: R*sin(Φk[n])/(r**2*((1-(Nm*(1/(1+9*(r/R)**γ))/Nc)-((R**2)/(r**2))*(math.sin(Φk[n]))**2)**0.5)), ak[n], R)
    Θk.append(Θ[0]) #вписываем в массив Θk только н1-е корни каждого
Θk
# %% codecell

# %% codecell
k
# %% codecell
θ2k=[]
for n in range(len(k)) :
    θ2=θk[n]-np.abs(Θk[n]*2)
    θ2k.append(θ2)
θ2k
# %% codecell

R1=[]
if 1-Nc/Nm>0:
    R1=R*(((Nc-Nm)/(9*Nm))**(1/γ))
elif (1-Nc/Nm)<=0:
    R1=0
R1=abs(R1)
R1

# %% codecell
k[0]
# %% codecell
from scipy.integrate import odeint

number=[0,1,2,3,4,5,6,7,8,9,10,11]

for k in range(len(number)) :
    def model(r,θ):
        drdθ =(r**2)*((1-(Nm*(1/(1+9*(r/R)**γ))/Nc)-((R**2)/(r**2))*(math.sin(Φk[number[k]]))**2)**0.5)/(R*sin(Φk[number[k]]))
        return drdθ
    r0 = R
    θ = np.linspace(θk[number[k]],θ2k[number[k]],1000)
    θ0=[]
    for n in range(len(θ)) :
            if θ[n]>=θk[number[k]]-(np.abs(Θk[number[k]]))+0.1*pi/180:
                θ0.append(θ[n])
                r= odeint(model,r0,θ0)
           # else:
                #r=0
    r = [i[0] for i in r]
    combined1=np.vstack((θ0,r[::])).T

    def model(r,θ):
        drdθ =(r**2)*((1-(Nm*(1/(1+9*(r/R)**γ))/Nc)-((R**2)/(r**2))*(math.sin(Φk[number[k]]))**2)**0.5)/(R*sin(Φk[number[k]]))
        return drdθ
    r0 = R
    θ = np.linspace(θk[number[k]],θ2k[number[k]],1000)
    θ0=[]
    for n in range(len(θ)) :
            if θ[n]<=θ2k[number[k]]+(np.abs(Θk[number[k]]))-0.1*pi/180:
                θ0.append(θ[n])
                r= odeint(model,r0,θ0)
            #else:
                #r=0
    r = [i[0] for i in r]
    combined2=np.vstack((θ0,r[::-1])).T

    r=np.concatenate((combined1, combined2 ))
#θ=np.concatenate((θ01, θ02))

    #plt.plot(r[:,0],r[:,1])
    α = np.linspace(0,pi*2,1000)
    x=r[:,1]*cos(r[:,0])
    y=r[:,1]*sin(r[:,0])

    #polar coordinate metod1
    # Calculating radius
    #RR = (x**2 + y**2)**0.5
    # Calculating angle (theta) in radian
    #theta = arctan(y/x)
    # Converting theta from radian to degree
    #theta = 180 * theta/math.pi

    #polar coordinate metod2
    #RR = r[:,1]
    #theta = r[:,0]
    #fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    #ax.plot(theta, RR)
    #np.savetxt("r polar"+str(number[k])+".txt",r)

    x2=R*sin(α)
    y2=R*cos(α)
    x3=R1*sin(α)
    y3=R1*cos(α)
    for n in range(1000) :
        if R1!=0:#не равно нулю
            xCENTER = np.linspace(-R,-R1,1000)
            yCENTER=α*0
        else:
            xCENTER = np.linspace(-R,R,1000)
            yCENTER=α*0
    plt.plot(x,y)
    plt.plot(x,-y)
    plt.plot(xCENTER,yCENTER)

    combined3=np.vstack((x,y)).T
    combined4=np.vstack((x,-y)).T
    circles=np.vstack((x2,y2,x3,y3,xCENTER,yCENTER)).T
    np.savetxt("r"+str(number[k])+".txt",combined3)
    np.savetxt("-r"+str(number[k])+".txt",combined4)
np.savetxt("circles.txt",circles)
plt.plot(R*sin(α),R*cos(α))
plt.plot(R1*sin(α),R1*cos(α))








    #new = pd.concat([bbb, bbb], axis=1)

    #aaa= pd.read_csv("r"+str(number[k+1]), delimiter=" ")
    #aaa.columns = ["x"+str(number[k+1]), "y"+str(number[k+1])]
    #new = pd.concat([bbb, aaa], axis=1)
    #print(bbb.columns = ["x"+str(number[k]), "y"+str(number[k])])
    #fig = px.line(new, new["x10"],new["y10"])
    #fig.show()

#np.savetxt("circles",circles)
#plt.plot(R*sin(α),R*cos(α))
#plt.plot(R1*sin(α),R1*cos(α))


#fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
#ax.plot(r[:,0], r[:,1])
#plt.show()
k=[0,1,2,3,4,5,6,7,8,9,10,11]


# %% codecell
R1
# %% codecell
#ВЫВОД РЕЗУЛЬТАТОВ В 1 файл
bbb0 = pd.read_csv("r0.txt", delimiter=" ")
bbb1 = pd.read_csv("r1.txt", delimiter=" ")
bbb2 = pd.read_csv("r2.txt", delimiter=" ")
bbb3 = pd.read_csv("r3.txt", delimiter=" ")
bbb4 = pd.read_csv("r4.txt", delimiter=" ")
bbb5 = pd.read_csv("r5.txt", delimiter=" ")
bbb6 = pd.read_csv("r6.txt", delimiter=" ")
bbb7 = pd.read_csv("r7.txt", delimiter=" ")
bbb8 = pd.read_csv("r8.txt", delimiter=" ")
bbb9 = pd.read_csv("r9.txt", delimiter=" ")
bbb10 = pd.read_csv("r10.txt", delimiter=" ")
bbb11 = pd.read_csv("r11.txt", delimiter=" ")
#bbb12 = pd.read_csv("r12.txt", delimiter=" ")
bbb13 = pd.read_csv("-r0.txt", delimiter=" ")
bbb14 = pd.read_csv("-r1.txt", delimiter=" ")
bbb15 = pd.read_csv("-r2.txt", delimiter=" ")
bbb16 = pd.read_csv("-r3.txt", delimiter=" ")
bbb17 = pd.read_csv("-r4.txt", delimiter=" ")
bbb18 = pd.read_csv("-r5.txt", delimiter=" ")
bbb19 = pd.read_csv("-r6.txt", delimiter=" ")
bbb20 = pd.read_csv("-r7.txt", delimiter=" ")
bbb21 = pd.read_csv("-r8.txt", delimiter=" ")
bbb22 = pd.read_csv("-r9.txt", delimiter=" ")
bbb23 = pd.read_csv("-r10.txt", delimiter=" ")
bbb24 = pd.read_csv("-r11.txt", delimiter=" ")

circles = pd.read_csv("circles.txt", delimiter=" ")

bbb = pd.concat([circles, bbb0, bbb1, bbb2, bbb3, bbb4, bbb5, bbb6, bbb7, bbb8, bbb9, bbb10, bbb11, bbb13, bbb14, bbb15, bbb16, bbb17, bbb18, bbb19, bbb20, bbb21, bbb22, bbb23, bbb24], axis=1)
#bbb.columns = ["Rx","Ry","R1x","R1y","x0", "y0","x1", "y1","x2", "y2","x3", "y3","x4", "y4","x5", "y5","x6", "y6","x7", "y7", "x8","y8", "x9","y9", "x10","y10", "x11","y11"]
bbb.to_excel('raw_data1.xls', index=False)

# %% markdown
# from scipy.integrate import odeint
#
# number=[0,1,2,3,4,5,6,7,8,9,10,11]
#
# for k in range(len(number)) :
#     def model(r,θ):
#         drdθ =(r**2)*((1-(Nm*(1-(r/R)**γ)/Nc)-((R**2)/(r**2))*(math.sin(Φk[number[k]]))**2)**0.5)/(R*sin(Φk[number[k]]))
#         return drdθ
#     r0 = R
#     θ = np.linspace(θk[number[k]],θ2k[number[k]],1000)
#     θ0=[]
#     for n in range(len(θ)) :
#             if θ[n]>=θk[number[k]]-(np.abs(Θk[number[k]]))+0.1*pi/180:
#                 θ0.append(θ[n])
#                 r= odeint(model,r0,θ0)
#            # else:
#                 #r=0
#     r = [i[0] for i in r]
#     combined1=np.vstack((θ0,r[::])).T
#
#     def model(r,θ):
#         drdθ =(r**2)*((1-(Nm*(1-(r/R)**γ)/Nc)-((R**2)/(r**2))*(math.sin(Φk[number[k]]))**2)**0.5)/(R*sin(Φk[number[k]]))
#         return drdθ
#     r0 = R
#     θ = np.linspace(θk[number[k]],θ2k[number[k]],1000)
#     θ0=[]
#     for n in range(len(θ)) :
#             if θ[n]<=θ2k[number[k]]+(np.abs(Θk[number[k]]))-0.1*pi/180:
#                 θ0.append(θ[n])
#                 r= odeint(model,r0,θ0)
#             #else:
#                 #r=0
#     r = [i[0] for i in r]
#     combined2=np.vstack((θ0,r[::-1])).T
#
#     r=np.concatenate((combined1, combined2 ))
# #θ=np.concatenate((θ01, θ02))
#
#     #plt.plot(r[:,0],r[:,1])
#     α = np.linspace(0,pi*2,1000)
#     x=r[:,1]*cos(r[:,0])
#     y=r[:,1]*sin(r[:,0])
#
#     #polar coordinate metod1
#     # Calculating radius
#     #RR = (x**2 + y**2)**0.5
#     # Calculating angle (theta) in radian
#     #theta = arctan(y/x)
#     # Converting theta from radian to degree
#     #theta = 180 * theta/math.pi
#
#     #polar coordinate metod2
#     RR = r[:,1]
#     theta = r[:,0]
#     fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
#     ax.plot(theta, RR)
#     np.savetxt("r polar"+str(number[k])+".txt",r)
#
#     x2=R*sin(α)
#     y2=R*cos(α)
#     x3=R1*sin(α)
#     y3=R1*cos(α)
#     #plt.plot(x,y)
#     combined3=np.vstack((x,y)).T
#     circles=np.vstack((x2,y2,x3,y3)).T
#     np.savetxt("r"+str(number[k])+".txt",combined3)
# np.savetxt("circles.txt",circles)
# #plt.plot(R*sin(α),R*cos(α))
# #plt.plot(R1*sin(α),R1*cos(α))
#
#
#
#
#
#
#
#
#     #new = pd.concat([bbb, bbb], axis=1)
#
#     #aaa= pd.read_csv("r"+str(number[k+1]), delimiter=" ")
#     #aaa.columns = ["x"+str(number[k+1]), "y"+str(number[k+1])]
#     #new = pd.concat([bbb, aaa], axis=1)
#     #print(bbb.columns = ["x"+str(number[k]), "y"+str(number[k])])
#     #fig = px.line(new, new["x10"],new["y10"])
#     #fig.show()
#
# #np.savetxt("circles",circles)
# #plt.plot(R*sin(α),R*cos(α))
# #plt.plot(R1*sin(α),R1*cos(α))
#
#
# #fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# #ax.plot(r[:,0], r[:,1])
# #plt.show()
# k=[0,1,2,3,4,5,6,7,8,9,10,11]
# %% codecell

# %% codecell

# %% codecell

# %% codecell


# %% markdown
#
# %% codecell

# %% codecell

# %% codecell

# %% codecell

# %% codecell

# %% codecell
