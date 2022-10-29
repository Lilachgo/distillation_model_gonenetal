import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit



#readind the data:
foram_data=pd.read_excel(r"C:/Users/gonen/Documents/Academy/research/EM_ISOTOPES/cores/ps009pc_newages.xlsx")
sst_data=pd.read_excel(r"C:/Users/gonen/Documents/Academy/research/EM_ISOTOPES/analysis/sst_from_marriner_2022.xlsx",sheet_name='Sheet2')
jeita_data=pd.read_excel(r"C:/Users/gonen/Documents/Academy/research/EM_ISOTOPES/spleothems/Jeita/Jeita_data_cheng_2015.xlsx", sheet_name='Sheet2')
soreq_data=pd.read_excel(r"C:/Users/gonen/Documents/Academy/research/EM_ISOTOPES/spleothems/Soreq/soreq_data_noaa.xlsx", sheet_name='Bar mathews 2003')
katal_data=pd.read_excel(r"C:/Users/gonen/Documents/Academy/research/EM_ISOTOPES/spleothems/Katalekhor/katal_data.xlsx")
kunaba_data=pd.read_excel(r"C:/Users/gonen/Documents/Academy/research/EM_ISOTOPES/spleothems/Kuna Ba/aax6656_table_s1.xlsx", sheet_name='data')
stack=pd.read_excel(r"C:\Users\gonen\Documents\Academy\research\EM_ISOTOPES\cores\_.xlsx")

#comprehensive model function:
def f_calc(cave_data,core_data, sst_data,temp):
    #linear intepolation of each of the 3 records:
    x_core=core_data["age"].to_numpy()
    y_core=core_data["d18o"].to_numpy()
    x_cave=cave_data["age"].to_numpy()
    y_cave=cave_data["d18o"].to_numpy()
    x_sst=sst_data["age"].to_numpy()
    y_sst=sst_data["sst"].to_numpy()
    f1=interp1d(x_core,y_core,bounds_error=False)
    f2=interp1d(x_cave,y_cave,bounds_error=False)
    f3=interp1d(x_sst,y_sst, bounds_error=False) 
    xnew = np.arange(0,10000,10)
    a=f1(xnew)
    b=f2(xnew)
    c=f3(xnew)
    interpolated_data=pd.DataFrame(columns=('age','core_d18o','cave_d18o'))
    interpolated_data['age']=np.arange(0,10000,10)
    interpolated_data['core_d18o']=a
    interpolated_data['cave_d18o']=b
    interpolated_data['sst']=c
    interpolated_data=interpolated_data.to_numpy() #0=age 1=core_d18o 2=cave_d18o 3=sst
    
    results=np.empty((0,10))
    for i in range(len(interpolated_data)):
        with_data=np.zeros((1,10)) 
        with_data[0,0]=interpolated_data[i][0]
        sst=(interpolated_data[i][3]+273.15) #sst in age in Kelvin
        with_data[0,3]=interpolated_data[i][3]
        t_difference=295-sst
        t_cave=temp-t_difference #modern t in the cave (K) -suitable for sensitivity A.
        t=t_cave-6
        alfa_caco3=math.exp(((18.03*(1000/t_cave))-32.42)/1000) #after kim and oneil 1997
        alfa_sst=math.exp((1137/(sst**2))-(0.4156/sst)-0.00207)
        alfa_eq=math.exp((1137/(t**2))-(0.4156/t)-0.00207)
        d_foram=interpolated_data[i][1]
        with_data[0,1]=interpolated_data[i][1]
        d_cave=interpolated_data[i][2]
        with_data[0,2]=interpolated_data[i][2]
        d_vapor=((((d_foram+((sst-273.15)-16.5)/4.8)+0.27)+1000)/alfa_sst)-1000-8.74*(1-0.51)
        f=((((1.03086*d_cave+30.86)+1000)/alfa_caco3)/(alfa_eq*((d_foram+(((((sst-273.15)-16.5)/4.8)+0.27)+1000)/alfa_sst)-8.74*(1-0.51))))**(1/(alfa_eq-1))
        with_data[0,4]=f    
        with_data[0,5]=d_vapor
        with_data[0,6]=t_cave-273.15
        with_data[0,7]=t-273.15
        with_data[0,8]=t_difference
        with_data[0,9]=interpolated_data[i][2]-interpolated_data[i][1]
        results=np.append(results,with_data, axis=0)
    results=pd.DataFrame(results,columns=('age','core_d18o','cave_d18o','sst','f','d_vapor','cave_t','eq_t','t_dif','Delta'))  
    return results
    
#applying the model to selected caves:
jeita=f_calc(jeita_data,foram_data,sst_data,290)
soreq=f_calc(soreq_data,foram_data,sst_data,286) 
katal=f_calc(katal_data,foram_data,sst_data,284) 
kunaba=f_calc(kunaba_data,foram_data,sst_data,286) 

#putting togther results of all the caves and applying gaussian smoothing:
all_caves=pd.DataFrame()
all_caves['age']=jeita['age']
all_caves['jeita_comp']=jeita['f']
all_caves['jeita_simp']=jeita['Delta']
all_caves['soreq_comp']=soreq['f']
all_caves['soreq_simp']=soreq['Delta']
all_caves['katalekhor_comp']=katal['f']
all_caves['katalekhor_simp']=katal['Delta']
all_caves['kuna_ba_comp']=kunaba['f']
all_caves['kuna_ba_simp']=kunaba['Delta']

smoothed=all_caves.copy()
smoothed['Jeita_comp']=smoothed['jeita_comp'].rolling(5,win_type='gaussian',center=True).mean(std=0.01)
smoothed['Jeita_simp']=smoothed['jeita_simp'].rolling(5,win_type='gaussian',center=True).mean(std=0.01)
smoothed['Soreq_comp']=smoothed['soreq_comp'].rolling(5,win_type='gaussian',center=True).mean(std=0.01)
smoothed['Soreq_simp']=smoothed['soreq_simp'].rolling(5,win_type='gaussian',center=True).mean(std=0.01)
smoothed['Katalekhor_comp']=smoothed['katalekhor_comp'].rolling(5,win_type='gaussian',center=True).mean(std=0.01)
smoothed['Katalekhor_simp']=smoothed['katalekhor_simp'].rolling(5,win_type='gaussian',center=True).mean(std=0.01)
smoothed['Kuna_ba_comp']=smoothed['kuna_ba_comp'].rolling(5,win_type='gaussian',center=True).mean(std=0.01)
smoothed['Kuna_ba_simp']=smoothed['kuna_ba_simp'].rolling(5,win_type='gaussian',center=True).mean(std=0.01)
d=np.arange(0,10000,50)
d=pd.Series(d)
smoothed['t']=smoothed['age'].round().isin(d.round())
smoothed=smoothed.drop(['jeita_comp','jeita_simp','soreq_comp','soreq_simp','katalekhor_comp','katalekhor_simp','kuna_ba_comp','kuna_ba_simp','t'],axis=1)

##plots:
#results plots (uncertainties are claulted seperatly):

fig, ax = plt.subplots(4, 1, sharex=True)
params={'font.size':24,'figure.figsize':(15,18)}
fig.subplots_adjust(hspace=0)
plt.rcParams.update(params)
ax[0].plot(smoothed['age'],smoothed['Jeita_comp'], color='forestgreen',label='comprehensive model')
#ax[0].invert_yaxis()
ax[0].set_ylim(0.89,0.67)
error1=0.04
ax[0].fill_between(smoothed['age'],smoothed['Jeita_comp']-error1, smoothed['Jeita_comp']+error1, color='forestgreen', alpha=0.2)
ax0 = ax[0].twinx()
ax0.plot(smoothed['age'], smoothed['Jeita_simp'],color='midnightblue',label='simplistic model')
#ax0.invert_yaxis()
ax0.set_ylim(-3.5,-7)
ax0.set_yticks([-7,-6,-5,-4])
error2=0.38
ax0.fill_between(smoothed['age'], smoothed['Jeita_simp']-error2, smoothed['Jeita_simp']+error2, color='midnightblue', alpha=0.2)
ax[0].set_ylabel('f')
ax0.set_ylabel(r'$\Delta^{18}O$ ‰')
ax[0].legend(frameon=False,loc='lower left')
ax0.legend(frameon=False, loc='upper right')

ax[1].plot(foram_data['age'],foram_data['d18o'],color='steelblue',label='ps009pc')
ax[1].invert_yaxis()
ax[1].plot(stack['age'],stack['stack'],color='navy',label='stack')
ax[1].legend(frameon=False,loc='lower left')
ax[1].set_ylabel('$\delta^{18}O$ \n (‰ VPDB)')
#ax1= ax[1].twinx()

#ax1.invert_yaxis()
#ax1.legend(frameon=False,loc='upper right')

ax[2].plot(jeita_data['age'],jeita_data['d18o'],color='black',label='Jeita')
ax[2].invert_yaxis()
ax[2].legend(frameon=False,loc='lower left')
ax[2].set_ylabel('$\delta^{18}O$ \n (‰ VPDB)')

ax[3].set_xlim(0,13000)
ax[3].set_xticks([500,1000,1500,2500,3000.3500,4500,5000,5500,6500,7000,7500,8500,9000,9500,10500,11000,11500,12500],minor=True)
ax[3].set_xlabel('Years BP')
ax[3].plot(sst_data['age'],sst_data['East_Med_T'],color='goldenrod',label='sst (Marriner et al. 2022)')
ax[3].legend(frameon=False,loc='lower left')
ax[3].set_ylabel(' deviation from average \n T\N{degree sign}C')
ax[3].invert_yaxis()
plt.show()    
    
#plotting sptial results of the comprehensive model:
fig, ax = plt.subplots()
ax.set_xticks([500,1000,1500,2500,3000.3500,4500,5000,5500,6500,7000,7500,8500,9000,9500],minor=True)
ax.set_xlabel('Years BP',fontsize=24)
params={'font.size':24,'figure.figsize':(18,10)}
plt.rcParams.update(params)
ax.plot(smoothed['age'], smoothed['Jeita_comp'], color='midnightblue',label='Jeita',linewidth=2) 
ax.plot(smoothed['age'], smoothed['Soreq_comp'], color='mediumblue',label='Soreq',linewidth=2)                    
ax.plot(smoothed['age'], smoothed['Kuna_ba_comp'], color='royalblue',label='Kuna Ba',linewidth=2)
ax.plot(smoothed['age'], smoothed['Katalekhor_comp'], color='dodgerblue',label='Katalekhor',linewidth=2)
ax.set(xlim=(0,10000))
ax.invert_yaxis()
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylabel('f')

#plotting sptial results of the simplistic model:
fig, ax = plt.subplots()
ax.set_xticks([500,1000,1500,2500,3000.3500,4500,5000,5500,6500,7000,7500,8500,9000,9500],minor=True)
ax.set_xlabel('Years BP',fontsize=24)
params={'font.size':24,'figure.figsize':(18,10)}
plt.rcParams.update(params)
ax.plot(smoothed['age'], smoothed['Jeita_simp'], color='darkolivegreen',label='Jeita',linewidth=2) 
ax.plot(smoothed['age'], smoothed['Soreq_simp'], color='green',label='Soreq',linewidth=2)                    
ax.plot(smoothed['age'], smoothed['Kuna_ba_simp'], color='mediumseagreen',label='Kuna Ba',linewidth=2)
ax.plot(smoothed['age'], smoothed['Katalekhor_simp'], color='springgreen',label='Katalekhor',linewidth=2)
ax.set(xlim=(0,10000))
ax.invert_yaxis()
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylabel(r'$\Delta^{18}O$ ‰')





