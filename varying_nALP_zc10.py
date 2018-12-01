
# coding: utf-8

# In[ ]:

# import necessary modules
#get_ipython().magic(u'matplotlib inline')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
from matplotlib import rc
import matplotlib.patches as patches

from scipy.interpolate import interp1d
from matplotlib.ticker import FixedLocator
from math import floor

# In[ ]:

# esthetic definitions for the plots

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
# matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]

l_TT,err_TT= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/PlanckerrorsbinnedTT.dat",unpack=True)
l_EE,err_EE= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/PlanckerrorsbinnedEE.dat",unpack=True)
lmin_phiphi,lmax_phiphi,cl_phiphi,err_phiphi= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/agressive_lensing.csv",unpack=True)

# In[ ]:

############################################
#
# Varying parameter (others fixed to default)
#
var_name1 = 'n_pheno_axion'
# var_name2 = 'Omega_fld_ac'
var_name2 = 'Omega_many_fld'
# var1_min = 1
# var1_max = 3
# var2_min = 1e14
# n  = [1,2,3]
n  = [1,1,1]
# n = [2,2]
# Omega_axion_ac = [1e14,5e14,1e15]
# Omega_axion_ac = [50,50,50]
# Omega_axion_ac = [3,4.5,7]
# Omega_many_fld = [3.1e-3,1.25e-5,0.4e-6]
# Omega_many_fld = [6.167e-3,11.95e-4,4.3e-4]
# Omega_many_fld = [0.02,4.289e-06,15.8e-9]######95%CL ac=1em5
# Omega_many_fld = [4.289e-06,4.289e-06]######95%CL ac=1em5
# Omega_many_fld = [3.1e-3,1.25e-5,0.3993e-06]######95%CL ac=1em3
# Omega_many_fld = [6.15e-3,1.2e-3,4.3e-4]#######95%CL ac=10
Omega_many_fld = [6.15e-3,3.1e-3,0.02]
# Omega_axion_ac = [2.278e+13,4.278e+14,1.8e+15]
a_c = [1e-1,1e-3,1e-5]
color = ['r','b--','g-.']
var_num = 3
# var_legend = r'$n_{\rm a}$'
legend_name= [r'$a_c =0.1$',r'$a_c =10^{-3}$',r'$a_c =10^{-5}$']
var_figname = 'nalp'
dashes_array = [[10,2,10,2],[3,1,3,1],[5,2,1,2]]
#
#############################################
#
# Fixed settings
#
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'format':'camb',
                   # LambdaCDM parameters
                   # 'h':0.67556,
                   '100*theta_s':1.04077,
                   'omega_b':0.02225,
                   # 'omega_cdm':0.12038,
                   # 'Omega_cdm':0.266,
                   'Omega_cdm':0.266,
                   # 'A_s':2.215e-9,
                   'ln10^{10}A_s':3.094,
                   'n_s':0.9645,
                   'tau_reio':0.079,
                   'N_ur':3.046,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   #'YHe':0.246,
                   # other output and precision parameters
                   'P_k_max_1/Mpc':3.0,
                   'l_switch_limber':9}
                   #,'temperature contributions':'tsw'}
                   #'background_verbose':1}
#
# arrays for output
#
kvec = np.logspace(-4,np.log10(3),1000)
legarray = []
twopi = 2.*math.pi
#
# Create figures
#
# fig_Pk, ax_Pk = plt.subplots()
# fig_TT, ax_TT = plt.subplots()
# fig_EE, ax_EE = plt.subplots()
# fig_PP, ax_PP = plt.subplots()

fig_TT, (ax_TT, ax_EE, ax_PP, ax_Pk) = plt.subplots(4,1, sharex=False,gridspec_kw=dict(height_ratios=[1,1,1,1]),figsize=(35,40))
# fig_Pk, ax_Pk = plt.subplots()
fig_TT.subplots_adjust(hspace=.3)
figxe, (axe,axT) = plt.subplots(2,1, sharex=True,gridspec_kw=dict(height_ratios=[1,1]),figsize=(15,13))
figxe.subplots_adjust(hspace=0)

ax_TT.tick_params(axis='x',labelsize=60,length=20,which='both',width=2,direction='inout')
ax_TT.tick_params(axis='y',labelsize=60,length=20,which='both',width=2,direction='inout')
ax_EE.tick_params(axis='x',labelsize=60,length=20,which='both',width=2,direction='inout')
ax_EE.tick_params(axis='y',labelsize=60,length=20,which='both',width=2,direction='inout')
axe.tick_params(axis='x',labelsize=60,length=20,which='both',width=2,direction='inout')
axe.tick_params(axis='y',labelsize=60,length=20,which='both',width=2,direction='inout')
axT.tick_params(axis='x',labelsize=60,length=20,which='both',width=2,direction='inout')
axT.tick_params(axis='y',labelsize=60,length=20,which='both',width=2,direction='inout')
ax_Pk.tick_params(axis='x',labelsize=60,length=20,which='both',width=2,direction='inout')
ax_Pk.tick_params(axis='y',labelsize=60,length=20,which='both',width=2,direction='inout')
ax_PP.tick_params(axis='x',labelsize=60,length=20,which='both',width=2,direction='inout')
ax_PP.tick_params(axis='y',labelsize=60,length=20,which='both',width=2,direction='inout')
ax_TT.tick_params(axis='x',labelsize=60,length=10,which='minor',direction='inout')
ax_TT.tick_params(axis='y',labelsize=60,length=10,which='minor',direction='inout')
ax_EE.tick_params(axis='x',labelsize=60,length=10,which='minor',direction='inout')
ax_EE.tick_params(axis='y',labelsize=60,length=10,which='minor',direction='inout')
axe.tick_params(axis='x',labelsize=60,length=10,which='minor',direction='inout')
axe.tick_params(axis='y',labelsize=60,length=10,which='minor',direction='inout')
axT.tick_params(axis='x',labelsize=60,length=10,which='minor',direction='inout')
axT.tick_params(axis='y',labelsize=60,length=10,which='minor',direction='inout')
ax_Pk.tick_params(axis='x',labelsize=60,length=10,which='minor',direction='inout')
ax_Pk.tick_params(axis='y',labelsize=60,length=10,which='minor',direction='inout')
ax_PP.tick_params(axis='x',labelsize=60,length=10,which='minor',direction='inout')
ax_PP.tick_params(axis='y',labelsize=60,length=10,which='minor',direction='inout')

# ax_TT.axis([2,2500,-0.06,0.06])
# ax_PP.axis([2,2500,-0.06,0.06])
# ax_EE.axis([2,2500,-0.06,0.06])
change = 50.0
# ax_PP.xaxis.set_major_locator(
#     FixedLocator(
#         np.concatenate((np.array([2, 10, change]),
#                         np.arange(500, 2500, 500)))))
# ax_PP.xaxis.set_minor_locator(
#     FixedLocator(
#         np.concatenate((np.arange(2, 10),
#                         np.arange(10, 50, 10),
#                         np.arange(floor(change/100), 2500, 100)))))
# ax_TT.xaxis.set_major_locator(
#     FixedLocator(
#         np.concatenate((np.array([2, 10, change]),
#                         np.arange(500, 2500, 500)))))
# ax_TT.xaxis.set_minor_locator(
#     FixedLocator(
#         np.concatenate((np.arange(2, 10),
#                         np.arange(10, 50, 10),
#                         np.arange(floor(change/100), 2500, 100)))))
# ax_EE.xaxis.set_major_locator(
#     FixedLocator(
#         np.concatenate((np.array([2, 10, change]),
#                         np.arange(500, 2500, 500)))))
# ax_EE.xaxis.set_minor_locator(
#     FixedLocator(
#         np.concatenate((np.arange(2, 10),
#                         np.arange(10, 50, 10),
#                         np.arange(floor(change/100), 2500, 100)))))

#
# loop over varying parameter values
#

M = Class()
M.set(common_settings)
M.compute()
#
# get Cls and Pk for LCDM
#
# clM = M.raw_cl(2500)

all_k = M.get_perturbations()  # this potentially constains scalars/tensors and all k values
conversion = pow(2.7255*1.e6,2)
# one_k = all_k['scalar'][0]     # this contains only the scalar perturbations for the requested k values
th_LCDM = M.get_thermodynamics()
z_LCDM = th_LCDM['z']
xe_LCDM = th_LCDM['x_e']
T_LCDM = th_LCDM['Tb [K]']
axT.loglog(z_LCDM,T_LCDM,'k',lw=5)

clM = M.lensed_cl(2500)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
clPP_LCDM = clM['pp'][2:]
# ax_TT.semilogx(ll_LCDM,(ll_LCDM)*(ll_LCDM+1)/(2*np.pi)*(clTT_LCDM),'k',lw=5)
# ax_EE.semilogx(ll_LCDM,(ll_LCDM)*(ll_LCDM+1)/(2*np.pi)*(clEE_LCDM),'k',lw=5)
# ax_PP.semilogx(ll_LCDM,(ll_LCDM)*(ll_LCDM+1)/(2*np.pi)*(clPP_LCDM),'k',lw=5)

pkM_LCDM = []

for k in kvec:
    pkM_LCDM.append(M.pk(k,0.))

M = Class()
M.set(common_settings)
M.set({'Omega_cdm':0.266+0.02})
M.compute()
clM2 = M.lensed_cl(2500)
ll_moreCDM = clM2['ell'][2:]
clTT_moreCDM = clM2['tt'][2:]
clEE_moreCDM = clM2['ee'][2:]
clPP_moreCDM = clM2['pp'][2:]
ax_TT.semilogx(ll_LCDM,(clTT_moreCDM-clTT_LCDM)/clTT_LCDM,'k',lw=5)
ax_EE.semilogx(ll_LCDM,(clEE_moreCDM-clEE_LCDM)/clEE_LCDM,'k',lw=5)
ax_PP.semilogx(ll_LCDM,(clPP_moreCDM-clPP_LCDM)/clPP_LCDM,'k',lw=5)
pkM = []
for k in kvec:
    pkM.append(M.pk(k,0.))

ax_Pk.semilogx(kvec,(np.array(pkM)-np.array(pkM_LCDM))/np.array(pkM_LCDM),'k',lw=5)
#
# M = Class()
# M.set(common_settings)
# M.set({'N_ur':3.5})
# M.compute()
# clM2 = M.lensed_cl(2500)
# ll_moreUr = clM2['ell'][2:]
# clTT_moreUr = clM2['tt'][2:]
# clEE_moreUr = clM2['ee'][2:]
# clPP_moreUr = clM2['pp'][2:]
# ax_TT.semilogx(ll_LCDM,(clTT_moreUr-clTT_LCDM)/clTT_LCDM,'k',lw=5)
# ax_EE.semilogx(ll_LCDM,(clEE_moreUr-clEE_LCDM)/clEE_LCDM,'k',lw=5)
# ax_PP.semilogx(ll_LCDM,(clPP_moreUr-clPP_LCDM)/clPP_LCDM,'k',lw=5)
# pkM = []
#
# for k in kvec:
#     pkM.append(M.pk(k,0.))
#
# ax_Pk.semilogx(kvec,(np.array(pkM)-np.array(pkM_LCDM))/np.array(pkM_LCDM),'k',lw=5,label='more Ur')



# dashes=([10,2,10,2])
### To plot the binned cosmic variance ####
l_min = 2.;
l_max = 2500.;
n_step = 15.;
j=0.
# def binned_cosmic_variance (result,l_ini,l_tot):
#     for i in range(0,int(l_tot)):
#      result = result + 2/(2*(l_ini+float(i))+1)
#     return np.sqrt(result/l_tot**2)
#
#
# while j < n_step:
#         result = 0.0
#         step = l_min*(l_max/l_min)**(j/n_step)
#         step_plus_1 = l_min*(l_max/l_min)**((j+1)/n_step)
#         print int(step), int(step_plus_1)
#         width = int(step_plus_1) - int(step)
#         ax_TT.add_patch(
#             patches.Rectangle(
#                 (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
#                 width,          # width
#                 2*binned_cosmic_variance(result,int(step),width),          # height
#                 color='r',
#                 alpha=0.1
#             )
#         )
#         ax_EE.add_patch(
#             patches.Rectangle(
#                 (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
#                 width,          # width
#                 2*binned_cosmic_variance(result,int(step),width),          # height
#                 color='r',
#                 alpha=0.1
#             )
#         )
#         # ax_PP.add_patch(
#         #     patches.Rectangle(
#         #         (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
#         #         width,          # width
#         #         2*binned_cosmic_variance(result,int(step),width),          # height
#         #         color='r',
#         #         alpha=0.1
#         #     )
#         # )
#         j = j+1
#         print j ,width
i=0
# print lmin_phiphi,lmax_phiphi,cl_phiphi,err_phiphi
# while i < len(lmin_phiphi):
#     ax_PP.add_patch(
#         patches.Rectangle(
#             (lmin_phiphi[i], -err_phiphi[i]),   # (x,y)
#             lmax_phiphi[i]-lmin_phiphi[i],          # width
#             err_phiphi[i]*2,          # height
#             color='grey',
#             alpha=0.2
#         )
#     )
#     i=i+1
fTT = interp1d(ll_LCDM,clTT_LCDM)
fEE = interp1d(ll_LCDM,clEE_LCDM)
# ax_TT.fill_between(l_TT,-err_TT/conversion/(l_TT*(l_TT+1)/2./np.pi)/fTT(l_TT),err_TT/conversion/(l_TT*(l_TT+1)/2./np.pi)/fTT(l_TT),facecolor='grey',alpha=0.2)
# ax_EE.fill_between(l_EE,-err_EE/conversion/(l_EE*(l_EE+1)/2./np.pi)/fEE(l_EE),err_EE/conversion/(l_EE*(l_EE+1)/2./np.pi)/fEE(l_EE),facecolor='grey',alpha=0.2)

# if i == 0:
#     var_color = 'k'
#     var_alpha = 1.
#     legarray.append(r'ref. $\Lambda CDM$')



for i in range(var_num):
    #
    # deal with varying parameters:
    #
    var_n = n[i]
    var_Omac = Omega_many_fld[i]
    dashes=dashes_array[i]
    # var_Omac = Omega_axion_ac[i]
    #
    print ' * Compute with %s=%d, %s=%e'%(var_name1,var_n,var_name2,var_Omac)
    #
    # deal with colors and legends
    #

    # var_legend = r'$n_{\rm a}=%d$'%(var_n)
    var_legend = legend_name[i]
    var_color = color[i]
    # var_alpha = 1.*i/(var_num-1.)
    legarray.append(var_legend)
    # if i == var_num-1:
    #     legarray.append(var_legend)
    #
    # call CLASS
    #
    M = Class()
    M.set(common_settings)
    M.set({var_name1:var_n,var_name2:var_Omac,
            'format':'camb',
            'w_fld_parametrization':'pheno_axion',
            'cs2_is_1':'no',
            'a_c':a_c[i],
            'Theta_initial_fld':3,
            'use_ppf':'no',
            'use_big_theta_fld':'yes',
            'fld_has_perturbations':'yes'})

    M.compute()
    #
    # get Cls
    #
    # clM = M.raw_cl(2500)
    th = M.get_thermodynamics()
    z = th['z']
    xe = th['x_e']
    T = th['Tb [K]']
    clM = M.lensed_cl(2500)
    ll = clM['ell'][2:]
    clTT = clM['tt'][2:]
    clEE = clM['ee'][2:]
    clPP = clM['pp'][2:]
    #
    # store P(k) for common k values
    #
    pkM = []
    for k in kvec:
        pkM.append(M.pk(k,0.))
    #
    # plot P(k)
    # ax_Pk.axis([1.e-4,1.,-0.21,0.21])##1e-5
    ax_Pk.axis([1.e-4,1.,-0.11,0.11])
    # ax_Pk.axis([1.e-4,1.,1e3,1e5])

    ax_Pk.semilogx(kvec,(np.array(pkM)-np.array(pkM_LCDM))/np.array(pkM_LCDM),var_color,lw=5,dashes=dashes,label=var_legend)
    # ax_Pk.loglog(kvec,np.array(pkM),var_color,lw=1)
    # ax_Pk.loglog(kvec,np.array(pkM_LCDM),'k',lw=1)
    #
    # plot C_l^TT
    #
    #ax_TT.semilogx(ll,clTT*ll*(ll+1)/twopi,var_color,lw=5)
    # ax1.plot(x, y, color="r") # regular plot (red)
    # ax1.set_xlabel('x')


    # ax_TT.axis([2,2500,-0.06,0.06])
    # ax_TT.semilogx(ll,(ll)*(ll+1)/(2*np.pi)*(clTT),var_color,lw=5)
    # ax_TT.loglog(ll,(clTT_LCDM),var_color,lw=5)
    ax_TT.semilogx(ll,(clTT-clTT_LCDM)/clTT_LCDM,var_color,lw=5,dashes=dashes,label=var_legend)
    # ax_TT.plot(ll,(clTT-clTT_LCDM)/clTT_LCDM,var_color,lw=5)
    # ax_TT.fill_between(l_TT,-err_TT/conversion/(l_TT*(l_TT+1)/2./np.pi)/fTT(l_TT),err_TT/conversion/(l_TT*(l_TT+1)/2./np.pi)/fTT(l_TT),facecolor='grey',alpha=0.2)


    # plot Cl phiphi
    #
    # ax_PP.loglog(ll,clPP*ll*(ll+1)*ll*(ll+1)/twopi,var_color,linestyle='-')
    #
    # ax_PP.axis([2,2500,-0.21,0.21])##ac1e-5
    ax_PP.axis([2,2500,-0.11,0.11])##ac1e-1
    # ax_PP.semilogx(ll,(ll)*(ll+1)/(2*np.pi)*(clPP),var_color,linestyle='-')
    # ax_PP.loglog(ll,(clPP_LCDM),var_color,linestyle='-')
    ax_PP.semilogx(ll,(clPP-clPP_LCDM)/clPP_LCDM,var_color,lw=5,dashes=dashes,label=var_legend)

    # ax_PP.fill_between((lmin_phiphi+lmax_phiphi)/2,-(err_phiphi)/cl_phiphi,(err_phiphi)/cl_phiphi,facecolor='grey',alpha=0.2)

    # ax_PP.axis([50,2500,-0.5,0.5])
    # ax_PP.plot(ll,(clPP-clPP_LCDM)/clPP_LCDM,var_color,linestyle='-')
    # ax_PP.fill_between((lmin_phiphi+lmax_phiphi)/2,-(err_phiphi)/cl_phiphi,(err_phiphi)/cl_phiphi,facecolor='grey',alpha=0.2)



    ax_EE.axis([2,2500,-0.11,0.11])
    # ax_EE.semilogx(ll,(ll)*(ll+1)/(2*np.pi)*(clEE),var_color,lw=5)
    # ax_EE.loglog(ll,(clEE_LCDM),var_color,lw=5)
    # ax_EE.semilogx(ll,(clEE-clEE_LCDM)/clEE_LCDM,var_color,lw=5)
    # fxe=interp1d(z_LCDM,xe_LCDM)
    # axe.semilogx(z,(xe-fxe(z))/fxe(z),var_color,lw=5)
    # axT.loglog(z,T,var_color,lw=5)
    # ax_EE.axis([50,2500,-0.06,0.06])
    ax_EE.semilogx(ll,(clEE-clEE_LCDM)/clEE_LCDM,var_color,lw=5,dashes=dashes,label=var_legend)

    # ax_EE.fill_between(l_EE,-err_EE/conversion/(l_EE*(l_EE+1)/2./np.pi)/fEE(l_EE),err_EE/conversion/(l_EE*(l_EE+1)/2./np.pi)/fEE(l_EE),facecolor='grey',alpha=0.2)

    #
    # reset CLASS
    #
    M.struct_cleanup()

    #

ax_Pk.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$',fontsize=62)
# ax_Pk.set_ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
# ax_Pk.text(0.007,0.11,r'$\frac{\Delta P(k)}{P(k)^{\Lambda{\rm CDM}}}$',fontsize=65)###ac1em5
ax_Pk.text(0.007,0.05,r'$\frac{\Delta P(k)}{P(k)^{\Lambda{\rm CDM}}}$',fontsize=65)
ax_PP.set_xlabel(r'$\ell$',fontsize=65)
ax_PP.text(50,0.03,r'$\frac{\Delta C_\ell^{\phi\phi}}{C_\ell^{\phi\phi}({\Lambda{\rm CDM}})}$',fontsize=65)
# ax_PP.text(50,0.06,r'$\frac{\Delta C_\ell^{\phi\phi}}{C_\ell^{\phi\phi}({\Lambda{\rm CDM}})}$',fontsize=65)
ax_EE.fill_between(l_EE,-err_EE/conversion/(l_EE*(l_EE+1)/2./np.pi)/fEE(l_EE),err_EE/conversion/(l_EE*(l_EE+1)/2./np.pi)/fEE(l_EE),facecolor='grey',alpha=0.2)
ax_PP.fill_between((lmin_phiphi+lmax_phiphi)/2,-(err_phiphi)/cl_phiphi,(err_phiphi)/cl_phiphi,facecolor='grey',alpha=0.2)
ax_TT.fill_between(l_TT,-err_TT/conversion/(l_TT*(l_TT+1)/2./np.pi)/fTT(l_TT),err_TT/conversion/(l_TT*(l_TT+1)/2./np.pi)/fTT(l_TT),facecolor='grey',alpha=0.2)

ax_TT.axis([2,2500,-0.06,0.06])
ax_TT.set_xlabel(r'$\ell$',fontsize=65)
ax_TT.text(50,0.02,r'$\frac{\Delta C_\ell^\mathrm{TT}}{C_\ell^\mathrm{TT}(\Lambda{\rm CDM})}$',fontsize=65)

ax_Pk.legend(frameon=False,prop={'size':45},loc='upper left',borderaxespad=0.)

ax_EE.set_xlabel(r'$\ell$',fontsize=65)

ax_EE.text(50,0.04,r'$\frac{\Delta C_\ell^\mathrm{EE}}{C_\ell^\mathrm{EE}(\Lambda{\rm CDM})}$',fontsize=65)


# fig_TT.savefig('spectra_nalp_zc1e3_ErrorPlanck.pdf')
# fig_TT.savefig('spectra_nalp_zc1e5_ErrorPlanck.pdf', bbox_inches='tight')
fig_TT.savefig('spectra_nalp_n1_ErrorPlanck.pdf', bbox_inches='tight')
# figxe.savefig('xe_nalp_zc1e3_ErrorPlanck.pdf')



# In[ ]:




# In[ ]:
