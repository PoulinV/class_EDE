
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

# In[ ]:

# esthetic definitions for the plots

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
# matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


# In[ ]:

############################################
#
# Varying parameter (others fixed to default)
#
var_name1 = 'n_pheno_axion'
var_name2 = 'Omega_fld_ac'
# var1_min = 1
# var1_max = 3
# var2_min = 1e14
n  = [1,2,3]
# Omega_axion_ac = [1e14,5e14,1e15]
Omega_axion_ac = [2.278e+5,4.278e+5,1.8e+5]
color = ['r','b','g']
var_num = 3
var_legend = r'$n_{\rm ALP}$'

var_figname = 'nalp'
#
#############################################
#
# Fixed settings
#
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   # LambdaCDM parameters
                   #'h':0.67556,
                   '100*theta_s':1.043,
                   'omega_b':0.022032,
                   'omega_cdm':0.12038,
                   'A_s':2.215e-9,
                   'n_s':0.9619,
                   'tau_reio':0.0925}
                   # ,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   #'YHe':0.246,
                   # other output and precision parameters
                   # 'P_k_max_1/Mpc':3.0,
                   # 'l_switch_limber':9}
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

fig_perts, (k1,k2) = plt.subplots(2,1, sharex=True,gridspec_kw=dict(height_ratios=[1,1]),figsize=(10,8))
# fig_Pk, ax_Pk = plt.subplots()
fig_perts.subplots_adjust(hspace=0)

k1.tick_params(axis='x',labelsize=25)
k1.tick_params(axis='y',labelsize=25)
k2.tick_params(axis='x',labelsize=25)
k2.tick_params(axis='y',labelsize=25)
# ax_Pk.tick_params(axis='x',labelsize=30)
# ax_Pk.tick_params(axis='y',labelsize=30)

# ax_TT.axis([2,2500,-0.7,0.7])
# ax_EE.axis([2,2500,-0.7,0.7])

#
# loop over varying parameter values
#

M = Class()
M.set(common_settings)
M.set({'k_output_values':'1,0.1,0.01,0.0001'})
M.compute()


all_k = M.get_perturbations()  # this potentially constains scalars/tensors and all k values

one_k = all_k['scalar'][1]     # this contains only the scalar perturbations for the requested k values
# print one_k
a = one_k['a']
delta_cdm = one_k['delta_cdm']
k1.loglog(a,np.abs(delta_cdm),'k',label=r'$\delta_{\rm cdm} (\Lambda{\rm CDM})$',lw=1)
k1.text(6e-5,5e2,r'$k=10^{-2}$Mpc$^{-1}$',fontsize=23)
var_legend = r'$\delta_{\rm cdm} (\Lambda{\rm CDM})$'
legarray.append(var_legend)
one_k = all_k['scalar'][3]     # this contains only the scalar perturbations for the requested k values
# print one_k

a = one_k['a']
delta_cdm = one_k['delta_cdm']
# store P(k) for common k values
k2.loglog(a,np.abs(delta_cdm),'k',lw=1)
k2.text(6e-5,1e3,r'$k=1$Mpc$^{-1}$',fontsize=23)
for i in range(var_num):
    #
    # deal with varying parameters:
    #
    var_n = n[i]
    var_Omac = Omega_axion_ac[i]
    #
    print ' * Compute with %s=%d, %s=%e'%(var_name1,var_n,var_name2,var_Omac)
    #
    # deal with colors and legends
    #

    var_legend = r'$\delta_{\rm ALP}, n_{\rm ALP}=%d$'%(var_n)
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
    M.set({'k_output_values':'0.1,0.01'})
    M.set({var_name1:var_n,var_name2:var_Omac,
            'w_fld_parametrization':'pheno_axion',
            'cs2_and_ca2_switch':'no',
            'cs2_is_1':'no',
            'a_c':1e-3,
            'Theta_initial_fld':3,
            'use_ppf':'no',
            'use_big_theta_fld':'yes',
            'fld_has_perturbations':'yes'})
    M.compute()
    #
    # get Cls
    #
    # clM = M.raw_cl(2500)

    all_k = M.get_perturbations()  # this potentially constains scalars/tensors and all k values

    one_k = all_k['scalar'][0]     # this contains only the scalar perturbations for the requested k values

    a = one_k['a']
    delta_fld = one_k['delta_fld[0]']
    # store P(k) for common k values
    print color[i]
    k1.loglog(a,np.abs(delta_fld),color[i],label=var_legend,lw=1)
    one_k = all_k['scalar'][1]     # this contains only the scalar perturbations for the requested k values

    a = one_k['a']
    delta_fld = one_k['delta_fld[0]']
    # store P(k) for common k values
    k2.loglog(a,np.abs(delta_fld),color[i],label=var_legend,lw=1)

    #
    # plot Cl phiphi
    #
    #ax_PP.loglog(ll,clPP*ll*(ll+1)*ll*(ll+1)/twopi,color=var_color,alpha=var_alpha,linestyle='-')
    #
    # reset CLASS
    #
    M.struct_cleanup()

k2.set_xlabel(r'$a$',fontsize=30)
k1.set_xlim([5e-5,1])
k2.set_xlim([5e-5,1])
# k1.set_xlabel(r'$\ell$',fontsize=40)
# ax_EE.set_ylabel(r'$[\ell(\ell+1)/2\pi]  C_\ell^\mathrm{EE}$')
# ax_EE.set_ylabel(r'$\frac{\Delta C_\ell^\mathrm{EE}}{C_\ell^\mathrm{EE}(\Lambda{\rm CDM})}$',fontsize=40)
fig_perts.tight_layout()
k1.legend(legarray,frameon=False,prop={'size':23})
fig_perts.savefig('perts_nalp.pdf')


# ax_EE.legend(legarray,frameon=False,prop={'size':16})
# fig_EE.tight_layout()
# fig_EE.savefig('spectra_%s_clee.pdf' % var_figname)
#
# output of C_l^pp figure
#
# ax_PP.set_xlim([10,2500])
# ax_PP.set_xlabel(r'$\ell$')
# ax_PP.set_ylabel(r'$[\ell^2(\ell+1)^2/2\pi]  C_\ell^\mathrm{\phi \phi}$')
# ax_PP.legend(legarray)
# fig_PP.tight_layout()
# fig_PP.savefig('spectra_%s_clpp.pdf' % var_figname)


# In[ ]:




# In[ ]:
