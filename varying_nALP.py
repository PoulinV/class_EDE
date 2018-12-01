
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
Omega_axion_ac = [50,50,50]
color = ['r','b--','g-.']
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
                   'tau_reio':0.0925,
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

fig_TT, (ax_TT, ax_EE, ax_Pk) = plt.subplots(3,1, sharex=False,gridspec_kw=dict(height_ratios=[1,1,1]),figsize=(15,15))
# fig_Pk, ax_Pk = plt.subplots()
fig_TT.subplots_adjust(hspace=0)

ax_TT.tick_params(axis='x',labelsize=30)
ax_TT.tick_params(axis='y',labelsize=30)
ax_EE.tick_params(axis='x',labelsize=30)
ax_EE.tick_params(axis='y',labelsize=30)
ax_Pk.tick_params(axis='x',labelsize=30)
ax_Pk.tick_params(axis='y',labelsize=30)

ax_TT.axis([2,2500,-0.5,0.5])
ax_EE.axis([2,2500,-0.5,0.5])

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
clM = M.lensed_cl(2500)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
clPP_LCDM = clM['pp'][2:]
pkM_LCDM = []
### To plot the binned cosmic variance ####
l_min = 2.;
l_max = 2500.;
n_step = 15.;
j=0.
def binned_cosmic_variance (result,l_ini,l_tot):
    for i in range(0,int(l_tot)):
     result = result + 2/(2*(l_ini+float(i))+1)
    return np.sqrt(result/l_tot**2)


while j < n_step:
        result = 0.0
        step = l_min*(l_max/l_min)**(j/n_step)
        step_plus_1 = l_min*(l_max/l_min)**((j+1)/n_step)
        print int(step), int(step_plus_1)
        width = int(step_plus_1) - int(step)
        ax_TT.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='r',
                alpha=0.1
            )
        )
        ax_EE.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='r',
                alpha=0.1
            )
        )
        j = j+1
        print j ,width
# if i == 0:
#     var_color = 'k'
#     var_alpha = 1.
#     legarray.append(r'ref. $\Lambda CDM$')
for k in kvec:
    pkM_LCDM.append(M.pk(k,0.))



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

    var_legend = r'$n_{\rm ALP}=%d$'%(var_n)
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
            'w_fld_parametrization':'pheno_axion',
            'cs2_and_ca2_switch':'no',
            'cs2_is_1':'yes',
            'a_c':0.1,
            'Theta_initial_fld':3,
            'use_ppf':'no',
            'use_big_theta_fld':'yes',
            'fld_has_perturbations':'yes'})
    M.compute()
    #
    # get Cls
    #
    # clM = M.raw_cl(2500)
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
    #
    ax_Pk.semilogx(kvec,(np.array(pkM)-np.array(pkM_LCDM))/np.array(pkM_LCDM),var_color,lw=1)
    #ax_Pk.loglog(kvec,np.array(pkM),var_color,lw=1)
    #
    # plot C_l^TT
    #
    #ax_TT.semilogx(ll,clTT*ll*(ll+1)/twopi,var_color,lw=1)
    ax_TT.semilogx(ll,(clTT-clTT_LCDM)/clTT_LCDM,var_color,lw=1)
    #
    # plot Cl EE
    #
    #ax_EE.loglog(ll,clEE*ll*(ll+1)/twopi,var_color,lw=1)
    # EEinterp = interp1d(ll,clEE)
    ax_EE.semilogx(ll_LCDM,(clEE-clEE_LCDM)/clEE_LCDM,var_color,lw=1)
    #
    # plot Cl phiphi
    #
    #ax_PP.loglog(ll,clPP*ll*(ll+1)*ll*(ll+1)/twopi,color=var_color,alpha=var_alpha,linestyle='-')
    #
    # reset CLASS
    #
    M.struct_cleanup()

#
# output of P(k) figure
#
ax_Pk.set_xlim([1.e-4,3.])
ax_Pk.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$',fontsize=40)
# ax_Pk.set_ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')
ax_Pk.text(0.01,0.1,r'$\frac{\Delta P(k)}{P(k)^{\Lambda{\rm CDM}}}$',fontsize=40)
# ax_Pk.set_ylabel(r'$\frac{\Delta P(k)}{P(k)^{\Lambda{\rm CDM}}}$',fontsize=40)
ax_Pk.legend(legarray,frameon=False,prop={'size':30})
# fig_Pk.tight_layout()
# fig_Pk.savefig('spectra_%s_Pk.pdf' % var_figname)
#
# output of C_l^TT figure
#
# ax_TT.set_xlim([2,2500])
# ax_TT.set_ylim([-1,1])
ax_TT.set_xlabel(r'$\ell$',fontsize=40)
# ax_TT.set_ylabel(r'$[\ell(\ell+1)/2\pi]  C_\ell^\mathrm{TT}$')
ax_TT.text(50,0.15,r'$\frac{\Delta C_\ell^\mathrm{TT}}{C_\ell^\mathrm{TT}(\Lambda{\rm CDM})}$',fontsize=40)
# ax_TT.set_ylabel(r'$\frac{\Delta C_\ell^\mathrm{TT}}{C_\ell^\mathrm{TT}(\Lambda{\rm CDM})}$',fontsize=40)
# ax_TT.legend(legarray,frameon=False,prop={'size':30})

#
# output of C_l^EE figure
#
# ax_EE.set_xlim([2,2500])
ax_EE.set_xlabel(r'$\ell$',fontsize=40)
# ax_EE.set_ylabel(r'$[\ell(\ell+1)/2\pi]  C_\ell^\mathrm{EE}$')
ax_EE.text(50,0.15,r'$\frac{\Delta C_\ell^\mathrm{EE}}{C_\ell^\mathrm{EE}(\Lambda{\rm CDM})}$',fontsize=40)
# ax_EE.set_ylabel(r'$\frac{\Delta C_\ell^\mathrm{EE}}{C_\ell^\mathrm{EE}(\Lambda{\rm CDM})}$',fontsize=40)
fig_TT.tight_layout()
fig_TT.savefig('spectra_%s_cltt.pdf' % var_figname)


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
