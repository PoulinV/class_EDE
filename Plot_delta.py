import numpy as np
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages




# gdm = np.loadtxt('output/Added_GDM09_background.dat')
# cdm = np.loadtxt('output/Added_GDM05_background.dat')
# ceff_cvis_0p3 = np.loadtxt('output/Added_GDM03_perturbations_k0_s.dat')

gdm = np.loadtxt('output/Added_GDM09_cl_lensed.dat')
cdm = np.loadtxt('output/Added_GDM05_cl_lensed.dat')

# gdm = np.loadtxt('output/Added_GDM09_perturbations_k0_s.dat')
# cdm = np.loadtxt('output/Added_GDM05_perturbations_k0_s.dat')


z_gdm = gdm[:,0]
# z = np.linspace(0.,2e3,num=2001,endpoint=True)
# rho_gdm = interp1d(gdm[:,0],gdm[:,11])
# H_gdm = interp1d(gdm[:,0],gdm[:,3])

# rho_cdm_gdm = interp1d(gdm[:,0],gdm[:,10])


z_cdm = cdm[:,0]
# rho_cdm = interp1d(cdm[:,0],cdm[:,10])
# H_cdm = interp1d(cdm[:,0],cdm[:,3])


l_gdm = gdm[:,0]
TT_gdm = interp1d( gdm[:,0] , gdm[:,1] )

l_cdm = cdm[:,0]
TT_cdm = interp1d( cdm[:,0] , cdm[:,1] )



# delta_gdm = interp1d(gdm[:,0],gdm[:,15])
# theta_gdm = interp1d(gdm[:,0],gdm[:,16])
# delta_cdm_gdm = interp1d(gdm[:,0],gdm[:,18])
# theta_cdm_gdm = interp1d(gdm[:,0],gdm[:,19])

# # z_cdm = cdm[:,0]
# delta_cdm = interp1d(cdm[:,0],cdm[:,15])
# theta_cdm = interp1d(cdm[:,0],cdm[:,16])


# plt.plot(z_cdm, np.abs(np.abs(rho_gdm[10:]/rho_cdm)))

# plt.plot(z, rho_gdm(z)/rho_cdm(z)-1.)
# plt.plot(z, H_gdm(z))
# plt.plot(z, H_cdm(z))


# plt.plot(z_gdm, (delta_gdm(z_gdm)-delta_cdm(z_cdm))/delta_cdm(z_cdm))
# plt.plot(z_gdm, delta_cdm_gdm(z_gdm))
# plt.plot(z_cdm, delta_cdm(z_cdm))

# plt.plot(z_cdm, theta_cdm)
# plt.plot(z_gdm, theta_cdm_gdm)
# plt.plot(z_gdm, theta_gdm)


# plt.plot(l_gdm, TT_gdm(l_gdm))
# plt.plot(l_cdm, TT_cdm(l_cdm))
plt.plot(l_gdm, ( TT_gdm(l_gdm) - TT_cdm(l_gdm))/TT_cdm(l_gdm))
# plt.xlabel(r'$\ell$', fontsize = 18)
# plt.ylabel(r'$\left(\frac{C_{\ell}^{GDM} - C_{\ell}^{CDM}}{C_{\ell}^{GDM}}\right)^{TT, lensed}$', fontsize = 18)
# plt.tight_layout()
# plt.title('Ratio of GDM to CDM in the TT poswer spectrum')

# plot_name = 'GDM_CDM_ratio_TT_unlensed.pdf'

# with PdfPages(plot_name) as pdf:
# 	pdf.savefig(bbox_inches='tight')


# plt.plot(z_gdm, rho_gdm)
# plt.plot(z_cdm, rho_cdm)

# plt.xscale("log")
# plt.ylim(-1e-12, 1e-12)
# plt.yscale("log")



plt.show()