from classy import Class
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from classy import CosmoComputationError
from time import time
from matplotlib.backends.backend_pdf import PdfPages
from sys import exit
import matplotlib
from matplotlib import rc
from scipy.interpolate import interp1d

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
# matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


start = time()

'''
	This is likely the only section you will need to modify
	All parameters defined here
'''

# Functionality for reading in .ini file isn't built in yet.
# Will build that if needed.

# Parameters we won't be changing
params = {'w_fld_parametrization' : 'pheno_axion',
            'cs2_is_1' : 'no',
			'output':'tCl,pCl,lCl,mPk',
            'lensing':'yes',
            # 'a_c' : 1e-3,
            'Theta_initial_fld' : 3,
            #, 'Omega_many_fld' : 1e-10
			'compute damping scale':'yes'
            }
n = 1
wn = (n-1.)/(n+1.)
params['n_pheno_axion'] = n
params['cs2_fld'] = 1
params['use_ppf'] = 'no'
params['use_big_theta_fld'] = 'yes'
# params['cs2_is_w'] = 'yes'
params['fld_has_perturbations'] = 'yes'

params['omega_b'] = 0.022383
params['omega_cdm'] = 0.12011
params['tau_reio'] = 0.0543
params['ln10^{10}A_s'] = 3.0448
params['n_s'] = 0.96605


H0_or_theta_s = 'H0'
v1_or_v2 = 'v2'
if H0_or_theta_s == 'theta_s':
	theta_s = 0.01040909
	params['100*theta_s'] = 100*theta_s
	param_list=['HP','rs_rec','rs_rec_over_rd_rec','da_rec','H0']
	param_legend = ['Peak Height',r'$r_s$',r'$r_s/r_d$',r'$d_a$',r'$H_0$']
	filename = "log_derivative_n1_vThetas_fiducial5percent.txt"
	filename2 = "log_derivative_n1_vThetas_fiducial5percent_alaTristan.txt"
elif H0_or_theta_s == 'H0':
	params['H0'] = 67
	param_list=['HP','rs_rec','rs_rec_over_rd_rec','da_rec','theta_s']
	param_legend = ['Peak Height',r'$r_s$',r'$r_s/r_d$',r'$d_a$',r'$100~\theta_s$']
	filename = "log_derivative_n1_vH0_fiducial5percent.txt"
	filename2 = "log_derivative_n1_vH0_fiducial5percent_alaTristan.txt"

# params['h'] = 0.67
# params['back_integration_stepsize'] = 5e-4
# params['reio_parametrization'] = 'reio_none'

N = 50 # Number of ac values
color_list = ['blue','red','green','purple','orange']
# params = ['omega_cdm','Om_ac','N_ur']
ac_values = np.logspace(-5, -2, N, endpoint = True)
der_param = np.zeros((N,len(param_list)))
Omega_ac_fid = np.zeros((N))

percentage_difference = 0.05
fraction_fiducial = 0.02

write_file = False
load_from_file = True

###############################################################################################################################################
# DEFINE FUNCTION TO CALCULATE THINGS
###############################################################################################################################################
def calculate_derivative_param(params,param_to_test,param_fiducial,percentage_difference,param_list,der_param,write_file,i):
	cosmo = Class() # Create an instance of the CLASS wrapper
	param = np.zeros((N,3,len(param_list)))
	for j in range(3):
		# going over a_c and Om_fld values
		# if j==0:
		# 	params['Omega_fld_ac'] = Omega_ac_fid[i]
		# if j==1:
		# 	params['Omega_fld_ac'] = Omega_ac_fid[i]+percentage_difference*Omega_ac_fid[i]
		# if j==2:
		# 	params['Omega_fld_ac'] = Omega_ac_fid[i]-percentage_difference*Omega_ac_fid[i]
		if j==0:
			params[param_to_test] = param_fiducial
		if j==1:
			params[param_to_test] = param_fiducial+percentage_difference*param_fiducial
		if j==2:
			params[param_to_test] = param_fiducial-percentage_difference*param_fiducial
		# if j==0:
		# 	params['Omega_many_fld'] = Omega0_fld
		# if j==1:
		# 	params['Omega_many_fld'] = Omega0_fld+percentage_difference*Omega0_fld
		# if j==2:
		# 	params['Omega_many_fld'] = Omega0_fld-percentage_difference*Omega0_fld

		# try to solve with a certain cosmology, no worries if it cannot
		# l_theta_s = np.pi/(theta_s)
		# print "here:",l_theta_s
		# try:
		cosmo.set(params) # Set the parameters to the cosmological code
		cosmo.compute() # solve physics
		cl = cosmo.lensed_cl(2500)
		ell = cl['ell'][2:]
		tt = cl['tt'][2:]
		fTT = interp1d(ell,tt*ell*(ell+1))

		for k in range(len(param_list)):
			#calculate height peak difference
			if(param_list[k] == 'HP'):
				# param[i][j][k] = max(tt)-fTT(10)
				param[i][j][k] = max(tt*ell*(ell+1))
			elif(param_list[k] == 'rs_rec'):
				param[i][j][k] = cosmo.rs_rec()
			elif(param_list[k] == 'rd_rec'):
				param[i][j][k] = cosmo.rd_rec()
			elif(param_list[k] == 'rs_rec_over_rd_rec'):
				print cosmo.rs_rec()/cosmo.rd_rec()
				param[i][j][k] = cosmo.rs_rec()/cosmo.rd_rec()
			elif(param_list[k] == 'da_rec'):
				param[i][j][k] = cosmo.da_rec()
			elif(param_list[k] == 'theta_s'):
				param[i][j][k] = cosmo.theta_s()
			elif(param_list[k] == 'H0'):
				param[i][j][k] = cosmo.Hubble(0)

		print max(tt*ell*(ell+1)), cosmo.theta_s(), cosmo.da_rec(), cosmo.rs_rec(), cosmo.z_rec()
		# for l in range(len(tt)):
		# 	# print l, tt[l], max(tt*ell*(ell+1))
		# 	if tt[l]*ell[l]*(ell[l]+1) == max(tt*ell*(ell+1)):
		# 		print "l:", l,max(tt*ell*(ell+1))

			# except CosmoComputationError: # this happens when CLASS fails
			# 	print CosmoComputationError
			# 	pass # eh, don't do anything

		#calculate derivative
		# der_HP[i]= (HP[i][1]-HP[i][2])/(Omega_ac_fid[i]+percentage_difference*Omega_ac_fid[i]-(Omega_ac_fid[i]-percentage_difference*Omega_ac_fid[i]))*Omega_ac_fid[0]/HP[i][0]

		cosmo.empty()
		cosmo.struct_cleanup()
	if write_file == True:
		f.write(str(ac_values[i])+'\t\t')
	print "calculating derivative"
	for k in range(len(param_list)):
	# der_HP[i]= (HP[i][1]-HP[i][2])/(fraction_fiducial+percentage_difference*fraction_fiducial-(fraction_fiducial-percentage_difference*fraction_fiducial))*fraction_fiducial/HP[i][0]
		der_param[i][k] = (param[i][1][k]-param[i][2][k])/(param_fiducial+percentage_difference*param_fiducial-(param_fiducial-percentage_difference*param_fiducial))*param_fiducial/param[i][0][k]
		if write_file == True:
			f.write(str(der_param[i][k])+'\t\t') # info on format
		print param_list[k],der_param[i][k]
	if write_file == True:
		f.write('\n')
	return
	# der_rs_drag[i] =
	# der_HP[i]= (HP[i][1]-HP[i][2])/(Omega0_fld+percentage_difference*Omega0_fld-(Omega0_fld-percentage_difference*Omega0_fld))*Omega0_fld/HP[i][0]



########### v1: CALCULATE LOG DERIVATIVE WRT f(a_c) AS A FUNCTION OF a_c  ###########
if v1_or_v2=='v1':
	param_to_test = 'fraction_axion_ac'
	param_fiducial = fraction_fiducial
	if load_from_file == False:
		if write_file == True:
			f = open(filename, 'w') # creates a new file to write in / overwrites
			f.write('# a_c, HP,rs_rec,rs_rec_over_rd_rec,da_rec,theta_s\n') # info on format
		for i in range(N):
			params['a_c'] = ac_values[i]
			print "a_c = ",ac_values[i]
			calculate_derivative_param(params,param_to_test,param_fiducial,percentage_difference,param_list,der_param,write_file,i)
		if write_file == True:
			f.close()

	elif load_from_file == True:
		file = np.loadtxt(filename, comments='#')
		if write_file == True:
			f2 = open(filename2, 'w') # creates a new file to write in / overwrites
	###############################################################################################################################################
	# PLOT THINGS
	###############################################################################################################################################

	end = time()
	print('\n\nTime taken for everything but plot = %.2f seconds' %(end - start))



	fig, ax = plt.subplots() # new plot
	ax.tick_params(axis='x',labelsize=23)
	ax.tick_params(axis='y',labelsize=23)
	if write_file == True:
		f2.write('# a_c, HP,rs_rec,rs_rec_over_rd_rec,da_rec,theta_s\n') # info on format
		for i in range(len(file[:,0])):
			f2.write(str(file[i,0])+'\t\t')
			for k in range(len(param_list)):
				f2.write(str(file[i,k+1]/fraction_fiducial)+'\t\t')
			f2.write('\n')

	for k in range(len(param_list)):
		if load_from_file == True:
			ax.plot(file[:,0],file[:,k+1],label=param_legend[k],color=color_list[k])
			ax.set_xlim((file[:,0].min(), file[:,0].max()))
			# ax.set_
		# ax.set_ylim((file[:,k+1].min(),der_param.max()))
		else:
			ax.plot(ac_values,der_param[:,k],label=param_legend[k],color=color_list[k])
			ax.set_xlim((ac_values.min(), ac_values.max()))
			ax.set_ylim((der_param.min(),der_param.max()))
	if H0_or_theta_s == 'H0':
		ax.set_title(r'$n = $ %d, $H_0 = 67$' %params['n_pheno_axion'], fontsize = 22)
	elif H0_or_theta_s == 'theta_s':
		ax.set_title(r'$n = $ %d, $100~\theta_s = 1.040909$' %params['n_pheno_axion'], fontsize = 22)

	if load_from_file == True and write_file == True:
		f2.close()
	ax.set_xlabel(r'$a_c$', fontsize = 23)
	ax.set_ylabel(r'dln Param/dln$f_{\rm EDE}(a_c)$', fontsize = 23)
	# ax.set_ylabel(r'dlog(X)/dlog(Y)', fontsize = 23)
	#
	# # ax.set_xscale("log")
	# # ax.set_yscale("log")
	plt.xscale('log')
	# plt.yscale('log')



	ax.legend(frameon=False,prop={'size':15},loc='best',borderaxespad=0.)

		# for k in range(len(param_list)):
		# 	ax.plot((10**-4)*(1+np.zeros((N,len(param_list)))),der_param[:,k],color=color_list[k],label=param_legend[k])
	# plt.savefig(output_file + '.png',bbox_inches='tight')

	plt.show()






########## v2: CALCULATE LOG DERIVATIVE WRT DIFFERENT COSMOLOGICAL PARAMETERS ##########

elif v1_or_v2 == 'v2':
	params_new = {'output':'tCl,pCl,lCl,mPk',
	        'lensing':'yes',
			'compute damping scale':'yes'}
	params_new['omega_b'] = 0.022383
	params_new['omega_cdm'] = 0.12011
	params_new['tau_reio'] = 0.0543
	params_new['ln10^{10}A_s'] = 3.0448
	params_new['n_s'] = 0.96605
	if H0_or_theta_s == 'H0':
		params_new['H0'] = 67
		param_to_test = ['N_ur','omega_cdm','H0','fraction_axion_ac']
		param_fiducial = [3.046,0.12,67,0.1]
		param_list=['HP','rs_rec','rs_rec_over_rd_rec','da_rec','theta_s']
		labels = ['blank','Peak Height',r'$r_s$',r'$r_d/r_d$',r'$d_a$',r'$100~\theta_s$']
		param_legend = [r'$N_{\rm eff}$',r'$\omega_{\rm cdm}$',r'$H_0$',r'$f_{\rm EDE}(a_c)$']
	if H0_or_theta_s == 'theta_s':
		params_new['100*theta_s'] = 1.040909
		param_to_test = ['N_ur','omega_cdm','100*theta_s','fraction_axion_ac']
		param_fiducial = [3.046,0.12,1.040909,0.1]
		param_list=['HP','rs_rec','rs_rec_over_rd_rec','da_rec','H0']
		labels = ['blank','Peak Height',r'$r_s$',r'$r_s/r_d$',r'$d_a$',r'$H_0$']
		param_legend = [r'$N_{\rm eff}$',r'$\omega_{\rm cdm}$',r'$100~\theta_s$',r'$f_{\rm EDE}(a_c)$']

	fig, ax = plt.subplots() # new plot
	ax.tick_params(axis='x',labelsize=23)
	ax.tick_params(axis='y',labelsize=23)
	list = [0,1,2,3,4]
	for i in range(len(param_to_test)):
			print "param = ",param_to_test[i]
			if param_to_test[i] == 'fraction_axion_ac':
				params_new['w_fld_parametrization']='pheno_axion'
				params_new['cs2_is_1']='no'
				params_new['Theta_initial_fld']= 3
				n = 3.
				wn = (n-1.)/(n+1.)
				params_new['n_pheno_axion'] = n
				params_new['cs2_fld'] = 1
				params_new['use_ppf'] = 'no'
				params_new['use_big_theta_fld'] = 'yes'
				# params_new['cs2_is_w'] = 'yes'
				params_new['fld_has_perturbations'] = 'yes'
				params_new['a_c']=0.00016
			calculate_derivative_param(params_new,param_to_test[i],param_fiducial[i],percentage_difference,param_list,der_param,False,i)
			ax.plot(list,der_param[i,:],color=color_list[i],label=param_legend[i])

	end = time()
	print('\n\nTime taken for everything = %.2f seconds' %(end - start))


	ax.set_xlabel(r'cosmological parameter', fontsize = 23)
	ax.set_xticklabels(labels)
	# ax.set_ylabel(r'dln Param/dln$f_{\rm EDE}(a_c)$', fontsize = 23)
	# ax.set_ylabel(r'dlog(X)/dlog$(f_{\rm EDE}(a_c))$', fontsize = 23)
	ax.set_ylabel(r'dlog(X)/dlog(Y)', fontsize = 23)

	# # ax.set_xscale("log")
	# # ax.set_yscale("log")
	# plt.xscale('log')
	# plt.yscale('log')



	ax.legend(frameon=False,prop={'size':20},loc='best',borderaxespad=0.)

		# for k in range(len(param_list)):
		# 	ax.plot((10**-4)*(1+np.zeros((N,len(param_list)))),der_param[:,k],color=color_list[k],label=param_legend[k])
	# plt.savefig(output_file + '.png',bbox_inches='tight')

	plt.show()
