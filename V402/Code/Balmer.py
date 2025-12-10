
from MyPyLib_v3 import *


# --------------------------- Zuordnung der Balmer-Linien ---------------------------



def degr( skalenteile ):
    return 0.15 * skalenteile

def rad( degree ):
    return degree/360 * 2 * np.pi

def exp_omega_G( lambdas, omegaB, n ):
    return ( np.pi + np.arcsin( n* lambdas/a - np.sin(omegaB)) ) / (2*np.pi)*360 - omegaB


    # verwendeter Wert Für Gitterkonstante: 24000 Striche / cm ---> [a] = nm
a_theo          = 1 / 24000 *10**7 
print(f"Hersteller: a_theo = {a_theo:.1f} nm")

    # finaler Wert für a 
a, a_err       = 423.46, 18.27 
a1, a1_err      = 517.6, 17.7
a2, a2_err      = 419.3, 33.3
print(f"Experiment: a_exp1  = ({a1:.1f}+-{a1_err:.1f}) nm")
print(f"Experiment: a_exp2  = ({a2:.1f}+-{a2_err:.1f}) nm")
print(f"Experiment: a_avg  = ({a:.1f}+-{a_err:.1f}) nm")



omegaB, omegaB_err = 140, 1  

colors          = ["rot, türkis, dunkelblau, violett"]    
omegaG_oku      = np.array([78, 59.5, 55, 53, 21.5 ])
omegaG_oku_err  = np.full(len(omegaG_oku), 1)

omegaG_oku2      = np.array([77.5, 63, 59 + degr(1.7), 59, np.NaN ])
omegaG_oku2_err  = np.full(len(omegaG_oku2), 1)

omegaG_CCD      = np.array([73.5, 63.5, np.NaN, np.NaN, 24])
omegaG_CCD_err  = np.full(len(omegaG_CCD), 1)

    # Wellenlängen in Luft alpha, beta, gamma, delta (NIST: https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=H&output_type=0&low_w=&upp_w=&unit=1&submit=Retrieve+Data&de=0&plot_out=0&I_scale_type=1&format=0&line_out=0&en_unit=0&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&unc_out=1&order_out=0&max_low_enrg=&show_av=2&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on)
lambdas_theo    = np.array([656.279, 486.1350, 434.0472, 410.1734, np.NaN])


def lambdas( omegaG, omegaG_err, omegaB, omegaB_err, a, a_err):
    alpha           = omegaG
    alpha_err       = omegaG_err
    beta            = omegaG + omegaB - 180 
    beta_err        = np.sqrt(omegaG_err**2 + omegaB_err**2)
    lambdas_exp     = a* (np.sin(rad(alpha)) + np.sin(rad(beta)))
    lambdas_exp_err = np.sqrt( (a*rad(alpha_err)*np.cos(rad(alpha)))**2 + (a*rad(beta_err)*np.cos(rad(beta)) + (lambdas_exp*a_err/a)**2 )**2 )
    return lambdas_exp, lambdas_exp_err

lambdas_oku, lambdas_oku_err = lambdas(omegaG_oku, omegaG_oku_err, omegaB, omegaB_err, a, a_err)
lambdas_oku2, lambdas_oku2_err = lambdas(omegaG_oku2, omegaG_oku2_err, omegaB, omegaB_err, a, a_err)
lambdas_CCD, lambdas_CCD_err = lambdas(omegaG_CCD, omegaG_CCD_err, omegaB, omegaB_err, a, a_err)


print("Erwartete lambda: ", lambdas_theo)
print("Okuar lambda:     ", lambdas_oku)
print("err:              ", lambdas_oku_err)
print("Okuar2 lambda:    ", lambdas_oku2)
print("err:              ", lambdas_oku2_err)
print("CCD lambda:       ", lambdas_CCD)
print("err:              ", lambdas_CCD_err)

header_Balmer   = [ r"$\omega_G\up{oku}$ / \unit{\degree}",
                    r"$\omega_G\up{oku2}$ / \unit{\degree}",
                    r"$\omega_G\up{CCD}$ / \unit{\degree}",
                    r"$\lambda\up{oku}$ / \unit{nm}", 
                    r"$\lambda\up{oku2}$ / \unit{nm}", 
                    r"$\lambda\up{CCD}$ / \unit{nm}",
                    r"$\lambda\sub{theo}$ / \unit{nm}",  ]
data_Balmer     = [omegaG_oku, omegaG_oku2, omegaG_CCD, lambdas_oku, lambdas_oku2, lambdas_CCD, lambdas_theo]
data_Balmer_err = [omegaG_oku_err, omegaG_oku2_err, omegaG_CCD_err, lambdas_oku_err, lambdas_oku2_err, lambdas_CCD_err, np.full(len(lambdas_theo), 0)]


array_to_latex_2( "../Data/Params_Balmerlinien.txt", data_Balmer, data_Balmer_err, header_Balmer, ".1f")



# --------------------------- Bestimmung der Rydbergkonstante ---------------------------

lambdas_exp     = lambdas_oku[0:4]
lambdas_exp_err = lambdas_oku_err[0:4]

pro_lambda      = 1/lambdas_exp
pro_lambda_err  = lambdas_exp_err/lambdas_exp**2  # vernachlässige Fehler

orders_exp      = np.array([3,4,5,6])
func_order_exp  = 1/2**2 - 1/orders_exp**2

dataset_Ryd     = [func_order_exp, None, pro_lambda, None]


fitbook_Ryd = {
    "data_set"          : dataset_Ryd,
    "region_of_interest": None,
    "fit_function"      : fitfunc_linear,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitdata_Ryd, params_Ryd, params_Ryd_err, chi_Ryd = ultimate_fit( fitbook_Ryd )

R, b            = params_Ryd
R_err, b_err    = params_Ryd_err

print("Params Rydberg-Plot:")
print(f"R      = ({R*10**9:.4e} +- {R_err*10**9:.4e}) / m")
print(f"b      = ({b:.4e} +- {b_err:.4e}) ")


m_e     = 9.1093837139 * 10**(-31)      # in kg
q_e     = 1.602176634  * 10**(-19)      # in C
eps_0   = 8.8541878188 * 10**(-12)      # in F / m
c_0     = 299792458                     # in m / s
h_theo  =  6.62607015 * 10**(-34)       # in J * s

R_theo     = m_e * q_e**4 / ( 8 * c_0 * eps_0**2 * h_theo **3 )

print(f"R_theo = {R_theo:.4e} / m")
print(f"h_theo = {h_theo:.4e} * Js")


h       = ( m_e * q_e**4 / ( 8 * c_0 * eps_0**2 * R *10**9 ) )**(1/3)
h_err   = 1/3 * h * R_err *10**9 / R

print(f"h      = ({h:.4e}+-{h_err:.4e}) J * s")

def Plot_Rydberg():
    
    sample_format_dict_1 = {
        "label"      : r"Messwerte ",          
        "fmt"        : 'x', 
        "color"      : "black",                               
        "markersize" : 10, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    sample_format_dict_2 = {
        "label"      : f"Linearer Fit ",                       
        "fmt"        : '--', 
        "color"      : sns.color_palette("bright")[0]   ,        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Rydberg-Formel $\left( \frac{1}{2^2} - \frac{1}{n^2} \right)$ / 1",
        "y_ax_label"  : r"Inverse Wellenlänge  $\frac{1}{\lambda} $ / nm$^{-1}$ "
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
    
    all_data                = [ dataset_Ryd, fitdata_Ryd ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/Rydberg_Fit.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

Plot_Rydberg()
Plot_Rydberg()



# --------------------------- Okular Isotopie-Aufspaltung ---------------------------

    # Aufspaltung konnte nur im Aufbau Oku2 beoabachtet werden, hier sind die Winkel jedoch völlig falsch 
    # ordne gemessene Aufspaltung vermuteten Linien zu

lambdas_theo

f_obj               = 300 # in mm

delta_skalent_oku2  = np.array([ degr(0.15), degr(0.05), np.NaN, np.NaN, np.NaN ]) #  in mm 
delta_skalent_oku2_err = np.full(5, 0.05 )

betas_rad_oku2      = (omegaB +  omegaG_oku2 - 180) / 360 * 2*np.pi 
betas_rad_oku2_err  = np.sqrt(omegaB_err**2 +  omegaG_oku2_err**2) / 360 * 2*np.pi 

delta_lambda_oku2   = a * np.cos(betas_rad_oku2) * delta_skalent_oku2 / f_obj
delta_lambda_oku2_err = delta_lambda_oku2 * np.sqrt( (a_err/a)**2 + (delta_skalent_oku2_err/delta_skalent_oku2)**2 + (betas_rad_oku2_err *np.tan(betas_rad_oku2))**2 )


print("Oku2 Delta lambdas:     ", delta_lambda_oku2 )
print("Oku2 Delta lambdas_err: ", delta_lambda_oku2_err )



# --------------------------- CCD Isotopie-Aufspaltung ---------------------------

sequence_replacement_txt("../Data/aufspaltung_24grad.txt", ",", ".")
sequence_replacement_txt("../Data/aufspaltung_63-5grad.txt", ",", ".")
sequence_replacement_txt("../Data/aufspaltung_73-5grad.txt", ",", ".")


data_null  = np.loadtxt("../Data/aufspaltung_24grad_improved.txt", skiprows=1)
data_beta  = np.loadtxt("../Data/aufspaltung_63-5grad_improved.txt", skiprows=1)
data_alpha = np.loadtxt("../Data/aufspaltung_73-5grad_improved.txt", skiprows=1)

betas_rad_CCD       = (omegaB +  omegaG_CCD - 180) / 360 * 2*np.pi 
betas_rad_CCD_err   = np.sqrt(omegaG_CCD_err**2 + omegaB_err**2) /360 *2*np.pi


def fitfunc_double_gauss(x, mu1, sigma1, A1, mu2, sigma2, A2, b):
    return fitfunc_gauss_amplitude(x, mu1, sigma1, A1) + fitfunc_gauss_amplitude(x, mu2, sigma2, A2) + b

# ---------- alpha

angle_alpha         = data_alpha[:, 0]  # in degree
intens_alpha        = data_alpha[:, 1]

angle_alpha_err     = np.full( len(angle_alpha), 1 )
intens_alpha_err    = intens_alpha * 0.05
dataset_alpha       = [angle_alpha, None, intens_alpha, intens_alpha_err]

fitbook_alpha = {
    "data_set"          : dataset_alpha,
    "region_of_interest": [-0.03, 0.13],
    "fit_function"      : fitfunc_double_gauss,
    "params_guess"      : [0.01, 0.04, 3, 0.09, 0.02, 1.2, 0.4],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata_alpha, params_alpha, params_alpha_err, chi_alpha = ultimate_fit(fitbook_alpha)

mu1, mu2            = params_alpha[0], params_alpha[3]
mu1_err, mu2_err    = params_alpha_err[0], params_alpha_err[3]

delta_winkel_alpha      = np.abs(mu1-mu2) 
delta_winkel_alpha_err  = np.sqrt( mu1_err**2 + mu2_err**2 )
delta_lambda_alpha      = a * np.cos(betas_rad_CCD[0]) * delta_winkel_alpha/360 *2*np.pi
delta_lambda_alpha_err  = delta_lambda_alpha* np.sqrt( (a_err/a)**2 + (delta_winkel_alpha_err/delta_winkel_alpha)**2 + (np.tan(betas_rad_CCD[0])*betas_rad_CCD_err[0])**2 )

# ---------- beta

angle_beta         = data_beta[:, 0]  # in degree
intens_beta        = data_beta[:, 1]

angle_beta_err     = np.full( len(angle_beta), 1 )
intens_beta_err    = intens_beta * 0.05
dataset_beta       = [angle_beta, None, intens_beta, intens_beta_err]

fitbook_beta = {
    "data_set"          : dataset_beta,
    "region_of_interest": [-0.333,-0.225],
    "fit_function"      : fitfunc_double_gauss,
    "params_guess"      : [-0.325, 0.01, 0.8, -0.265, 0.02, 2.2, 0.4],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata_beta, params_beta, params_beta_err, chi_beta = ultimate_fit(fitbook_beta)

mu1, mu2            = params_beta[0], params_beta[3]
mu1_err, mu2_err    = params_beta_err[0], params_beta_err[3]

delta_winkel_beta      = np.abs(mu1-mu2) 
delta_winkel_beta_err  = np.sqrt( mu1_err**2 + mu2_err**2 )
delta_lambda_beta      = a * np.cos(betas_rad_CCD[1]) * delta_winkel_beta/360 *2*np.pi
delta_lambda_beta_err  = delta_lambda_beta* np.sqrt( (a_err/a)**2 + (delta_winkel_beta_err/delta_winkel_beta)**2 + (np.tan(betas_rad_CCD[1])*betas_rad_CCD_err[1])**2 )


# ---------- nullte Ordnung

angle_null         = data_null[:, 0]  # in degree
intens_null        = data_null[:, 1]

angle_null_err     = np.full( len(angle_null), 1 )
intens_null_err    = intens_null * 0.05
dataset_null       = [angle_null, None, intens_null, intens_null_err]


fitbook_null = {
    "data_set"          : dataset_null,
    "region_of_interest": [-1.45,-1.15],
    "fit_function"      : fitfunc_double_gauss,
    "params_guess"      : [-1.41, 0.05, 5, -1.2, 0.1, 60, 0.4],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata_null, params_null, params_null_err, chi_null = ultimate_fit(fitbook_null)


mu1, mu2            = params_null[0], params_null[3]
mu1_err, mu2_err    = params_null_err[0], params_null_err[3]

delta_winkel_null      = np.abs(mu1-mu2) 
delta_winkel_null_err  = np.sqrt( mu1_err**2 + mu2_err**2 )
delta_lambda_null      = a * np.cos(betas_rad_CCD[4]) * delta_winkel_null/360 *2*np.pi 
delta_lambda_null_err  = delta_lambda_null* np.sqrt( (a_err/a)**2 + (delta_winkel_null_err/delta_winkel_null)**2 + (np.tan(betas_rad_CCD[4])*betas_rad_CCD_err[4])**2 )


# ---------- Fits


def plot_aufspaltung( id, dataset, range, fitdata=None, fitparams=None, chi=None ):
        
    if( not isinstance(fitdata, list) ):
        chi = 0


    sample_format_dict_exp = {
        "label"      : r"Messwerte ",          
        "fmt"        : 'o', 
        "color"      : "black",                               
        "markersize" : 2, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_fit = {
        "label"      : f"Gaußfit $(\\chi^2 = {chi:.1f})$",                       
        "fmt"        : '-', 
        "color"      : sns.color_palette("bright")[0],        
        "markersize" : 4, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }  
    sample_format_dict_vertical_1 = {
        "label"      : r"$\mu_1$",                       
        "fmt"        : '-.', 
        "color"      : sns.color_palette("bright")[4],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  
    sample_format_dict_vertical_2 = {
        "label"      : r"$\mu_2$",                       
        "fmt"        : '-.', 
        "color"      : sns.color_palette("bright")[6],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  
    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Winkel $\beta$ / Grad",
        "y_ax_label"  : r"Intensität $I$ / 1"
    }
    
    general_format_dict = standard_format_dict.copy()
    general_format_dict["custom_x_range"] = [True, range[0], range[1]]
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()

    all_data                = [ dataset ]                               
    all_sample_format_dicts = [ sample_format_dict_exp ]
    
    if( isinstance(fitdata, list) ):
        all_data.append(fitdata)
        all_sample_format_dicts.append(sample_format_dict_fit)

        data_vertical_1 = [ np.full(2, fitparams[0]), None, np.array([0, max(dataset[2])]), None ]
        data_vertical_2 = [ np.full(2, fitparams[3]), None, np.array([0, max(dataset[2])]), None ]
        
        all_data.append(data_vertical_1)
        all_data.append(data_vertical_2)

        all_sample_format_dicts.append(sample_format_dict_vertical_1)
        all_sample_format_dicts.append(sample_format_dict_vertical_2)

    save_plot = True, "../Figures/Aufspaltung_CCD_"+id+".jpg" 
   
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

plot_aufspaltung("alpha", dataset_alpha, [-0.1,0.15], fitdata_alpha, params_alpha, chi_alpha)
plot_aufspaltung("beta", dataset_beta, [-0.38,-0.18], fitdata_beta, params_beta, chi_beta)
plot_aufspaltung("null", dataset_null, [-1.5,-1.05], fitdata_null, params_null, chi_null)


# ---------- Tabellierung der Aufspaltungen

delta_winkel_CCD        = [delta_winkel_alpha, delta_winkel_beta, delta_winkel_null]
delta_winkel_CCD_err    = [delta_winkel_alpha_err, delta_winkel_beta_err, delta_winkel_null_err]
delta_lambda_CCD        = [delta_lambda_alpha, delta_lambda_beta, delta_lambda_null]
delta_lambda_CCD_err    = [delta_lambda_alpha_err, delta_lambda_beta_err, delta_lambda_null_err]


print("Winkeldifferenzen:  ", delta_winkel_CCD)
print("Winkeldifferenzen:  ", delta_winkel_CCD_err)
print("lambda-differenzen: ", delta_lambda_CCD)
print("lambda-differenzen: ", delta_lambda_CCD_err)


header_aufspalt     = [
    r"$\mu_1$ / \unit{\degree}",
    r"$\sigma_1$ / \unit{\degree}",
    r"$A_1$ / 1",
    r"$\mu_2$ / \unit{\degree}",
    r"$\sigma_2$ / \unit{\degree}",
    r"$A_2$ / 1",
    r"$b$ / 1",
    r"$\Delta \beta$ / \unit{\degree}",
    r"$\Delta \lambda$ / \unit{nm}"
]
params_aufspalt     = np.array([ params_alpha, params_beta, params_null ]).T 
params_aufspalt_err = np.array([ params_alpha_err, params_beta_err, params_null_err ]).T 

params_aufspalt     = np.concatenate(( params_aufspalt, [delta_winkel_CCD], [delta_lambda_CCD] ))
params_aufspalt_err = np.concatenate(( params_aufspalt_err, [delta_winkel_CCD_err], [delta_lambda_CCD_err] ))

array_to_latex_2("../Data/Params_Aufspaltung_CCD.txt", params_aufspalt, params_aufspalt_err, header_aufspalt, ".4f")


# ---------- Linienbreiten




    # betrachte nur das sigma2 vom größeren Peak
sigmas_alpha, sigma_beta, null              = params_aufspalt[4]
sigmas_alpha_err, sigma_beta_err, null_err  = params_aufspalt_err[4]

FWHM_alpha              = 2*np.sqrt(2*np.log(2)) * sigmas_alpha
FWHM_alpha_err          = 2*np.sqrt(2*np.log(2)) * sigmas_alpha_err

FWHM_beta               = 2*np.sqrt(2*np.log(2)) * sigma_beta
FWHM_beta_err           = 2*np.sqrt(2*np.log(2)) * sigma_beta_err

lambda_FWHM_alpha       = a * np.cos(betas_rad_CCD[0]) * FWHM_alpha/360 *2*np.pi 
lambda_FWHM_alpha_err   = delta_lambda_null* np.sqrt( (a_err/a)**2 + (FWHM_alpha_err/FWHM_alpha)**2 + (np.tan(betas_rad_CCD[0])*betas_rad_CCD_err[0])**2 )

lambda_FWHM_beta        = a * np.cos(betas_rad_CCD[1]) * FWHM_beta/360 *2*np.pi 
lambda_FWHM_beta_err    = delta_lambda_null* np.sqrt( (a_err/a)**2 + (FWHM_beta_err/FWHM_beta)**2 + (np.tan(betas_rad_CCD[1])*betas_rad_CCD_err[1])**2 )



    # Doppler-Formel: https://juser.fz-juelich.de/record/874600/files/Energie_Umwelt_489.pdf

k_B         = 1.380649*10**(-23) 
m_p         = 1.67262192595*10**(-27) 
factor_dopp = np.sqrt(8*np.log(2)* 1000*k_B/(18*m_p*c_0**2) ) # Wassermolekül

lambda_dopp_alpha, lambda_dopp_alpha_err   = lambdas_theo[0] * factor_dopp, 0
lambda_dopp_beta, lambda_dopp_alpha_err    = lambdas_theo[1] * factor_dopp, 0

tau_alpha          = 1 /(4.410 *10**(-1) *10**8)
tau_beta           = 1 /(8.419 *10**(-2) *10**8) 

print(f"tau_alpha = {tau_alpha*10**9:.1e} ns")
print(f"tau_beta  = {tau_beta*10**9:.1e} ns")


lambda_nat_alpha   = 10**(-18) * lambdas_theo[0]**2 /(2*np.pi *tau_alpha *c_0)
lambda_nat_beta    = 10**(-18) * lambdas_theo[1]**2 /(2*np.pi *tau_beta *c_0)

print(f"lambda_nat_alpha = {lambda_nat_alpha*10**9:.1e} nm")
print(f"lambda_nat_beta  = {lambda_nat_beta*10**9:.1e} nm")
print(lambdas_theo[1])

column1     = [lambdas_theo[0], lambdas_theo[1]]
column1_err = [0,0]
column2     = [sigmas_alpha, sigma_beta]
column2_err = [sigmas_alpha_err, sigma_beta_err]
column3     = [FWHM_alpha, FWHM_beta]
column3_err = [FWHM_alpha_err, FWHM_beta_err]
column4     = [lambda_FWHM_alpha,lambda_FWHM_beta]
column4_err = [lambda_FWHM_alpha_err,lambda_FWHM_beta_err]
column5     = [lambda_dopp_alpha, lambda_dopp_beta]
column5_err = [0,0]
column6     = [tau_alpha*10**9, tau_beta*10**9]
column6_err = [0,0]
column7     = [lambda_nat_alpha*10**9, lambda_nat_beta*10**9]
column7_err = [0,0]

header_linbreite = [
    r"$\lambda \up{theo}$ / \unit{nm}",
    r"$\sigma_2$ / \unit{\degree}",
    r"$FWHM_2$ / \unit{\degree}",
    r"$\Delta \lambda_2 $ / \unit{nm}",
    r"$\Delta \lambda \sub{dop} $ / \unit{nm}",
    r"$\tau $ / \unit{ns}",
    r"$\Delta \lambda \sub{nat} $ / \unit{nm}",
]

data_linbreite = [column1,column2, column3, column4, column5, column6, column7]
data_linbreite_err = [column1_err,column2_err, column3_err, column4_err, column5_err, column6_err, column7_err]

array_to_latex_2("../Data/Params_Linbreite.txt", data_linbreite, data_linbreite_err, header_linbreite, ".8f")




def lambda_theo_formel(n):
    factor_H = 1/(1 + m_e/m_p)
    factor_D = 1/(1 + m_e/(2*m_p))

    lambda_H = (factor_H*R_theo * (1/4 - 1/n**2) )**(-1)
    lambda_D = (factor_D*R_theo * (1/4 - 1/n**2) )**(-1)
    print(f"lambda_theo_formel H (n = {n:.0f})     = {lambda_H*10**9:.4f} nm")
    print(f"lambda_theo_formel D (n = {n:.0f})     = {lambda_D*10**9:.4f} nm")
    print(f"delta_lambda_theo_formel (n = {n:.0f}) = {(lambda_H-lambda_D)*10**9:.4f} nm")



lambda_theo_formel(3)
lambda_theo_formel(4)
lambda_theo_formel(5)

A_alpha     = lambdas_theo[0]/delta_lambda_alpha
A_alpha_err = lambdas_theo[0]/delta_lambda_alpha**2 * delta_lambda_alpha_err

A_beta     = lambdas_theo[0]/delta_lambda_beta
A_beta_err = lambdas_theo[0]/delta_lambda_beta**2 * delta_lambda_beta_err

print(f"Auflösungsvermögen der H_alpha-Linie: {A_alpha:.2f}+-{A_alpha_err:.2f}")
print(f"Auflösungsvermögen der H_beta-Linie:  {A_beta:.2f}+-{A_beta_err:.2f}")