
from MyPyLib_v3 import *


# --------------------------- Bestimmung Gitterkonstante Hg-Lampe ---------------------------


def degr( skalenteile ):
    return 0.15 * skalenteile

def rad( degree ):
    return degree/360 * 2 * np.pi


omegaG_Hg       = np.array([51, 51.5, 52, 54, 54 + degr(1), 54 + degr(2), 59, 68, 68 + degr(1.4), 71, 71.5,73 ])
omegaG_Hg_err   = np.full(len(omegaG_Hg), 1)

lambdas_Hg      = np.array([404.656, 407.783, 410.805, 433.922, 434.749, 435.833, 546.074, 576.960, 579.066, 623.440, 671.643, 690.752 ])

omegaB, omegaB_err = 140, 1

exp_Intens_Hg   = np.array([1800, 150, 40, 250, 400, 4000, 1100, 240, 280, 30, 160, 250])
pers_confidence = np.array([100, 70, 5, 20, 20, 5, 100, 10, 10, 20, 20, 30 ])

alpha_Hg        = omegaG_Hg
alpha_Hg_err    = omegaG_Hg_err

beta_Hg         = omegaG_Hg + omegaB - 180 
beta_Hg_err     = np.sqrt(omegaG_Hg_err**2 + omegaB_err**2)

func_of_lambda_Hg       = np.sin(rad(alpha_Hg)) + np.sin(rad(beta_Hg))
func_of_lambda_Hg_err   = np.sqrt( (rad(alpha_Hg_err)*np.cos(rad(alpha_Hg)))**2 + (rad(beta_Hg_err)*np.cos(rad(beta_Hg)))**2 )

mask_reliance_2     = (exp_Intens_Hg > 500)
mask_reliance       = (pers_confidence > -1)
mask_unreli         = (mask_reliance == False)

dataset_Hg      = [ lambdas_Hg[mask_reliance], None, func_of_lambda_Hg[mask_reliance], func_of_lambda_Hg_err[mask_reliance] ]
dataset_Hg_ign  = [ lambdas_Hg[mask_unreli], None, func_of_lambda_Hg[mask_unreli], func_of_lambda_Hg_err[mask_unreli] ]

def fitfunc_Gitter(lambdas, a, b):
    return lambdas / a + b  

fitbook_Hg = {
    "data_set"          : dataset_Hg,
    "region_of_interest": None,
    "fit_function"      : fitfunc_Gitter,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitdata_Hg, params_Hg, params_Hg_err, chi_Hg = ultimate_fit( fitbook_Hg )


def Plot_Hg_Gitterkonst():
    
    sample_format_dict_1 = {
        "label"      : r"Messwerte (n=1)",          
        "fmt"        : '.', 
        "color"      : "black",                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    sample_format_dict_2 = {
        "label"      : f"Linearer Fit ($\\chi^2 = {chi_Hg:.1f}$)",                       
        "fmt"        : '--', 
        "color"      : sns.color_palette("bright")[0]   ,        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"$\lambda$ / nm",
        "y_ax_label"  : r"$\sin{(\alpha)} + \sin{(\beta)}$ / 1"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
    
    all_data                = [ dataset_Hg, fitdata_Hg ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/Hg_Gitterkonstante_Fit.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

Plot_Hg_Gitterkonst()
Plot_Hg_Gitterkonst()

print(f"Gitterkonstante: a = ({params_Hg[0]} +- {params_Hg_err[0]}) nm")
print(f"Offset:          b = ({params_Hg[1]} +- {params_Hg_err[1]}) ")

a_iterative     = lambdas_Hg/func_of_lambda_Hg
a_iterative_err = np.abs( a_iterative * func_of_lambda_Hg_err/func_of_lambda_Hg ) 

#print("Naive Berechnug von a: \n", a_iterative)



# --------------------------- Bestimmung Gitterkonstante Hg-Lampe (2) ---------------------------


omegaG_Hg2       = np.array([-27.5, -24.5, -19, -14 ])
omegaG_Hg2_err   = np.full(len(omegaG_Hg2), 1)

lambdas_Hg2      = np.array([576.960, 546.074, 491.607, 435.833])
omegaB2, omegaB2_err = 140, 1

exp_Intens_Hg2   = np.array([240, 1100, 80, 4000])

alpha_Hg2        = omegaG_Hg2
alpha_Hg2_err    = omegaG_Hg2_err

beta_Hg2         = omegaG_Hg2 + omegaB2 - 180 
beta_Hg2_err     = np.sqrt(omegaG_Hg2_err**2 + omegaB2_err**2)

func_of_lambda_Hg2       = np.sin(rad(alpha_Hg2)) + np.sin(rad(beta_Hg2))
func_of_lambda_Hg2_err   = np.sqrt( (rad(alpha_Hg2_err)*np.cos(rad(alpha_Hg2)))**2 + (rad(beta_Hg2_err)*np.cos(rad(beta_Hg2)))**2 )

dataset_Hg2      = [ lambdas_Hg2, None, func_of_lambda_Hg2, func_of_lambda_Hg2_err ]


fitbook_Hg2 = {
    "data_set"          : dataset_Hg2,
    "region_of_interest": None,
    "fit_function"      : fitfunc_Gitter,
    "params_guess"      : [-400,0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitdata_Hg2, params_Hg2, params_Hg2_err, chi_Hg2 = ultimate_fit( fitbook_Hg2 )

def Plot_Hg_Gitterkonst2():
    
    sample_format_dict_1 = {
        "label"      : r"Messwerte (n = -1)",          
        "fmt"        : '.', 
        "color"      : "black",                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    sample_format_dict_2 = {
        "label"      : f"Linearer Fit ($\\chi^2 = {chi_Hg2:.1f}$)",                       
        "fmt"        : '--', 
        "color"      : sns.color_palette("bright")[0]   ,        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"$\lambda$ / nm",
        "y_ax_label"  : r"$\sin{(\alpha)} + \sin{(\beta)}$ / 1"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
    
    all_data                = [ dataset_Hg2, fitdata_Hg2 ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/Hg_Gitterkonstante_Fit2.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

Plot_Hg_Gitterkonst2()
Plot_Hg_Gitterkonst2()

print(f"Gitterkonstante: a2 = ({params_Hg2[0]} +- {params_Hg2_err[0]}) nm")
print(f"Offset:          b2 = ({params_Hg2[1]} +- {params_Hg2_err[1]}) ")

a_iterative2 = -lambdas_Hg2/func_of_lambda_Hg2
a_iterative2_err = np.abs( a_iterative2 * func_of_lambda_Hg2_err/func_of_lambda_Hg2 ) 

#print("Naive Berechnug von a2: \n", a_iterative2)


# --------------------------- Bestimmung Gitterkonstante Hg-Lampe joint venture ---------------------------

a_clomun        = np.concatenate((a_iterative2, a_iterative))
a_clomun_err    = np.concatenate((a_iterative2_err, a_iterative_err ))


a_average,a_average_err     = stichproben_varianz(a_clomun)

print(f"a_average           = ({a_average:.2f}+-{a_average_err:.2f}) nm")

#print("a_iteratives 2,1:")
#print( a_clomun, a_clomun_err )