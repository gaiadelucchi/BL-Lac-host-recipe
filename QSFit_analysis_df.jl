
## Testing spectra on QSFit - Julia 

using GModelFit
using QSFit; using QSFit.QSORecipes; using GModelFitViewer

using DataFrames
using DataStructures
using Statistics

using CSV 
using FITSIO

# --------------------------------------------------------------------------------------------------------------------------

import QSFit.QSORecipes.add_qso_continuum!

import QSFit: set_lines_dict!
import QSFit.QSORecipes: fit!

recipe = CRecipe{Type1}()

abstract type MyRecipe <: Type1 end
function set_lines_dict!(recipe::CRecipe{<: MyRecipe})
    recipe.lines = OrderedDict{Symbol, QSFit.SpectralLine}()
    return get_lines_dict(recipe)
end

function fit!(recipe::CRecipe{<: MyRecipe}, resid::GModelFit.Residuals)
    GModelFit.update!(resid.meval)
    if GModelFit.nfree(resid) == 0
        return nothing
    end
    return @invoke fit!(recipe::CRecipe{<: Type1}, resid)
end
#-----------------------------------------------------------------------------------------------------------
function QSFit_analysis(spectrum_path, host_template, z)

    spec = Spectrum(Val(:ASCII), spectrum_path, label="$(spectrum_path)_FITTED_WITH_$(host_template)")

    recipe = CRecipe{MyRecipe}(redshift=z, Av=0.0) #extintion coeff Av = 0, synthetic spectra are not attenuated

    push!(recipe.host_template, :library => "swire", :template => host_template)
    recipe.n_nuisance = 0 # QSFit does not consider nuisance lines in the fit
    recipe.use_ironuv = false # same with optical iron 
    recipe.use_ironopt = false # same with uv iron 
    recipe.use_balmer = false # it does not consider Balmer lines

    res = analyze(recipe, spec)
    show(res.bestfit)

    #viewer(res)
  
    return res     
end 

# --------------- Data frame creation with QSFit results and parameters ------------------------------------
function retrieve_df(name_file, redshift, host_in, host_temp, res, shift) 
    
    fit_stat = res.fitstats.fitstat
    
    # QSO values
    PL = res.bestfit.params[:QSOcont].norm.val
    sigma_PL = res.bestfit.params[:QSOcont].norm.unc
    alpha = res.bestfit.params[:QSOcont].alpha.val
    sigma_alpha= res.bestfit.params[:QSOcont].alpha.val
            
    #Galaxy values
    val_host = res.bestfit.params[:Galaxy].norm.val
    sigma_host = res.bestfit.params[:Galaxy].norm.unc
    rel_err_host = (res.bestfit.params[:Galaxy].norm.unc)/(res.bestfit.params[:Galaxy].norm.val)

    # Creare un DataFrame con i dati estratti
    df = DataFrame(
        filename = [name_file],
        redshift = [redshift],
        input_host = [host_in],
        QSFit_host = [host_temp],
        fit_stat = [fit_stat],
        PL = [PL],
        sigma_PL = [sigma_PL],
        alpha = [alpha],
        sigma_alpha = [sigma_alpha],
        val_host = [val_host],
        sigma_host = [sigma_host],
        rel_err_host = [rel_err_host],
        sed = [shift]
    )

    return df
end

#--------------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------
function retrieve_input(file)
   
    #file name structure: spec_HOST TYPE_LUM HOST_BIN LUM BLLAC_RATIO HOST/BLAZAR_REDSHIFT.txt
    #split(a, '_')[2] = HOST TYPE
    #split(a, '_')[3] = LUMINOSITY OF HOST
    #split(a, '_')[4] = GAMMA LUMINOSITY BIN OF BL LAC
    #split(a, '_')[5] = RATIO HOST/BLAZAR
    #split(a, '_')[6] = REDSHIFT
    # If the 7th part exist: split(a, '_')[7] = SHIFT
    
    a = replace(file, ".txt" => "") # delete  .txt
    parts = split(a, "_")
    
    host_in = parts[2]
    host_L = parts[3]
    gamma_bin = parts[4]
    ratio_host_blazar = parts[5]
    z = parse(Float64, parts[6])  #parse transforms a string type in a float64
    println(host_in , " " , host_L , " " , gamma_bin,  " ",  ratio_host_blazar ,  " " , z)

    parts = split(a, "_")
    if length(parts) >= 7
        shift = parts[7]
    else
        shift = "ORIGINAL"
    end
        
    
        
    return host_in, host_L, gamma_bin, ratio_host_blazar, z, shift
end

function retrieve_files(path)
    files = readdir(path)
    # If there are different types of files the following row filters out .txt files
    text_files = [file for file in files if endswith(file, ".txt")]
    return text_files
end
#-----------------------------------------------------------------------------------------------------------------------------

#dir_path = "./spectra_shiftedSED/" 
#dir_path = "./noisy_fluxes_10_plusshiftedSED/" 

dir_path = "./noisy_fluxes_5/" 


files = retrieve_files(dir_path)

host_temp = ["Ell5", "Sb"]  # intermediate kinds of host galaxies to analyse synthetic spectra


merge_df = DataFrame()

for file in files
    host_input, host_L, Bl_gamma_bin, ratio_h_b, redshift, shift = retrieve_input(file)
    spec_name = file  
    spec_path = dir_path*file 
    
    for host in host_temp 
        QSFit_res = QSFit_analysis(spec_path, host, redshift)
        df = retrieve_df(file, redshift, host_input, host, QSFit_res, shift)      
        append!(merge_df, df)
            
    end
       
end


name_csv = "./df_noise_5.csv"
merge_df[!,:ratio_host_bl] = merge_df.val_host ./ merge_df.PL

CSV.write(name_csv, merge_df)

#------------------------------------------------------------------
