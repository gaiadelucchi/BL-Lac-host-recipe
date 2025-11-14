using BLLacHostRecipe
using DataFrames
using QSFit
using CSV 

using DataStructures
using Statistics
using FITSIO

#-----------------------------------------------------------------------------------------------------------
function QSFit_analysis(spectrum_path, specobjid, host_template, z, ext_V)

    spec = Spectrum(Val(:SDSS_DR10), spectrum_path, label="$(specobjid)_$(host_template)")

    recipe = CRecipe{MyRecipe}(redshift=z, Av=ext_V) #coeff di estinzione Av 0 perché spettro simulato non è attenuato

    push!(recipe.host_template, :library => "swire", :template => host_template)
    recipe.n_nuisance = 0 #rimuove le nuisance lines
    recipe.use_ironuv = false #rimuove il ferro nell’ottico
    recipe.use_ironopt = false #rimuove il ferro nell’uv
    recipe.use_balmer = false #rimuove Balmer

    res = analyze(recipe, spec)
    show(res.bestfit)

    #viewer(res)
  
    return res     
end 

# --------------- Data frame creation with QSFit results and parameters ------------------------------------
function safe_param(params, comp_sym, param_sym)
    try
        p = params[comp_sym][param_sym]
        return p
    catch e
        return nothing
    end
end

function retrieve_df(name_file, specobjid, redshift, host_temp, res)
    fit_stat = res.fitstats.fitstat
    params = res.bestfit.params

    # QSO values
    qso_norm = safe_param(params, :QSOcont, :norm)
    PL = qso_norm !== nothing ? qso_norm.val : NaN
    sigma_PL = qso_norm !== nothing ? qso_norm.unc : NaN

    qso_alpha = safe_param(params, :QSOcont, :alpha)
    alpha = qso_alpha !== nothing ? qso_alpha.val : NaN
    sigma_alpha = qso_alpha !== nothing ? qso_alpha.unc : NaN

    # Galaxy values
    gal_norm = safe_param(params, :Galaxy, :norm)
    if gal_norm === nothing
        val_host = 0.0
        sigma_host = 0.0
        rel_err_host = 0.0
    else
        val_host = gal_norm.val
        sigma_host = gal_norm.unc
        rel_err_host = (sigma_host != 0 && val_host != 0) ? sigma_host / val_host : 0.0
    end

    df = DataFrame(
        filename = [name_file],
        specobj_id = [specobjid],
        redshift = [redshift],
        QSFit_host = [host_temp],
        fit_stat = [fit_stat],
        PL = [PL],
        sigma_PL = [sigma_PL],
        alpha = [alpha],
        sigma_alpha = [sigma_alpha],
        val_host = [val_host],
        sigma_host = [sigma_host],
        rel_err_host = [rel_err_host]
    )

    return df
end

#--------------------------------------------------------------------------------------------------------------------------------
function extract_Av(namefile)
    f = FITS(namefile)
    df = DataFrame(read(f[1]))
    close(f)
    return df
end

# Carica tabella galactic extinction
function load_extinction_table(txtfile::String)
    df = CSV.read(txtfile, DataFrame; delim=' ', ignorerepeated=true,
                  header=["ra", "dec", "extinction"])
    return df
end

# Estrae plate, mjd, fiberid dal nome file
function parse_spec_filename(path::String)
    fname = split(basename(path), '.')[1]
    parts = split(fname, '-')
    if length(parts) != 4 || parts[1] != "spec"
        error("Nome file non valido: $path")
    end
    plate, mjd, fiberid = parse.(Int, parts[2:4])
    return plate, mjd, fiberid
end

# Trova riga corrispondente per plate, mjd, fiberid
function find_ra_dec(sdss::DataFrame, plate::Int, mjd::Int, fiberid::Int)
    row = filter(r -> r["#plate"] == plate && r["mjd"] == mjd && r["fiberid"] == fiberid, sdss)
    if nrow(row) == 0
        error("Nessuna corrispondenza trovata per $plate-$mjd-$fiberid")
    end
    return row.ra[1], row.dec[1]
end

# Trova estinzione galattica con tolleranza
function find_extinction(extdf::DataFrame, ra::Float64, dec::Float64; tol=1e-5)
    row = filter(r -> isapprox(r.ra, ra; atol=tol) && isapprox(r.dec, dec; atol=tol), extdf)
    if nrow(row) == 0
        error("Nessuna estinzione trovata per sorgente RA=$ra DEC=$dec")
    end
    return row.extinction[1]
end

# Funzione principale
function extract_metadata(specfile::String, sdss::DataFrame, extdf::DataFrame)
    plate, mjd, fiberid = parse_spec_filename(specfile)
    ra, dec = find_ra_dec(sdss, plate, mjd, fiberid)
    extinction = find_extinction(extdf, ra, dec)
    return (plate=plate, mjd=mjd, fiberid=fiberid, ra=ra, dec=dec, extinction=extinction)
end


#-----------------------------------------------------------------------------------------------------------------------------
function retrieve_files(path)
    files = readdir(path)
    SDSS_files = [file for file in files if endswith(file, ".fits")]
    return SDSS_files
end
#-----------------------------------------------------------------------------------------------------------------------------

host_temp = ["Ell5", "Sb"] 

#dir_path_SDSS = "./SDSS_DR17_spectra/" 
#file = "spec-0266-51630-0529.fits"
#spec_name = file  
#spec_path = dir_path_SDSS*file 
#obj_hdu = f[3]
#z = read(obj_hdu, "z")[1]
#redshift = Float64(z)  
#specobjid = read(obj_hdu, "specobjid")[1]

#println("specobjid = ", specobjid)
#println("Z = ", z)

#host = "Sb"
#QSFit_res = QSFit_analysis(spec_path, specobjid, host, redshift)
#df = retrieve_df(file, specobjid, redshift, host, QSFit_res)          


dir_path_SDSS = "./SDSS_DR17_spectra/" 

#dir_path_SDSS = "./" 

files = retrieve_files(dir_path_SDSS)
#files = "spec-7752-58072-0116.fits"
merge_df = DataFrame()

for file in files
    
    spec_name = file  
    spec_path = dir_path_SDSS*file 
    

    f = FITS(spec_path)
    obj_hdu = f[3]
    z = read(obj_hdu, "z")[1]
    redshift = Float64(z)  
    specobjid = read(obj_hdu, "specobjid")[1]
    specobj_id = string(specobjid)

    println("specobjid = ", specobjid)
    println("Z = ", z)

    extintionV = 0.0
    
    for host in host_temp 
        QSFit_res = QSFit_analysis(spec_path, specobjid, host, redshift, extintionV)
        df = retrieve_df(file, specobjid, redshift, host, QSFit_res)          
        append!(merge_df, df; promote=true)      
    end
       
end


name_csv = "./df_analysis_SDSSDR17.csv"
merge_df[!,:ratio_host_bl] = merge_df.val_host ./ merge_df.PL

CSV.write(name_csv, merge_df)

#------------------------------------------------------------------
