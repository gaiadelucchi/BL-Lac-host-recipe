module BLLacHostRecipe

using GModelFit
using QSFit; using QSFit.QSORecipes; using GModelFitViewer
using DataStructures

# --------------------------------------------------------------------------------------------------------------------------

import QSFit.QSORecipes.add_qso_continuum!

import QSFit: set_lines_dict!
import QSFit.QSORecipes: fit!

export MyRecipe

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




end # module BLLacHostRecipe
