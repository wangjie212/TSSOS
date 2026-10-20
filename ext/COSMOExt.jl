module COSMOExt
using COSMO
using TSSOS
using JuMP

function cosmo_optimizer()
    #https://github.com/oxfordcontrol/COSMO.jl/blob/master/src/settings.jl
    optimizer_with_attributes(COSMO.Optimizer, "eps_abs"=>1e-8, "eps_rel"=>1e-8, "eps_prim_inf"=>1e-8, "eps_dual_inf"=>1e-8)
end

end