module COSMOExt
using COSMO
using TSSOS
using JuMP

#mosek_setting is left here temporarily before removing it from interfaces
function cosmo_optimizer(mosek_setting::TSSOS.MosekParameters)
    optimizer_with_attributes(COSMO.Optimizer, "eps_abs"=>1e-8, "eps_rel"=>mosek_setting.tol_relgap, "eps_prim_inf"=>mosek_setting.tol_pfeas, "eps_dual_inf"=>mosek_setting.tol_dfeas, "verbose"=>true)
end

function cosmo_optimizer(::Nothing)
    #https://github.com/oxfordcontrol/COSMO.jl/blob/master/src/settings.jl
    optimizer_with_attributes(COSMO.Optimizer, "eps_abs"=>1e-8, "eps_rel"=>1e-8, "eps_prim_inf"=>1e-8, "eps_dual_inf"=>1e-8)
end

end