module COSMOExt
using COSMO
using TSSOS
using JuMP

#mosek_setting is left here temporarily before removing it from interfaces
function cosmo_optimizer(mosek_setting::TSSOS.MosekParameters)
    optimizer_with_attributes(COSMO.Optimizer)
end

function cosmo_optimizer(::Nothing)
    optimizer_with_attributes(COSMO.Optimizer)
end

end