module MosekExt
using MosekTools
using TSSOS
using JuMP

function mosek_optimizer(mosek_setting::TSSOS.MosekParameters)
    optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas, 
                "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap, "MSK_DPAR_OPTIMIZER_MAX_TIME" => mosek_setting.time_limit, "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads)
end

function mosek_optimizer(::Nothing)
    optimizer_with_attributes(Mosek.Optimizer)
end

end