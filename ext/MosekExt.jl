module MosekExt
using MosekTools
using TSSOS
using JuMP

function mosek_optimizer()
    optimizer_with_attributes(Mosek.Optimizer, 
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-8,  # primal feasibility tolerance
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-8,  # dual feasibility tolerance
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-8,  # relative primal-dual gap tolerance
    "MSK_DPAR_OPTIMIZER_MAX_TIME" => 1e4,  # limit of running time
    "MSK_IPAR_NUM_THREADS" => 0 # number of threads available for Mosek
)
end

end