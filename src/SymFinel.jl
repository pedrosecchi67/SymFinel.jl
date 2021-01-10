module SymFinel

    using Calculus

    include("CoordExpr.jl")
    include("Interpolator.jl")
    include("Regressor.jl")
    include("Domain.jl")
    include("Finel.jl")
    include("Residuals.jl")

end # module

module SymFinelTests

    using Calculus

    include("CoordExpr.jl")
    include("Interpolator.jl")
    include("Regressor.jl")
    include("Domain.jl")
    include("Finel.jl")
    include("Residuals.jl")

    export get_basefun, get_basefun_derivative, get_linpar_contribution, get_interp_linmat,
        _get_gradfunc, regress, finel_get_regr, get_abs_coords!, from_finel_jac!, get_punctual_resfunc, 
            vals_to_polys, polys_to_vars

end # module
