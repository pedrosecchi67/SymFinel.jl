"""
Struct defining a regressor

* `interp`: interpolation definition
* `points`: vector of vectors with as many coordinates as it's necessary to univocally determining
    the linar coefficients at `interp.linpars`
* `interpmat`: matrix that defines coefficients for variable polynomials given
    evaluations at points
"""
mutable struct Regressor
    interp::InterpolationFunction
    pts::Vector{
        Vector{Float64}
    }
    interpmat::Matrix{Float64}
end

"""
Constructor for Regressor

* `interp`: interpolation definition
* `points`: vector of vectors with as many coordinates as it's necessary to univocally determining
    the linar coefficients at `interp.linpars`
"""
function Regressor(interp::InterpolationFunction, pts::Vector{Vector{Float64}})
    interpmat=inv(
        get_interp_linmat(interp, pts)
    )

    Regressor(interp, pts, interpmat)
end

"""
Given a regressor, find polynomial coefficients to interpolate a variable, with variables provided as values at points
"""
function regress(regr::Regressor, values::Vector{Float64})
    return regr.interpmat*values
end
