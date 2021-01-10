export InterpolationFunction

"""
Struct defining an interpolation basis

* `ex`: expression defining the basis function format
* `coords`: vector of symbols for coordinates in 3D space
* `linpars`: vector of symbols for linear parameters, defining the interpolation
    as a linear transformation
* `DOFs`: `Vector{Vector{Int64}}`. Each inner vector defines a set of derivative orders
    for a degree of freedom of the basis function, in relationship to the variables at `coords`.

    Example:
        a degree of freedom defined by `[0, 1, 0]` characterizes the first derivative, in the second
        cartesian axis at a node in the three-dimensional space.
"""
mutable struct InterpolationFunction
    ex::Expr
    coords::Vector{Symbol}
    linpars::Vector{Symbol}
    DOFs::Vector{
        Vector{Int64}
    }
    _DOF_analysis_functions::Matrix{Function}
    _funcs::Array{Function}
    _npts::Int64
end

"""
Constructor for `InterpolationFunction` struct.

* `ex`: expression defining the basis function format
* `coords`: vector of symbols for coordinates in 3D space
* `linpars`: vector of symbols for linear parameters, defining the interpolation
    as a linear transformation
* `DOFs`: `Vector{Vector{Int64}}`. Each inner vector defines a set of derivative orders
    for a degree of freedom of the basis function, in relationship to the variables at `coords`.

    Example:
        a degree of freedom defined by `[0, 1, 0]` characterizes the first derivative, in the second
        cartesian axis at a node in the three-dimensional space.

* `maximum_derivative`: maximum derivative order to be kept in function form
"""
function InterpolationFunction(
    ex::Expr,
    coords::Vector{Symbol},
    linpars::Vector{Symbol},
    DOFs::Vector{
        Vector{Int64}
    };
    maximum_derivative::Int64=1
)
    _DOF_analysis_functions=Matrix{Function}(undef, length(DOFs), length(linpars))

    for (i, DOF) in enumerate(DOFs)
        for (j, lp) in enumerate(linpars)
            _DOF_analysis_functions[i, j]=get_linpar_contribution(ex, coords, linpars, lp, DOF)
        end
    end

    nd=length(DOFs[1])

    _funcs=Array{Function, nd}(undef, [maximum_derivative+1 for i=1:nd]...)

    for fcinds in CartesianIndices(_funcs)
        _funcs[fcinds]=get_basefun_derivative(
            ex,
            coords,
            linpars,
            collect(Tuple(fcinds)).-1
        )
    end

    InterpolationFunction(
        ex,
        coords,
        linpars,
        DOFs,
        _DOF_analysis_functions,
        _funcs,
        length(linpars)/length(DOFs)
    )
end

"""
Gets a matrix defining DOFs at given sample points as a linear transformation of `linpars`
"""
function get_interp_linmat(interp::InterpolationFunction, coord_values::Vector{Vector{Float64}})
    z=zeros(Float64, length(interp.linpars))

    vcat(
        [
            [interp._DOF_analysis_functions[i, j](pt..., z...) 
                for i=1:size(interp._DOF_analysis_functions, 1), j=1:size(interp._DOF_analysis_functions, 2)]
                    for pt in coord_values
        ]...
    )
end

"""
Get interpolation function of given derivative order from interpolation object
"""
function get_interpfun(interp::InterpolationFunction; deriv_order::Union{Nothing, Vector{Int64}}=nothing)
    return interp._funcs[
        (isnothing(deriv_order) ? 1 : (deriv_order.+1))...
    ]
end
