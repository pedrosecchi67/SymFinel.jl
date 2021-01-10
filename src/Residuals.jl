export Residual, get_volume_residual_function, get_surface_residual_function

using FastGaussQuadrature

"""
Function to obtain Jacobian from `Finel` struct
"""
function from_finel_jac!(Jac::Matrix{Float64}, fin::Finel, x...)
    deriv_order=zeros(Int64, length(x))
    for j=1:length(x)
        deriv_order[j]=1

        f=get_interpfun(fin.domain.interp; deriv_order=deriv_order)

        for (i, coefs) in enumerate(fin._coord_linpars)
            Jac[i, j]=f(x..., coefs...)
        end

        deriv_order[j]=0
    end
end

"""
Function to return polynomial coefficients from array notation
"""
function get_polycoefs(vars::Tuple{Vector{Float64}, InterpolationFunction}...)
    [
        regress(finel_get_regr(fin, intp), vals) for (vals, intp) in vars
    ]
end

"""
Definition for volume residual obtention configurations struct

* `variables`: vector of tuples. First element indicates a symbol for a variable; the second, an `InterpolationFunction`
    struct to serve as its interpolator. The third should be a vector of tuples, each with a symbol and a vecor of ints
    indicating derivative order in each axis
    
    Example:
        ```[
            (
                :u, interp1, [
                    (:du!dx, [1, 0]), (du!dy, [0, 1])
                ]
            ),
            (
                :v, interp2, []
            )
        ]```
* `ex`: expression evaluating the residuals.

    Examples:
        * `:(u*du!dx+v-1.0)`
        * `begin w=sqrt(u); w*v-du!dx end`
* `extra_args`: symbols for any additional arguments to be passed
"""
mutable struct Residual
    variables::Vector{
        Tuple{
            Symbol, InterpolationFunction, Vector{
                Tuple{Symbol, Vector{Int64}}
            }
        }
    }
    ex::Expr
    _nvr::Int64
    _weights::Vector{Float64}
    _abscissas::Vector{Float64}
    extra_args::Vector{Symbol}
end

"""
Constructor for Residual taking only public values

* `variables`: vector of tuples. First element indicates a symbol for a variable; the second, an `InterpolationFunction`
    struct to serve as its interpolator. The third should be a vector of tuples, each with a symbol and a vecor of ints
    indicating derivative order in each axis
    
    Example:
        ```[
            (
                :u, interp1, [
                    (:du!dx, [1, 0]), (du!dy, [0, 1])
                ]
            ),
            (
                :v, interp2, []
            )
        ]```
* `ex`: expression evaluating the residuals.

    Examples:
        * `:(u*du!dx+v-1.0)`
        * `begin w=sqrt(u); w*v-du!dx end`

* `nquad`: Gaussian quadrature order (Gauss Legendre)
* `extra_args`: symbols for any extra arguments to be passed for residual computation
"""
function Residual(
    variables,
    ex::Expr,
    nquad::Int64;
    extra_args=[]
)
    _nvr=length(variables)

    for (_, _, vv) in variables
        _nvr+=length(vv)
    end

    x, w=gausslegendre(nquad)
    x=(1.0.+x)./2
    w./=2

    Residual(
        variables,
        ex,
        _nvr,
        w,
        x,
        extra_args
    )
end

"""
Function to obtain function for punctual evaluation volume residual

* `resd`: an `Residual` instance

* return: a function (in the format `resfun(v1, dv1!dx, dv2!dy, v2...)` evaluating the residual given other variables)
"""
function get_punctual_resfunc(resd::Residual)
    pronames=Vector{Symbol}()

    for v in resd.variables
        push!(pronames, v[1])
        
        for ds in v[3]
            push!(pronames, ds[1])
        end
    end

    q=quote
        bfun=($(pronames...), $(resd.extra_args...)) -> ( $(resd.ex) )
    end

    try
        eval(q)

        return bfun
    catch
        throw(
            error(
                "get_punctual_resfunc:ERROR:unable to compile punctual residual function for expression $(resd.ex)"
            )
        )
    end
end

"""
Transform vectors of values into vectors of absolute coordinate polynomial coefficients
"""
function vals_to_polys(fin::Finel, resd::Residual, vals::Vector{Vector{Float64}})
    return [
        regress(finel_get_regr(fin, vintp[2]), vs) for (vs, vintp) in zip(vals, resd.variables)
    ]
end

"""
Obtains variables of interest given absolute coordinates and polynomial coefficients
"""
function polys_to_vars(resd::Residual, pols::Vector{Vector{Float64}}, pos::Vector{Float64})
    vs=Vector{Float64}(undef, resd._nvr)

    nv=1
    for (v, p) in zip(resd.variables, pols)
        f=get_interpfun(v[2])
        
        vs[nv]=f(pos..., p...)
        nv+=1

        for dv in v[3]
            f=get_interpfun(v[2]; deriv_order=dv[2])

            vs[nv]=f(pos..., p...)
            nv+=1
        end
    end

    vs
end

"""
Function to produce local coordinates given percentages in the space between restrictions
"""
function loc_coordinate_inrestr!(dmn::Domain, percs::Vector{Float64}, coords::Vector{Float64}; exempt::Union{Int64, Nothing}=nothing)
    w=1.0

    for (i, p) in enumerate(percs)
        br, fr=dmn._restriction_funcs[i][1](coords...), dmn._restriction_funcs[i][2](coords...)

        if i!=exempt
            w*=(fr-br)
        end

        coords[i]=(1.0-p)*br+p*fr
    end

    w
end

"""
Get volume evaluation of residual
"""
function get_volume_residual_function(resd::Residual)
    locfunc=get_punctual_resfunc(resd)
    nquad=length(resd._abscissas)

    volres=(vvals::Vector{Vector{Float64}}, fin::Finel, args...) -> begin
        naux=length(fin.domain.interp.coords)

        pols=vals_to_polys(fin, resd, vvals)
        rngs=[1:nquad for i=1:naux]
        locpercs=Vector{Float64}(undef, naux)
        locpos=Vector{Float64}(undef, naux)

        nd=length(fin._coord_linpars)

        Jac=Matrix{Float64}(undef, nd, naux)
        pos=Vector{Float64}(undef, nd)

        R=0.0

        for cinds in Iterators.product(rngs...)
            w=1.0

            for (i, c) in enumerate(cinds)
                w*=resd._weights[c]

                locpercs[i]=resd._abscissas[c]
            end

            w*=loc_coordinate_inrestr!(fin.domain, locpercs, locpos)

            get_abs_coords!(pos, fin, locpos)

            from_finel_jac!(Jac, fin, locpos...)

            w*=orthonormal_det!(Jac)

            R+=locfunc(polys_to_vars(resd, pols, pos)..., args...)*w
        end

        R
    end

    volres
end

"""
Get surface evaluation of residual
"""
function get_surface_residual_function(resd::Residual) # havent changed anything yet
    locfunc=get_punctual_resfunc(resd)
    nquad=length(resd._abscissas)

    surfres=(vvals::Vector{Vector{Float64}}, fin::Finel, args...) -> begin
        naux=length(fin.domain.interp.coords)

        pols=vals_to_polys(fin, resd, vvals)
        rngs=[1:nquad for i=1:(naux-1)]
        locpercs=Vector{Float64}(undef, naux)
        locpos=Vector{Float64}(undef, naux)

        nd=length(fin._coord_linpars)

        Jac=Matrix{Float64}(undef, nd, naux)
        pos=Vector{Float64}(undef, nd)

        R=0.0

        for dim=1:length(fin.domain.interp.coords)
            for nb=1:2
                for cinds in Iterators.product(rngs...)
                    w=1.0

                    j=1
                    for i=1:length(fin.domain.interp.coords)
                        if i!=dim
                            c=cinds[j]

                            w*=resd._weights[c]

                            locpercs[i]=resd._abscissas[c]

                            j+=1
                        else
                            locpercs[i]=(
                                nb==1 ?
                                    0.0 : 1.0
                            )
                        end
                    end

                    w*=loc_coordinate_inrestr!(fin.domain, locpercs, locpos; exempt=dim)

                    get_abs_coords!(pos, fin, locpos)

                    from_finel_jac!(Jac, fin, locpos...)

                    restr_grad=fin.domain._restriction_derivatives[dim][nb](locpos...)

                    locvals=locfunc(polys_to_vars(resd, pols, pos)..., args...)
                    shat=normal_volspace!(Jac, dim, restr_grad, nb==1)

                    R+=dot(locvals, 
                        shat)*w
                end
            end
        end

        R
    end

    surfres
end
