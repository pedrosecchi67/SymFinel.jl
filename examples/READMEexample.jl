"""
Background for the example in the README
"""

loc_coords=InterpolationFunction(
    :(a0+a1*xi+a2*eta), # expression
    [:xi, :eta], # symbols for coordinates
    [:a0, :a1, :a2], # symbols for interpolation coefficients
    [ # degrees of freedom
        [0, 0] # values at corners (see explanation below)
    ];
    maximum_derivative=1 # we'll only really need first order derivatives
)

abs_coords=InterpolationFunction(
    :(a0+a1*x+a2*y),
    [:x, :y],
    [:a0, :a1, :a2],
    [
        [0, 0]
    ];
    maximum_derivative=1
)

dmn=Domain(
    loc_coords, # local coordinates InterpolationFunction struct
    [ # frontier of the domain in each coordinate of loc_coords
        ( # xi
            0.0, 1.0
        ),
        ( # eta
            0.0, :(1.0-xi) # linear triangular element
        )
    ],
    [ # points at which the interpolation from local to absolute coordinates
        [0.0, 0.0], # is imposed (triangle corners)
        [1.0, 0.0],
        [0.0, 1.0] # (xi, eta)
    ]
)

fin=Finel(
    [ # corners
        [1.0, 1.0],
        [2.0, 1.0],
        [2.0, 2.0]
    ],
    dmn # domain of integration
)

resd=Residual(
    [
        ( # first variable: weight function
            :w, abs_coords, [] # no derivatives calculated
        ),
        ( # second variable: u
            :u, abs_coords, [
                (:du!dx, [1, 0]) # first order derivative in x axis calculated
            ]
        )
    ],
    :(w*(du!dx-v)), # expression for the residual
    3; # Gauss-Legendre quadrature order
    extra_args=[:v] # v, a constant
)

resfunc=get_volume_residual_function(resd)

ws=[0.0, 1.0, 1.0] # format compatible with abs_coords InterpolationFunction,
us=[0.5, 1.0, 1.0] # as detailed before
v=-1.0

R=resfunc([ws, us], fin, v) # and any other extra_args symbols

# Test it!! :)
@assert isapprox(R, 0.5)

surface_resd=Residual(
    [
        ( # first variable: weight function
            :w, abs_coords, [] # no derivatives calculated
        ),
        ( # second variable: u
            :u, abs_coords, [
                (:du!dx, [1, 0]) # first order derivative in x axis calculated
            ]
        )
    ],
    :(w.*[du!dx, v]), # expression for the residual
    3; # Gauss-Legendre quadrature order, now at domain faces
    extra_args=[:v] # v, a constant
)

surface_resfunc=get_surface_residual_function(surface_resd)

R=surface_resfunc([ws, us], fin, v)

# test it!! :)
@assert isapprox(R, 0.25)
