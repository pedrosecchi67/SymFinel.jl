"""
Example modelling a linear beam element

Formulation:

K=-E*A*[dw_i!dx*du_j!dx]
"""

# local coordinates:
local_interp=InterpolationFunction(
    :(a0+a1*xi), # interpolation function expression
    [:xi], # local coordinate system degrees of freedom
    [:a0, :a1], # linear coefficients for regression
    [
        [0] # conditions to impose with regression: absolute value (order 0 derivative in the first axis)
    ];
    maximum_derivative=1
)

# absolute coordinates
abs_interp=InterpolationFunction(
    :(a0+a1*x), # interpolation function expression
    [:x], # absolute coordinate system degrees of freedom
    [:a0, :a1], # linear coefficients for regression
    [
        [0] # conditions to impose with regression: absolute value (order 0 derivative)
    ];
    maximum_derivative=1
)

# defining an integration domain:
dmn=Domain(
    local_interp, # local coordinate system
    [
        [0.0, 1.0] # boundaries for local coordinates: xi between 0 and 1
    ],
    [
        [0.0], # positions at which to impose values for interpolation:
        [1.0] # xi=0, xi=1
    ]
)

# defining a residual for within the beam:
domain_resd=Residual(
    [
        ( # variable w (weight), interpolated by InterpolationFunction abs_interp,
            :w, abs_interp, [
                (:dw!dx, [1]) # accompanied by its first order derivative
            ]
        )
        ( # variable u (displacement), interpolated by InterpolationFunction abs_interp,
            :u, abs_interp, [
                (:du!dx, [1]) # accompanied by its first order derivative
            ]
        )
    ],
    :(-E*(dw!dx*du!dx)),
    3; # Gauss-Legendre quadrature order
    extra_args=[:E, :A] # additional argument: elasticity
)

dmn_resd_func=get_volume_residual_function(domain_resd)

# defining a specific finite element:
fin=Finel(
    [
        [0.0],
        [1.0] # defined by corner points at x=0 and x=1 (corresponding to xi=0 and xi=1 in local coordinates: see Domain)
    ],
    dmn
)

# obtaining elasticity matrix:
K=zeros(Float64, 2, 2)

ws=zeros(Float64, 2)
us=zeros(Float64, 2)

E=1.0
A=1.0

for i=1:2
    ws[i]=1.0

    for j=1:2
        us[j]=1.0

        variable_vals=[ws, us]
        K[i, j]=dmn_resd_func(variable_vals, fin, E, A)

        us[j]=0.0
    end

    ws[i]=0.0
end

@assert isapprox(
    K,
    [
        -1.0 1.0;
        1.0 -1.0
    ]
)
