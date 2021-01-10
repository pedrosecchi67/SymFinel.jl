# SymFinel.jl

`SymFinel.jl` is a package written to ease the programming of versatile Generalized Finite Element Method formulations for any element geometry.

It allows you to generate integration functions for equations in finite elements of any geometry, number of dimensions, order and interpolation function, without having to rewrite your code to change these details once the residual definition is made.

**Use SymFinel if you want to enjoy Julia's metaprogramming capabilities to generate finite element residuals for a great variaty of formulations, with a small amout of changes to your code.**

## An Example

To obtain the residual of $\partial u/\partial x-v=0$, in a triangular, simplex element using `SymFinel`, you can follow the steps below.

### Interpolation Functions

We'll need an interpolation function in the local coordinate system and one in the global coordinate system. To define them, we can use:

```
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
```

* The first argument is an expression for the interpolation function;
* The second argument lists coordinate symbols;
* The third argument lists symbols for coefficients used to interpolate desired variables over a domain;
* The fourth and last variable indicates degrees of freedom used for obtaining an interpolation. For example:

```
DOFs=[
    [0, 0],
    [1, 0],
    [0, 1]
]
```

would incurr that input such as the following was used for interpolation for variable `u`:

```
us=[
    u_1, du!dx_1, du!dy_1, # first point
    u_2, du!dx_2, du!dy_2, # second point
    # ... goes on until enough impositions to 
    # determine all coefficients are provided
]
```

Since, for this example, only first order elements will be used, `[[0, 0]]`should be enough.

### Defining a Domain

We can define our integration domain as reproduced below:

```
dmn=Domain(
    loc_coords, # local coordinates' InterpolationFunction struct
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
```

* The first argument is an `InterpolationFunction` struct with a definition for our local coordinate system;
* The second argument is a list of tuples defining boundaries for each coordinate, as expressions or literals;
* The last argument indicates points at which coordinates (in a format compatible with the one described in the "Interpolation Functions" section) will be imposed for interpolation - in this case, the positions of element corners in local coordinates.

### Creating Finite Elements

To create a finite element with corners $(x, y)=(1, 1)$, $(2, 1)$ and $(2, 2)$, we can use:

```
fin=Finel(
    [ # corners
        [1.0, 1.0],
        [2.0, 1.0],
        [2.0, 2.0]
    ],
    dmn # domain of integration
)
```

### Defining Residuals

To obtain the following integral:

$$
\int w(\partial u/\partial x-v) d\Omega
$$

In which $w$ is a weight function, we can use function `resfun`, defined below:

```
resd=Residual(
    [
        ( # first variable: weight function
            :w, abs_interp, [] # no derivatives calculated
        ),
        ( # second variable: u
            :u, abs_interp, [
                (:du!dx, [1, 0]) # first order derivative in x axis calculated
            ]
        )
    ],
    :(w*(du!dx-v)), # expression for the residual
    3; # Gauss-Legendre quadrature order
    extra_args=[:v] # v, a constant
)

resfunc=get_volume_residual_function(resd)
```

The argument structure is reasonably self-explanatory. For any reference, however, check the docstring for `Residual`! :)

Function `resfunc` can be evaluated with:

```
ws=[0.0, 1.0, 1.0] # format compatible with abs_coords InterpolationFunction,
us=[0.5, 1.0, 1.0] # as detailed before
v=-1.0

R=resfunc([ws, us], fin, v) # and any other extra_args symbols

# Test it!! :)
@assert isapprox(R, 0.5)
```

For some formulations, one may need to obtain integrals such as:

$$
\int_{\partial \Omega} w (\partial u/\partial x, v) \cdot \hat{n} d\Gamma
$$

Evaluated at the frontier of the domain, $\hat{n}$ being an outward facing vector. One can do this by defining a `Residual` function with a vectorial expression:

```
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
    :(w.*[du!dx, v]), # expression for the residual (now a vector!!)
    3; # Gauss-Legendre quadrature order, now at domain faces
    extra_args=[:v] # v, a constant
)

surface_resfunc=get_surface_residual_function(surface_resd)

R=surface_resfunc([ws, us], fin, v)

# test it!! :)
@assert isapprox(R, 0.25)
```

## Installation

SymFinel should be installed from git using:

```
] add https://github.com/pedrosecchi67/SymFinel.jl
```

## Troubleshooting

Make sure all expressions passed to `Domain` and `InterpolationFunction` are friendly to function `Calculus.differentiate`, from Calculus.jl.

If any other problem is experienced, please leave an issue statement at our repository page.
