@assert begin
    bf=get_basefun(
        :(a+x*b+y*c),
        [
            :x, :y
        ],
        [
            :a,
            :b,
            :c
        ]
    )

    isapprox(
        bf(1.0, 2.0, -1.0, 1.0, 2.0), 4.0
    )
end

@assert begin
    bfd=get_basefun_derivative(
        :(a+x*b+y*c),
        [
            :x, :y
        ],
        [
            :a,
            :b,
            :c
        ],
        [1, 0]
    )

    isapprox(
        bfd(
            1.0, -2.0, 4.0, 1.0, -3.0
        ),
        1.0
    )
end

@assert begin
    blpc=get_linpar_contribution(
        :(a+b*x+c*y),
        [
            :x,
            :y
        ],
        [
            :a,
            :b,
            :c
        ],
        :c,
        [0, 0]
    )

    isapprox(
        blpc(
            1.0, -2.0, 4.0, 1.0, -3.0
        ),
        -2.0
    )
end
