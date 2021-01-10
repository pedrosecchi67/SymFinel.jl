include("../src/SymFinel.jl")

using .SymFinelTests

include("ExprTests.jl")
include("InterpTests.jl")
include("DomainTests.jl")
include("FinelTests.jl")
include("ResidualTest.jl")

include("../examples/AxialBeamElement.jl")
include("../examples/READMEexample.jl")
