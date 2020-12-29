module IntegralUtils

using LinearAlgebra
using Decimals

export parse_function, parse_num, Φ


"""
    ∂(f::Function, var::Symbol, P₀::Array{T, 1}; Δ::Number = 1e-3)::Union{Number, Array{Number, 1}} where T <: Number

Compute a partial derivative of `f` with respect to `var` at a point `P₀`.

# Examples
```
julia> ∂((x, y, z) -> cos(x)sin(y)exp(z), :y, [π, 0, log(2)])
-1.9999999999996667

julia> ∂((x, y, z) -> cos(x)sin(y)exp(z), :y, [π, 0, log(2)]; Δ = 1e-9)
-2.0

julia> ∂((x, y, z) -> [x, y, z], :z, [1, 1, 1])
3-element Array{Float64,1}:
 0.0
 0.0
 0.9999999999177334

julia> ∂((u, v) -> [cos(u), sin(u), v], :u, [0, 1])
3-element Array{Float64,1}:
 -5.000444502911705e-7
  0.9999999999998334
  0.0
```
"""
function ∂(f::Function, var::Symbol, P₀::Array{T1, 1}; Δ::Number = 1e-6)::Union{Number, Array{T2, 1} where T2 <: Number} where T1 <: Number
    if length(P₀) == 3
        limits = Dict(sym=>vec for (sym, vec) in
                zip((:x, :y, :z), [Δ*Diagonal(ones((3, 3)))[3i-2:3i] for i in 1:3]))
    elseif length(P₀) == 2
        limits = Dict(sym=>vec for (sym, vec) in
                zip((:u, :v), [Δ*Diagonal(ones((2, 2)))[2i-1:2i] for i in 1:2]))
    end
    return (f((P₀ + limits[var])...) - f(P₀...)) / Δ
end


"""
    divergence(F::Function, P₀::Array{T, 1})::Number where T <: Number

Calculate the divergence of a vector field `F` at a point `P₀`.

# Examples
```
julia> divergence((x, y, z) -> [z * cos(x), z * sin(y), z - cos(x+y)], [π, π, 1])
4.998225722196992e-7

julia> divergence((x, y, z) -> [z * cos(x), z * sin(y), z - cos(x+y)], [0, 0, 1])
1.9999994998731165
```
"""
function divergence(F::Function, P₀::Array{T, 1})::Number where T <: Number
   return sum(hcat([∂(F, sym, P₀) for sym in (:x, :y, :z)]...)' .* Diagonal(ones(3, 3)))
end


"""
    normal(r::Function, P₀::Array{T, 1})::Array{Number, 1} where T <: Number

Compute a normal vector to the surface parametrized by `r` at point `P₀`.

# Examples
```
julia> normal((u, v) -> [cos(u)cos(v), sin(u)cos(v), sin(v)], [π/2, 0])
3-element Array{Number,1}:
 -5.000444502910872e-7
  0.9999999999174001
  5.000444502499502e-7

julia> normal((u, v) -> [cos(u), sin(u), v], [π/2, 1])
3-element Array{Number,1}:
 -5.000444502500335e-7
 0.9999999998353001
 0.0
```
"""
function normal(r::Function, P₀::Array{T, 1})::Array{Number, 1} where T <: Number
    rᵤ = ∂(r, :u, P₀)
    rᵥ = ∂(r, :v, P₀)
    return cross(rᵤ, rᵥ)
end


"""
    transform(F::Function, r::Function, P₀::Array{T, 1})::Number where T <: Number

Convert ``F(P₀)⋅dS`` to ``F(r(P₀))dudv``.

# Examples
```
julia> transform((x, y, z) -> [x, y, z], (u, v) -> [cos(u)cos(v), sin(u)cos(v), sin(v)], [π/2, 0])
0.9999999999174001

julia> transform((x, y, z) -> [-y, x, 0], (u, v) -> [(cos(u) + 2)cos(v), (cos(u) + 2)sin(v), sin(u)], [0, 0])
-4.500400052619784e-6
```
"""
function transform(F::Function, r::Function, P₀::Array{T, 1})::Number where T <: Number
   return dot(F(r(P₀...)...), normal(r, P₀))
end


"""
    create_weights(X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function, ξ::Int, υ::Int, ζ::Int)::Array{Float64, 3}

Build a 3d matrix according to composite Simpson's rule for triple integrals.

# Examples
```
julia> create_weights((0, 1), (x, y) -> 0, (x, y) -> 1 - y - x, x -> 0, x -> 1 - x, 2, 2, 2)
3×3×3 Array{Float64,3}:
[:, :, 1] =
 1.0  1.0  0.0
 2.0  2.0  0.0
 0.0  0.0  0.0

[:, :, 2] =
 4.0  4.0  0.0
 8.0  8.0  0.0
 0.0  0.0  0.0

[:, :, 3] =
 1.0  1.0  0.0
 2.0  2.0  0.0
 0.0  0.0  0.0

 julia> create_weights((-1, 1), (x, y) -> sqrt(2-y^2-x^2), (x, y) -> 1, x -> 0, x -> sqrt(1-x^2), 4, 4, 2)
 5×5×3 Array{Float64,3}:
 [:, :, 1] =
  0.0  -1.11847  -0.828427  -1.11847  0.0
  0.0  -4.22673  -3.13553   -4.22673  0.0
  0.0  -1.73205  -1.2915    -1.73205  0.0
  0.0  -2.11231  -1.59166   -2.11231  0.0
  0.0   0.0       0.0        0.0      0.0

 [:, :, 2] =
  0.0   -4.4739    -3.31371   -4.4739   0.0
  0.0  -16.9069   -12.5421   -16.9069   0.0
  0.0   -6.9282    -5.16601   -6.9282   0.0
  0.0   -8.44925   -6.36665   -8.44925  0.0
  0.0    0.0        0.0        0.0      0.0

 [:, :, 3] =
  0.0  -1.11847  -0.828427  -1.11847  0.0
  0.0  -4.22673  -3.13553   -4.22673  0.0
  0.0  -1.73205  -1.2915    -1.73205  0.0
  0.0  -2.11231  -1.59166   -2.11231  0.0
  0.0   0.0       0.0        0.0      0.0
```
"""
function create_weights(X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function,
    ξ::Int, υ::Int, ζ::Int)::Array{Float64, 3}
    Xᵢ = LinRange(X... , ξ + 1)
    lengths = ψ.(Xᵢ) - ϕ.(Xᵢ)
    areas = reshape([(η(x, y) - ρ(x, y)) * (ψ(x) - ϕ(x)) for x in Xᵢ for y in LinRange(ϕ(x), ψ(x), υ + 1)], (ξ + 1, υ + 1))
    return vcat([1], 4 * ones(ξ - 2) - repeat([0, 2], (ξ - 2) ÷ 2), [4, 1]) .*
            vcat([1], 4 * ones(υ - 2) - repeat([0, 2], (υ - 2) ÷ 2), [4, 1])' .* areas .*
            reshape(vcat([1], 4 * ones(ζ - 2) - repeat([0, 2], (ζ - 2) ÷ 2), [4, 1]), (1, 1, ζ + 1))
end


"""
    create_weights(U::Tuple{Number, Number}, ϕ::Function, ψ::Function, μ::Int, ν::Int)::Array{Float64, 2}

Create a 2d matrix of weights according to composite simpson rule for double integrals.

# Examples
```
julia> create_weights((0, 1), x -> 0, x -> x, 6, 4)
7×5 Array{Float64,2}:
 0.0        0.0      0.0       0.0      0.0
 0.666667   2.66667  1.33333   2.66667  0.666667
 0.666667   2.66667  1.33333   2.66667  0.666667
 2.0        8.0      4.0       8.0      2.0
 1.33333    5.33333  2.66667   5.33333  1.33333
 3.33333   13.3333   6.66667  13.3333   3.33333
 1.0        4.0      2.0       4.0      1.0

julia> create_weights((0, π/4), x -> sin(x), x -> cos(x), 8, 4)
9×5 Array{Float64,2}:
 1.0           4.0          2.0           4.0          1.0
 3.58867      14.3547       7.17734      14.3547       3.58867
 1.57139       6.28556      3.14278       6.28556      1.57139
 2.66662      10.6665       5.33325      10.6665       2.66662
 1.08239       4.32957      2.16478       4.32957      1.08239
 1.6421        6.56839      3.2842        6.56839      1.6421
 0.551799      2.2072       1.1036        2.2072       0.551799
 0.554469      2.21787      1.10894       2.21787      0.554469
 1.11022e-16   4.44089e-16  2.22045e-16   4.44089e-16  1.11022e-16
```
"""
function create_weights(U::Tuple{Number, Number}, ϕ::Function, ψ::Function, μ::Int, ν::Int)::Array{Float64, 2}
    Uᵢ = LinRange(U..., μ + 1)
    lengths = ψ.(Uᵢ) - ϕ.(Uᵢ)
    return vcat([1], 4 * ones(μ - 2) - repeat([0, 2], (μ - 2) ÷ 2), [4, 1]) .*
            vcat([1], 4 * ones(ν - 2) - repeat([0, 2], (ν - 2) ÷ 2), [4, 1])' .* lengths
end


"""
    split_region(X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function, ξ::Int, υ::Int, ζ::Int)::Array{Array{Float64, 1}, 3}

Divide a 3d region into points on space.

# Examples
```
julia> split_region((0, 1), (x, y) -> 0, (x, y) -> 1, x -> 0, x -> 1, 4, 2, 2)
5×3×3 Array{Array{Float64,1},3}:
[:, :, 1] =
 [0.0, 0.0, 0.0]  [0.0, 0.5, 1.0]   [0.25, 0.0, 0.5]
 [0.0, 0.0, 0.5]  [0.0, 1.0, 0.0]   [0.25, 0.0, 1.0]
 [0.0, 0.0, 1.0]  [0.0, 1.0, 0.5]   [0.25, 0.5, 0.0]
 [0.0, 0.5, 0.0]  [0.0, 1.0, 1.0]   [0.25, 0.5, 0.5]
 [0.0, 0.5, 0.5]  [0.25, 0.0, 0.0]  [0.25, 0.5, 1.0]

[:, :, 2] =
 [0.25, 1.0, 0.0]  [0.5, 0.0, 1.0]  [0.5, 1.0, 0.5]
 [0.25, 1.0, 0.5]  [0.5, 0.5, 0.0]  [0.5, 1.0, 1.0]
 [0.25, 1.0, 1.0]  [0.5, 0.5, 0.5]  [0.75, 0.0, 0.0]
 [0.5, 0.0, 0.0]   [0.5, 0.5, 1.0]  [0.75, 0.0, 0.5]
 [0.5, 0.0, 0.5]   [0.5, 1.0, 0.0]  [0.75, 0.0, 1.0]

[:, :, 3] =
 [0.75, 0.5, 0.0]  [0.75, 1.0, 1.0]  [1.0, 0.5, 0.5]
 [0.75, 0.5, 0.5]  [1.0, 0.0, 0.0]   [1.0, 0.5, 1.0]
 [0.75, 0.5, 1.0]  [1.0, 0.0, 0.5]   [1.0, 1.0, 0.0]
 [0.75, 1.0, 0.0]  [1.0, 0.0, 1.0]   [1.0, 1.0, 0.5]
 [0.75, 1.0, 0.5]  [1.0, 0.5, 0.0]   [1.0, 1.0, 1.0]

julia> split_region((-1, 1), (x, y) -> 0, (x, y) -> sqrt(2-x^2-y^2), x -> -sqrt(1-x^2), x -> sqrt(1-x^2), 4, 4, 2)
5×5×3 Array{Array{Float64,1},3}:
[:, :, 1] =
 [-1.0, 0.0, 0.0]  [-1.0, 0.0, 1.0]  [-1.0, 0.0, 0.5]  [-0.5, -0.866025, 0.0]    [-0.5, -0.433013, 1.25]
 [-1.0, 0.0, 0.5]  [-1.0, 0.0, 0.0]  [-1.0, 0.0, 1.0]  [-0.5, -0.866025, 0.5]    [-0.5, 0.0, 0.0]
 [-1.0, 0.0, 1.0]  [-1.0, 0.0, 0.5]  [-1.0, 0.0, 0.0]  [-0.5, -0.866025, 1.0]    [-0.5, 0.0, 0.661438]
 [-1.0, 0.0, 0.0]  [-1.0, 0.0, 1.0]  [-1.0, 0.0, 0.5]  [-0.5, -0.433013, 0.0]    [-0.5, 0.0, 1.32288]
 [-1.0, 0.0, 0.5]  [-1.0, 0.0, 0.0]  [-1.0, 0.0, 1.0]  [-0.5, -0.433013, 0.625]  [-0.5, 0.433013, 0.0]

[:, :, 2] =
 [-0.5, 0.433013, 0.625]  [0.0, -1.0, 0.0]       [0.0, -0.5, 1.32288]  [0.0, 0.5, 0.661438]  [0.5, -0.866025, 0.0]
 [-0.5, 0.433013, 1.25]   [0.0, -1.0, 0.5]       [0.0, 0.0, 0.0]       [0.0, 0.5, 1.32288]   [0.5, -0.866025, 0.5]
 [-0.5, 0.866025, 0.0]    [0.0, -1.0, 1.0]       [0.0, 0.0, 0.707107]  [0.0, 1.0, 0.0]       [0.5, -0.866025, 1.0]
 [-0.5, 0.866025, 0.5]    [0.0, -0.5, 0.0]       [0.0, 0.0, 1.41421]   [0.0, 1.0, 0.5]       [0.5, -0.433013, 0.0]
 [-0.5, 0.866025, 1.0]    [0.0, -0.5, 0.661438]  [0.0, 0.5, 0.0]       [0.0, 1.0, 1.0]       [0.5, -0.433013, 0.625]

[:, :, 3] =
 [0.5, -0.433013, 1.25]  [0.5, 0.433013, 0.625]  [1.0, 0.0, 0.0]  [1.0, 0.0, 1.0]  [1.0, 0.0, 0.5]
 [0.5, 0.0, 0.0]         [0.5, 0.433013, 1.25]   [1.0, 0.0, 0.5]  [1.0, 0.0, 0.0]  [1.0, 0.0, 1.0]
 [0.5, 0.0, 0.661438]    [0.5, 0.866025, 0.0]    [1.0, 0.0, 1.0]  [1.0, 0.0, 0.5]  [1.0, 0.0, 0.0]
 [0.5, 0.0, 1.32288]     [0.5, 0.866025, 0.5]    [1.0, 0.0, 0.0]  [1.0, 0.0, 1.0]  [1.0, 0.0, 0.5]
 [0.5, 0.433013, 0.0]    [0.5, 0.866025, 1.0]    [1.0, 0.0, 0.5]  [1.0, 0.0, 0.0]  [1.0, 0.0, 1.0]
```
"""
function split_region(X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function,
    ξ::Int, υ::Int, ζ::Int)::Array{Array{Float64, 1}, 3}
    Xᵢ = LinRange(X..., ξ + 1)
    return reshape([[x, y, z] for x in Xᵢ for y in LinRange(ϕ(x), ψ(x), υ + 1) for z in LinRange(ρ(x, y), η(x, y), ζ + 1)], ξ + 1, υ + 1, :)
end


"""
    split_region(Uᵢ::Tuple{Number, Number}, ϕ::Function, ψ::Function, μ::Int, ν::Int)::Array{Array{Float64, 1}, 2}

Divide a 2d region into points on plane.

# Examples
```
julia> split_region((0, 1), x -> 0, x -> x, 6, 4)
7×5 Array{Array{Float64,1},2}:
 [0.0, 0.0]             [0.166667, 0.0833333]  [0.333333, 0.333333]  [0.666667, 0.166667]  [0.833333, 0.625]
 [0.0, 0.0]             [0.166667, 0.125]      [0.5, 0.0]            [0.666667, 0.333333]  [0.833333, 0.833333]
 [0.0, 0.0]             [0.166667, 0.166667]   [0.5, 0.125]          [0.666667, 0.5]       [1.0, 0.0]
 [0.0, 0.0]             [0.333333, 0.0]        [0.5, 0.25]           [0.666667, 0.666667]  [1.0, 0.25]
 [0.0, 0.0]             [0.333333, 0.0833333]  [0.5, 0.375]          [0.833333, 0.0]       [1.0, 0.5]
 [0.166667, 0.0]        [0.333333, 0.166667]   [0.5, 0.5]            [0.833333, 0.208333]  [1.0, 0.75]
 [0.166667, 0.0416667]  [0.333333, 0.25]       [0.666667, 0.0]       [0.833333, 0.416667]  [1.0, 1.0]

julia> split_region((0, 2), x -> -sqrt(4-x^2), x -> sqrt(4-x^2), 6, 6)
7×7 Array{Array{Float64,1},2}:
 [0.0, -2.0]       [0.333333, -1.97203]   [0.666667, -1.88562]   [1.0, -1.73205]  [1.33333, -1.49071]   [1.66667, -1.10554]   [2.0, 0.0]
 [0.0, -1.33333]   [0.333333, -1.31468]   [0.666667, -1.25708]   [1.0, -1.1547]   [1.33333, -0.993808]  [1.66667, -0.737028]  [2.0, 0.0]
 [0.0, -0.666667]  [0.333333, -0.657342]  [0.666667, -0.628539]  [1.0, -0.57735]  [1.33333, -0.496904]  [1.66667, -0.368514]  [2.0, 0.0]
 [0.0, 0.0]        [0.333333, 0.0]        [0.666667, 0.0]        [1.0, 0.0]       [1.33333, 0.0]        [1.66667, 0.0]        [2.0, 0.0]
 [0.0, 0.666667]   [0.333333, 0.657342]   [0.666667, 0.628539]   [1.0, 0.57735]   [1.33333, 0.496904]   [1.66667, 0.368514]   [2.0, 0.0]
 [0.0, 1.33333]    [0.333333, 1.31468]    [0.666667, 1.25708]    [1.0, 1.1547]    [1.33333, 0.993808]   [1.66667, 0.737028]   [2.0, 0.0]
 [0.0, 2.0]        [0.333333, 1.97203]    [0.666667, 1.88562]    [1.0, 1.73205]   [1.33333, 1.49071]    [1.66667, 1.10554]    [2.0, 0.0]
```
"""
function split_region(U::Tuple{Number, Number}, ϕ::Function, ψ::Function, μ::Int, ν::Int)::Array{Array{Float64, 1}, 2}
    Uᵢ = LinRange(U..., μ + 1)
    return reshape([[u, v] for u in Uᵢ for v in LinRange(ϕ(u), ψ(u), ν + 1)], μ + 1, :)
end


"""
    coeff(A::Tuple{Number, Number}, n::Int)::Float64

Compute a coefficient for composite Simpson's rule with `n` midpoints and interval `A`.

# Examples
```
julia> coeff((0, 1), 2)
0.16666666666666666

julia> coeff((0, 4), 4)
0.3333333333333333
```
"""
coeff(A::Tuple{Number, Number}, n::Number)::Float64 = (A[2] - A[1]) / (3 * n)


"""
    ∯(F::Function, X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function, ξ::Int, υ::Int, ζ::Int)::Float64

Determine flux of vector field `F` through a 3d region using Gauss-Ostrogradski theorem and Simpson's rule for ``ξ⋅υ⋅ζ`` nodes.

# Examples
```
julia> ∯((x, y, z) -> [x, y, z], (0, 1), (x, y) -> 0, (x, y) -> 1, x -> 0, x -> 1, 2, 2, 2)
3.0000000000163776

julia> ∯((x, y, z) -> [x^2, y^2, z^2], (-1, 1), (x, y) -> 0, (x, y) -> x^2+y^2, x -> -1, x -> 1, 6, 4, 2)
2.5488967009810803
```
"""
function ∯(F::Function, X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function, ξ::Int, υ::Int, ζ::Int)::Float64

    weights = create_weights(X, ρ, η, ϕ, ψ, ξ, υ, ζ)
    points = split_region(X, ρ, η, ϕ, ψ, ξ, υ, ζ)
    return sum(divergence.(F, points) .* weights) * prod([coeff(interval, steps)
            for (interval, steps) in zip((X, (0, 1), (0, 1)), (ξ, υ, ζ))])
end


"""
    ∯(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function, μ::Int, ν::Int)::Float64

Determine flux of vector field `F` through a 2d region using Jacobian and Simpson's rule for ``μ⋅ν`` nodes.

# Examples
```
julia> ∯((x, y, z) -> [x, y, z], (u, v) -> [cos(u), sin(u), v], (0, 2π), u -> 0, u -> u, 2, 2)
19.739208804825942

julia> ∯((x, y, z) -> [x^2, y^2, z^2], (u, v) -> [u, v, 4-u-v], (0, 4), u -> 0, u -> 4-u, 8, 10)
64.00000000410644
```
"""
function ∯(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function, μ::Int, ν::Int)::Float64
    points = split_region(U, ϕ, ψ, μ, ν)
    weights = create_weights(U, ϕ, ψ, μ, ν)
    return sum(transform.(F, r, points) .* weights') * prod(coeff(interval, steps)
            for (interval, steps) in zip((U, (0, 1)), (μ, ν)))
end


"""
    draw_points(Xᵢ::LinRange{Float64}, ϕ::Function, ψ::Function, N::Int)::Array{Array{Float64, 1}, 1}

Choose randomly 2d points from a rectangle that bounds a region.

# Examples
```
julia> draw_points(LinRange(0, 1, 11), x -> 0, x -> 1, 11)
11-element Array{Array{Float64,1},1}:
 [0.2, 0.4]
 [0.3, 0.7]
 [0.5, 1.0]
 [0.2, 0.8]
 [0.1, 0.5]
 [0.1, 0.3]
 [1.0, 0.3]
 [0.9, 0.0]
 [0.6, 0.5]
 [0.7, 0.8]
 [0.1, 1.0]

julia> draw_points(LinRange(1, 2, 5), x -> 0, x -> log(x), 5)
5-element Array{Array{Float64,1},1}:
 [2.0, 0.0]
 [1.25, 0.34657359027997264]
 [1.5, 0.17328679513998632]
 [1.75, 0.17328679513998632]
 [1.5, 0.0]
```
"""
function draw_points(Xᵢ::LinRange{Float64}, ϕ::Function, ψ::Function, N::Int)::Array{Array{Float64, 1}, 1}
    Yᵢ = LinRange(min(ϕ.(Xᵢ)...), max(ψ.(Xᵢ)...), N)
    return [[rand(Xᵢ), rand(Yᵢ)] for i in 1:N]
end


"""
    draw_points(Xᵢ::LinRange{Float64}, Yᵢ::LinRange{Float64}, ρ::Function, η::Function, N::Int)::Array{Array{Float64, 1}, 1}

Choose randomly 3d points from a cuboid that bounds a region.

# Examples
```
julia> draw_points(LinRange(0, 1, 5), LinRange(0, 1, 5), (x, y) -> 0, (x, y) -> 1, 11)
11-element Array{Array{Float64,1},1}:
 [0.0, 0.5, 0.0]
 [0.0, 0.75, 0.2]
 [0.5, 0.0, 0.8]
 [0.5, 0.25, 0.1]
 [0.75, 0.25, 0.3]
 [0.25, 0.5, 0.5]
 [0.25, 0.25, 0.4]
 [0.5, 0.25, 0.4]
 [0.5, 0.75, 0.4]
 [0.25, 0.0, 1.0]
 [0.75, 1.0, 0.1]

julia> draw_points(LinRange(0, 1, 5), LinRange(0, 1, 5), (x, y) -> 0, (x, y) -> sqrt(x+y), 11)
11-element Array{Array{Float64,1},1}:
 [0.75, 0.0, 0.28284271247461906]
 [0.75, 0.25, 0.4242640687119285]
 [0.5, 0.5, 0.5656854249492381]
 [1.0, 0.5, 1.2727922061357857]
 [0.5, 0.0, 1.1313708498984762]
 [0.5, 1.0, 0.5656854249492381]
 [0.75, 0.0, 1.1313708498984762]
 [0.75, 0.0, 0.848528137423857]
 [1.0, 0.5, 0.7071067811865476]
 [0.75, 1.0, 1.4142135623730951]
 [0.5, 1.0, 0.848528137423857]
```
"""
function draw_points(Xᵢ::LinRange{Float64}, Yᵢ::LinRange{Float64}, ρ::Function, η::Function, N::Int)::Array{Array{Float64, 1}, 1}
    unzip(f, P) = f(P...)
    Zᵢ = LinRange(min(unzip.(ρ, zip(Xᵢ, Yᵢ))...), max(unzip.(η, zip(Xᵢ, Yᵢ))...), N)
    return [[rand(Xᵢ), rand(Yᵢ), rand(Zᵢ)] for i in 1:N]
end


"""
    ∯(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function, N::Int)::Float64

Determine flux of vector field `F` through a 3d region using Gauss-Ostrogradski theorem and Monte Carlo method for `N` random points.

# Examples
```
julia> ∯((x, y, z) -> [x, y, z], (0, 1), (x, y) -> 0, (x, y) -> 1, x -> 0, x -> 1, 10000)
3.00000000002332

julia> ∯((x, y, z) -> [x^2, y^2, z^2], (-1, 1), (x, y) -> 0, (x, y) -> x^2+y^2, x -> -1, x -> 1, 1000000)
2.4837039282914986
```
"""
function ∯(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function, N::Int)::Float64
    Uᵢ = LinRange(U..., N)
    points = draw_points(Uᵢ, ϕ, ψ, N)
    validate(P::Array{Float64, 1}) = ϕ(P[1]) <= P[2] <= ψ(P[1]) ? transform(F, r, P) : 0.0
    return sum(validate.(points)) * (U[2] - U[1]) * (max(ψ.(Uᵢ)...) - min(ϕ.(Uᵢ)...)) / N
end


"""
    ∯(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function, N::Int)::Float64

Determine flux of vector field `F` through a 2d region using Jacobian and Monte Carlo method for `N` random points.

# Examples
```
julia> ∯((x, y, z) -> [x, y, z], (u, v) -> [cos(u), sin(u), v], (0, 2π), u -> 0, u -> u, 10000)
19.755000172505373

julia> ∯((x, y, z) -> [x^2, y^2, z^2], (u, v) -> [u, v, 4-u-v], (0, 4), u -> 0, u -> 4-u, 1000000)
64.06677011718612
```
"""
function ∯(F::Function, X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function, N::Int)::Float64
    Xᵢ = LinRange(X..., N)
    Yᵢ = LinRange(min(ϕ.(Xᵢ)...), max(ψ.(Xᵢ)...), N)
    unzip(f, P) = f(P...)
    points = draw_points(Xᵢ, Yᵢ, ρ, η, N)
    validate(P::Array{Float64, 1}) = ϕ(P[1]) <= P[2] <= ψ(P[1]) && ρ(P[1:2]...) <= P[3] <= η(P[1:2]...) ? divergence(F, P) : 0.0
    return sum(validate.(points)) * (X[2] - X[1]) * (max(Yᵢ...) - min(Yᵢ...)) * (max(unzip.(η, zip(Xᵢ, Yᵢ))...) - min(unzip.(ρ, zip(Xᵢ, Yᵢ))...)) / N
end


"""
    round_float(α::Float64, ϵ::Float64)::Union{Float64, Int}

Round `α` to the same number of digits as `ϵ` or nearest integer with `ϵ` accuracy.

# Examples
```
julia> round_float(π, 1e-12)
3.14159265359

julia> round_float(exp(1), 1e-3)
2.718

julia> round_float(1.999, 1e-3)
2
```
"""
round_float(α::Float64, ϵ::Float64)::Union{Int, Float64} = abs(α - round(α)) < ϵ ? Int(round(α)) : round(α, digits = abs(Decimal(ϵ).q))


"""
    arguments(args::Symbol...)

Parse arguments' symbols to a string.

# Examples
```
julia> arguments(:u, :v)
"(u, v)"

julia> arguments(:x, :y, :z)
"(x, y, z)"
```
"""
arguments(args::Symbol...)::String = '(' * String(args[1]) * reduce(*, ", " .* String.(args[2:end])) * ')'


"""
    parse_function(name::Symbol, body::String, args::Symbol...)::Function

Parse a body of a function and arguments' and function name's symbols to the actual function.

# Examples
```
julia> parse_function(:f, "cos(x)sin(y)+exp(z)", :x, :y, :z)
f (generic function with 1 method)

julia> parse_function(:g, "u^2-cos(exp(1/v))", :u, :v)
g (generic function with 1 method)
```
"""
parse_function(body::String, args::Symbol...)::Function = eval(Meta.parse(arguments(args...) * "->" * body))


"""
    parse_function(name::Symbol, body::String, args::Symbol)::Function

Parse a body of a function and argument and function name's symbols to the actual function.

# Examples
```
julia> parse_function(:f, "2cos(x)", :x)
f (generic function with 1 method)

julia> parse_function(:g, "u^2-7u", :u)
g (generic function with 1 method)
```
"""
parse_function(body::String, args::Symbol)::Function = eval(Meta.parse(string(args) * "->" * body))


"""
    parse_num(value::String)::Float64

Parse a number or an expression from string to float.

# Examples
```
julia> parse_num("2pi")
6.283185307179586

julia> parse_num("2pi - 6/3 + 7.2")
11.483185307179586
```
"""
parse_num(value::String)::Float64 = eval(Meta.parse(value))


"""
    Φ(F::Function, X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function; ϵ::Number = 1e-3, technique::String = "Simpson")::Union{Int, Float64}

Determine a flux of a vector field `F` through a surface.

# Examples
```
julia> Φ((x, y, z) -> [x, y, z], (0, 1), (x, y) -> 0, (x, y) -> 1, x -> 0, x -> 1)
3

julia> Φ((x, y, z) -> [-y, x, 0], (-1, 1), (x, y) -> 0, (x, y) -> sqrt(2-y^2-x^2), x -> -sqrt(1-x^2), x -> sqrt(1-x^2); technique = "Monte Carlo", ϵ = 0.1)
0
```
"""
function Φ(F::Function, X::Tuple{Number, Number}, ρ::Function, η::Function, ϕ::Function, ψ::Function;
     ϵ::Number = 1e-3, technique::String = "Simpson")::Union{Int, Float64}
    if technique in ("Simpson", "Monte Carlo")
        surface_integral(n) = technique == "Simpson" ? ∯(F, X, ρ, η, ϕ, ψ, 2n, 2n, 2n) : ∯(F, X, ρ, η, ϕ, ψ, 100n)
    else
        throw(ArgumentError("Invalid technique name."))
    end
    n = 1
    Φₙ = surface_integral(n)
    n += 1
    Φₙ₊₁ = surface_integral(n)
    while abs(Φₙ₊₁ - Φₙ) > ϵ
        n += 1
        Φₙ, Φₙ₊₁ = Φₙ₊₁, surface_integral(n)
    end
    return round_float(Φₙ₊₁, ϵ)
end


"""
    Φ(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function; ϵ::Number = 1e-3, technique::String = "Simpson")::Union{Int, Float64}

Determine a flux of a vector field `F` through a surface parametrized by `r`.

# Examples
```
julia> Φ((x, y, z) -> (x^2, y^2, z^2), (u, v) -> [u, v, 4-u-v], (0, 4), u -> 0, u -> 4-u)
64

julia> Φ((x, y, z) -> [x-y, y-z, z-x], (u, v) -> [(cos(u)+2)cos(v), (cos(u)+2)sin(v), sin(u)], (0, 2π), u -> 0, u -> 2π; technique = "Monte Carlo", ϵ = 0.01)
-118.25
```
"""
function Φ(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function;
     ϵ::Number = 1e-3, technique::String = "Simpson")::Union{Int, Float64}
    if technique in ("Simpson", "Monte Carlo")
        surface_integral(n::Int) = technique == "Simpson" ? ∯(F, r, U, ϕ, ψ, 2n, 2n) : ∯(F, r, U, ϕ, ψ, 100n)
    else
        throw(ArgumentError("Invalid technique name."))
    end
    n = 1
    Φₙ = surface_integral(n)
    n += 1
    Φₙ₊₁ = surface_integral(n)
    while abs(Φₙ₊₁ - Φₙ) > ϵ
        n += 1
        Φₙ, Φₙ₊₁ = Φₙ₊₁, surface_integral(n)
    end
    return round_float(Φₙ₊₁, ϵ)
end


end
