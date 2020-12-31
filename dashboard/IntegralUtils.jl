module IntegralUtils

using LinearAlgebra
using Decimals

export parse_function, parse_num, Φ


"""
    ∂(f::Function, var::Symbol, P₀::Array{T, 1}; Δ::Real = 1e-3)::Union{Real, Array{Real, 1}} where T <: Real

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
function ∂(f::Function, var::Symbol, P₀::Array{T1, 1}; Δ::Real = 1e-6)::Union{Real, Array{T2, 1} where T2 <: Real} where T1 <: Real
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
    divergence(F::Function, P₀::Array{T, 1})::Real where T <: Real

Calculate the divergence of a vector field `F` at a point `P₀`.

# Examples
```
julia> divergence((x, y, z) -> [z * cos(x), z * sin(y), z - cos(x+y)], [π, π, 1])
4.998225722196992e-7

julia> divergence((x, y, z) -> [z * cos(x), z * sin(y), z - cos(x+y)], [0, 0, 1])
1.9999994998731165
```
"""
function divergence(F::Function, P₀::Array{T, 1})::Real where T <: Real
   return sum(hcat([∂(F, sym, P₀) for sym in (:x, :y, :z)]...)' .* Diagonal(ones(3, 3)))
end


"""
    normal(r::Function, P₀::Array{T, 1})::Array{Real, 1} where T <: Real

Compute a normal vector to the surface parametrized by `r` at point `P₀`.

# Examples
```
julia> normal((u, v) -> [cos(u)cos(v), sin(u)cos(v), sin(v)], [π/2, 0])
3-element Array{Real,1}:
 -5.000444502910872e-7
  0.9999999999174001
  5.000444502499502e-7

julia> normal((u, v) -> [cos(u), sin(u), v], [π/2, 1])
3-element Array{Real,1}:
 -5.000444502500335e-7
 0.9999999998353001
 0.0
```
"""
function normal(r::Function, P₀::Array{T, 1})::Array{Real, 1} where T <: Real
    rᵤ = ∂(r, :u, P₀)
    rᵥ = ∂(r, :v, P₀)
    return cross(rᵤ, rᵥ)
end


"""
    transform(F::Function, r::Function, P₀::Array{T, 1})::Real where T <: Real

Convert ``F(P₀)⋅dS`` to ``F(r(P₀))dudv``.

# Examples
```
julia> transform((x, y, z) -> [x, y, z], (u, v) -> [cos(u)cos(v), sin(u)cos(v), sin(v)], [π/2, 0])
0.9999999999174001

julia> transform((x, y, z) -> [-y, x, 0], (u, v) -> [(cos(u) + 2)cos(v), (cos(u) + 2)sin(v), sin(u)], [0, 0])
-4.500400052619784e-6
```
"""
function transform(F::Function, r::Function, P₀::Array{T, 1})::Real where T <: Real
   return dot(F(r(P₀...)...), normal(r, P₀))
end


function try_eval(f::Function, args::Float64...)::Float64
    try
        return f(args...)
    catch e
        return 0.0
    end
end


"""
    create_weights(X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function, ξ::Int64, υ::Int64, ζ::Int64)::Array{Float64, 3}

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
function create_weights(X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function,
    ξ::Int64, υ::Int64, ζ::Int64)::Array{Float64, 3}
    Xᵢ = LinRange(X... , ξ + 1)
    areas = reshape([try_eval((a, b) -> η(a, b) - ρ(a, b), x, y) * try_eval(a -> ψ(a) - ϕ(a), x) for x in Xᵢ for y in LinRange(ϕ(x), ψ(x), υ + 1)], (1, υ + 1, ξ + 1))
    areas[areas .< 0] .= 0.0
    return vcat([1], 4 * ones(ξ - 2) - repeat([0, 2], (ξ - 2) ÷ 2), [4, 1])  .*
            vcat([1], 4 * ones(υ - 2) - repeat([0, 2], (υ - 2) ÷ 2), [4, 1])' .*
            reshape(vcat([1], 4 * ones(ζ - 2) - repeat([0, 2], (ζ - 2) ÷ 2), [4, 1]), (1, 1, ζ + 1)) .* areas
end


"""
    create_weights(U::Tuple{Real, Real}, ϕ::Function, ψ::Function, μ::Int64, ν::Int64)::Array{Float64, 2}

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
function create_weights(U::Tuple{Real, Real}, ϕ::Function, ψ::Function, μ::Int64, ν::Int64)::Array{Float64, 2}
    Uᵢ = LinRange(U..., μ + 1)
    lengths = ψ.(Uᵢ) - ϕ.(Uᵢ)
    return vcat([1], 4 * ones(μ - 2) - repeat([0, 2], (μ - 2) ÷ 2), [4, 1]) .*
            vcat([1], 4 * ones(ν - 2) - repeat([0, 2], (ν - 2) ÷ 2), [4, 1])' .* lengths
end


"""
    split_region(X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function, ξ::Int64, υ::Int64, ζ::Int64)::Array{Array{Float64, 1}, 3}

Divide a 3d region Int64o points on space.

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
function split_region(X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function,
    ξ::Int64, υ::Int64, ζ::Int64)::Array{Array{Float64, 1}, 3}
    Xᵢ = LinRange(X..., ξ + 1)
    return reshape([[x, y, z] for x in Xᵢ for y in LinRange(ϕ(x), ψ(x), υ + 1) for z in LinRange(try_eval(ρ, x, y), try_eval(η, x, y), ζ + 1)], :, υ + 1, ξ + 1)
end


"""
    split_region(Uᵢ::Tuple{Real, Real}, ϕ::Function, ψ::Function, μ::Int64, ν::Int64)::Array{Array{Float64, 1}, 2}

Divide a 2d region Int64o points on plane.

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
function split_region(U::Tuple{Real, Real}, ϕ::Function, ψ::Function, μ::Int64, ν::Int64)::Array{Array{Float64, 1}, 2}
    Uᵢ = LinRange(U..., μ + 1)
    return reshape([[u, v] for u in Uᵢ for v in LinRange(ϕ(u), ψ(u), ν + 1)], μ + 1, :)
end


"""
    coeff(A::Tuple{Real, Real}, n::Int64)::Float64

Compute a coefficient for composite Simpson's rule with `n` midpoints and interval `A`.

# Examples
```
julia> coeff((0, 1), 2)
0.16666666666666666

julia> coeff((0, 4), 4)
0.3333333333333333
```
"""
coeff(A::Tuple{Real, Real}, n::Real)::Float64 = (A[2] - A[1]) / (3 * n)


"""
    ∯(F::Function, X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function, ξ::Int64, υ::Int64, ζ::Int64)::Float64

Determine flux of vector field `F` through a 3d region using Gauss-Ostrogradski theorem and Simpson's rule for ``ξ⋅υ⋅ζ`` nodes.

# Examples
```
julia> ∯((x, y, z) -> [x, y, z], (0, 1), (x, y) -> 0, (x, y) -> 1, x -> 0, x -> 1, 2, 2, 2)
3.0000000000163776

julia> ∯((x, y, z) -> [x^2, y^2, z^2], (-1, 1), (x, y) -> 0, (x, y) -> x^2+y^2, x -> -1, x -> 1, 6, 4, 2)
2.5488967009810803
```
"""
function ∯(F::Function, X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function, ξ::Int64, υ::Int64, ζ::Int64)::Float64

    weights = create_weights(X, ρ, η, ϕ, ψ, ξ, υ, ζ)
    points = split_region(X, ρ, η, ϕ, ψ, ξ, υ, ζ)
    return sum(divergence.(F, points) .* weights) * prod([coeff(interval, steps)
            for (interval, steps) in zip((X, (0, 1), (0, 1)), (ξ, υ, ζ))])
end


"""
    ∯(F::Function, r::Function, U::Tuple{Real, Real}, ϕ::Function, ψ::Function, μ::Int64, ν::Int64)::Float64

Determine flux of vector field `F` through a 2d region using Jacobian and Simpson's rule for ``μ⋅ν`` nodes.

# Examples
```
julia> ∯((x, y, z) -> [x, y, z], (u, v) -> [cos(u), sin(u), v], (0, 2π), u -> 0, u -> u, 2, 2)
19.739208804825942

julia> ∯((x, y, z) -> [x^2, y^2, z^2], (u, v) -> [u, v, 4-u-v], (0, 4), u -> 0, u -> 4-u, 8, 10)
64.00000000410644
```
"""
function ∯(F::Function, r::Function, U::Tuple{Real, Real}, ϕ::Function, ψ::Function, μ::Int64, ν::Int64)::Float64
    points = split_region(U, ϕ, ψ, μ, ν)
    weights = create_weights(U, ϕ, ψ, μ, ν)
    return sum(transform.(F, r, points) .* weights') * prod(coeff(interval, steps)
            for (interval, steps) in zip((U, (0, 1)), (μ, ν)))
end


"""
    draw_points(Xᵢ::LinRange{Float64}, Yᵢ::LinRange{Float64}, Zᵢ::LinRange{Float64}, N::Int64)::Array{Array{Float64, 1}, 1}

Choose randomly 3d points from a cuboid.

# Examples
```
julia> draw_points(LinRange(0, 1, 5), LinRange(0, 1, 5), LinRange(0, 1, 5), 11)
11-element Array{Array{Float64,1},1}:
 [0.75, 0.0, 0.5]
 [1.0, 0.25, 1.0]
 [1.0, 0.0, 1.0]
 [1.0, 0.5, 0.25]
 [0.5, 0.25, 1.0]
 [0.5, 0.0, 0.0]
 [0.5, 0.5, 0.75]
 [1.0, 0.75, 1.0]
 [0.5, 1.0, 0.0]
 [0.5, 1.0, 1.0]
 [0.75, 0.75, 1.0]

julia> draw_points(LinRange(-1, 1, 5), LinRange(0, π, 19), LinRange(2, 2, 2), 3)
3-element Array{Array{Float64,1},1}:
 [-0.5, 2.0943951023931953, 2.0]
 [0.0, 2.96705972839036, 2.0]
 [0.0, 2.443460952792061, 2.0]
```
"""
function draw_points(Xᵢ::LinRange{Float64}, Yᵢ::LinRange{Float64}, Zᵢ::LinRange{Float64}, N::Int64)::Array{Array{Float64, 1}, 1}
    return [[rand(Xᵢ), rand(Yᵢ), rand(Zᵢ)] for i in 1:N]
end


"""
    draw_points(Xᵢ::LinRange{Float64}, Yᵢ::LinRange{Float64}, N::Int64)::Array{Array{Float64, 1}, 1}

Choose randomly 2d points from a rectangle.

# Examples
```
julia> draw_points(LinRange(0, 1, 5), LinRange(0, 1, 5), 11)
11-element Array{Array{Float64,1},1}:
 [0.5, 1.0]
 [0.75, 0.75]
 [0.5, 0.75]
 [0.0, 0.25]
 [1.0, 1.0]
 [0.75, 0.25]
 [0.25, 0.25]
 [0.5, 0.75]
 [1.0, 0.5]
 [0.5, 1.0]
 [0.0, 0.25]

julia> draw_points(LinRange(-1, 1, 5), LinRange(2, 2, 19), 3)
3-element Array{Array{Float64,1},1}:
 [0.5, 2.0]
 [-1.0, 2.0]
 [-0.5, 2.0]
```
"""
function draw_points(Xᵢ::LinRange{Float64}, Yᵢ::LinRange{Float64}, N::Int64)::Array{Array{Float64, 1}, 1}
    return [[rand(Xᵢ), rand(Yᵢ)] for i in 1:N]
end


"""
    surface_part(F::Function, r::Function, S::Array{LinRange{Float64}, 1}, N::Int64, ϕ::Function, ψ::Function)::Float64

Compute flux over one of parts of a rectangle using Monte Carlo method.

# Examples
```
julia> surface_part((x, y, z) -> [x, y, z], (u, v) -> [u, v, 0], [LinRange(0, 1, 100), LinRange(0, 1, 100)], 1000, u -> 0, u -> u)
0.0

julia> surface_part((x, y, z) -> [x^2, y^2, z^2], (u, v) -> [u, v, u-v], [LinRange(0, 1, 100), LinRange(-1, 1, 100)], 100000, u -> -u, u -> u)
0.33893245180194237
```
"""
function surface_part(F::Function, r::Function, S::Array{LinRange{Float64}, 1}, N::Int64, ϕ::Function, ψ::Function)::Float64
    Uᵢ, Vᵢ = S
    points = draw_points(Uᵢ, Vᵢ, N)
    validate(P::Array{Float64, 1}) = ϕ(P[1]) <= P[2] <= ψ(P[1]) ? transform(F, r, P) : 0.0
    return sum(validate.(points)) * (max(Uᵢ...) - min(Uᵢ...)) * (max(Vᵢ...) - min(Vᵢ...)) / N
end


"""
    surface_part(F::Function, S::Array{LinRange{Float64}, 1}, N::Int64, ρ::Function, η::Function, ϕ::Function, ψ::Function)::Float64

Compute flux over one of parts of a cuboid using Monte Carlo method.

# Examples
```
julia> surface_part((x, y, z) -> [x, y, z], [LinRange(0, 1, 100), LinRange(0, 1, 100), LinRange(0, 1, 100)], 1000, (x, y) -> 0, (x, y) -> 1, x -> 0, x -> x)
1.4940000000104536

julia> surface_part((x, y, z) -> [x^2, y^2, z^2], [LinRange(0, 1, 100), LinRange(0, 1, 100), LinRange(-1, 1, 100)], 100000, (x, y) -> 0, (x, y) -> x+y, x -> -x, x -> x)
1.301682061621842
```
"""
function surface_part(F::Function, S::Array{LinRange{Float64}, 1}, N::Int64, ρ::Function, η::Function, ϕ::Function, ψ::Function)::Float64
    Xᵢ, Yᵢ, Zᵢ = S
    points = draw_points(Xᵢ, Yᵢ, Zᵢ, N)
    validate(P::Array{Float64, 1}) = ϕ(P[1]) <= P[2] <= ψ(P[1]) && ρ(P[1:2]...) <= P[3] <= η(P[1:2]...) ? divergence(F, P) : 0.0
    return sum(validate.(points)) * (max(Xᵢ...) - min(Xᵢ...)) * (max(Yᵢ...) - min(Yᵢ...)) * (max(Zᵢ...) - min(Zᵢ...)) / N
end


"""
    ∯(F::Function, r::Function, U::Tuple{Real, Real}, ϕ::Function, ψ::Function, N::Int64, parts::Float64)::Float64

Determine flux of vector field `F` through a 3d region using Gauss-Ostrogradski theorem and Monte Carlo method for `N` random points.

# Examples
```
julia> ∯((x, y, z) -> [x, y, z], (0, 1), (x, y) -> 0, (x, y) -> 1, x -> 0, x -> 1, 10000, 16.0)
2.999400020649781

julia> ∯((x, y, z) -> [x^2, y^2, z^2], (-1, 1), (x, y) -> 0, (x, y) -> x^2+y^2, x -> -1, x -> 1, 100, 256.0)
2.4227719867790514
```
"""
function ∯(F::Function, r::Function, U::Tuple{Real, Real}, ϕ::Function, ψ::Function, N::Int64, parts::Float64)::Float64
    x_part = Int64(parts / 2)
    y_part = Int64(parts / x_part)
    Uᵢ = LinRange(U..., x_part*N)
    Vᵢ = LinRange(min(ϕ.(Uᵢ)...), max(ψ.(Uᵢ)...), y_part*N)
    S = [[Uᵢ[N*i+1:(i+1)N], Vᵢ[N*j+1:(j+1)N]] for i in 0:(x_part-1), j in 0:(y_part-1)]
    return sum(surface_part.(F, r, S, N, ϕ, ψ))
end


"""
    ∯(F::Function, X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function, N::Int64, parts::Float64)::Float64

Determine flux of vector field `F` through a 2d region using Jacobian and Monte Carlo method for `N` random points.

# Examples
```
julia> ∯((x, y, z) -> [x, y, z], (u, v) -> [cos(u), sin(u), v], (0, 2π), u -> 0, u -> u, 10000, 16.0)
19.776707787586332

julia> ∯((x, y, z) -> [x^2, y^2, z^2], (u, v) -> [u, v, 4-u-v], (0, 4), u -> 0, u -> 4-u, 100000, 32.0)
63.9916073475161
```
"""
function ∯(F::Function, X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function, N::Int64, parts::Float64)::Float64
    x_part = y_part = Int64(parts / 4)
    z_part = Int64(parts / 2x_part)
    Xᵢ = LinRange(X..., x_part*N)
    Yᵢ = LinRange(min(ϕ.(Xᵢ)...), max(ψ.(Xᵢ)...), y_part*N)
    unzip(f, P) = try_eval(f, P...)
    Zᵢ = LinRange(min(unzip.(ρ, zip(Xᵢ, Yᵢ))...), max(unzip.(η, zip(Xᵢ, Yᵢ))...), z_part*N)
    S = [[Xᵢ[N*i+1:(i+1)N], Yᵢ[N*j+1:(j+1)N], Zᵢ[N*k+1:(k+1)N]] for i in 0:(x_part-1), j in 0:(y_part-1), k in 0:(z_part-1)]
    return sum(surface_part.(F, S, N, ρ, η, ϕ, ψ))
end


"""
    round_float(α::Float64, ϵ::Float64)::Union{Float64, Int64}

Round `α` to the same Real of digits as `ϵ` or nearest integer with `ϵ` accuracy.

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
round_float(α::Float64, ϵ::Float64)::Union{Int64, Float64} = abs(α - round(α)) < ϵ ? Int64(round(α)) : round(α, digits = abs(Decimal(ϵ).q))


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

Parse a Real or an expression from string to float.

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
    Φ(F::Function, X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function; N::Int64 = 100, technique::String = "Simpson")::Union{Int64, Float64}

Determine a flux of a vector field `F` through a surface.

# Examples
```
julia> Φ((x, y, z) -> [x, y, z], (0, 1), (x, y) -> 0, (x, y) -> 1, x -> 0, x -> 1)
3

julia> Φ((x, y, z) -> [-y, x, 0], (-1, 1), (x, y) -> 0, (x, y) -> sqrt(2-y^2-x^2), x -> -sqrt(1-x^2), x -> sqrt(1-x^2); technique = "Monte Carlo", ϵ = 0.1)
0
```
"""
function Φ(F::Function, X::Tuple{Real, Real}, ρ::Function, η::Function, ϕ::Function, ψ::Function;
     N::Int64 = 100, technique::String = "Simpson")::Union{Int64, Float64}
    if technique in ("Simpson", "Monte Carlo")
        surface_integral(n) = technique == "Simpson" ? ∯(F, X, ρ, η, ϕ, ψ, 2n, 2n, 2n) : ∯(F, X, ρ, η, ϕ, ψ, 100n, 4.0n)
    else
        throw(ArgumentError("Invalid technique name."))
    end
    return round_float(surface_integral(N), 1e-3)
end


"""
    Φ(F::Function, r::Function, U::Tuple{Real, Real}, ϕ::Function, ψ::Function; N::Int64 = 100, technique::String = "Simpson")::Union{Int64, Float64}

Determine a flux of a vector field `F` through a surface parametrized by `r`.

# Examples
```
julia> Φ((x, y, z) -> (x^2, y^2, z^2), (u, v) -> [u, v, 4-u-v], (0, 4), u -> 0, u -> 4-u)
64

julia> Φ((x, y, z) -> [x-y, y-z, z-x], (u, v) -> [(cos(u)+2)cos(v), (cos(u)+2)sin(v), sin(u)], (0, 2π), u -> 0, u -> 2π; technique = "Monte Carlo", ϵ = 0.01)
-118.25
```
"""
function Φ(F::Function, r::Function, U::Tuple{Real, Real}, ϕ::Function, ψ::Function;
     N::Int64 = 100, technique::String = "Simpson")::Union{Int64, Float64}
    if technique in ("Simpson", "Monte Carlo")
        surface_integral(n::Int64) = technique == "Simpson" ? ∯(F, r, U, ϕ, ψ, 2n, 2n) : ∯(F, r, U, ϕ, ψ, 100n, 4.0n)
    else
        throw(ArgumentError("Invalid technique name."))
    end
    return round_float(surface_integral(N), 1e-3)
end


end
