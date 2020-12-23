#module IntegralUtils

using LinearAlgebra
using Decimals

#export divergence, transform, create_weights,
# split_region, coeff, surface_integral, round_float, parse_function

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
    create_weights(ξ::Int, υ::Int, ζ::Int)::Array{Float64, 3}

Create a 3d matrix of weights according to composite simpson rule for triple integrals.

# Examples
```
julia> create_weights(2, 2, 2)
3×3×3 Array{Float64,3}:
[:, :, 1] =
 1.0   4.0  1.0
 4.0  16.0  4.0
 1.0   4.0  1.0

[:, :, 2] =
  4.0  16.0   4.0
 16.0  64.0  16.0
  4.0  16.0   4.0

[:, :, 3] =
 1.0   4.0  1.0
 4.0  16.0  4.0
 1.0   4.0  1.0

 julia> create_weights(6, 4, 2)
 7×5×3 Array{Float64,3}:
 [:, :, 1] =
  1.0   4.0  2.0   4.0  1.0
  4.0  16.0  8.0  16.0  4.0
  2.0   8.0  4.0   8.0  2.0
  4.0  16.0  8.0  16.0  4.0
  2.0   8.0  4.0   8.0  2.0
  4.0  16.0  8.0  16.0  4.0
  1.0   4.0  2.0   4.0  1.0

 [:, :, 2] =
   4.0  16.0   8.0  16.0   4.0
  16.0  64.0  32.0  64.0  16.0
   8.0  32.0  16.0  32.0   8.0
  16.0  64.0  32.0  64.0  16.0
   8.0  32.0  16.0  32.0   8.0
  16.0  64.0  32.0  64.0  16.0
   4.0  16.0   8.0  16.0   4.0

 [:, :, 3] =
  1.0   4.0  2.0   4.0  1.0
  4.0  16.0  8.0  16.0  4.0
  2.0   8.0  4.0   8.0  2.0
  4.0  16.0  8.0  16.0  4.0
  2.0   8.0  4.0   8.0  2.0
  4.0  16.0  8.0  16.0  4.0
  1.0   4.0  2.0   4.0  1.0
```
"""
function create_weights(ξ::Int, υ::Int, ζ::Int)::Array{Float64, 3}
    return vcat([1], 4 * ones(ξ - 2) - repeat([0, 2], (ξ - 2) ÷ 2), [4, 1]) .*
            vcat([1], 4 * ones(υ - 2) - repeat([0, 2], (υ - 2) ÷ 2), [4, 1])' .*
            reshape(vcat([1], 4 * ones(ζ - 2) - repeat([0, 2], (ζ - 2) ÷ 2), [4, 1]), (1, 1, ζ + 1))
end


"""
    create_weights(μ::Int, ν::Int)::Array{Float64, 2}

Create a 2d matrix of weights according to composite simpson rule for double integrals.

# Examples
```
julia> create_weights(2, 2)
3×3 Array{Float64,2}:
 1.0   4.0  1.0
 4.0  16.0  4.0
 1.0   4.0  1.0

julia> create_weights(6, 4)
7×5 Array{Float64,2}:
 1.0   4.0  2.0   4.0  1.0
 4.0  16.0  8.0  16.0  4.0
 2.0   8.0  4.0   8.0  2.0
 4.0  16.0  8.0  16.0  4.0
 2.0   8.0  4.0   8.0  2.0
 4.0  16.0  8.0  16.0  4.0
 1.0   4.0  2.0   4.0  1.0
```
"""
function create_weights(μ::Int, ν::Int)::Array{Float64, 2}
    return vcat([1], 4 * ones(μ - 2) - repeat([0, 2], (μ - 2) ÷ 2), [4, 1]) .*
            vcat([1], 4 * ones(ν - 2) - repeat([0, 2], (ν - 2) ÷ 2), [4, 1])'
end


"""
    step(k::Int, N::Int, A::Tuple{Number, Number})::Float64

Calculate a `k`th split point of interval `A`.

# Examples
```
julia> step(2, 4, (0, 4))
2.0

julia> step(3, 6, (1, 4))
2.5
```
"""
step(k::Int, N::Int, A::Tuple{Number, Number})::Float64 = k * A[2] / N + A[1] * (1 - k / N)


"""
    split_region(X::Tuple{Number, Number}, Y::Tuple{Number, Number}, Z::Tuple{Number, Number}, ξ::Int, υ::Int, ζ::Int)::Array{Array{Float64, 1}, 3}

Divide a cuboid into points in 3d space.

# Examples
```
julia> split_region((0, 1), (0, 1), (0, 1), 2, 2, 2)
3×3×3 Array{Array{Float64,1},3}:
[:, :, 1] =
 [0.0, 0.0, 0.0]  [0.0, 0.5, 0.0]  [0.0, 1.0, 0.0]
 [0.5, 0.0, 0.0]  [0.5, 0.5, 0.0]  [0.5, 1.0, 0.0]
 [1.0, 0.0, 0.0]  [1.0, 0.5, 0.0]  [1.0, 1.0, 0.0]

[:, :, 2] =
 [0.0, 0.0, 0.5]  [0.0, 0.5, 0.5]  [0.0, 1.0, 0.5]
 [0.5, 0.0, 0.5]  [0.5, 0.5, 0.5]  [0.5, 1.0, 0.5]
 [1.0, 0.0, 0.5]  [1.0, 0.5, 0.5]  [1.0, 1.0, 0.5]

[:, :, 3] =
 [0.0, 0.0, 1.0]  [0.0, 0.5, 1.0]  [0.0, 1.0, 1.0]
 [0.5, 0.0, 1.0]  [0.5, 0.5, 1.0]  [0.5, 1.0, 1.0]
 [1.0, 0.0, 1.0]  [1.0, 0.5, 1.0]  [1.0, 1.0, 1.0]

 julia> split_region((0, 1), (0, 1), (0, 1), 8, 4, 2)
 9×5×3 Array{Array{Float64,1},3}:
 [:, :, 1] =
  [0.0, 0.0, 0.0]    [0.0, 0.25, 0.0]    [0.0, 0.5, 0.0]    [0.0, 0.75, 0.0]    [0.0, 1.0, 0.0]
  [0.125, 0.0, 0.0]  [0.125, 0.25, 0.0]  [0.125, 0.5, 0.0]  [0.125, 0.75, 0.0]  [0.125, 1.0, 0.0]
  [0.25, 0.0, 0.0]   [0.25, 0.25, 0.0]   [0.25, 0.5, 0.0]   [0.25, 0.75, 0.0]   [0.25, 1.0, 0.0]
  [0.375, 0.0, 0.0]  [0.375, 0.25, 0.0]  [0.375, 0.5, 0.0]  [0.375, 0.75, 0.0]  [0.375, 1.0, 0.0]
  [0.5, 0.0, 0.0]    [0.5, 0.25, 0.0]    [0.5, 0.5, 0.0]    [0.5, 0.75, 0.0]    [0.5, 1.0, 0.0]
  [0.625, 0.0, 0.0]  [0.625, 0.25, 0.0]  [0.625, 0.5, 0.0]  [0.625, 0.75, 0.0]  [0.625, 1.0, 0.0]
  [0.75, 0.0, 0.0]   [0.75, 0.25, 0.0]   [0.75, 0.5, 0.0]   [0.75, 0.75, 0.0]   [0.75, 1.0, 0.0]
  [0.875, 0.0, 0.0]  [0.875, 0.25, 0.0]  [0.875, 0.5, 0.0]  [0.875, 0.75, 0.0]  [0.875, 1.0, 0.0]
  [1.0, 0.0, 0.0]    [1.0, 0.25, 0.0]    [1.0, 0.5, 0.0]    [1.0, 0.75, 0.0]    [1.0, 1.0, 0.0]

 [:, :, 2] =
  [0.0, 0.0, 0.5]    [0.0, 0.25, 0.5]    [0.0, 0.5, 0.5]    [0.0, 0.75, 0.5]    [0.0, 1.0, 0.5]
  [0.125, 0.0, 0.5]  [0.125, 0.25, 0.5]  [0.125, 0.5, 0.5]  [0.125, 0.75, 0.5]  [0.125, 1.0, 0.5]
  [0.25, 0.0, 0.5]   [0.25, 0.25, 0.5]   [0.25, 0.5, 0.5]   [0.25, 0.75, 0.5]   [0.25, 1.0, 0.5]
  [0.375, 0.0, 0.5]  [0.375, 0.25, 0.5]  [0.375, 0.5, 0.5]  [0.375, 0.75, 0.5]  [0.375, 1.0, 0.5]
  [0.5, 0.0, 0.5]    [0.5, 0.25, 0.5]    [0.5, 0.5, 0.5]    [0.5, 0.75, 0.5]    [0.5, 1.0, 0.5]
  [0.625, 0.0, 0.5]  [0.625, 0.25, 0.5]  [0.625, 0.5, 0.5]  [0.625, 0.75, 0.5]  [0.625, 1.0, 0.5]
  [0.75, 0.0, 0.5]   [0.75, 0.25, 0.5]   [0.75, 0.5, 0.5]   [0.75, 0.75, 0.5]   [0.75, 1.0, 0.5]
  [0.875, 0.0, 0.5]  [0.875, 0.25, 0.5]  [0.875, 0.5, 0.5]  [0.875, 0.75, 0.5]  [0.875, 1.0, 0.5]
  [1.0, 0.0, 0.5]    [1.0, 0.25, 0.5]    [1.0, 0.5, 0.5]    [1.0, 0.75, 0.5]    [1.0, 1.0, 0.5]

 [:, :, 3] =
  [0.0, 0.0, 1.0]    [0.0, 0.25, 1.0]    [0.0, 0.5, 1.0]    [0.0, 0.75, 1.0]    [0.0, 1.0, 1.0]
  [0.125, 0.0, 1.0]  [0.125, 0.25, 1.0]  [0.125, 0.5, 1.0]  [0.125, 0.75, 1.0]  [0.125, 1.0, 1.0]
  [0.25, 0.0, 1.0]   [0.25, 0.25, 1.0]   [0.25, 0.5, 1.0]   [0.25, 0.75, 1.0]   [0.25, 1.0, 1.0]
  [0.375, 0.0, 1.0]  [0.375, 0.25, 1.0]  [0.375, 0.5, 1.0]  [0.375, 0.75, 1.0]  [0.375, 1.0, 1.0]
  [0.5, 0.0, 1.0]    [0.5, 0.25, 1.0]    [0.5, 0.5, 1.0]    [0.5, 0.75, 1.0]    [0.5, 1.0, 1.0]
  [0.625, 0.0, 1.0]  [0.625, 0.25, 1.0]  [0.625, 0.5, 1.0]  [0.625, 0.75, 1.0]  [0.625, 1.0, 1.0]
  [0.75, 0.0, 1.0]   [0.75, 0.25, 1.0]   [0.75, 0.5, 1.0]   [0.75, 0.75, 1.0]   [0.75, 1.0, 1.0]
  [0.875, 0.0, 1.0]  [0.875, 0.25, 1.0]  [0.875, 0.5, 1.0]  [0.875, 0.75, 1.0]  [0.875, 1.0, 1.0]
  [1.0, 0.0, 1.0]    [1.0, 0.25, 1.0]    [1.0, 0.5, 1.0]    [1.0, 0.75, 1.0]    [1.0, 1.0, 1.0]
```
"""
function split_region(X::Tuple{Number, Number}, Y::Tuple{Number, Number}, Z::Tuple{Number, Number},
        ξ::Int, υ::Int, ζ::Int)::Array{Array{Float64, 1}, 3}

    return [[step(x, ξ, X), step(y, υ, Y), step(z, ζ, Z)] for x in 0:ξ, y in 0:υ, z in 0:ζ]
end


"""
    split_region(U::Tuple{Number, Number}, V::Tuple{Number, Number}, μ::Int, ν::Int)::Array{Array{Float64, 1}, 2}

Divide a rectangle into points in 2d space.

# Examples
```
julia> split_region((0, 1), (0, 1), 2, 2)
3×3 Array{Array{Float64,1},2}:
 [0.0, 0.0]  [0.0, 0.5]  [0.0, 1.0]
 [0.5, 0.0]  [0.5, 0.5]  [0.5, 1.0]
 [1.0, 0.0]  [1.0, 0.5]  [1.0, 1.0]

julia> split_region((0, 1), (0, 1), 4, 2)
5×3 Array{Array{Float64,1},2}:
 [0.0, 0.0]   [0.0, 0.5]   [0.0, 1.0]
 [0.25, 0.0]  [0.25, 0.5]  [0.25, 1.0]
 [0.5, 0.0]   [0.5, 0.5]   [0.5, 1.0]
 [0.75, 0.0]  [0.75, 0.5]  [0.75, 1.0]
 [1.0, 0.0]   [1.0, 0.5]   [1.0, 1.0]
```
"""
function split_region(U::Tuple{Number, Number}, V::Tuple{Number, Number},
        μ::Int, ν::Int)::Array{Array{Float64, 1}, 2}

    return [[step(u, μ, U), step(v, ν, V)] for u in 0:μ, v in 0:ν]
end


validate(A::Array{Float64, 1}, ϕ::Function, ψ::Function)::Float64 = ϕ(A[1]) <= A[2] <= ψ(A[1]) ? 1.0 : NaN


function split_region(U::Tuple{Number, Number}, ϕ::Function, ψ::Function, μ::Int, ν::Int)::Array{Array{Float64, 1}, 2}
    Uᵢ = LinRange(U..., μ + 1)
    Vᵢ = LinRange(min(ϕ.(Uᵢ)...), max(ψ.(Uᵢ)...), ν + 1)
    R = [[u, v] for u in Uᵢ, v in Vᵢ]
    return validate.(R, ϕ, ψ) .* R
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
coeff(A::Tuple{Number, Number}, n::Int)::Float64 = (A[2] - A[1]) / (3 * n)


"""
    ∯(F::Function, X::Tuple{Number, Number}, Y::Tuple{Number, Number}, Z::Tuple{Number, Number}; ξ::Int = 2, υ::Int = 2, ζ::Int = 2)::Float64

Determine flux of vector field `F` through a cuboid using Gauss-Ostrogradski theorem and Simpson's rule for ``ξ⋅υ⋅ζ`` nodes.

# Examples
```
julia> surface_integral((x, y, z) -> [x, y, z], (0, 1), (0, 1), (0, 1))
3.0000000000163776

julia> surface_integral((x, y, z) -> [-y, x, 0], (-1, 1), (-1, 1), (-1, 1))
0.0
```
"""
function ∯(F::Function, X::Tuple{Number, Number}, Y::Tuple{Number, Number}, Z::Tuple{Number, Number};
    ξ::Int = 2, υ::Int = 2, ζ::Int = 2)::Float64

    weights = create_weights(ξ, υ, ζ)
    points = split_region(X, Y, Z, ξ, υ, ζ)
    return sum(divergence.(F, points) .* weights) * prod([coeff(interval, steps)
            for (interval, steps) in zip((X, Y, Z), (ξ, υ, ζ))])
end


"""
    ∯(F::Function, r::Function, U::Tuple{Number, Number}, V::Tuple{Number, Number}; μ::Int = 2, ν::Int = 2)::Float64

Determine flux of vector field `F` through a surface parametrized by `r` using Simpson's rule for ``μ⋅ν`` nodes.

# Examples
```
julia> surface_integral((x, y, z) -> [-y, x, 0], (u, v) -> [(cos(u) + 2)cos(v), (cos(u) + 2)sin(v), sin(u)],  (0, 2π), (0, π/2))
-1.1513612148289776e-5

julia> surface_integral((x, y, z) -> [x, y, z], (u, v) -> [cos(u)cos(v), sin(u)cos(v), sin(v)], (0, 2π), (-π/2, π/2), 2, 8)
12.56806186014683
```
"""
function ∯(F::Function, r::Function, U::Tuple{Number, Number}, V::Tuple{Number, Number};
        μ::Int = 2, ν::Int = 2)::Float64

    weights = create_weights(μ, ν)
    points = split_region(U, V, μ, ν)
    return sum(transform.(F, r, points) .* weights) * prod([coeff(interval, steps)
            for (interval, steps) in zip((U, V), (μ, ν))])
end


function ∯(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function;
        μ::Int = 2, ν::Int = 2)::Float64

    weights = create_weights(μ, ν)
    points = split_region(U, ϕ, ψ, μ, ν)
    Uᵢ = LinRange(U..., μ)
    return sum(replace!(transform.(F, r, points), NaN => 0) .* weights) * prod([coeff(interval, steps)
            for (interval, steps) in zip((U, (min(ϕ.(Uᵢ)...), max(ψ.(Uᵢ)...))), (μ, ν))])
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
round_float(α::Float64, ϵ::Float64)::Union{Int, Float64} = abs(α - round(α)) < ϵ ? Int(round(α)) : 2 + round(α, digits = abs(Decimal(ϵ).q))


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

Parse a body of a function and arguments' and function name's symbols to a actual function.

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
    Φ(F::Function, X::Tuple{Number, Number}, Y::Tuple{Number, Number}, Z::Tuple{Number, Number}; ϵ::Number = 1e-3, n::Int = 1)::Union{Int, Float64}

Determine a flux of a vector field `F` through a cuboid.

# Examples
```
julia> Φ((x, y, z) -> [x, y, z], (0, 1), (0, 1), (0, 1))
2.9999999999958114

julia> Φ((x, y, z) -> [-y, x, 0], (-1, 1), (-1, 1), (-1, 1))
0.0
```
"""
function Φ(F::Function, X::Tuple{Number, Number}, Y::Tuple{Number, Number}, Z::Tuple{Number, Number};
     ϵ::Number = 1e-3, n::Int = 1)::Union{Int, Float64}
    Φₙ = ∯(F, X, Y, Z; ξ = 2n, υ = 2n, ζ = 2n)
    n += 1
    Φₙ₊₁ = ∯(F, X, Y, Z; ξ = 2n, υ = 2n, ζ = 2n)
    while abs(Φₙ₊₁ - Φₙ) > ϵ
        n += 1
        Φₙ, Φₙ₊₁ = Φₙ₊₁, ∯(F, X, Y, Z; ξ = 2n, υ = 2n, ζ = 2n)
    end
    return round_float(Φₙ₊₁, ϵ)
end



"""
    Φ(F::Function, r::Function, U::Tuple{Number, Number}, V::Tuple{Number, Number}; ϵ::Number = 1e-3, n::Int = 1)::Union{Int, Float64}

Determine a flux of a vector field `F` through a surface parametrized by `r`.

# Examples
```
julia> Φ((x, y, z) -> [-y, x, 0], (u, v) -> [(cos(u) + 2)cos(v), (cos(u) + 2)sin(v), sin(u)],  (0, 2π), (0, π/2))
-9.869203248214035e-6

julia> Φ((x, y, z) -> [x, y, z], (u, v) -> [cos(u)cos(v), sin(u)cos(v), sin(v)], (0, 2π), (-π/2, π/2))
12.566373879290955
```
"""
function Φ(F::Function, r::Function, U::Tuple{Number, Number}, V::Tuple{Number, Number};
     ϵ::Number = 1e-3, n::Int = 1)::Union{Int, Float64}
    Φₙ = ∯(F, r, U, V; μ = 2n, ν = 2n)
    n += 1
    Φₙ₊₁ = ∯(F, r, U, V; μ = 2n, ν = 2n)
    while abs(Φₙ₊₁ - Φₙ) > ϵ
        n += 1
        Φₙ, Φₙ₊₁ = Φₙ₊₁, ∯(F, r, U, V; μ = 2n, ν = 2n)
    end
    return round_float(Φₙ₊₁, ϵ)
end


function Φ(F::Function, r::Function, U::Tuple{Number, Number}, ϕ::Function, ψ::Function;
     ϵ::Number = 1e-3, n::Int = 1)::Union{Int, Float64}
    Φₙ = ∯(F, r, U, ϕ, ψ; μ = 2n, ν = 2n)
    n += 1
    Φₙ₊₁ = ∯(F, r, U, ϕ, ψ; μ = 2n, ν = 2n)
    while abs(Φₙ₊₁ - Φₙ) > ϵ
        n += 1
        println(Φₙ)
        Φₙ, Φₙ₊₁ = Φₙ₊₁, ∯(F, r, U, ϕ, ψ; μ = 2n, ν = 2n)
    end
    return round_float(Φₙ₊₁, ϵ)
end

#end
