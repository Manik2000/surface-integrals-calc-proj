module IntegralUtils

using LinearAlgebra

export divergence, transform, create_weights, split_region, coeff

function ∂(f::Function, var::Symbol, P₀::Array{T1, 1}; Δ::Number = 1e-3)::Union{Number, Array{T2, 1} where T2 <: Number} where T1 <: Number
    if length(P₀) == 3
        limits = Dict(sym=>vec for (sym, vec) in
                zip((:x, :y, :z), [Δ*Diagonal(ones((3, 3)))[3i-2:3i] for i in 1:3]))
    elseif length(P₀) == 2
        limits = Dict(sym=>vec for (sym, vec) in
                zip((:u, :v), [Δ*Diagonal(ones((2, 2)))[2i-1:2i] for i in 1:2]))
    end
    return (f((P₀ + limits[var])...) - f(P₀...)) / Δ
end


function divergence(𝐅::Function, P₀::Array{T, 1})::Number where T <: Number
   return sum(hcat([∂(𝐅, sym, P₀) for sym in (:x, :y, :z)]...)' .* Diagonal(ones(3, 3)))
end


function 𝐍(𝐫::Function, P₀::Array{T, 1})::Array{Number, 1} where T <: Number
    𝐫ᵤ = ∂(𝐫, :u, P₀)
    𝐫ᵥ = ∂(𝐫, :v, P₀)
    return cross(𝐫ᵤ, 𝐫ᵥ)
end


function transform(𝐅::Function, 𝐫::Function, P₀::Array{T, 1})::Number where T <: Number
   return dot(𝐅(𝐫(P₀...)...), 𝐍(𝐫, P₀))
end


function create_weights(ξ::Int, υ::Int, ζ::Int)::Array{Float64, 3}
    return vcat([1], 4 * ones(ξ - 2) - repeat([0, 2], (ξ - 2) ÷ 2), [4, 1]) .*
            vcat([1], 4 * ones(υ - 2) - repeat([0, 2], (υ - 2) ÷ 2), [4, 1])' .*
            reshape(vcat([1], 4 * ones(ζ - 2) - repeat([0, 2], (ζ - 2) ÷ 2), [4, 1]), (1, 1, ζ + 1))
end


function create_weights(μ::Int, ν::Int)::Array{Float64, 2}
    return vcat([1], 4 * ones(μ - 2) - repeat([0, 2], (μ - 2) ÷ 2), [4, 1]) .*
            vcat([1], 4 * ones(ν - 2) - repeat([0, 2], (ν - 2) ÷ 2), [4, 1])'
end


step(k, N, 𝒜) = k * 𝒜[2] / N + 𝒜[1] * (1 - k / N)


function split_region(𝒳::Tuple{Number, Number}, 𝒴::Tuple{Number, Number}, 𝒵::Tuple{Number, Number},
        ξ::Int, υ::Int, ζ::Int)::Array{Array{Float64, 1}, 3}

    return [[step(x, ξ, 𝒳), step(y, υ, 𝒴), step(z, ζ, 𝒵)] for x in 0:ξ, y in 0:υ, z in 0:ζ]
end


function split_region(𝒰::Tuple{Number, Number}, 𝒱::Tuple{Number, Number},
        μ::Int, ν::Int)::Array{Array{Float64, 1}, 2}

    return [[step(u, μ, 𝒰), step(v, ν, 𝒱)] for u in 0:μ, v in 0:ν]
end


coeff(𝒜, n) = (𝒜[2] - 𝒜[1]) / (3 * n)

end
