module GraphingUtils
using PlotlyJS

export graph_all


function graph_all(x, y, z, Fx::Function, Fy::Function, Fz::Function,
                x_min::Number, x_max::Number,
                y_min::Number, y_max::Number,
                z_min::Number, z_max::Number,
                density::Int64)::Plot

    x_span = range(x_min, stop=x_max, length=density)
    y_span = range(y_min, stop=y_max, length=density)
    z_span = range(z_min, stop=z_max, length=density)
    𝓍 = repeat(x_span, density^3)
    𝓎 = repeat(y_span, inner=density, outer=density^2)
    𝓏 = repeat(z_span, inner=density^3)

    trace₁ = surface(;x=x, y=y, z=z)
    trace₂ = cone(;x=𝓍, y=𝓎, z=𝓏, u=Fx.(𝓍), v=Fy.(𝓎), w=Fz.(𝓏), showscale=false)
    layout = Layout(autosize=false, width=600, height=600)
    return Plot([trace₁, trace₂], layout)
end

# check(u, v, f) = v_min(u) <= v && v_max(u) >= v ? f(u, v) : NaN

end
