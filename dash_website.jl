using Dash
using DashCoreComponents
using DashHtmlComponents
using PlotlyJS
include("IntegralUtils.jl")
include("GraphingUtils.jl")
using .IntegralUtils
using .GraphingUtils


mathjax = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML"
app = dash(external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"],
            assets_folder="assets",
            external_scripts=[mathjax])

options_ = ["S: r(u, v)", "S: z = f(x, y), y = g(x)"]


app.layout = html_div() do
    html_div(children = [
    html_h1("\$\$\\text{Surface integrals calculator}\$\$", style=Dict("textAlign" => "center", "margin" => "3%")),
    html_div(id="header",
        children=[
        html_h3(id="header_h3", "\$\$ \\text{Hello, let the fun begin!}\$\$")
                ]
            ),
    html_hr(style=Dict("borderTop" => "1px dashed #d69c2f")),
    html_div(id="surface", children=[
        html_label(id="surface_label", "\$\$\\text{Choose the way of defining the surface:}\$\$"),
        dcc_dropdown(
            id = "dropdown",
            options = [(label = i, value = i) for i in options_],
            value = options_[1],
            ),
    ]),
    html_hr(style=Dict("borderTop" => "1px dashed #d69c2f")),
    html_div(id="main", children=[
        html_div(id="inputs", children=[
            html_label(id="vector_field_l", "\$\$\\text{Enter your vector field formula}\$\$"),
            html_div(id="field", children=[
            html_label(id="mFx", "\$\$\\vec{F}\\cdot\\hat{i}\$\$"),
            dcc_input(id="Fx", type="text", value="x"),
            html_label(id="mFy", "\$\$\\vec{F}\\cdot\\hat{j}\$\$"),
            dcc_input(id="Fy", type="text", value="y"),
            html_label(id="mFz", "\$\$\\vec{F}\\cdot\\hat{k}\$\$"),
            dcc_input(id="Fz", type="text", value="z"),
            ]
        ),
        html_label(id="method_dropdown", "\$\$\\text{Choose method of numeric integration:}\$\$"),
        dcc_dropdown(
            id = "method",
            options = [(label = i, value = i) for i in ("Simpson", "Monte Carlo")],
            value = "Simpson"
            ),
        html_label(id="accurary", "\$\$\\text{Choose the accuracy of numeric integration:}\$\$"),
        dcc_slider(
            id = "integral_accuracy",
            min = 10,
            max = 100,
            marks = Dict([Symbol(v) => Symbol(v) for v in 10:10:100]),
            value = 10,
            step = 10,
            ),
            html_div(id="parametric_bounds", children=[
            html_label(id="r1lab", "\$\$\\vec{r}(u, v)\\cdot\\hat{i} =\$\$"),
            dcc_input(id="r1", type="text", value="sqrt(1/4 + u^2) * cos(v)"),
            html_label(id="r2lab", "\$\$\\vec{r}(u, v)\\cdot\\hat{j} =\$\$"),
            dcc_input(id="r2", type="text", value="sqrt(1/4 + u^2) * sin(v)"),
            html_label(id="r3lab", "\$\$\\vec{r}(u, v)\\cdot\\hat{k} =\$\$"),
            dcc_input(id="r3", type="text", value="u"),
            html_label(id="u_range1", "\$\$u_{min} =\$\$"),
            dcc_input(id="u_range1i", type="text", value="-1"),
            html_label(id="u_range2", "\$\$u_{max} =\$\$"),
            dcc_input(id="u_range2i", type="text", value="1"),
            html_label(id="v_range1", "\$\$\\phi(u) =\$\$"),
            dcc_input(id="v_range1i", type="text", value="0"),
            html_label(id="v_range2", "\$\$\\theta(u) =\$\$"),
            dcc_input(id="v_range2i", type="text", value="2pi")
            ], style=Dict("display" => "grid")
            ),
        html_div(id="fxy_bounds", children=[
            html_label(id="x_range1", "\$\$x_{\\text{min}} =\$\$"),
            dcc_input(id="x_range1i", type="text", value="-2"),
            html_label(id="x_range2", "\$\$ x_{\\text{max}} =\$\$"),
            dcc_input(id="x_range2i", type="text", value="3"),
            html_label(id="y_range1", "\$\$ \\zeta(x) =\$\$"),
            dcc_input(id="y_range1i", type="text", value="-2"),
            html_label(id="y_range2", "\$\$ \\eta(x) =\$\$"),
            dcc_input(id="y_range2i", type="text", value="2"),
            html_label(id="z_range1", "\$\$ f(x, y) =\$\$"),
            dcc_input(id="z_range1i", type="text", value="-(x^2 + y^2)"),
            html_label(id="z_range2", "\$\$ g(x, y) =\$\$"),
            dcc_input(id="z_range2i", type="text", value="x^2+y^2"),
        ], style=Dict("display" => "none")
            ), html_button(id = "submit-button-state", children = "submit", n_clicks = 0),
        html_div(id="output_")
        ]),
    html_div(id="graphing", children=[
            dcc_graph(
                id = "surface_plot",
                figure = (
                    data = [
                        (x = [0],
                         y = [0],
                         z = [[0], [0]],
                        type = "surface")
                    ],
                    layout = (
                                autosize=false,
                                width=600,
                                height=600
                            )
                )
            ),
        ]
    ),
]),
html_div(id="plot_attributes",
    children=[
        html_label(id="field_density_lab", "\$\$\\text{Choose the density of the vector field.}\$\$"),
        dcc_slider(
            id = "field_density",
            min = 0,
            max = ,10
            marks = Dict([Symbol(v) => Symbol(v) for v in 0:1:10]),
            value = 5,
            step = 1,
            ),
            html_label(id="graph_accuracy_lab", "\$\$\\text{Choose the accuracy of plotting.}\$\$"),
            dcc_slider(
            id = "graph_accuracy",
            min = 10,
            max = 100,
            marks = Dict([Symbol(v) => Symbol(v) for v in 10:10:100]),
            value = 20,
            step = 10,
            )
    ]
), html_footer("\$\$\\text{MIT License}. \\ \\text{MM, ML, MK.}\$\$", style=Dict("marginTop" => "3em", "textAlign" => "center"))
])
end

callback!(app,
        Output("parametric_bounds", "style"),
        Output("fxy_bounds", "style"),
        Input("dropdown", "value")
        ) do user_choice
        if user_choice == options_[1]
            return Dict("display" => "grid"), Dict("display" => "none")
        else
            return Dict("display" => "none"), Dict("display" => "grid")
        end
end


callback!(app,
        Output("output_", "children"),
        Input("submit-button-state", "n_clicks"),
        State("dropdown", "value"),
        State("r1", "value"),
        State("r2", "value"),
        State("r3", "value"),
        State("Fx", "value"),
        State("Fy", "value"),
        State("Fz", "value"),
        State("method", "value"),
        State("integral_accuracy", "value"),
        State("u_range1i", "value"),
        State("u_range2i", "value"),
        State("v_range1i", "value"),
        State("v_range2i", "value"),
        State("x_range1i", "value"),
        State("x_range2i", "value"),
        State("y_range1i", "value"),
        State("y_range2i", "value"),
        State("z_range1i", "value"),
        State("z_range2i", "value")
        ) do n_clicks, dropdown_value, r1, r2, r3, Fx, Fy, Fz, technique, n, u_min, u_max,
            v_min, v_max, x_min, x_max, y_min, y_max, z_min, z_max
        if dropdown_value == options_[1]
            try
                u_min, u_max = parse_num.([u_min, u_max])
                p = parse_function(join(["[" * r1, r2, r3 *"]" ], ", "), :u, :v)
                q = parse_function(join(["[" * Fx, Fy, Fz * "]"], ", "), :x, :y, :z)
                ϕ = parse_function(v_min, :u)
                ψ = parse_function(v_max, :u)
                value(f, g, a, b) = Φ(f, g, (u_min, u_max), a, b; N = Int(n ÷ 2), technique = technique)
                value_ = @eval abs(($value($q, $p, $ϕ, $ψ)))
                return html_h5("|Φ| = $(value_).", style=Dict("textAlign" => "center"))
            catch e
                return html_h5("\$\$\\text{Wrong input. Try again.}\$\$", style=Dict("textAlign" => "center"))
            end
        else
            try
                x_min, x_max = parse_num.([x_min, x_max])
                q = parse_function(join(["[" * Fx, Fy, Fz * "]"], ", "), :x, :y, :z)
                ϕ = parse_function(y_min, :x)
                ψ = parse_function(y_max, :x)
                ρ = parse_function(z_min, :x, :y)
                η = parse_function(z_max, :x, :y)
                value(f, a, b, c, d) = Φ(f, (x_min, x_max), a, b, c, d; N = Int(n ÷ 2), technique = technique)
                value_ = @eval abs(($value($q, $ρ, $η, $ϕ, $ψ)))
                return html_h5("|Φ| = $(value_).", style=Dict("textAlign" => "center"))
            catch e
                return html_h5("\$\$ \\text{Wrong input. Try again.} \$\$", style=Dict("textAlign" => "center"))
            end
        end
end


callback!(app,
        Output("surface_plot", "figure"),
        Input("submit-button-state", "n_clicks"),
        State("r1", "value"),
        State("r2", "value"),
        State("r3", "value"),
        State("dropdown", "value"),
        State("graph_accuracy", "value"),
        State("u_range1i", "value"),
        State("u_range2i", "value"),
        State("v_range1i", "value"),
        State("v_range2i", "value"),
        State("x_range1i", "value"),
        State("x_range2i", "value"),
        State("y_range1i", "value"),
        State("y_range2i", "value"),
        State("z_range1i", "value"),
        State("z_range2i", "value"),
        State("Fx", "value"),
        State("Fy", "value"),
        State("Fz","value"),
        State("field_density", "value")
        ) do n_clicks, r1, r2, r3, dropdown_value, N, u_min, u_max, v_min, v_max,
            x_min, x_max, y_min, y_max, z_min, z_max, F1, F2, F3, density

            if dropdown_value == options_[1] # _______________parametric

                    try
                        u_min = parse_num(u_min)
                        u_max = parse_num(u_max)
                        v_min = parse_function(v_min, :u)
                        v_max = parse_function(v_max, :u)
                    catch e
                        u_min = -1
                        u_max = 1
                        v_min = parse_function("0", :u)
                        v_max = parse_function("2pi", :u)
                    end

                    if false
                        return Plot(surface(;
                                x = [0],
                                y = [0],
                                z = [[0], [0]]))
                    else
                        Fx = Fy = Fz = "_"
                        X(u, v) = 0
                        Y(u, v) = 0
                        Z(u, v) = 0

                        try
                            X = parse_function(r1, :u, :v)
                            Y = parse_function(r2, :u, :v)
                            Z = parse_function(r3, :u, :v)
                            Fx = parse_function(F1, :x, :y, :z)
                            Fy = parse_function(F2, :x, :y, :z)
                            Fz = parse_function(F3, :x, :y, :z)
                            @eval ($X(-2, 3), $Y(-2, 3), $Z(-2, 3))
                            @eval (isa($X(-2, 3), Array{Float64, 1}))
                            @eval (isa($Y(-2, 3), Array{Float64, 1}))
                            @eval (isa($Z(-2, 3), Array{Float64, 1}))
                            @eval (isa($Fx(-2, -2, -2), Number))
                            @eval (isa($Fy(-2, -2, -2), Number))
                            @eval (isa($Fz(-2, -2, -2), Number))
                        catch e
                            X = parse_function("0", :u, :v)
                            Y = parse_function("0", :u, :v)
                            Z = parse_function("0", :u, :v)
                            Fx = parse_function("0", :x, :y, :z)
                            Fy = parse_function("0", :x, :y, :z)
                            Fz = parse_function("0", :x, :y, :z)
                        end


                        us = LinRange(u_min, u_max, N)
                        value_(g, u::LinRange) = g.(u)
                        value_(g, u::Number) = g(u)
                        vs = @eval (LinRange(minimum($value_($v_min, $us)), maximum($value_($v_max, $us)), $N))

                        check(u, v, f) = @eval($value_($v_min, $u) <= $v && $value_($v_max, $u) >= $v ? $f($u, $v) : NaN)

                        value(g) = check.(us', vs, g)
                        x1 = @eval ($value($X))
                        y1 = @eval ($value($Y))
                        z1 = @eval ($value($Z))

                        return @eval($graph_all($x1, $y1, $z1, $Fx, $Fy, $Fz,
                                minimum(filter(!isnan, $x1)), maximum(filter(!isnan, $x1)),
                                minimum(filter(!isnan, $y1)), maximum(filter(!isnan, $y1)),
                                minimum(filter(!isnan, $z1)), maximum(filter(!isnan, $z1)), $density))
                    end
            else  # ______________________________________________f(x,y)

                Fx = Fy = Fz = "_"

                try
                    x_min = parse_num(x_min)
                    x_max = parse_num(x_max)
                    y_min = parse_function(y_min, :x)
                    y_max = parse_function(y_max, :x)
                    z_min = parse_function(z_min, :x, :y)
                    z_max = parse_function(z_max, :x, :y)
                    Fx = parse_function(F1, :x, :y, :z)
                    Fy = parse_function(F2, :x, :y, :z)
                    Fz = parse_function(F3, :x, :y, :z)
                    @eval (isa($Fx(-2, -2, -2), Number))
                    @eval (isa($Fy(-2, -2, -2), Number))
                    @eval (isa($Fz(-2, -2, -2), Number))
                catch e
                    x_min = -2
                    y_min = parse_function("-2", :x)
                    z_min = parse_function("-2", :x, :y)
                    x_max = 2
                    y_max = parse_function("2", :x)
                    z_max = parse_function("2", :x, :y)
                    Fx = parse_function("0", :x, :y, :z)
                    Fy = parse_function("0", :x, :y, :z)
                    Fz = parse_function("0", :x, :y, :z)
                end

                if false  #(x_min > x_max) | (y_min > y_max) | (z_min > z_max)
                    return Plot(surface(;
                            x = [0],
                            y = [0],
                            z = [[0], [0]]))
                else
                    xs = LinRange(x_min, x_max, N)
                    value2_(g::Function, x::LinRange{Float64}) = g.(x)
                    value2_(g::Function, x::Number) = g(x)
                    ys = @eval (LinRange(minimum($value2_($y_min, $xs)), maximum($value2_($y_max, $xs)), $N))

                    check2(x::Number, y::Number, h::Function)::Float64 = @eval($value2_($y_min, $x) <= $y && $value2_($y_max, $x) >= $y ? $h($x, $y) : NaN)

                    value2(g) = check2.(xs', ys, g)

                    zs = @eval($value2($z_max))
                    z₀ = @eval($value2($z_min))
                    temp = deepcopy(zs)
                    zs[zs .< z₀] .= NaN
                    z₀[z₀ .> temp] .= NaN

                    return @eval($graph_all($xs, $ys, $zs, $Fx, $Fy, $Fz,
                            minimum($xs), maximum($xs),
                            minimum(filter(!isnan, $ys)), maximum(filter(!isnan, $ys)),
                            minimum(filter(!isnan, $z₀)), maximum(filter(!isnan, $zs)), $density, $z₀))
                end
            end
        end


run_server(app, "0.0.0.0")
