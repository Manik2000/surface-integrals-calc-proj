using Dash
using DashCoreComponents
using DashHtmlComponents
using PlotlyJS
using .IntegralUtils


app = dash(external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"])
markdown_text = """
                We are **extremaly** happy that you've visited our website.
                Below we've prepared a calulator of surface integrals.
                1. Firstly choose a way for describing the surface,
                2. Then insert the data for bla bla bla
                """
options_ = ["parametric representation", "z = f(x, y)"]


app.layout = html_div() do
    html_div(style=Dict("textAlign" => "center"),
    children=[
        html_h1("Surface integrals calculator", style=Dict("textAlign" => "center")),
        dcc_markdown(markdown_text),
            html_div(
                children=[
                    dcc_dropdown(
                        id = "dropdown",
                        options = [(label = i, value = i) for i in options_],
                        value = "parametric representation"
                        ),
                    html_hr(),
                    html_label(id="representation", "Enter your input"),
                    dcc_input(id="input", type="text", style=Dict("width" => "40%"), value="[sqrt(1/4 + u^2) * cos(v), sqrt(1/4 + u^2) * sin(v),  u]"),
                    html_div(id="parametric_bounds", children=[
                        html_label(id="u_range1", "Insert the lower bound for u"),
                        dcc_input(id="u_range1i", type="text", value="-1"),
                        html_label(id="u_range2", "Insert the upper bound for u"),
                        dcc_input(id="u_range2i", type="text", value="1"),
                        html_label(id="v_range1", "Insert the lower bound for v"),
                        dcc_input(id="v_range1i", type="text", value="0"),
                        html_label(id="v_range2", "Insert the upper bound for v"),
                        dcc_input(id="v_range2i", type="text", value="2pi")
                        ], style=Dict("display" => "block")
                        ),
                    html_div(id="f(x,y)_bounds", children=[
                        html_label(id="x_range1", "Insert the lower bound for x"),
                        dcc_input(id="x_range1i", type="text", value="0"),
                        html_label(id="x_range2", "Insert the upper bound for x"),
                        dcc_input(id="x_range2i", type="text", value="4"),
                        html_label(id="y_range1", "Insert the lower bound for y"),
                        dcc_input(id="y_range1i", type="text", value="0"),
                        html_label(id="y_range2", "Insert the upper bound for y"),
                        dcc_input(id="y_range2i", type="text", value="4"),
                        html_label(id="z_range1", "Insert the lower bound for z"),
                        dcc_input(id="z_range1i", type="text", value="0"),
                        html_label(id="z_range2", "Insert the upper bound for z"),
                        dcc_input(id="z_range2i", type="text", value="4")
                    ], style=Dict("display" => "none")
                    ),
                    html_div(id="output_")
                    ], style = (width = "48%", display = "inline-block", float = "left")),
            html_div(
                children=[
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
                    )
                ], style = (width = "48%", height="80%", display = "inline-block")
            )
        ]
    )
end


# callback!(app,
#         Output("input", "value"),
#         Input("radio", "value")
#         ) do user_choice
#         return options_[user_choice]
# end
#

callback!(app,
        Output("parametric_bounds", "style"),
        Output("f(x,y)_bounds", "style"),
        Input("dropdown", "value")
        ) do user_choice
        if user_choice == "parametric representation"
            return Dict("display" => "block"), Dict("display" => "none")
        else
            return Dict("display" => "none"), Dict("display" => "block")
        end
end


callback!(app,
        Output("output_", "children"),
        Input("input", "value"),
        Input("dropdown", "value")
        ) do return_value, dropdown_value
        if dropdown_value == options_[1]
            try
                p = parse_function(return_value, :u, :v)
                value(g) = Φ((x, y, z) -> [x, y, z], g, (0, 1), (0, 1))
                return @eval ($value($p))
            catch e
                return "wrong input"
            end
        end
end


callback!(app,
        Output("surface_plot", "figure"),
        Input("input", "value"),
        Input("dropdown", "value"),
        Input("u_range1i", "value"),
        Input("u_range2i", "value"),
        Input("v_range1i", "value"),
        Input("v_range2i", "value"),
        Input("x_range1i", "value"),
        Input("x_range2i", "value"),
        Input("y_range1i", "value"),
        Input("y_range2i", "value"),
        Input("z_range1i", "value"),
        Input("z_range2i", "value")
        ) do return_value, dropdown_value, u_min, u_max, v_min, v_max, x_min, x_max, y_min, y_max, z_min, z_max
            if dropdown_value == options_[1]

                    try
                        u_min = eval(Meta.parse(u_min))
                        u_max = eval(Meta.parse(u_max))
                        v_min = eval(Meta.parse(v_min))
                        v_max = eval(Meta.parse(v_max))
                    catch e
                        u_min = -1
                        u_max = 1
                        v_min = 0
                        v_max = 2π
                    end

                    if (u_min > u_max) | (v_min > v_max)
                        return Plot(surface(;
                                x = [0],
                                y = [0],
                                z = [[0], [0]]))
                    else
                        N = 100  # interpolation parameter

                        xyz = "_"
                        X(u, v) = 0
                        Y(u, v) = 0
                        Z(u, v) = 0

                        try
                            xyz = split((strip(return_value, ['[', ']'])), ",")
                            X = parse_function(string(xyz[1]), :u, :v)
                            Y = parse_function(string(xyz[2]), :u, :v)
                            Z = parse_function(string(xyz[3]), :u, :v)
                            @eval ($X(-2, 3), $Y(-2, 3), $Z(-2, 3))
                            @eval (isa($X(-2, 3), Array{Float64, 1}))
                            @eval (isa($Y(-2, 3), Array{Float64, 1}))
                            @eval (isa($Z(-2, 3), Array{Float64, 1}))
                        catch e
                            X = parse_function("0", :u, :v)
                            Y = parse_function("0", :u, :v)
                            Z = parse_function("0", :u, :v)
                        end

                        # println(X(2, 3), Y(2, 3), Z(2, 3))
                        vs = range(v_min, v_max, length=N)
                        us = range(u_min, u_max, length=N)
                        value(g) = g.(us', vs)
                        x1 = @eval ($value($X))
                        y1 = @eval ($value($Y))
                        z1 = @eval ($value($Z))

                        return Plot(surface(;
                                x = x1,
                                y = y1,
                                z = z1))
                    end
            else
                return Plot(surface(;
                        x = [0],
                        y = [0],
                        z = [[0], [0]]))
            end
        end


run_server(app, "0.0.0.0")
