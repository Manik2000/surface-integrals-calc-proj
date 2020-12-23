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
                    dcc_input(id="input", type="text", style=Dict("width" => "20em"), value="[u + v, v - u, v^2 + u^2]"),
                    html_div(id="ranges"),
                    html_div(id="output_")
                    ])
        ]
    )
end


# callback!(app,
#         Output("input", "value"),
#         Input("radio", "value")
#         ) do user_choice
#         return options_[user_choice]
# end


callback!(app,
        Output("ranges", "children"),
        Input("dropdown", "value")
        ) do user_choice
        if user_choice == "parametric representation"
            return html_div(children=[
                html_label(id="u_range1", "Insert the lower bound for u"),
                dcc_input(id="u_range1i", type="text", value="0"),
                html_label(id="u_range2", "Insert the upper bound for u"),
                dcc_input(id="u_range2i", type="text", value="4"),
                html_label(id="v_range1", "Insert the lower bound for v"),
                dcc_input(id="v_range1i", type="text", value="0"),
                html_label(id="v_range2", "Insert the upper bound for v"),
                dcc_input(id="v_range2i", type="text", value="4")
                ]
                )
        else
            return html_div(children=[
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
            ])
        end
end


callback!(app,
        Output("output_", "children"),
        Input("input", "value"),
        Input("dropdown", "value")
        ) do return_value, dropdown_value
        if dropdown_value == options_[1]
            f = Φ((x, y, z) -> [x, y, z], parse_function(return_value, :u, :v), (0, 1), (0, 1))
            return f
        end
        # @eval foo() = $(Meta.parse("x->x"))
end

run_server(app, "0.0.0.0")
