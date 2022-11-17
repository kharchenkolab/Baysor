import JSON

function vega_header(title::String="Plot", other::String="")
    return """
      <head>
        <title>$title</title>
        <meta charset="UTF-8">
        <script>$(read(VL.asset("vega.min.js"), String))</script>
        <script>$(read(VL.asset("vega-lite.min.js"), String))</script>
        <script>$(read(VL.asset("vega-embed.min.js"), String))</script>
        $other
      </head>
    """
end

function vega_style()
    return """
    <style media="screen">
      .vega-actions a {
        margin-right: 10px;
        font-family: sans-serif;
        font-size: x-small;
        font-style: italic;
      }
    </style>
    """
end

function vega_plot_html(specs::Dict{String})
    res = """
      <script type="text/javascript">
        var opt = {
          mode: "vega-lite",
          renderer: "$(VL.Vega.RENDERER)",
          actions: $(VL.Vega.ACTIONSLINKS)
        }
    """
    for (i, (divid,spec)) in enumerate(specs)
        res *= "var spec$i = $(JSON.json(VL.add_encoding_types(VL.Vega.getparams(spec))))\nvegaEmbed('#$divid', spec$i, opt);\n\n"
    end

    return res * "</script>"
end
