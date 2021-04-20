import JSON
using VegaLite

function vega_header(title::String="Plot")
    return """
      <head>
        <title>$title</title>
        <meta charset="UTF-8">
        <script>$(read(VegaLite.asset("vega.min.js"), String))</script>
        <script>$(read(VegaLite.asset("vega-lite.min.js"), String))</script>
        <script>$(read(VegaLite.asset("vega-embed.min.js"), String))</script>
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

function vega_plot_html(specs::Dict{String, VegaLite.VLSpec})
    res = """
      <script type="text/javascript">
        var opt = {
          mode: "vega-lite",
          renderer: "$(VegaLite.Vega.RENDERER)",
          actions: $(VegaLite.Vega.ACTIONSLINKS)
        }
    """
    for (i, (divid,spec)) in enumerate(specs)
        res *= "var spec$i = $(JSON.json(VegaLite.add_encoding_types(VegaLite.Vega.getparams(spec))))\nvegaEmbed('#$divid', spec$i, opt);\n\n"
    end

    return res * "</script>"
end
