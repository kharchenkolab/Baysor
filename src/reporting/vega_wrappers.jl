import JSON
using Pkg.Artifacts

asset(url...) = normpath(joinpath(artifact"vegalite_app", "minified", url...))

function vega_header(title::String="Plot", other::String="")
    return """
      <head>
        <title>$title</title>
        <meta charset="UTF-8">
        <script>$(read(asset("vega.min.js"), String))</script>
        <script>$(read(asset("vega-lite.min.js"), String))</script>
        <script>$(read(asset("vega-embed.min.js"), String))</script>
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

function vega_plot_html(specs::Dict{String}) # The dict is deliberately abstract, as these specs are created as Dict{String, Any}
    res =
    """
      <script type="text/javascript">
        var opt = {mode: "vega-lite"}
        function showError(el, error){
          el.innerHTML = ('<div class="error" style="color:red;">'
                          + '<p>JavaScript Error: ' + error.message + '</p>'
                          + "<p>This usually means there's a typo in your chart specification. "
                          + "See the javascript console for the full traceback.</p>"
                          + '</div>');
          throw error;
        }
    """
    for (i, (divid,spec)) in enumerate(specs)
        res *= "var spec$i = $(JSON.json(spec))\nvegaEmbed('#$divid', spec$i, opt).catch(error => showError(el, error));\n\n"
    end

    return res * "</script>"
end
