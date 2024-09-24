using Documenter
using Baysor

makedocs(
    sitename="Baysor Documentation",
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Running Baysor" => [
            "Summary" => "run.md",
            "Cell segmentation" => "segmentation.md",
            "Preview" => "preview.md",
            "Segmentation-free analysis" => "segfree.md",
            "Advanced configuration" => "configuration.md",
        ],
        "Citation info" => "citation.md",
        "Examples" => "examples.md"
        # TODO:
        # API
    ]
)

deploydocs(repo = "github.com/kharchenkolab/Baysor.git")