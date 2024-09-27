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
    ],
    format=Documenter.HTML(;
        assets=[
            asset(
                "https://cloud.umami.is/script.js", class=:js,
                attributes=Dict(Symbol("data-website-id") => "b9c0c023-9682-4a74-8d3e-81a1d9a19e98", :defer => "")
            )
	    ],
    )
)

if "deploy" in ARGS
    ENV["TRAVIS_REPO_SLUG"] = "github.com/kharchenkolab/Baysor.git"
    ENV["TRAVIS_PULL_REQUEST"] = "false"
    ENV["TRAVIS_OS_NAME"] = "linux"
    ENV["TRAVIS_JULIA_VERSION"] = "1.1"
    ENV["TRAVIS_TAG"] = ""
    ENV["TRAVIS_BRANCH"] = "master"

    deploydocs(
        repo="github.com/kharchenkolab/Baysor.git",
        target="build",
        branch="gh-pages",
        devbranch="master",
    )
end