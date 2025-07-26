using Documenter

# For reproducibility
mkpath(joinpath(@__DIR__, "src", "assets"))
cp(
    joinpath(@__DIR__, "Manifest.toml"),
    joinpath(@__DIR__, "src", "assets", "Manifest.toml");
    force=true,
)
cp(
    joinpath(@__DIR__, "Project.toml"),
    joinpath(@__DIR__, "src", "assets", "Project.toml");
    force=true,
)

# For binder
using Literate
OUTPUT_MD = joinpath(@__DIR__, "src")
OUTPUT_NB = joinpath(@__DIR__, "src", "notebooks")
OUTPUT_JL = joinpath(@__DIR__, "src", "scripts")
files = ["index.jl", "bistable.jl", "oscillator.jl"]
function markdown_postprocess(content)
    content = replace(content, "gh-pages?filepath=" => "gh-pages?urlpath=%2Fdoc%2Ftree%2F") # change jupyter notebook to jupyter lab
    return content
end
for file in files
    INPUT = joinpath(@__DIR__, "src-literate", file)
    Literate.markdown(INPUT, OUTPUT_MD; postprocess=markdown_postprocess)
    Literate.notebook(INPUT, OUTPUT_NB; execute=true)
    Literate.script(INPUT, OUTPUT_JL)
end

# Files to copy
function copy_file(file, dir_source, dir_destination)
    cp(
        joinpath(@__DIR__, dir_source, file),
        joinpath(@__DIR__, dir_destination, file);
        force=true,
    )
    return nothing
end
files = ["bistable.png", "openloop.png"]
for file in files
    copy_file(file, "src-literate", "src")
    copy_file(file, "src-literate", joinpath("src", "notebooks"))
end

# 
repo_url = joinpath("github.com", "agustinyabo", "PWLdynamics.jl")

makedocs(;
    draft=false, # if draft is true, then the julia code from .md is not executed
    # to disable the draft mode in a specific markdown file, use the following:
    # ```@meta
    # Draft = false
    # ```
    remotes=nothing,
    warnonly=:cross_references,
    sitename="PWLdynamics",
    format=Documenter.HTML(;
        repolink="https://" * repo_url,
        prettyurls=false,
        size_threshold_ignore=["index.md", "bistable.md", "oscillator.md"],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=[
        "Introduction" => "index.md",
        "Bistable toggle switch" => "bistable.md",
        "Damped genetic oscillator" => "oscillator.md",
    ],
)

deploydocs(; push_preview=true, repo=repo_url * ".git", devbranch="main")
