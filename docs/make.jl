using Documenter
# Pkg.activate("E:\\Programs\\blockpop\\TSSOS")
using TSSOS

makedocs(
# root = "E:\\Programs\\blockpop\\TSSOS\\docs",
sitename = "TSSOS",
pages = ["Home" => "index.md",
         "Polynomial Optimization" => "spop.md",
         "Noncommutative Polynomial Optimization" => "ncpop.md"],
modules = [TSSOS],
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/wangjie212/TSSOS.git"
)
