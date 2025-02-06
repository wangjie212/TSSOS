using Documenter
using TSSOS

makedocs(
sitename = "TSSOS",
pages = ["Home" => "index.md",
         "Polynomial Optimization" => "pop.md",
         "Structures" => "structure.md",
         "Techniques" => "technique.md",
         "Sum-Of-Squares Optimization" => "sos.md",
         "Polynomial Matrix Optimization" => "pmo.md",
         "Complex Polynomial Optimization" => "cpop.md",
         "Examples" => ["AC Optimal Power Flow" => "opf.md",
         "Sum-Of-Rational-Functions Optimization" => "sorf.md",
         ],
         ],
modules = [TSSOS],
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/wangjie212/TSSOS.git"
)
