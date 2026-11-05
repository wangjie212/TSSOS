using Documenter
using TSSOS

makedocs(
sitename = "TSSOS",
pages = ["Home" => "index.md",
         "Polynomial Optimization" => "pop.md",
         "Structures" => "structure.md",
         "Techniques" => "technique.md",
         "Symmetry and Atomic Extraction" => "symmetry_extraction.md",
         "Sum-Of-Squares Optimization" => "sos.md",
         "Polynomial Matrix Optimization" => "pmo.md",
         "Complex Polynomial Optimization" => "cpop.md",
         "Quaternion Polynomial Optimization" => "quaternion.md",
         "Examples" => ["AC Optimal Power Flow" => "opf.md",
         "Sum-Of-Rational-Functions Optimization" => "sorf.md",
         "Inclusion Constants for Free Spectrahedra" => "inclusion.md",          
         ],
         ],
modules = [TSSOS],
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/wangjie212/TSSOS.git"
)
