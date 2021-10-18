import Pkg

Pkg.activate("../output/") ##Add packages to this project (https://docs.julialang.org/en/v1/stdlib/Pkg/index.html)

for package in ["CSV", "LaTeXStrings", "Interpolations", "Plots", "Roots"]
  Pkg.add(package)
end
