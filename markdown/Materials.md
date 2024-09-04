

# Assigning materials and colors per triangle of a mesh

Alejandro Morales

Centre for Crop Systems Analysis - Wageningen University

In this example we show how to assign different materials and colors to each triangle of a
mesh. This assignment is performed when creating the mesh and storing it in the turtle. The
use is responsible for storing the colors and materials in an object they can access later
(typically in a node inside the graph).

As usual we start with importing the necessary libraries

```julia
using VirtualPlantLab
import ColorTypes: RGB
using Plots
import GLMakie
```

We will create a simple scene with a tiled ellipse. We control the number of triangles of
each tile with the argument `n`, which is then used to generate random materials and colors
for each triangle (of course in a real application these would not be random but depend on
some properties or states of the system). Note that we also keep track of the areas of each
triangle as we will use them later to calculate the irradiance on each triangle.

The colors are zeroed out as later we will scale them by the irradiance absorbed by each
triangle. The materials are randomly generated with different reflectance and transmittance
per triangle.

```julia
mutable struct Tile <: Node
    length::Float64
    width::Float64
    n::Int64
    materials::Vector{Lambertian{1}}
    colors::Vector{RGB}
    areas::Vector{Float64}
end
function Tile(length, width, n)
    Tile(length, width, n, [Lambertian(τ = rand()*0.3, ρ = rand()*0.3) for _ in 1:n],
         zeros(RGB, n), zeros(n))
end
```

In the feed method we first create the mesh with the right dimensions and number of triangles.
We then extract the areas of each triangle in the order of storage and only then assign them
to the turtle alongside the colors and materials using `Mesh!`. Note that if we had directly
used `Ellipse!` to create and add the mesh to the turtle we would not have been able to
access the areas of each triangle (technically we could calculate them by reconstructing the
triangles from the vertices stored in the turtle but this approach is more cumbersome).

```julia
function VirtualPlantLab.feed!(turtle::Turtle, t::Tile, data)
    e = Ellipse(length = t.length, width = t.width, n = t.n)
    t.areas = areas(e) ## Note that we use "areas" and not "area"
    Mesh!(turtle, e, materials = t.materials, colors = t.colors)
end
```

We can now create a simple scene with a tiled floor and a directional light source to test
the calculations of absorbed irradiance.

```julia
g = Graph(axiom = RA(90.0) +  Tile(1.0, 1.0, 40));
sc = Scene(g);
source = DirectionalSource(sc, θ = π/4, Φ = π/2, radiosity = 1.0, nrays = 5_000_000);
nothing #hide
```

Let's run the raytracer and extract the power per triangle and calculate the irradiance.

```julia
rt = RayTracer(sc, source, settings = RTSettings(parallel = true));
trace!(rt)
tile = apply(g, Query(Tile))[1];
pow = [power(m)[1] for m in tile.materials]
irradiance = pow./tile.areas
```

Let's compute the absorptance of each pizza slice. The amount of irradiance absorbed per tile
should scale with the absorptance of the tile.

```julia
τ = [m.τ[1] for m in tile.materials]
ρ = [m.ρ[1] for m in tile.materials]
α = 1.0 .- τ .- ρ
pl = scatter(α, irradiance, ylabel = "Irradiance", xlabel = "Absorptance", legend = false)
Plots.abline!(pl, 1.0, 0.0)
```

Finally, let's scale the colors of each triangle by the irradiance absorbed by each triangle
assuming different shades of green and let's visualize them by turning of the shader in the
renderer to get the exact colors we specify.

```julia
tile.colors = RGB.(0.0, α, 0.0)
sc = Scene(g)
render(sc, wireframe = true, shading = GLMakie.NoShading)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

