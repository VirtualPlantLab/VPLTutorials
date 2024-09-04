using VirtualPlantLab
import ColorTypes: RGB
using Plots
import GLMakie

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

function VirtualPlantLab.feed!(turtle::Turtle, t::Tile, data)
    e = Ellipse(length = t.length, width = t.width, n = t.n)
    t.areas = areas(e) ## Note that we use "areas" and not "area"
    Mesh!(turtle, e, materials = t.materials, colors = t.colors)
end

g = Graph(axiom = RA(90.0) +  Tile(1.0, 1.0, 40));
sc = Scene(g);
source = DirectionalSource(sc, θ = π/4, Φ = π/2, radiosity = 1.0, nrays = 5_000_000);

rt = RayTracer(sc, source, settings = RTSettings(parallel = true));
trace!(rt)
tile = apply(g, Query(Tile))[1];
pow = [power(m)[1] for m in tile.materials]
irradiance = pow./tile.areas

τ = [m.τ[1] for m in tile.materials]
ρ = [m.ρ[1] for m in tile.materials]
α = 1.0 .- τ .- ρ
pl = scatter(α, irradiance, ylabel = "Irradiance", xlabel = "Absorptance", legend = false)
Plots.abline!(pl, 1.0, 0.0)

tile.colors = RGB.(0.0, α, 0.0)
sc = Scene(g)
render(sc, wireframe = true, shading = GLMakie.NoShading)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
