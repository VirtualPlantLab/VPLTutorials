using VirtualPlantLab
import GLMakie # Import rather than "using" to avoid masking Scene
using ColorTypes # To define colors for the rendering
module sn
    import VirtualPlantLab
    struct E <: VirtualPlantLab.Node
        length::Float64
    end
end
import .sn

const L = 1.0
axiom = sn.E(L) + VirtualPlantLab.RU(120.0) + sn.E(L) + VirtualPlantLab.RU(120.0) + sn.E(L)

function Kochsnowflake(x)
    L = data(x).length
    sn.E(L/3) + RU(-60.0) + sn.E(L/3) + RU(120.0) + sn.E(L/3) + RU(-60.0) + sn.E(L/3)
end
rule = Rule(sn.E, rhs = Kochsnowflake)

Koch = Graph(axiom = axiom, rules = Tuple(rule))

function VirtualPlantLab.feed!(turtle::Turtle, e::sn.E, vars)
    HollowCylinder!(turtle, length = e.length, width = e.length/10,
                    height = e.length/10, move = true,
                    colors = RGB(rand(), rand(), rand()))
    return nothing
end

sc = Scene(Koch)
render(sc, axes = false)

rewrite!(Koch)
render(Scene(Koch), axes = false)

for i in 1:3
    rewrite!(Koch)
end
render(Scene(Koch), axes = false)

function Kochsnowflake2(x)
   L = data(x).length
   sn.E(L/3) + RU(60.0) + sn.E(L/3) + RU(-120.0) + sn.E(L/3) + RU(60.0) + sn.E(L/3)
end
rule2 = Rule(sn.E, rhs = Kochsnowflake2)
Koch2 = Graph(axiom = axiom, rules = Tuple(rule2))

rewrite!(Koch2)
render(Scene(Koch2), axes = false)

rewrite!(Koch2)
render(Scene(Koch2), axes = false)

rewrite!(Koch2)
render(Scene(Koch2), axes = false)

axiomCesaro = sn.E(L) + RU(90.0) + sn.E(L) + RU(90.0) + sn.E(L) + RU(90.0) + sn.E(L)
Cesaro = Graph(axiom = axiomCesaro, rules = (rule2,))
render(Scene(Cesaro), axes = false)

rewrite!(Cesaro)
render(Scene(Cesaro), axes = false)

rewrite!(Cesaro)
render(Scene(Cesaro), axes = false)

rewrite!(Cesaro)
render(Scene(Cesaro), axes = false)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
