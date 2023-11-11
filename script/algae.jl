using VirtualPlantLab

module algae
    import VirtualPlantLab: Node
    struct A <: Node end
    struct B <: Node end
end
import .algae

axiom = algae.A()

rule1 = Rule(algae.A, rhs = x -> algae.A() + algae.B())
rule2 = Rule(algae.B, rhs = x -> algae.A())

function rule_1(x)
    algae.A() + algae.B()
end

organism = Graph(axiom = axiom, rules = (rule1, rule2))

rewrite!(organism)

import GLMakie
draw(organism)

for i in 1:4
    rewrite!(organism)
end

draw(organism)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
