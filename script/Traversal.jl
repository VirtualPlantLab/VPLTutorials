using VirtualPlantLab
import GLMakie
module TravTypes
    import PlantGraphs: Node
    export A, B
    @kwdef mutable struct A <: Node
        value::Float64 = 0.0
        state::Float64 = 0.0
        level::Int     = 0
    end
    @kwdef struct B <: Node
        state::Float64 = 0.0
        level::Int     = 0
    end
end

using .TravTypes

# Function that generates the arguments to create a node
fill(level) = (level = level, state = rand())
# Motif of the graph (`;fill(l)...` is equivalent to typing `level = l, state = rand()`)
motif(l) = A(;fill(l)...) + (B(;fill(l)...), B(;fill(l)...) + (B(;fill(l)...), B(;fill(l)...)))
# Create the graph by repeating the motif 10 times
axiom = motif(1)
for i in 2:10
    axiom += motif(i)
end

PlantGraphs.node_label(node::A, id) = "A: $(node.level) - $(round(node.value, digits = 2))"
PlantGraphs.node_label(node::B, id) = "B: $(node.level)"
g = Graph(axiom = axiom)
draw(g)

function accumulate(node)
    # Extract the first child and put it into a stack
    stack = [children(node)...]
    state = 0.0
    while !isempty(stack)
        # Pop the last child from the stack
        child = pop!(stack)
        # Accumulate the state
        state += data(child).state
        # Add the children of the child to the stack
        if has_children(child)
            for child in children(child)
                push!(stack, child)
            end
        end
    end
    data(node).value = state
    return true
end

# We not construct the query object and apply it
q = Query(A, condition = accumulate)
apply(g, q)
draw(g)

qA = Query(A)
qB = Query(B)
nodesA = apply(g, qA)
nodesB = apply(g, qB)

accum = zeros(length(nodesA))

for node in nodesB
    level = node.level
    accum[1:level] .+= node.state ## Add node.state to elements 1:level
end
for node in nodesA
    level = node.level
    level > 1 && (accum[1:level-1] .+= node.state)
end
for node in nodesA
    node.value = accum[node.level]
end
draw(g)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
