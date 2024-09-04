

# Traversing a graph

In this guide we will see how to address one possible situation where there is need to
traverse a graph to accumulate information in a relational manner. The key idea in this
example is that a given node of a graph requires accumulated information from all the
descendent nodes of a particular type. This pattern shows up in many different contexts
such as (i) calculating the probability of bud break as affected by apical dominance (without
explict simulation of hormone transport) or (ii) the pipe model (where the cross-sectional
area of a stem at a given point is proportional to the leaf area above that point).

This problem can be solved in a number of ways, which differ in terms of complexity and
performance:

1. Relational query: For every target node, query the graph for all the descendent nodes of
a particular type and accumulate the information. This does not require any change to the
rewriting rules and allows for complex relationships (see tutorials for examples).

2. Annotated topology: Assign each target node and descendents with a level tag such that
we can reconstruct later all the descendents of a particular node (i.e., all nodes of the
same level or higher). This requires a change to the rewriting rules to add the level tag
but the traversal is simple as we simply extract all the relevant nodes from the graph and
perform the logic outside of the graph by filtering for different tag levels.

3. Ad-hoc traversal: Use the `traverse_dfs()` or `traverse_bfs()` functions to write an
ad-hoc traversal algorithm. Since the accumulation must happen in reverse (from the leaf
nodes towards the root node), the ad-hoc function needs to be written in a recursive manner,
most likely keeping its own stack of nodes. This is the most complex solution but also
(potentially) the most efficient one as no nodes need to be extracted and each graph is
visited exactly once.

To illustrate all these approaches, let's create a simple graph that represents the essence
of the problem. We are going to assume to types of nodes (`A` and `B`), where the target
nodes are of type `A` and the descendents of interest include both types. The information to
be accumulated is a simple integer that we will call `state`. The result of accumulating the
states will be stored in the `value`. We also add the `level` tag for second approach. The
types are defined in a module as usual:

```julia
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
```

We specialize `darw()` to include the level tag and the state of each node when we draw the
graph. Remember that this can be achieved by defining methods for `node_label()` for each
type of node.

```julia
PlantGraphs.node_label(node::A, id) = "A: $(node.level) - $(round(node.value, digits = 2))"
PlantGraphs.node_label(node::B, id) = "B: $(node.level)"
g = Graph(axiom = axiom)
draw(g)
```

## Relational queries

The first approach is to use relational queries to accumulate the information as a
side-effect. Note that we do not actually return the nodes as that is not required, we
simple use the query mechanism to access the context of each node so that we can traverse
the graph in a relational manner.

We could also implement this logic inside a rule if the computation performed for each node
would affect whether a particular rule needs to be triggered or not (for example, if we
were calculating the probability of a lateral bud breaking).

First, we construct the function that will be used to accumulate the information and modify
the node in-place as a side-effect. We then return the nodes.

```julia
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
```

## Annotated topology

In the second approach, we assign a level tag to each node so that we can reconstruct the
descendents of a particular node. We already incorporated this when constructing the axiom,
so here we simply extract the nodes and perform the accumulation based on on those level
tags. The queries are very simply, as we just extract all nodes of a given type.

```julia
qA = Query(A)
qB = Query(B)
nodesA = apply(g, qA)
nodesB = apply(g, qB)
```

Now we construct a table that will store the accumulated information for each level as a
dictionary. We know that the total number of levels is equal to the number of `A` nodes.

```julia
accum = zeros(length(nodesA))
```

The logic now is to add the `state` of each `A` or `B` node to all the levels that are equal
or lower (we treat the positions in the array `accum` as the levels). We can then assign
these values to the `value` field of the `A` nodes. Note how for `A` nodes we need to avoid
adding the `state` of the node that defines the level (hence the `1:level-1`) since only
descendents should be used.

```julia
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
```

## Ad-hoc traversal

*Work in progress*

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

