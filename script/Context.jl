using VirtualPlantLab
module types
    using VirtualPlantLab
    struct Cell <: Node
        state::Int64
    end
end
import .types: Cell
function transfer(context)
    if has_parent(context)
        return (true, (parent(context), ))
    else
        return (false, ())
    end
end
rule = Rule(Cell, lhs = transfer, rhs = (context, father) -> Cell(data(father).state), captures = true)
axiom = Cell(1) + Cell(0) + Cell(0)
pop = Graph(axiom = axiom, rules = rule)

getCell = Query(Cell)
apply(pop, getCell)

rewrite!(pop)
apply(pop, getCell)

rewrite!(pop)
apply(pop, getCell)

pop  = Graph(axiom = axiom, rules = rule)
states = Int64[]
traverse_dfs(pop, fun = node -> push!(states, node.state))
states

rewrite!(pop)
states = Int64[]
traverse_dfs(pop, fun = node -> push!(states, node.state))
states

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
