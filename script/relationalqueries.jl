using VirtualPlantLab
import GLMakie

module Queries
    using VirtualPlantLab
    struct A <: Node end

    struct C <: Node
        val::Float64
    end

    struct B <: Node
        ID::Int
    end
end
import .Queries as Q

motif(n, i = 0) = Q.A() + (Q.C(45.0) + Q.A() + (Q.C(45.0) +  Q.A() + Q.B(i + 1),
                                           Q.C(-45.0) + Q.A() + Q.B(i + 2),
                                                       Q.A() + Q.B(i + 3)),
                         Q.C(- 45.0) + sum(Q.A() for i in 1:n) + Q.B(i + 4))
axiom =  motif(3, 0) + motif(2, 4) + motif(1, 8) + Q.A() + Q.A() + Q.B(13)
graph = Graph(axiom = axiom);
draw(graph)

VirtualPlantLab.node_label(n::Q.B, id) = "B-$(n.ID)"
VirtualPlantLab.node_label(n::Q.A, id) = "A"
VirtualPlantLab.node_label(n::Q.C, id) = "C"
draw(graph)

Q1 = Query(Q.B)

A1 = apply(graph, Q1)

function Q2_fun(n)
    check, steps = has_ancestor(n, condition = x -> data(x) isa Q.C)
    !check
end

Q2 = Query(Q.B, condition = Q2_fun)
A2 = apply(graph, Q2)

function branch_motif(p)
    data(p) isa Q.A &&
    has_descendant(p, condition = x -> data(x) isa Q.C && data(x).val < 0.0)[1] &&
    has_ancestor(p, condition = x -> data(x) isa Q.C && data(x).val > 0.0)[1]
end

function Q3_fun(n, nsteps)
    # Condition 1
    check, steps = has_ancestor(n, condition = branch_motif)
    !check && return false
    # Condition 2
    p = parent(n, nsteps = steps)
    check, steps = has_ancestor(p, condition = is_root)
    steps != nsteps && return false
    return true
end

Q3 = Query(Q.B, condition = n -> Q3_fun(n, 2))
A3 = apply(graph, Q3)

function Q4_fun(n)
    !has_ancestor(n, condition = x -> is_root(x) && length(children(x)) > 1)[1]
end

Q4 = Query(Q.B, condition = Q4_fun)
A4 = apply(graph, Q4)

function Q5_fun(n)
    check, steps = has_ancestor(n, condition = is_root)
    steps == 4
end

Q5 = Query(Q.B, condition = Q5_fun)
A5 = apply(graph, Q5)

function Q6_fun(n, nA)
    check = Q3_fun(n, nA)
    !check && return false
    p2 = parent(n, nsteps = 2)
    data(p2) isa Q.A
end

Q6 = Query(Q.B, condition = n -> Q6_fun(n, 3))
A6 = apply(graph, Q6)

Q7 = Query(Q.B, condition = n -> Q6_fun(n, 4) || Q2_fun(n))
A7 = apply(graph, Q7)

function Q8_fun(n)
    p1 = parent(n)
    p2 = parent(n, nsteps = 2)
    p3 = parent(n, nsteps = 3)
    data(p1) isa Q.A && data(p2) isa Q.C && data(p2).val > 0.0 && data(p3) isa Q.A
end

Q8 = Query(Q.B, condition = Q8_fun)
A8 = apply(graph, Q8)

function Q9_fun(n)
    p1 = parent(n)
    p2 = parent(n, nsteps = 2)
    p3 = parent(n, nsteps = 3)
    p4 = parent(n, nsteps = 4)
    data(p1) isa Q.A && data(p2) isa Q.C && data(p2).val < 0.0 &&
       data(p3) isa Q.A && data(p4) isa Q.C
end

Q9 = Query(Q.B, condition = Q9_fun)
A9 = apply(graph, Q9)

function Q10_fun(n)
    Q6_fun(n, 3) && return true ## Check node 7
    Q9_fun(n) && has_ancestor(n, condition = is_root)[2] == 6 && return true ## Check node 6
    has_ancestor(n, condition = is_root)[2] == 5 && data(parent(n, nsteps = 3)) isa Q.C && return true ## Check node 8 (and not 4!)
end

Q10 = Query(Q.B, condition = Q10_fun)
A10 = apply(graph, Q10)

function Q11_fun(n)
    Q5_fun(n) && return true ## 3
    Q6_fun(n, 3) && return true ## 7
    Q6_fun(n, 4) && return true ## 11
    has_ancestor(n, condition = is_root)[2] == 5 && data(parent(n, nsteps = 2)) isa Q.C &&
        data(parent(n, nsteps = 4)) isa Q.A && return true ## 12
end

Q11 = Query(Q.B, condition = Q11_fun)
A11 = apply(graph, Q11)

function Q12_fun(n)
    Q6_fun(n, 3) && return true # 7
    has_ancestor(n, condition = is_root)[2] == 5 && data(parent(n, nsteps = 2)) isa Q.C &&
        data(parent(n, nsteps = 4)) isa Q.A && return true ## 12
end

Q12 = Query(Q.B, condition = Q12_fun)
A12 = apply(graph, Q12)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
