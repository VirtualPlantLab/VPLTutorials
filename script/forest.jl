using VirtualPlantLab
using Distributions, Plots, ColorTypes
import GLMakie
# Data types
module TreeTypes
    import VirtualPlantLab
    # Meristem
    struct Meristem <: VirtualPlantLab.Node end
    # Bud
    struct Bud <: VirtualPlantLab.Node end
    # Node
    struct Node <: VirtualPlantLab.Node end
    # BudNode
    struct BudNode <: VirtualPlantLab.Node end
    # Internode (needs to be mutable to allow for changes over time)
    Base.@kwdef mutable struct Internode <: VirtualPlantLab.Node
        length::Float64 = 0.10 # Internodes start at 10 cm
    end
    # Leaf
    Base.@kwdef struct Leaf <: VirtualPlantLab.Node
        length::Float64 = 0.20 # Leaves are 20 cm long
        width::Float64  = 0.1 # Leaves are 10 cm wide
    end
    # Graph-level variables
    Base.@kwdef struct treeparams
        growth::Float64 = 0.1
        budbreak::Float64 = 0.25
        phyllotaxis::Float64 = 140.0
        leaf_angle::Float64 = 30.0
        branch_angle::Float64 = 45.0
    end
end

import .TreeTypes

function VirtualPlantLab.feed!(turtle::Turtle, i::TreeTypes.Internode, data)
    # Rotate turtle around the head to implement elliptical phyllotaxis
    rh!(turtle, data.phyllotaxis)
    HollowCylinder!(turtle, length = i.length, height = i.length/15, width = i.length/15,
                move = true, color = RGB(0.5,0.4,0.0))
    return nothing
end

function VirtualPlantLab.feed!(turtle::Turtle, l::TreeTypes.Leaf, data)
    # Rotate turtle around the arm for insertion angle
    ra!(turtle, -data.leaf_angle)
    # Generate the leaf
    Ellipse!(turtle, length = l.length, width = l.width, move = false,
             color = RGB(0.2,0.6,0.2))
    # Rotate turtle back to original direction
    ra!(turtle, data.leaf_angle)
    return nothing
end

function VirtualPlantLab.feed!(turtle::Turtle, b::TreeTypes.BudNode, data)
    # Rotate turtle around the arm for insertion angle
    ra!(turtle, -data.branch_angle)
end

meristem_rule = Rule(TreeTypes.Meristem, rhs = mer -> TreeTypes.Node() +
                                              (TreeTypes.Bud(), TreeTypes.Leaf()) +
                                         TreeTypes.Internode() + TreeTypes.Meristem())

function prob_break(bud)
    # We move to parent node in the branch where the bud was created
    node =  parent(bud)
    # We count the number of internodes between node and the first Meristem
    # moving down the graph
    check, steps = has_descendant(node, condition = n -> data(n) isa TreeTypes.Meristem)
    steps = Int(ceil(steps/2)) # Because it will count both the nodes and the internodes
    # Compute probability of bud break and determine whether it happens
    if check
        prob =  min(1.0, steps*graph_data(bud).budbreak)
        return rand() < prob
    # If there is no meristem, an error happened since the model does not allow
    # for this
    else
        error("No meristem found in branch")
    end
end
branch_rule = Rule(TreeTypes.Bud,
            lhs = prob_break,
            rhs = bud -> TreeTypes.BudNode() + TreeTypes.Internode() + TreeTypes.Meristem())

function create_tree(origin, growth, budbreak, orientation)
    axiom = T(origin) + RH(orientation) + TreeTypes.Internode() + TreeTypes.Meristem()
    tree =  Graph(axiom = axiom, rules = (meristem_rule, branch_rule),
                  data = TreeTypes.treeparams(growth = growth, budbreak = budbreak))
    return tree
end

getInternode = Query(TreeTypes.Internode)

function elongate!(tree, query)
    for x in apply(tree, query)
        x.length = x.length*(1.0 + data(tree).growth)
    end
end

function growth!(tree, query)
    elongate!(tree, query)
    rewrite!(tree)
end

function simulate(tree, query, nsteps)
    new_tree = deepcopy(tree)
    for i in 1:nsteps
        growth!(new_tree, query)
    end
    return new_tree
end

origins = [Vec(i,j,0) for i = 1:2.0:20.0, j = 1:2.0:20.0]

orientations = [rand()*360.0 for i = 1:2.0:20.0, j = 1:2.0:20.0]

growths = rand(LogNormal(-2, 0.3), 10, 10)
histogram(vec(growths))

budbreaks = rand(Beta(2.0, 10), 10, 10)
histogram(vec(budbreaks))

forest = vec(create_tree.(origins, growths, budbreaks, orientations));

newforest = [simulate(tree, getInternode, 2) for tree in forest];

render(Scene(newforest))

newforest = [simulate(tree, getInternode, 15) for tree in newforest];
render(Scene(newforest))

using Base.Threads
newforest = deepcopy(forest)
@threads for i in eachindex(forest)
    newforest[i] = simulate(forest[i], getInternode, 6)
end
render(Scene(newforest), parallel = true)

newforest = deepcopy(forest)
for step in 1:15
    @threads for i in eachindex(newforest)
        newforest[i] = simulate(newforest[i], getInternode, 1)
    end
end
render(Scene(newforest), parallel = true)

scene = Scene(newforest);

soil = Rectangle(length = 21.0, width = 21.0)
rotatey!(soil, pi/2)
VirtualPlantLab.translate!(soil, Vec(0.0, 10.5, 0.0))

VirtualPlantLab.add!(scene, mesh = soil, color = RGB(1,1,0))

render(scene, axes = false)

res = calculate_resolution(width = 16.0, height = 16.0, dpi = 1_000)
output = render(scene, axes = false, resolution = res)
export_scene(scene = output, filename = "nice_trees.png")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
