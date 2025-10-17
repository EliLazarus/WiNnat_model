using Plots

nodesute = Dict(
    "Utility" => (.5, 6),
    "Final Consumption (σ₁)" => (0.5, 5.5),
    "Transportation (σ₂)" => (-1.5, 5),
    "Non-Transportation (σ₃)" => (2.2, 5),
    "Consumer goods\n and services (σ₅)" => (3.7, 4.4),
    "Non-personal\n transportation" => (-2.2,4.4),
    "Housing \nexpenditures (σ₆)" => (1.7, 4.4),
    "Personal \ntransportation (σ₄)" => (-.6,4.4),
    "d1"=> (2.8,3.7),
    "d2"=> (3.4,3.7),
    "d3..."=> (4.,3.7),
    "..dn"=> (4.6,3.7),
    "Home energy \nexpenditures (σ₈)" => (2.6, 3.1),
    "Owner-occupied and \nrental expenditures" => (.9, 3.7),
    "Electricity" => (1.8, 2.4),
    "Vehicle and service \nexpenditures" => (-1.5,3.7),
    "Fuels (σ₇)" => (-0.1,3.1),
    "Home fuels (σ₉)" => (3.4, 2.4),
    "Natural gas" => (2.7, 1.7),
    "Fuel oil" => (4.3, 1.7),
    "Motor vehicle fuels" => (-0.9,2.4),
    "Veh \nelectricity" => (.7,2.4),
    )

    
# Define edges (connections)
edgesute = [
    ("Utility", "Final Consumption (σ₁)"),
    ("Final Consumption (σ₁)", "Transportation (σ₂)"),
    ("Final Consumption (σ₁)", "Non-Transportation (σ₃)"),
    ("Non-Transportation (σ₃)", "Consumer goods\n and services (σ₅)"),
    ("Non-Transportation (σ₃)", "Housing \nexpenditures (σ₆)"),
    ("Housing \nexpenditures (σ₆)", "Home energy \nexpenditures (σ₈)"),
    ("Housing \nexpenditures (σ₆)", "Owner-occupied and \nrental expenditures"),
    ("Home energy \nexpenditures (σ₈)", "Electricity"),
    ("Home energy \nexpenditures (σ₈)", "Home fuels (σ₉)"),
    ("Home fuels (σ₉)","Natural gas"),
    ("Home fuels (σ₉)","Fuel oil"),
    ("Transportation (σ₂)", "Personal \ntransportation (σ₄)"),
    ("Transportation (σ₂)", "Non-personal\n transportation"),
    ("Personal \ntransportation (σ₄)", "Vehicle and service \nexpenditures"),
    ("Personal \ntransportation (σ₄)", "Fuels (σ₇)"),
    ("Fuels (σ₇)", "Motor vehicle fuels"),
    ("Fuels (σ₇)", "Veh \nelectricity"),
    ("Consumer goods\n and services (σ₅)", "d1"),
    ("Consumer goods\n and services (σ₅)", "d2"),
    ("Consumer goods\n and services (σ₅)", "d3..."),
    ("Consumer goods\n and services (σ₅)", "..dn"),
]

# Create scatter plot for nodes
x = [pos[1] for pos in values(nodesute)]
y = [pos[2] for pos in values(nodesute)]
labels = collect(keys(nodesute))
y2 = deepcopy(y) .+.25

plot(x, y, seriestype = :scatter, markersize = 4, markercolor=:grey, label = "", legend = false, yaxis=false, xaxis=false, grid=false)

# Draw arrows for edges
for (from, to) in edgesute
    x0, y0 = nodesute[from]
    x1, y1 = nodesute[to]
    plot!([x0, x1], [y0, y1], series = :line, lw = 1.7, color = :grey, label = "")
end

# title!("Nested CES Tree Structure (in σ terms)")
plot!()
annotate!([(x[i], y2[i], (labels[i],9,:steelblue, "Palatino Roman")) for i in 1:length(x)], xlim=(-2.5,5),ylim=(1.5,6.5)) 
# =(10,"Palatino Roman"))

# png("./Results/UtilityTree.png")











# ### For the Intermediate Production Nesting Structure

# nodes = Dict(
#     "Yj" => (0, 6),
#     "Top CES (σ)" => (0, 5.5),
#     "ID" => (-1.5, 5),
#     "VAE (σ_VAE)" => (1.5, 5),
#     "VA" => (.5, 4),
#     "Energy (σ_E)" => (2.5, 4),
#     "Fossil" => (1.5, 3),
#     "Electricity" => (3.5, 3),
#     "Rnw" => (4.4, 2),
#     "Coal" => (3.5, 2),    
#     "Gas" => (2.5, 2),
#     "Oil" => (.5, 2),
#     "y1" => (-2.4,4),
#     "y2" => (-1.8,4),
#     "y3..." => (-1.2,4),
#     "...y72" => (-.6,4)
# )

# # Define edges (connections)
# edges = [
#     ("Yj", "Top CES (σ)"),
#     ("Top CES (σ)", "ID"),
#     ("Top CES (σ)", "VAE (σ_VAE)"),
#     ("VAE (σ_VAE)", "VA"),
#     ("VAE (σ_VAE)", "Energy (σ_E)"),
#     ("Energy (σ_E)", "Fossil"),
#     ("Energy (σ_E)", "Electricity"),
#     ("Electricity","Rnw"),
#     ("Electricity","Coal"),
#     ("Electricity","Gas"),
#     ("Fossil", "Gas"),
#     ("Fossil", "Oil"),
#     ("ID", "y1"),
#     ("ID", "y2"),
#     ("ID", "y3..."),
#     ("ID", "...y72")
# ]

# # Create scatter plot for nodes
# x = [pos[1] for pos in values(nodes)]
# y = [pos[2] for pos in values(nodes)]
# labels = collect(keys(nodes))
# y2 = deepcopy(y) .+.25

# plot(x, y, seriestype = :scatter, markersize = 4, markercolor=:grey, label = "", legend = false, yaxis=false, xaxis=false, grid=false)
# annotate!([(x[i], y2[i], (labels[i],11,:steelblue, "Palatino Roman")) for i in 1:length(x)], xlim=(-2.5,5),ylim=(1.5,6.5)) 
# # =(10,"Palatino Roman"))

# # Draw arrows for edges
# for (from, to) in edges
#     x0, y0 = nodes[from]
#     x1, y1 = nodes[to]
#     plot!([x0, x1], [y0, y1], series = :line, lw = 1.5, color = :black, label = "")
# end

# # title!("Nested CES Tree Structure (in σ terms)")
# plot!()
# # png("./Results/IntDemandtree.png")