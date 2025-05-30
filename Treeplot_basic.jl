nodes = Dict(
    "Yj" => (0, 6),
    "Top CES (σ)" => (0, 5.5),
    "ID" => (-1.5, 5),
    "VAE (σ_VAE)" => (1.5, 5),
    "VA" => (.5, 4),
    "Energy (σ_E)" => (2.5, 4),
    "Fossil" => (1.5, 3),
    "Electricity" => (3.5, 3),
    "Rnw" => (4.4, 2),
    "Coal" => (3.5, 2),    
    "Gas" => (2.5, 2),
    "Oil" => (.5, 2),
    "y1" => (-2.4,4),
    "y2" => (-1.8,4),
    "y3..." => (-1.2,4),
    "...y72" => (-.6,4)
)

# Define edges (connections)
edges = [
    ("Yj", "Top CES (σ)"),
    ("Top CES (σ)", "ID"),
    ("Top CES (σ)", "VAE (σ_VAE)"),
    ("VAE (σ_VAE)", "VA"),
    ("VAE (σ_VAE)", "Energy (σ_E)"),
    ("Energy (σ_E)", "Fossil"),
    ("Energy (σ_E)", "Electricity"),
    ("Electricity","Rnw"),
    ("Electricity","Coal"),
    ("Electricity","Gas"),
    ("Fossil", "Gas"),
    ("Fossil", "Oil"),
    ("ID", "y1"),
    ("ID", "y2"),
    ("ID", "y3..."),
    ("ID", "...y72")
]

# Create scatter plot for nodes
x = [pos[1] for pos in values(nodes)]
y = [pos[2] for pos in values(nodes)]
labels = collect(keys(nodes))
y2 = deepcopy(y) .+.25

plot(x, y, seriestype = :scatter, markersize = 4, markercolor=:grey, label = "", legend = false, yaxis=false, xaxis=false, grid=false)
annotate!([(x[i], y2[i], (labels[i],11,:steelblue, "Palatino Roman")) for i in 1:length(x)], xlim=(-2.5,5),ylim=(1.5,6.5)) 
# =(10,"Palatino Roman"))

# Draw arrows for edges
for (from, to) in edges
    x0, y0 = nodes[from]
    x1, y1 = nodes[to]
    plot!([x0, x1], [y0, y1], series = :line, lw = 1.5, color = :black, label = "")
end

# title!("Nested CES Tree Structure (in σ terms)")
plot!()
# png("./Results/IntDemandtree.png")