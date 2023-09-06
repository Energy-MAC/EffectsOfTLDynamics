using Plots

# Create some sample data
x = 1:10
y = sin.(x)
y2 = cos.(x)

# Create a 3x3 grid of plots
plots = []

for i in 1:3
    for j in 1:3
        p = plot(x, y, label="Plot $i,$j")
        plot!(p, x, y2,label="Plot $i,$j")
        plot!(p, legend = true)
        push!(plots, p)
    end
end

# Combine the individual plots into a single layout
combined_plot = plot(plots..., layout=(3, 3))