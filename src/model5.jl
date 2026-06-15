using Catalyst, Latexify, DifferentialEquations, BifurcationKit, HomotopyContinuation, Plots, CairoMakie 

# Load parameters and plotting things
include("parameters_model5.jl")
include("plotting.jl")

# Define 'reaction network'
model5 = @reaction_network begin
	# modLDL
	σ_l, ∅ => l
	mm(l,ν_l*m,ξ_l), l => ∅
	d_l, l --> ∅
	# p cytokine
	mm(l,σ_p1,β_p), ∅ => p
	σ_p2*q, ∅ => p
	d_p, p --> ∅
	mm(l,ν_p*m,ξ_l), ∅ => p
	# q cytokine
	mm(l,ν_q*m,ξ_l), ∅ => q
	d_q, q --> ∅
	σ_q, q --> ∅
	# macrophage
	σ_m*(1+α*q)*(p), ∅ => m 
	mm(l,ν_m*m,ξ_l), m => ∅
	mm(h,θ*ν_N*N,ξ_h), ∅ => m
	d_m, m --> ∅
	# foam cells
	mm(l,ν_m*m,ξ_l), ∅ => N
	mm(h,ν_N*N,ξ_h), N => ∅
	d_N, N --> ∅
	# HDL
	σ_h, ∅ => h
	mm(h,ν_h*N,ξ_h), h => ∅
	d_h, h --> ∅
end

# Set initial conditions, time span and parameter values
u0 = [:l => 0.0, :p => 0.0, :q => 0.0, :m => 0.0, :N => 0.0, :h => 0.0]
tspan = (0.0, 20.0)
ps = [ # decay rates
	 :d_m => d_m, :d_l => d_l, :d_h => d_h, :d_N => d_N,
	 :d_p => d_p, :d_q => d_q,
	 # Hill scaling factor
	 :ν_l => ν_l, :ν_m => ν_m, :ν_p => ν_p,
	 :ν_q => ν_q, :ν_h => ν_h, :ν_N => ν_N,
	 # recruitment rates
	 :σ_m => σ_m, :σ_l => σ_l, :σ_p1 => σ_p1, :σ_p2 => σ_p2,
	 :σ_q => σ_q, :σ_h => σ_h,
	 # saturation parameters
     :ξ_l => ξ_l, :ξ_h => ξ_h, :β_p => β_p,   
	 # other parameters
	 :θ => θ, :α => α
	 ]

# Define ODE problem
oprob = ODEProblem(model5, u0, tspan, ps)
# Solve
sol = DifferentialEquations.solve(oprob)
fig = Figure(fontsize=25)
ax = Axis(fig[1,1], xlabel="Time", ylabel="Concentration") # these aren't actually concentrations in every case though...
lines!(ax,sol.t,sol[1,:], label="l", linewidth=2)
lines!(ax,sol.t,sol[2,:], label="p", linewidth=2)
lines!(ax,sol.t,sol[3,:], label="q", linewidth=2)
lines!(ax,sol.t,sol[4,:], label="m", linewidth=2)
lines!(ax,sol.t,sol[5,:], label="N", linewidth=2)
lines!(ax,sol.t,sol[6,:], label="h", linewidth=2)
axislegend(ax, position = :rt, framevisible = false, backgroundcolor = :transparent)
display(fig)

# Steady states
ss = hc_steady_states(model5, ps)

###### BIFURCATIONS #####
p_start = ps
bif_par = :σ_m
bif_par_str = String(bif_par)
# Sensitive to initial conditions!
u_guess = [:l => ss[1][1], :p => ss[1][2], :q => ss[1][3], :m => ss[1][4], :N => ss[1][5], :h => ss[1][6]]
# u_guess = [:l => 0.67, :p => 0.9, :q => 0.16, :m => 0.82, :N => 0.22, :h => 0.9]
p_start = ps
# Span for the bifurcation parameter
p_span = (-0.1, 1.0)

for plot_var in [:l, :m, :p]
	plot_var_str = String(plot_var)
	bprob = BifurcationProblem(model5, u_guess, p_start, bif_par; plot_var)

	# Modified options to find the bifurcations
	opts_br = ContinuationPar(
		p_min = p_span[1],
		p_max = p_span[2],
		detect_fold = true,
		ds = 1e-4,
		dsmax = 1e-3,
		dsmin = 1e-6,          
		detect_bifurcation = 3, 
		n_inversion = 8,        
		max_bisection_steps = 25, 
		tol_bisection_eigenvalue = 1e-12, 
		max_steps = 20000,
		newton_options = NewtonPar(tol = 1e-12)  
	)

	# Bifurcation diagram
	bif_dia = bifurcationdiagram(bprob, PALC(), 2, opts_br; bothside = true)

	# Plotting
	col = [stb ? "#313973" : "#cd205a" for stb in bif_dia.γ.branch.stable]
	ls = [stb ? :solid : :dash for stb in bif_dia.γ.branch.stable]
	scene = Plots.plot(bif_dia; xguide = bif_par, yguide = plot_var, branchlabel = "Steady state", 
		linewidthstable = 3.0, linewidthunstable = 2.0, linestylestable = :solid,     
		linestyleunstable = :dash, color=col, markercolor = "#c59420", markersize=6,
		xlabel=L"%$(bif_par_str)", ylabel=L"%$(plot_var_str)^*",
		background_color_legend = nothing, 
		foreground_color_legend = nothing,
		legend=:left) 

	# Make nicer axes limits
	Plots.xlims!(0.0, 0.3)
	default_ylims = Plots.ylims()
	# Extract index of plotting variable
	idx = findfirst(p -> Symbolics.getname(p) == plot_var, species(model5))
	Plots.ylims!(0.0, min(default_ylims[2],bif_dia.γ.specialpoint[end-1].x[idx]*2))
	

	# To make legend
	Plots.plot!(scene, [NaN], [NaN], 
		color = "#313973",          
		linestyle = :solid, 
		linewidth = 3.0,
		label = "Stable"
	)
	Plots.plot!(scene, [NaN], [NaN], 
		color = "#cd205a",          
		linestyle = :solid, 
		linewidth = 2.0,
		label = "Unstable"
	)
	Plots.scatter!(scene, [NaN], [NaN],
		# markershape = :circle,    # match whatever shape BifurcationKit uses for bp
		markerstrokewidth = 0,
		markercolor = "#c59420", 
		markersize = 6,          
		label = "Branch point"   
	)

	# Find and blank the auto-generated branch labels
	for s in scene.series_list
		if s[:label] == "Steady state"
			s[:label] = ""
		end
	end
	for s in scene.series_list
		if s[:label] == "bp"
			s[:label] = ""
		end
	end
	Plots.plot!(scene)

	savefig("model5_"*plot_var_str*"_sigma_m.png")
	savefig("model5_"*plot_var_str*"_sigma_m.svg")
end
