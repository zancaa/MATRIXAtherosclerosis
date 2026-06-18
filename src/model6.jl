using Catalyst, Latexify, DifferentialEquations, BifurcationKit, HomotopyContinuation, Plots, CairoMakie, ColorSchemes 

# Load parameters and plotting things
include("parameters_model6.jl")
include("plotting.jl")

# Define 'reaction network'
model6 = @reaction_network begin
	# modLDL
	σ_l, ∅ => l
	mm(l,ν_l*m,ξ_l), l => ∅
	d_l, l --> ∅
	# p cytokine
	mm(l,σ_p,β_p), ∅ => p
	d_p, p --> ∅
	mm(l,ν_p*m,ξ_l), ∅ => p
	# macrophage
	σ_m*(1+α*(m*l/(ξ_l + l)))*(p), ∅ => m
	mm(l,ν_m*m,ξ_l), m => ∅
	d_m, m --> ∅
end

# Set initial conditions, time span and parameter values
# To find multiple branches in bifurcation diagram, need to use multiple initial conditions
# u0 = [:l => 2.5, :p => -6.0, :m => -2.1] # weird stuff!!
u0_br2 = [:l => 1.5, :p => -10.0, :m => -10.0] # other branch for default parameters - but this is unstable, numerics won't find it (without some faffing)
u0_br3 = [:l => 5/4, :p => 1/18, :m => 0] # another branch
# u0 = [:l => 0.003, :p => 9e7, :m => 2e8]
u0 = [:l => 0.0, :p => 1.0, :m => 1.0]
tspan = (0.0, 1000.0)
ps = [ # decay rates
	 :d_m => d_m, :d_l => d_l,
	 :d_p => d_p,
	 # Hill scaling factor
	 :ν_l => ν_l, :ν_m => ν_m, :ν_p => ν_p,
	 # recruitment rates
	 :σ_m => σ_m, :σ_l => σ_l, :σ_p => σ_p, 
	 # saturation parameters
     :ξ_l => ξ_l, :β_p => β_p,  
	 # other parameters
	 :α => α
	 ]

# Define ODE problem
oprob = ODEProblem(model6, u0, tspan, ps)
# Solve
sol = DifferentialEquations.solve(oprob)
fig = Figure(fontsize=25)
ax = Axis(fig[1,1])
lines!(ax,sol.t,sol[1,:], label="l")
lines!(ax,sol.t,sol[2,:], label="p")
lines!(ax,sol.t,sol[3,:], label="m")
axislegend(ax, position = :rt, framevisible = false, backgroundcolor = :transparent)
display(fig)

# Steady states
ss = hc_steady_states(model6, ps)

###### BIFURCATIONS #####
#### Analytic results (to overlay numerics)) ####
## To plot sigma_m versus l* 
numerat_l(l_star) = ((d_m/ν_l)*((ξ_l + l_star)/l_star) + ν_m/ν_l) * (σ_l - d_l*l_star)
denom_l(l_star) = (1/d_p) * ((σ_p*l_star)/(β_p + l_star) + (ν_p/ν_l)*(σ_l - d_l*l_star)) * (1 + (α/ν_l)*(σ_l - d_l*l_star))
sigma_m_l(l_star) = numerat_l(l_star)/denom_l(l_star)
## To plot sigma_m versus p* 
p_star(l_star) = (1/d_p) * ((σ_p * l_star)/(β_p + l_star) + (ν_p/ν_l)*(σ_l - d_l*l_star))
## To plot sigma_m versus m* 
m_star(l_star) = (ξ_l + l_star)/(ν_l * l_star) * (σ_l - d_l * l_star)

l_star_vals = 0.0:0.01:2.5

# The below code plots 'nicer' curves than the OG
# Function to extract line segments (otherwise get weird behaviours
# -- i.e. random bits of lines that aren't coloured consistently
function get_contiguous_segments(br, want_stable::Bool)
    segments_p = Vector{Vector{Float64}}()
    segments_x = Vector{Vector{Float64}}()
    
    current_p = Float64[]
    current_x = Float64[]
    
    for pt in br.branch
        if pt.stable == want_stable
            push!(current_p, pt.param)
            push!(current_x, pt.x[1]) 
        else
            if !isempty(current_p)
                push!(segments_p, current_p)
                push!(segments_x, current_x)
                current_p = Float64[]
                current_x = Float64[]
            end
        end
    end
    
    if !isempty(current_p)
        push!(segments_p, current_p)
        push!(segments_x, current_x)
    end
    
    return segments_p, segments_x
end

# Plots!
σ_l_vals = [σ_l] # [0.9, 1.2]
σ_l_ps_idx = findfirst(p -> first(p) == :σ_l, ps)
par_order = Catalyst.parameters(model6)

bif_par = :σ_m
bif_par_str = String(bif_par)
u_guess = [:l => ss[1][1], :p => ss[1][2], :m => ss[1][3]] 
# u_guess = [:l => 0.67, :p => 0.9, :m => 0.82]
p_start = ps
# Span for the bifurcation parameter
p_span = (0.5, 1.5)
# Change some options to make sure bifurcations are 'seen'
opts_br = ContinuationPar(p_min = p_span[1], p_max = p_span[2], detect_fold = true, 
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

plot_var_vals = [:l, :m, :p]
plot_vec = Vector{Plots.Plot}(undef, length(plot_var_vals))

for (plt_idx, plot_var) in enumerate(plot_var_vals) 

	plot_var_str = String(plot_var)

	for (σ_l_idx, σ_l) in enumerate(σ_l_vals) # used a range of sigma_l's in places
		p_start[σ_l_ps_idx] = :σ_l => σ_l
		if length(σ_l_vals) == 1
			plot_vec[plt_idx] = Plots.plot()
			for (u_g_idx, u_g) in enumerate([u_guess, u0_br2, u0_br3])
				bprob = BifurcationProblem(model6, u_g, p_start, bif_par; plot_var)

				# global branch1 = continuation(bprob, PALC(), opts_br; verbosity = 1, bothside = true)
				bif_dia = bifurcationdiagram(bprob, PALC(), 2, opts_br; bothside = true)
			
				col = [stb ? "#313973" : "#cd205a" for stb in bif_dia.γ.branch.stable]
				ls = [stb ? :solid : :dash for stb in bif_dia.γ.branch.stable]

				# Plot each branch
				for branch in [bif_dia.γ]  # add child branches as needed
					segs_p_stable, segs_x_stable = get_contiguous_segments(branch, true)
					segs_p_unstable, segs_x_unstable = get_contiguous_segments(branch, false)
					
					for (i, (p_seg, x_seg)) in enumerate(zip(segs_p_stable, segs_x_stable))
						Plots.plot!(plot_vec[plt_idx], p_seg, x_seg,
							color = "#313973",
							linewidth = 3,
							linestyle = :solid,
							label = (i == 1 && u_g_idx == 1) ? "Stable" : "",
							background_color_legend = nothing, 
							foreground_color_legend = nothing,
							xlabel=L"%$(bif_par_str)", ylabel=L"%$(plot_var_str)^*",
						)
					end
					
					bp_p = Float64[]
					bp_x = Float64[]
					for (i, (p_seg, x_seg)) in enumerate(zip(segs_p_unstable, segs_x_unstable))
						Plots.plot!(plot_vec[plt_idx], p_seg, x_seg,
							color = "#cd205a", 
							linewidth = 3,
							linestyle = :solid,
							label = (i == 1 && u_g_idx == 1) ? "Unstable" : "",
						)
					end

					# Extract and plot branch points
					for pt in branch.specialpoint
						if pt.type == :bp
							push!(bp_p, pt.param)
							# handle both scalar and vector x
							if pt.x isa AbstractVector
								push!(bp_x, pt.printsol.x[1])
							else
								push!(bp_x, pt.printsol.x)
							end
						end
					end
					
					Plots.scatter!(plot_vec[plt_idx], bp_p, bp_x,
						markershape = :circle,
						markercolor = "#c59420",
						markerstrokewidth = 0,
						markersize = 6,
						label = (u_g_idx == 1) ? "Branch point" : "",
						background_color_legend = nothing, # Transparent background
						foreground_color_legend = nothing,
						xlabel=L"%$(bif_par_str)", ylabel=L"%$(plot_var_str)^*",
					)
				end

				if u_g_idx == 1
					# Make nicer axes limits
					Plots.xlims!(0.6, 1.0)
					default_ylims = Plots.ylims()
					# Extract index of plotting variable
					idx = findfirst(p -> Symbolics.getname(p) == plot_var, species(model6))
					Plots.ylims!(0.0, min(default_ylims[2],bif_dia.γ.specialpoint[end-1].x[idx]*2))
				end
			end
			# Overlay analytic solution (could have these the other way around)
			if plot_var == :l
				Plots.plot!(
					plot_vec[plt_idx],
					sigma_m_l.(l_star_vals),
					l_star_vals,
					lw = 2.5,
					color = "#d5acee",
					linestyle = :dot,
					label = "Analytic solution"
				)
			elseif plot_var == :m
				Plots.plot!(
					plot_vec[plt_idx],
					sigma_m_l.(l_star_vals),
					m_star.(l_star_vals),
					lw = 2.5,
					color = "#d5acee",
					linestyle = :dot,
					label = "Analytic solution"
				)
			elseif plot_var == :p
				Plots.plot!(
					plot_vec[plt_idx],
					sigma_m_l.(l_star_vals),
					p_star.(l_star_vals),
					lw = 2.5,
					color = "#d5acee",
					linestyle = :dot,
					label = "Analytic solution"
				)
			end

			# Save!
			savefig("model6_"*plot_var_str*"_sigma_m_sigma_l_"*replace(string(σ_l), "." => "p")*".png")
			savefig("model6_"*plot_var_str*"_sigma_m_sigma_l_"*replace(string(σ_l), "." => "p")*".svg")

			# Uncomment below if only want legend in first figure
			# if plt_idx > 1
			# 	plot!([],[], legend=false)
			# end
		elseif length(σ_l_vals) > 1
			bprob = BifurcationProblem(model6, u_guess, p_start, bif_par; plot_var)
			# global branch1 = continuation(bprob, PALC(), opts_br; verbosity = 1, bothside = true)
			bif_dia = bifurcationdiagram(bprob, PALC(), 2, opts_br; bothside = true)
			# Below code accounts for plotting when there are multiple branches (occurs for some values of σ_l)
			# To be used when more than one σ_l is being plotted simultaneously
			# Plot each branch
			if σ_l_idx == 1
				# Initialise plot
				global fig = Plots.plot()
				# Sample distinct colors
				n_colors = length(σ_l_vals)
				global color_array = [ColorSchemes.Accent_5[z] for z in range(0, stop=1, length=n_colors)]
			end
			for branch in [bif_dia.γ]  
				segs_p_stable, segs_x_stable = get_contiguous_segments(branch, true)
				segs_p_unstable, segs_x_unstable = get_contiguous_segments(branch, false)
				
				for (i, (p_seg, x_seg)) in enumerate(zip(segs_p_stable, segs_x_stable))
					Plots.plot!(fig, p_seg, x_seg,
						color = color_array[σ_l_idx],
						linewidth = 3,
						linestyle = :solid,
						label = i == 1 ? L"\sigma_l = %$σ_l" : "", 
						background_color_legend = nothing, 
						foreground_color_legend = nothing,
						xlabel=L"%$(bif_par_str)", ylabel=L"%$(plot_var_str)^*",
					)
				end
				
				bp_p = Float64[]
				bp_x = Float64[]
				for (i, (p_seg, x_seg)) in enumerate(zip(segs_p_unstable, segs_x_unstable))
					Plots.plot!(fig, p_seg, x_seg,
						color = color_array[σ_l_idx],
						linewidth = 2,
						linestyle = :dot,
						label = "" 
					)
				end

			end

			for s in fig.series_list
				if s[:label] == "Branch point"
					s[:label] = ""
				end
			end
			if σ_l_idx == length(σ_l_vals)
				savefig(fig,"model6_"*plot_var_str*"_vs_sigma_m_sigma_l_range.png")
				savefig(fig,"model6_"*plot_var_str*"_vs_sigma_m_sigma_l_range.svg")
			end
		end
	end
end

if length(σ_l_vals) == 1 
	# This plotting is currently hideous, will need to update...
	fig = Plots.plot(plot_vec..., layout=(length(plot_vec),1))
	savefig(fig,"model6_sigma_m_sigma_l_"*replace(string(σ_l), "." => "p")*".png")
	savefig(fig,"model6_sigma_m_sigma_l_"*replace(string(σ_l), "." => "p")*".svg")
elseif length(σ_l_vals) > 1 
	fig = Plots.plot(plot_vec..., layout=(1,length(plot_vec)))
	savefig(fig,"model6_sigma_m_sigma_l_range.png")
	savefig(fig,"model6_sigma_m_sigma_l_range.svg")
end

####### Codimension-2 (plotting parameters against each other) #######
# To be able to do this, things need to be initialised to actually find 
# branch points before the cusp (in this case) can be tracked
idx1 = findfirst(p -> Symbol(p) == bif_par, par_order)
lens1 = @optic _[idx1]
# Would need to change the below line if using different parameters
# Note: the index in the p_start vector does not match par_order (could make them match)
σ_m_ps_idx = findfirst(p -> first(p) == :σ_m, ps)
p_start[σ_m_ps_idx] = :σ_m => 0.8
plot_var = :l
plot_var_str = String(plot_var)
# Tomfoolery with finding the right starting values to find the branch point
ss = hc_steady_states(model6, p_start)
u_guess = [:l => ss[2][1], :p => ss[2][2], :m => ss[2][3]] 
bprob = BifurcationProblem(model6, u_guess, p_start, bif_par; plot_var)
bprob_codim2 = @set bprob.lens = lens1
# Make span small enough to find the branch points
opts_br = ContinuationPar(p_min = 0.5, p_max = 1.1, detect_fold = true, 
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
br = BifurcationKit.continuation(bprob_codim2, PALC(), opts_br)

# Find branch point
ind = findfirst(x -> x.type == :bp, br.specialpoint)
# Choose second parameter 
second_par = :σ_l

# Indexing
idx2 = findfirst(p -> Symbol(p) == second_par, par_order)
# Would need to change the below line if using different parameters
# Note: the index in the p_start vector does not match par_order (could make them match)
p_start[σ_l_ps_idx] = :σ_l => 1.2
lens2 = @optic _[idx2]

# Codimension-2 continuation
sn_codim2 = BifurcationKit.continuation(br, ind, lens2,
    BifurcationKit.ContinuationPar(opts_br,
        p_min = 0.0, p_max = 2.0,
        ds = 1e-4, dsmax = 1e-3,
        detect_bifurcation = 3,
        n_inversion = 8,
        max_steps = 20000);
    normC = norminf,
    bothside = true,
    detect_codim2_bifurcation = 2,
)

scene = Plots.plot(sn_codim2, vars = (:p1, :p9), branchlabel = "Fold", 
	xlabel=L"\sigma_l", ylabel=L"\sigma_m", linewidth=3, markersize=6, 
	markercolor = "#eec54a", color = "#8b73bd", background_color_legend = nothing, # Transparent background
    foreground_color_legend = nothing)
scene.series_list[3][:label] = "Cusp"
Plots.plot!(scene)
savefig("model6_cusp.png")
savefig("model6_cusp.svg")
