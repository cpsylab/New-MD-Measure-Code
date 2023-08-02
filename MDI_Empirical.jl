using QuadGK,MAT, StatsBase, Plots, DataFrames, Distributions, Random, Statistics, StatsKit, GLM, HypothesisTests, StatsPlots, MixedModels, LsqFit

function read_mat_file(file_path)
    file = matopen(file_path)
    variables = names(file)
    data = Dict()

    for var in variables
        data[var] = read(file, var)
    end

    close(file)

    return data
end

# This is the function that we will try to fit
@. logistic5(data,params) = params[4]+(params[1]-params[4])/((1+(data/params[3])^params[2])^params[5])

# Randomly initialize the initial parameters to help prevent non-convergence
function p_init()
    a = rand(Uniform(0, 0.1))
    b = rand(Uniform(1, 10))
    c = rand(Uniform(0.1, 0.7))
    d = rand(Uniform(0.9, 1))
    e = rand(Uniform(0.5, 1.5))
    return [a, b, c, d, e]
end

# Tries to fit the function until it succeeds
function fit_function(data_x, data_y; model=logistic5, lower=[0.,0,0,0,0], upper=[1.,20,1,1,10], kwargs...)
    # Try again with different initial parameter values until curve_fit returns
    while true
        try
            return curve_fit(model,data_x,data_y, p_init(); lower, upper, kwargs...)
        catch
            continue
        end
    end
end

# https://github.com/mdlee/mpt4mst/tree/main/data
data = read_mat_file("Data/msttDataWithin.mat")["d"]

truth = vec(data["truthOSN"]) # 1 = old 2 = new 3 = lure
lureBin = vec(data["lureBinOSN"]) # low value = similar, high value = distinct
decision = vec(data["decisionOSN"]) # 1 = old 2 = new 3 = similar

LDI = data["LDI"]
REC = data["REC"]

oldStim = truth .== 1
newStim = truth .== 2
lureStim = truth .== 3

# Recoding lurebin for old and new cases
lureBin[oldStim] .= 0
lureBin[newStim] .= 6

# Scale lureBin to be between 0 and 6
scale(data) = (data .- minimum(data)) ./ (maximum(data) - minimum(data))
lureBin .= scale(lureBin)

curve_fit_params = Matrix(undef, 21, 5)
participant_plot = Vector(undef, 21)
participant_auc_l5 = Vector(undef, 21)
participant_auc_l5_scale = Vector(undef, 21)
participant_l5_max = Vector(undef, 21)
participant_l5_min = Vector(undef, 21)

for i in 1:21 # for each participant

    idx = (192 * (i - 1)) + 1:(192 * i)

    # Exclude NaN decisions + associated lurebin
    cases_noNaN = .!isnan.(decision[idx])

    # set decision metric (represents pNew)
    decision_metric = [Dict(1 => 0, 2 => 1, 3 => 1)[x] for x in decision[idx][cases_noNaN]]

    # Fit participant data to logistic5 curve
    fit_logistic5 = fit_function(lureBin[idx][cases_noNaN], decision_metric)
    curve_fit_params[i,:] = fit_logistic5.param


    # Get max, calculate area between max line and l5 curve
    participant_l5_max[i] = logistic5(maximum(lureBin), curve_fit_params[i,:])
    participant_l5_min[i] = logistic5(0, curve_fit_params[i,:])
    area_diff(x) = participant_l5_max[i] - logistic5(x, curve_fit_params[i,:])
    participant_auc_l5[i], = quadgk(area_diff, 0, maximum(lureBin))
    participant_auc_l5_scale[i] = participant_auc_l5[i]/(participant_l5_max[i]-participant_l5_min[i])


    # Plotting
    x_values = range(0, stop=maximum(lureBin), length=500)
    plot_values = logistic5(x_values, fit_logistic5.param)
    participant_plot[i] = scatter(lureBin[idx][cases_noNaN], decision_metric, alpha = 0.1, xlabel = "Lure Bin", ylabel  = "New/Similar response")
    plot!(participant_plot[i], x_values, plot_values;
            plot_title="Participant #$(string(i))",
            annotation=[(0.5, 0.5, "Area above curve: $(round(participant_auc_l5[i], digits=2))"),
                        (0.5, 0.4, "Scaled area above curve: $(round(participant_auc_l5_scale[i], digits=2))")],
            legend = false,
            dpi=300,
            )

    hline!(participant_plot[i], [participant_l5_max[i]], linestyle=:dash, color=:red, label="Max")
end

LDI_l5_plot = scatter(participant_auc_l5, LDI[1,:];
                        legend = false,
                        xlabel = "Area",
                        ylabel = "LDI",
                        xlim = (-0.1,1),
                        ylim = (-0.1,1),
                        dpi=300)
for i in 1:length(participant_auc_l5)
    annotate!(LDI_l5_plot, participant_auc_l5[i], LDI[1,i] + 0.03, text(i, 8, :black))
end

REC_l5_plot = scatter(participant_auc_l5, REC[1,:];
                        legend = false,
                        xlabel = "Area",
                        ylabel = "REC",
                        xlim = (-0.1,1),
                        ylim = (-0.1,1),
                        dpi=300)
for i in 1:length(participant_auc_l5)
    annotate!(REC_l5_plot, participant_auc_l5[i], REC[1,i] + 0.03, text(i, 8, :black))
end

LDI_l5_scale_plot = scatter(participant_auc_l5_scale, LDI[1,:];
                             legend = false,
                             xlabel = "Scaled Area",
                             ylabel = "LDI",
                             xlim = (-0.1,1),
                             ylim = (-0.1,1),
                             dpi=300)
for i in 1:length(participant_auc_l5_scale)
    annotate!(LDI_l5_scale_plot, participant_auc_l5_scale[i], LDI[1,i] + 0.03, text(i, 8, :black))
end

REC_l5_scale_plot = scatter(participant_auc_l5_scale, REC[1,:];
                             legend = false,
                             xlabel = "Scaled Area",
                             ylabel = "REC",
                             xlim = (-0.1,1),
                             ylim = (-0.1,1),
                             dpi=300)
for i in 1:length(participant_auc_l5_scale)
    annotate!(REC_l5_scale_plot, participant_auc_l5_scale[i], REC[1,i] + 0.03, text(i, 8, :black))
end

REC_max_plot = scatter(participant_l5_max, REC[1,:];
                        legend = false,
                        xlabel = "Curve Maximum",
                        ylabel = "REC",
                        xlim = (-0.1,1),
                        ylim = (-0.1,1),
                        dpi=300)
for i in 1:length(participant_l5_max)
    annotate!(REC_max_plot, participant_l5_max[i], REC[1,i] + 0.03, text(i, 8, :black))
end

REC_scaled_plot = scatter(participant_l5_max .- participant_l5_min, REC[1,:];
                            legend = false,
                            xlabel = "Curve Maximum - Curve Minimum",
                            ylabel = "REC",
                            xlim = (-0.1,1),
                            ylim = (-0.1,1),
                            dpi=300)
for i in 1:length(participant_l5_max)
    annotate!(REC_scaled_plot, participant_l5_max[i].-participant_l5_min[i], REC[1,i] + 0.03, text(i, 8, :black))
end

# Save all the figures
savefig(LDI_l5_plot, "MDI_Figures/LDI_l5_plot.png")
savefig(LDI_l5_scale_plot, "MDI_Figures/LDI_l5_scale_plot.png")
savefig(REC_l5_plot, "MDI_Figures/REC_l5_plot.png")
savefig(REC_l5_scale_plot, "MDI_Figures/REC_l5_scale_plot.png")
savefig(REC_max_plot, "MDI_Figures/REC_max_plot.png")
savefig(REC_scaled_plot, "MDI_Figures/REC_scaled_plot.png")


for i in eachindex(participant_plot)
    savefig(participant_plot[i], "MDI_Figures/participant_plot_$i.png")
end

participant_plot[2]
