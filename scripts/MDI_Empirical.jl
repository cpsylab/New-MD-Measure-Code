using DataFrames
using MAT
using StatsPlots
using CSV
using MDI

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

# https://github.com/mdlee/mpt4mst/tree/main/data
data = read_mat_file("./Data/msttDataWithin.mat")["d"]

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

# Scale lureBin to be between 0 and 1
scale(data) = (data .- minimum(data)) ./ (maximum(data) - minimum(data))
lureBin .= scale(lureBin)

curve_fit_params = Matrix(undef, 21, 5)
participant_plot = Vector(undef, 21)
participant_barplot = Vector(undef, 21)
participant_pairplot = Vector(undef, 21)
participant_auc_l5 = Vector(undef, 21)
participant_auc_l5_scale = Vector(undef, 21)
participant_l5_max = Vector(undef, 21)
participant_l5_min = Vector(undef, 21)

seed = 1234

#= Code below verifies the provided LDI and REC measures are correct
i = 4
idx = (192 * (i - 1)) + 1:(192 * i)
decision[idx]
cases_noNaN = .!isnan.(decision[idx])

lure_decisions = (decision[idx] .== 3) .& .!ismissing.(decision[idx])
old_decisions = (decision[idx] .== 1) .& .!ismissing.(decision[idx])

p_lure_lure = sum(lure_decisions[truth[idx] .== 3]) / length(lure_decisions[truth[idx] .== 3])
p_lure_foil = sum(lure_decisions[truth[idx] .== 2]) / length(lure_decisions[truth[idx] .== 2])

p_old_old = sum(old_decisions[truth[idx] .== 1]) / length(old_decisions[truth[idx] .== 1])
p_old_foil = sum(old_decisions[truth[idx] .== 2]) / length(old_decisions[truth[idx] .== 2])

p_lure_lure - p_lure_foil
LDI[4]

p_old_old - p_old_foil
REC[4]
=#

for i in 1:21 # for each participant

    idx = (192 * (i - 1)) + 1:(192 * i)

    # Exclude NaN decisions + associated lurebin
    cases_noNaN = .!isnan.(decision[idx])

    # set decision metric (represents pNew)
    decision_metric = [Dict(1 => 0, 2 => 1, 3 => 1)[x] for x in decision[idx][cases_noNaN]]

    # Fit participant data to logistic5 curve
    fit_logistic5 = fit_model(lureBin[idx][cases_noNaN], decision_metric, seed = seed)
    curve_fit_params[i,:] = fit_logistic5.param


    # Get max, calculate area between max line and l5 curve
    # participant_l5_max[i] = logistic5(maximum(lureBin), curve_fit_params[i,:])
    # participant_l5_min[i] = logistic5(0, curve_fit_params[i,:])
    # area_diff(x) = participant_l5_max[i] - logistic5(x, curve_fit_params[i,:])
    # participant_auc_l5[i], = quadgk(area_diff, 0, maximum(lureBin))
    # participant_auc_l5_scale[i] = participant_auc_l5[i]/(participant_l5_max[i]-participant_l5_min[i])
    participant_auc_l5[i], participant_auc_l5_scale[i], participant_l5_min[i], participant_l5_max[i] = get_aucs(curve_fit_params[i,:])

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

    unique_lureBins = unique(lureBin[idx])
    stacked_data = zeros(Float64, length(unique_lureBins), 4) # Including an extra column for NaNs

    for j in 1:length(unique_lureBins)
        bin = unique_lureBins[j]
        for k in idx
            if lureBin[k] == bin
                for decision_value in 1:3
                    if decision[k] == decision_value
                        stacked_data[j, decision_value] += 1
                    end
                end
                if isnan(decision[k])
                    stacked_data[j, 4] += 1 # Counting NaNs in the fourth column
                end
            end
        end
    end
    participant_barplot[i] = groupedbar(unique_lureBins, stacked_data,bar_position = :stack, bar_width=0.1, labels=["Old" "New" "Similar" "NaN"], xlabel="Lure Bin", ylabel="Count", title="Participant #$(string(i))", legend = :outertopright)

    stacked_pairs = zeros(Int64, 3,3)
    for j in 1:length(truth[idx][cases_noNaN])
        stacked_pairs[Int(truth[idx][cases_noNaN][j]), Int(decision[idx][cases_noNaN][j])] += 1
    end

    participant_pairplot[i] = heatmap(stacked_pairs', xlabel = "Truth", ylabel = "Decision", xticks = ([1:1:3;],["Old", "New", "Similar"]), yticks = ([1:1:3;],["Old", "New", "Similar"]))
end

participant_barplot[1]
participant_plot[15]
participant_pairplot[15] # I think this doesnt work properly

df = DataFrame(participant_auc_l5 = participant_auc_l5,
                participant_auc_l5_scale = participant_auc_l5_scale,
                participant_l5_min = participant_l5_min,
                participant_l5_max = participant_l5_max,
                REC = REC[1,:],
                LDI = LDI[1,:])

#CSV.write("data_L5.csv", df)

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
#=
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
=#
participant_plot[1]

participant_barplot[1]

LDI_l5_scale_plot

