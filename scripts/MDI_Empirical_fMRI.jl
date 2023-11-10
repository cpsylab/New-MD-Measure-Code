using CSV
using DataFrames
using StatsPlots
using MDI

# https://osf.io/v7cb4/
df = DataFrame(CSV.File("./Data/mst_fmri_test.csv"))

count = 0
for i in 1:72
    if length(unique(df[(108 * (i - 1)) + 1:(108 * i), "Subject"])) == 1
        global count += 1
    end
end
count # each subject has exactly 108 rows

decision = vec(df[!, "TestObj.RESP"]) # v = old, b = similar, n = new, na = missing
truth = vec(df[!, "ItemType"])
lureBin = vec(df[!, "LureBin"])

#=
truth = vec(data["truthOSN"]) # 1 = old 2 = new 3 = lure
lureBin = vec(data["lureBinOSN"]) # low value = similar, high value = distinct
decision = vec(data["decisionOSN"]) # 1 = old 2 = new 3 = similar

LDI = data["LDI"]
REC = data["REC"]
=#

oldStim = truth .== "Repeated"
newStim = truth .== "Foil"
lureStim = truth .== "Lure"

unique(lureBin)

# Recoding lurebin for old and new cases
lureBin[oldStim] .= 0
lureBin[newStim] .= 6

# Scale lureBin to be between 0 and 1
scale(data) = (data .- minimum(data)) ./ (maximum(data) - minimum(data))
lureBin = scale(lureBin)

num_participants = 72

curve_fit_params = Matrix(undef, num_participants, 5)
participant_plot = Vector(undef, num_participants)
participant_barplot = Vector(undef, num_participants)
participant_pairplot = Vector(undef, num_participants)
participant_auc_l5 = Vector(undef, num_participants)
participant_auc_l5_scale = Vector(undef, num_participants)
participant_l5_max = Vector(undef, num_participants)
participant_l5_min = Vector(undef, num_participants)
LDI = Vector(undef, num_participants)
REC = Vector(undef, num_participants)

decision_recoded = [get(Dict("v" => 1, "b" => 3, "n" => 2), x, x) for x in decision]

#= diagnostics code
i = 19
idx = (108 * (i - 1)) + 1:(108 * i)
decision[idx]
cases_noNaN = decision[idx] .!== missing

decision[idx][cases_noNaN]


length(lure_decisions[truth[idx] .== "Foil"])

unique_lureBins = unique(lureBin[idx])
stacked_data = zeros(Float64, length(unique_lureBins), 4) # Including an extra column for NaNs

lure_decisions = (decision[idx] .== "b") .& .!ismissing.(decision[idx])
old_decisions = (decision[idx] .== "v") .& .!ismissing.(decision[idx])

p_lure_lure = sum(lure_decisions[truth[idx] .== "Lure"]) / length(lure_decisions[truth[idx] .== "Lure"])
p_lure_foil = sum(lure_decisions[truth[idx] .== "Foil"]) / length(lure_decisions[truth[idx] .== "Foil"])

p_old_old = sum(old_decisions[truth[idx] .== "Repeated"]) / length(old_decisions[truth[idx] .== "Repeated"])
p_old_foil = sum(old_decisions[truth[idx] .== "Foil"]) / length(old_decisions[truth[idx] .== "Foil"])

p_lure_lure - p_lure_foil
p_old_old - p_old_foil
=#


seed = 1234


for i in 1:num_participants # for each participant

    idx = (108 * (i - 1)) + 1:(108 * i)

    # Exclude NaN decisions + associated lurebin
    cases_noNaN = decision[idx] .!== missing

    # set decision metric (represents pNew)
    decision_metric = [Dict("v" => 0, "b" => 1, "n" => 1)[x] for x in decision[idx][cases_noNaN]]

    # Fit participant data to logistic5 curve
    fit_logistic5 = fit_model(lureBin[idx][cases_noNaN], decision_metric, seed = seed)
    curve_fit_params[i,:] = fit_logistic5.param


    # Get max, calculate area between max line and l5 curve
    participant_auc_l5[i], participant_auc_l5_scale[i], participant_l5_min[i], participant_l5_max[i] = get_aucs(curve_fit_params[i,:])

    lure_decisions = (decision[idx] .== "b") .& .!ismissing.(decision[idx])
    old_decisions = (decision[idx] .== "v") .& .!ismissing.(decision[idx])

    p_lure_lure = sum(lure_decisions[truth[idx] .== "Lure"]) / length(lure_decisions[truth[idx] .== "Lure"])
    p_lure_foil = sum(lure_decisions[truth[idx] .== "Foil"]) / length(lure_decisions[truth[idx] .== "Foil"])

    p_old_old = sum(old_decisions[truth[idx] .== "Repeated"]) / length(old_decisions[truth[idx] .== "Repeated"])
    p_old_foil = sum(old_decisions[truth[idx] .== "Foil"]) / length(old_decisions[truth[idx] .== "Foil"])

    LDI[i] = p_lure_lure - p_lure_foil
    REC[i] = p_old_old - p_old_foil

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
                    if decision_recoded[k] === decision_value
                        stacked_data[j, decision_value] += 1
                    end
                end
                if decision_recoded[k] .=== missing
                    stacked_data[j, 4] += 1 # Counting NaNs in the fourth column
                end
            end
        end
    end
    participant_barplot[i] = groupedbar(unique_lureBins, stacked_data,bar_position = :stack, bar_width=0.1, labels=["Old" "New" "Similar" "NaN"], xlabel="Lure Bin", ylabel="Count", title="Participant #$(string(i))", legend = :outertopright)
end

df = DataFrame(participant_auc_l5 = participant_auc_l5,
                participant_auc_l5_scale = participant_auc_l5_scale,
                participant_l5_min = participant_l5_min,
                participant_l5_max = participant_l5_max,
                REC = REC[:],
                LDI = LDI[:])

#CSV.write("data_L5_fMRI.csv", df)

# Save all the figures
#=

for i in eachindex(participant_plot)
    savefig(participant_plot[i], "MDI_Figures_Wahlheim/L5_plots/participant_plot_$i.png")
end

for i in eachindex(participant_plot)
    savefig(participant_barplot[i], "MDI_Figures_Wahlheim/Bar_plots/participant_barplot_$i.png")
end
=#
