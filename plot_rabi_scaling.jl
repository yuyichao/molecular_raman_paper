#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival, get_plot_range
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const freqs = [560 15 3.282 0.042
               560 6 1.182 0.033
               560 3 0.488 0.015

               503 15 2.073 0.051
               503 9 1.074 0.036
               503 6 0.680 0.022
               503 3 0.2544 0.0077

               605 15 6.39 0.15
               605 6 1.985 0.043
               605 3 0.845 0.042

               625 15 8.35 0.25
               625 6 2.513 0.036
               625 3 0.996 0.029
               ]

function model_omega(x, p)
    return p[1] .* x.^1.29
end

broacast_array(f, x) = f(x)
broacast_array(f, x::AbstractArray) = f.(x)

# x: freq, power
# p: offset, strength
function gen_model(fpa0)
    function model(xs, p)
        function real_model(x)
            freq, power = x
            offset, strength = p
            p_omega = (offset - strength / (freq - fpa0),)
            return model_omega(power, p_omega)
        end
        return broacast_array(real_model, xs)
    end
end

function gen_data(data_in)
    freq_ids = Dict{Float64,Int}()
    xs = Tuple{Float64,Float64,Int}[]
    ys = Float64[]
    uncs = Float64[]
    idmax = 0
    for i in 1:size(data_in, 1)
        freq = data_in[i, 1]
        power = data_in[i, 2]
        freq_id = get!(freq_ids, freq) do
            idmax += 1
            return idmax
        end
        push!(xs, (freq, power, freq_id))
        push!(ys, data_in[i, 3])
        push!(uncs, data_in[i, 4])
    end
    return (x=xs, y=ys, unc=uncs, freqs=collect(keys(freq_ids)))
end

# function fit_freq_model(model, data, _p0)
#     p0 = [_p0; zeros(length(data.freqs))]
#     return fit_data(model, data.x, data.y, data.unc, p0, plotx=false)
# end

function get_plot_data_freq(data, fit, model, freq)
    local freq_id
    xs = Float64[]
    ys = Float64[]
    uncs = Float64[]
    for i in 1:length(data.x)
        x = data.x[i]
        x[1] == freq || continue
        if !@isdefined(freq_id)
            freq_id = x[3]
        else
            @assert(freq_id == x[3])
        end
        push!(xs, x[2])
        push!(ys, data.y[i])
        push!(uncs, data.unc[i])
    end
    function plot_func(x)
        return model((freq, x, freq_id), fit.param)
    end
    plotx = get_plot_range(xs)
    return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
end

function get_plot_data_power(data, fit, model, power)
    xs = Float64[]
    xids = Float64[]
    ys = Float64[]
    uncs = Float64[]
    for i in 1:length(data.x)
        x = data.x[i]
        x[2] == power || continue
        push!(xs, x[1])
        push!(ys, data.y[i])
        push!(uncs, data.unc[i])
    end
    function plot_func(x)
        return model((x, power,), fit.param)
    end
    plotx = get_plot_range(data.freqs)
    return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
end

const model = gen_model(705)
const data = gen_data(freqs)
const fit1 = fit_data(model, data.x, data.y, [0.1, 0.1], plotx=false)
@show fit1.uncs

const prefix = joinpath(@__DIR__, "imgs", "rabi_scaling")

figure(figsize=[7.2, 5.4])
const powers_plot = [15, 6, 3]
# const plot_freqs = get_plot_range(data.freqs)
for i in 1:length(powers_plot)
    power = powers_plot[i]
    pd = get_plot_data_power(data, fit1, model, power)
    plot(pd.plotx .- 703.6, pd.ploty, "C$(i - 1)")
    errorbar(pd.x .- 703.6, pd.y, pd.unc, fmt="C$(i - 1)o", label="$(power) mW")
end
# xticks([288500, 288530, 288560, 288590, 288620])
# p: offset, strength
# text(500, 5.0, ("\$\\left(a-\\dfrac{b}{f-705 GHz}\\right)\\cdot P^{1.29}\$"))
# text(500, 2.8, ("\$a=$(fit1.uncs[1] * 1000)\$ Hz/mW\$^{1.29}\$\n" *
#               "\$b=$(fit1.uncs[2])\$ kHz\$\\cdot\$GHz/mW\$^{1.29}\$"), fontsize="small")
legend(loc=(0.74, 0.32), fontsize="small", handlelength=0.6, handletextpad=0.3)
grid()
yticks([0, 4, 8])
xticks([-200, -160, -120, -80])
ylim([0, 9.7])
xlabel("Detuning (GHz)")
ylabel("\$\\Omega_{R} (2\\pi\\cdot \\mathrm{kHz})\$")
NaCsPlot.maybe_save("$(prefix)")

figure(figsize=[3.6, 3.0])
freq = 560
plot2 = get_plot_data_freq(data, fit1, model, freq)
bg = matplotlib.patches.Rectangle((2.8, 0.42), 10.2, 3.78, facecolor="white", alpha=0.9)
gca().add_patch(bg)
bg = matplotlib.patches.Rectangle((13, 0.6), 4, 3.6, facecolor="white", alpha=0.9)
gca().add_patch(bg)
# errorbar(plot2.x, plot2.y, plot2.unc, fmt="C0o")
errorbar(plot2.x[1], plot2.y[1], plot2.unc[1], fmt="C0o")
errorbar(plot2.x[2], plot2.y[2], plot2.unc[2], fmt="C1o")
errorbar(plot2.x[3], plot2.y[3], plot2.unc[3], fmt="C2o")
plot(plot2.plotx, plot2.ploty, "k")
xscale("log")
yscale("log")
minorticks_off()
xticks([])
yticks([1, 4], ["1", "4"])
xlim([2.8, 17])
ylim([0.42, 4.2])
ax = gca()
ax.spines["top"].set_visible(false)
ax.spines["bottom"].set_visible(false)
ax.spines["right"].set_visible(false)
NaCsPlot.maybe_save("$(prefix)_560")

NaCsPlot.maybe_show()
