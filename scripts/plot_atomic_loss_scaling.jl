#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival, get_plot_range
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const rates = [15 32.5 3.1
               6 2.77 0.27
               3 0.573 0.078]

function gen_power_model(pwr)
    return function (x, p)
        return p[1] .* x.^pwr
    end
end
fit_gamma_a2 = fit_data(gen_power_model(2.58), rates[:, 1], rates[:, 2], rates[:, 3],
                        [1.0]; plot_lo=2.7)
@show fit_gamma_a2.uncs

const prefix = joinpath(@__DIR__, "../imgs", "atomic_loss_scaling")

# figure()
# errorbar(rates[:, 1], rates[:, 2], rates[:, 3], fmt="C0o")
# plot(fit_gamma_a2.plotx, fit_gamma_a2.ploty, "C0")
# xscale("log")
# yscale("log")
# grid()
# xticks([3, 4, 6, 8, 12], ["3", "4", "6", "8", "12"])
# yticks([1, 3, 10, 30], ["1", "3", "10", "30"])
# ax = gca()
# ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.set_xticks(3:16, minor=true)
# xlabel("Tweezer Power (mW)")
# ylabel("\$\\Gamma_{\\mathrm{atom}} (2\\pi\\times \\mathrm{Hz})\$")
# NaCsPlot.maybe_save("$(prefix)")

figure(figsize=[4.8, 3.6 * 0.7] .* 1.08792)
ax = gca()
bg = matplotlib.patches.Rectangle((0.625, 0.31), 3.625, 39.69, facecolor="white", alpha=0.9)
ax.add_patch(bg)
errorbar(rates[:, 1] ./ 4, rates[:, 2], rates[:, 3], fmt="C0o")
plot(fit_gamma_a2.plotx ./ 4, fit_gamma_a2.ploty, "C0")
xscale("log")
yscale("log")
grid()
xlim([0.625, 4.25])
ylim([0.31, 40])
xticks([1, 2, 4], ["1", "2", "4"])
yticks([1, 3, 10, 30], ["1", "3", "10", "30"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks(0.7:0.1:4, minor=true)
xlabel("Tweezer Power (mW)", fontsize=17)
ylabel("\$\\Gamma_{\\mathrm{atom}} (2\\pi\\times \\mathrm{Hz})\$", fontsize=17)
tight_layout()
NaCsPlot.maybe_save("$(prefix)_inset")

NaCsPlot.maybe_show()
