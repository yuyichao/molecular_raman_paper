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

const prefix = joinpath(@__DIR__, "imgs", "atomic_loss_scaling")

figure()
errorbar(rates[:, 1], rates[:, 2], rates[:, 3], fmt="C0.")
plot(fit_gamma_a2.plotx, fit_gamma_a2.ploty, "C0")
xscale("log")
yscale("log")
grid()
xticks([3, 4, 6, 8, 12], ["3", "4", "6", "8", "12"])
yticks([1, 3, 10, 30], ["1", "3", "10", "30"])
xlabel("Tweezer Power (mW)")
ylabel("\$\\Gamma_{atom} (2\\pi\\cdot \\mathrm{Hz})\$")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
