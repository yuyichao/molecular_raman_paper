#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
# using DataStructures
# using LsqFit
using LibArchive

using DelimitedFiles

const fname = joinpath(@__DIR__, "../damop-2020/20MHz_linewidth_c3SigmaOnly_3322_80kHzConfinement.csv.zst")
const data = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    readdlm(reader, ',', skipstart=1)
end
const prefix = joinpath(@__DIR__, "imgs", "raman")

fig = figure(figsize=[1, 1] * 4.8)

ax0 = fig.add_subplot(111)    # The big subplot
ax0.spines["top"].set_color("none")
ax0.spines["bottom"].set_color("none")
ax0.spines["left"].set_color("none")
ax0.spines["right"].set_color("none")
ax0.tick_params(labelcolor="w", top="off", bottom="off", left="off", right="off")
ylabel("\$\\Omega_{Raman}\$ and \$\\Gamma_{scatter}\$ (\$2\\pi\\cdot \\mathrm{kHz})\$", fontsize="small")
ax0.get_yaxis().set_label_coords(-0.13, 0.5)
ax0.tick_params(axis="x", length=0)
ax0.tick_params(axis="y", length=0)
tax0 = ax0.twinx()
tax0.tick_params(labelcolor="w", top="off", bottom="off", left="off", right="off")
tax0.tick_params(axis="x", length=0)
tax0.tick_params(axis="y", length=0)
ylabel("\$\\Omega_{Raman}/\\Gamma_{scatter}\$", fontsize="small", color="C2")

ax1 = fig.add_subplot(211)
plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ 2π / 1000), "C0",
     label="\$\\Omega_{Raman}\$")
plot(data[:, 1] .- 288625.081, abs.(data[:, 4] ./ 2π / 1000), "C1",
     label="\$\\Gamma_{scatter}\$")
text(-32, 2.5, "v=0", fontsize="small")
ylim([0, 30])
grid()
tax1 = ax1.twinx()
tax1.plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ data[:, 4]), "C2",
          label="\$\\frac{\\Omega_{Raman}}{\\Gamma_{scatter}}\$")
xlim([-35, 35])
ylim([0, 60])
tax1.tick_params(axis="y", labelcolor="C2")
setp(ax1.get_xticklabels(), visible=false)
setp(ax1.get_xticklabels(), visible=false)
ax1.tick_params(axis="x", length=0)

ax2 = fig.add_subplot(212)
subplots_adjust(hspace=0.0)
l1 = plot(data[:, 1] .- 351271.53, abs.(data[:, 2] ./ 2π / 1000), "C0",
          label="\$\\Omega_{Raman}\$")
l2 = plot(data[:, 1] .- 351271.53, abs.(data[:, 4] ./ 2π / 1000), "C1",
          label="\$\\Gamma_{scatter}\$")
text(-32, 130, "v=63", fontsize="small")
xlabel("One-Photon Detuning (GHz)")
ylim([0, 596])
grid()
ax = gca()
tax2 = ax.twinx()
l3 = tax2.plot(data[:, 1] .- 351271.53, abs.(data[:, 2] ./ data[:, 4]), "C2",
              label="\$\\frac{\\Omega_{Raman}}{\\Gamma_{scatter}}\$")
ls = [l1; l2; l3]
legend(ls, [l.get_label() for l in ls], fontsize="x-small", loc="upper right")
xlim([-35, 35])
ylim([0, 1.49])
tax2.tick_params(axis="y", labelcolor="C2")
NaCsPlot.maybe_save("$(prefix)_v0_v63")

NaCsPlot.maybe_show()
