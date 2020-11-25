#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

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

# const fname = joinpath(@__DIR__, "../../damop-2020/20MHz_linewidth_c3SigmaOnly_3322_80kHzConfinement.csv.zst")
const fname = joinpath(@__DIR__, "../data/50MHz_linewidth_8TotalExcitedStates_3322_80kHzConfinement_3.75mWPowerAtAtom.csv.zst")
const data = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    readdlm(reader, ',', skipstart=1)
end
const prefix = joinpath(@__DIR__, "../imgs", "raman")

fig = figure(figsize=[1.11, 1.5] * 4.8)

ax1 = fig.add_subplot(311)
plot(data[:, 1] .- 339724.57, abs.(data[:, 2] ./ 2π / 1000), "C1",
     label="\$v'=40\$")
plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ 2π / 1000), "C0",
     label="\$v'=0\$")
ylabel("\$\\Omega_{R}~(2\\pi\\!\\times\\!\\mathrm{kHz})\$")
xlim([-300, 300])
ylim([6, 20000])
yscale("log")
# yticks([10, 100, 1000, 10000], ["10", "100", "1000", "10000"])
ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
grid()
setp(ax1.get_xticklabels(), visible=false)
setp(ax1.get_xticklabels(), visible=false)
ax1.tick_params(axis="x", length=0)

ax2 = fig.add_subplot(312)
subplots_adjust(hspace=0.018)
plot(data[:, 1] .- 288625.081, abs.(data[:, 4] ./ 2π / 1000), "C0",
     label="\$v'=0\$")
plot(data[:, 1] .- 339724.57, abs.(data[:, 4] ./ 2π / 1000), "C1",
     label="\$v'=40\$")
ylabel("\$\\Gamma_{s}~(2\\pi\\!\\times\\!\\mathrm{kHz})\$")
legend(fontsize=13.88, loc="upper left", handlelength=1, handletextpad=0.3)
xlim([-300, 300])
ylim([0.12, 200000])
yscale("log")
# yticks([1, 10, 100, 1000, 10000], ["1", "10", "100", "1000", "10000"])
ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
grid()
setp(ax2.get_xticklabels(), visible=false)
setp(ax2.get_xticklabels(), visible=false)
ax2.tick_params(axis="x", length=0)

ax3 = fig.add_subplot(313)
subplots_adjust(hspace=0.018)
plot(data[:, 1] .- 288625.081, abs.(data[:, 2] ./ data[:, 4]), "C0",
     label="\$v'=0\$")
plot(data[:, 1] .- 339724.57, abs.(data[:, 2] ./ data[:, 4]), "C1",
     label="\$v'=40\$")
xlabel("One-Photon Detuning (GHz)")
ylabel("\$\\Omega_{R}/\\Gamma_{s}\$")
xlim([-300, 300])
ylim([0.02, 90])
yscale("log")
yticks([0.1, 1, 10], ["0.1", "1", "10"])
ax3.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
grid()
NaCsPlot.maybe_save("$(prefix)_v0_vhi")

NaCsPlot.maybe_show()
