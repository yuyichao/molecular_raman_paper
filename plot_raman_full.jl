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
const prefix = joinpath(@__DIR__, "imgs", "raman_theory_full")

figure(figsize=[12.8, 4.8])
l1 = plot(data[:, 1] ./ 1000, abs.(data[:, 2] ./ 2π / 1000), "C0",
          label="\$\\Omega_{R}\$")
l2 = plot(data[:, 1] ./ 1000, abs.(data[:, 4] ./ 2π / 1000), "C1",
          label="\$\\Gamma_{s}\$")
ylabel("\$2\\pi\\cdot \\mathrm{kHz}\$")
xlabel("Raman Single-Photon Frequency (THz)")
yscale("log")
ylim([0.003, 100])
yticks([0.01, 0.1, 1, 10, 100], ["0.01", "0.1", "1", "10", "100"])
grid()
legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, 0.95), frameon=false)
gca().yaxis.set_label_coords(-0.04, 0.55)
xlim([287.5, 351.5])
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
