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

read_zst_csv(fname) = LibArchive.Reader(fname) do reader
    LibArchive.support_format_raw(reader)
    LibArchive.support_filter_all(reader)
    LibArchive.next_header(reader)
    readdlm(reader, ',', skipstart=1)
end

# const fname = joinpath(@__DIR__, "../../damop-2020/20MHz_linewidth_c3SigmaOnly_3322_80kHzConfinement.csv.zst")
const fname = joinpath(@__DIR__, "../data/50MHz_linewidth_8TotalExcitedStates_3322_80kHzConfinement_3.75mWPowerAtAtom.csv.zst")
const fname2 = joinpath(@__DIR__, "../data/50MHz_linewidth_8TotalExcitedStates_3322_80kHzConfinement_3.75mWPowerAtAtom_ThresholdContribution.csv.zst")
const data = read_zst_csv(fname)
const data2 = read_zst_csv(fname2)
const prefix = joinpath(@__DIR__, "../imgs", "raman_theory_full")

figure(figsize=[13.8, 4.8])
plot(data[:, 1] ./ 1000, abs.(data[:, 2] ./ 2π / 1000), "C0", label="\$\\Omega_{R}\$")
plot(data2[:, 1] ./ 1000, abs.(data2[:, 4] ./ 2π / 1000), "r", linewidth=4, alpha=0.7)
plot(data[:, 1] ./ 1000, abs.(data[:, 4] ./ 2π / 1000), "C1", label="\$\\Gamma_{s}\$")
# plot(data[:, 1] ./ 1000, abs.(data[:, 2] ./ data[:, 4]), "C2")
# es = [2.8863, 2.9026, 2.9184, 2.9339, 2.9492, 2.9643, 2.9792, 2.9939, 3.0084, 3.0227,
#       3.0369, 3.0510, 3.0649, 3.0787, 3.0925, 3.1062, 3.1198, 3.1334, 3.1469, 3.1603,
#       3.1736, 3.1868, 3.1999, 3.2129, 3.2258, 3.2385, 3.2511, 3.2635, 3.2757, 3.2876,
#       3.2993, 3.3106, 3.3215, 3.3322, 3.3425, 3.3525, 3.3621, 3.3714, 3.3803, 3.3889,
#       3.3972, 3.4052, 3.4129, 3.4204, 3.4276, 3.4345, 3.4412, 3.4475, 3.4537, 3.4597,
#       3.4654, 3.4708, 3.4759, 3.4808, 3.4852, 3.4894, 3.4932, 3.4967, 3.4999, 3.5028,
#       3.5054, 3.5076, 3.5096, 3.5113, 3.5127, 3.5139, 3.5149, 3.5157, 3.5162, 3.5167,
#       3.5170, 3.5171, 3.5172, 3.5173]
# for v in es
#     axvline(v * 100, color="C3", linewidth=5)
# end
ylabel("\$2\\pi\\times \\mathrm{kHz}\$")
xlabel("Tweezer Frequency (THz)")
yscale("log")
ylim([0.12, 400])
yticks([1, 10, 100], ["1", "10", "100"])
grid()
legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, 0.95), frameon=false)
gca().yaxis.set_label_coords(-0.034, 0.55)
xlim([287, 340])
NaCsPlot.maybe_save("$(prefix)")

figure(figsize=[13.8, 2.6])
plot(data2[:, 1] ./ 1000, abs.(data2[:, 4] ./ 2π / 1000), "r", linewidth=4, alpha=0.7, ls=":")
plot(data[:, 1] ./ 1000, abs.(data[:, 4] ./ 2π / 1000), "r", alpha=0.4)
axvline(288.625081, color="C0", ls="--", linewidth=2)
axvline(339.724570, color="C0", ls="--", linewidth=2)
text(288.625081, 10, "v'=0", rotation=90, va="center", ha="right", color="C0", fontsize=16)
text(339.724570, 3, "v'=40", rotation=90, va="center", ha="right", color="C0", fontsize=16)
# plot(data[:, 1] ./ 1000, abs.(data[:, 2] ./ data[:, 4]), "C2")
# es = [2.8863, 2.9026, 2.9184, 2.9339, 2.9492, 2.9643, 2.9792, 2.9939, 3.0084, 3.0227,
#       3.0369, 3.0510, 3.0649, 3.0787, 3.0925, 3.1062, 3.1198, 3.1334, 3.1469, 3.1603,
#       3.1736, 3.1868, 3.1999, 3.2129, 3.2258, 3.2385, 3.2511, 3.2635, 3.2757, 3.2876,
#       3.2993, 3.3106, 3.3215, 3.3322, 3.3425, 3.3525, 3.3621, 3.3714, 3.3803, 3.3889,
#       3.3972, 3.4052, 3.4129, 3.4204, 3.4276, 3.4345, 3.4412, 3.4475, 3.4537, 3.4597,
#       3.4654, 3.4708, 3.4759, 3.4808, 3.4852, 3.4894, 3.4932, 3.4967, 3.4999, 3.5028,
#       3.5054, 3.5076, 3.5096, 3.5113, 3.5127, 3.5139, 3.5149, 3.5157, 3.5162, 3.5167,
#       3.5170, 3.5171, 3.5172, 3.5173]
# for v in es
#     axvline(v * 100, color="C3", linewidth=5)
# end
ylabel("\$\\Gamma_{s}~(2\\pi\\times \\mathrm{kHz})\$")
xlabel("Tweezer Frequency (THz)")
yscale("log")
ylim([0.12, 400])
yticks([1, 10, 100], ["1", "10", "100"])
grid()
gca().yaxis.set_label_coords(-0.04, 0.5)
xlim([287, 340])
NaCsPlot.maybe_save("$(prefix)_gamma")

NaCsPlot.maybe_show()
