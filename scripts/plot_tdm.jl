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

function read_csv_compressed(fname, args...; kwargs...)
    LibArchive.Reader(fname) do reader
        LibArchive.support_format_raw(reader)
        LibArchive.support_filter_all(reader)
        LibArchive.next_header(reader)
        readdlm(reader, args...; kwargs...)
    end
end

const fname_a3b3 = joinpath(@__DIR__, "../data/TDM_Rosario/DM_a3S_b3P_Full.csv.zst")
const fname_a3c3 = joinpath(@__DIR__, "../data/TDM_Rosario/DM_a3S_c3S_Full.csv.zst")
const fname_X1A1 = joinpath(@__DIR__, "../data/TDM_Rosario/DM_X1S_A1S_Full.csv.zst")
const fname_X1B1 = joinpath(@__DIR__, "../data/TDM_Rosario/DM_X1S_B1P_Full.csv.zst")
const data_a3b3 = read_csv_compressed(fname_a3b3, ',', Float64)
const data_a3c3 = read_csv_compressed(fname_a3c3, ',', Float64)
const data_X1A1 = read_csv_compressed(fname_X1A1, ',', Float64)
const data_X1B1 = read_csv_compressed(fname_X1B1, ',', Float64)

const prefix = joinpath(@__DIR__, "../imgs", "TDM")

figure(figsize=[13.8, 4.8])
plot(data_a3b3[:, 1] .* 0.529177210903, data_a3b3[:, 2], "C1",
     label="\$a^3\\Sigma\\rightarrow b^3\\Pi\$")
plot(data_a3c3[:, 1] .* 0.529177210903, data_a3c3[:, 2], "C2",
     label="\$a^3\\Sigma\\rightarrow c^3\\Sigma\$")
plot(data_X1A1[:, 1] .* 0.529177210903, data_X1A1[:, 2], "C3",
     label="\$X^1\\Sigma\\rightarrow A^1\\Sigma\$")
plot(data_X1B1[:, 1] .* 0.529177210903, data_X1B1[:, 2], "C4",
     label="\$X^1\\Sigma\\rightarrow B^1\\Pi\$")
grid()
xscale("log")
xlim([2.7, 25])
xticks([3, 4, 5, 6, 7, 8, 9, 10, 20], ["3", "4", "5", "6", "7", "8", "9", "10", "20"])
# xticks([4, 6, 8, 10, 20], ["4", "6", "8", "10", "20"])
# xticks([3, 10], ["3", "10"])
gca().set_xticks([3:10; 20], minor=true)
gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
xlabel("R (Angstroms)")
ylabel("Transition Dipole Moment (\$ea_0\$)")
legend(fontsize=13.88)
gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
NaCsPlot.maybe_save("$(prefix)_wide")

figure(figsize=[6.4, 4.8])
plot(data_a3b3[:, 1] .* 0.529177210903, data_a3b3[:, 2], "C1",
     label="\$a^3\\Sigma\\rightarrow b^3\\Pi\$")
plot(data_a3c3[:, 1] .* 0.529177210903, data_a3c3[:, 2], "C2",
     label="\$a^3\\Sigma\\rightarrow c^3\\Sigma\$")
plot(data_X1A1[:, 1] .* 0.529177210903, data_X1A1[:, 2], "C3",
     label="\$X^1\\Sigma\\rightarrow A^1\\Sigma\$")
plot(data_X1B1[:, 1] .* 0.529177210903, data_X1B1[:, 2], "C4",
     label="\$X^1\\Sigma\\rightarrow B^1\\Pi\$")
grid()
xscale("log")
xlim([2.7, 25])
xticks([3, 4, 5, 6, 7, 8, 9, 10, 20], ["3", "4", "5", "6", "7", "8", "9", "10", "20"])
# xticks([4, 6, 8, 10, 20], ["4", "6", "8", "10", "20"])
# xticks([3, 10], ["3", "10"])
gca().set_xticks([3:10; 20], minor=true)
gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
xlabel("R (Angstroms)")
ylabel("Transition Dipole Moment (\$ea_0\$)")
legend(fontsize=13.88)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
