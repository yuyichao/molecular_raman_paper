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

const fname = joinpath(@__DIR__, "../data/FCFtoc3Sigma.csv")
const data = readdlm(fname, ',', skipstart=1)
const prefix = joinpath(@__DIR__, "../imgs", "fcf")

figure(figsize=[13.8, 4.8])
plot(data, "C0o-")
xlabel("\$v'\$")
ylabel("FCF")
yscale("log")
xlim([0, 73])
grid()
NaCsPlot.maybe_save("$(prefix)_c3sigma")

NaCsPlot.maybe_show()
