#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("../../experiments/nacs_202003/molecular_raman_model.jl")

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["nacs_202010/data/data_20201001_203322.mat",
                "nacs_202010/data/data_20201002_203848.mat",
                "nacs_202010/data/data_20201002_085901.mat",
                "nacs_202010/data/data_20201003_102015.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "../../experiments/", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [[0; [0.02, 0.04, 0.07, 0.085, 0.10, 0.13, 0.145, 0.16,
                    0.175, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34, 0.37,
                    0.40, 0.43, 0.46, 0.49] .* 1.6 .- 0.01], # 15 mW, 770.5704 MHz
               ([0], # t=0
                [0], # pi pulse
                [0, 0.03, 0.06, 0.11, 0.15, 0.2, 0.3, 0.4, 0.8, 1.5]), # molecule lifetime
               ([0.0], # 15 mW, 0 ms
                570.4 .+ (-8:1.6:8), # 15 mW, 0.13 ms
                570.4 .+ (-11:2.2:11), # 15 mW, 0.26 ms
                ),
               ([0.0, 2, 4, 10],
                [0.0, 20, 40, 80, 120],
                [0.0, 20, 40, 80, 160] .* 2,
                [0.0, 25, 50, 100],
                [0.0, 50, 100, 200],
                [0.0, 50, 100, 200] .* 2),
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((-1, 2), (-3, 4,)), maxcnts, specs)
const datas_na = select_datas(datas, NaCsData.select_single((1, -2), (3, -4,)), maxcnts, specs)

const data_nacs_t = datas_nacs[1] # Survival 1
const data_nacs_00 = datas_nacs[3][1] # Survival 2
const data_nacs_12 = datas_nacs[3][2] # Survival 2
const data_nacs_25 = datas_nacs[3][3] # Survival 2

const data_pa_na = [datas_na[4][4]; datas_na[4][1]]
const data_pa_cs = datas_cs[4][4]
const data_pa_a = datas_nacs[4][4]
const data_pa_m = datas_nacs[4][1]

const data_nacs_m = datas_nacs[2][3]

const data_fit = [NaCsData.map_params((i, v) -> (1, 570.4, v, 2), data_nacs_00);
                  NaCsData.map_params((i, v) -> (1, v, 0.12, 2), data_nacs_12);
                  NaCsData.map_params((i, v) -> (1, v, 0.25, 2), data_nacs_25);
                  NaCsData.map_params((i, v) -> (1, 570.4, v, 1), data_nacs_t);
                  NaCsData.map_params((i, v) -> (2, v, 0.0, 3), data_pa_na);
                  NaCsData.map_params((i, v) -> (3, v, 0.0, 4), data_pa_cs);
                  NaCsData.map_params((i, v) -> (4, v, 0.0, 5), data_pa_a);
                  NaCsData.map_params((i, v) -> (5, v, 0.0, 6), data_pa_m);
                  NaCsData.map_params((i, v) -> (6, v, 0.0, 7), data_nacs_m)]

const prefix = joinpath(@__DIR__, "imgs", "exp_data")

function model_exp(x, p)
    if length(p) > 2
        p[1] .* exp.(.- x .* p[2]) .+ p[3]
    else
        p[1] .* exp.(.- x .* p[2])
    end
end

function model_ramsey(x, p)
    p[1] .* exp.(.- x .* p[2]) .* cos.(x .* p[4] .* 2π).^2 .+ p[3]
end

function get_model_param(p, idx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r = p
    p1 = p[9 + idx]
    p0 = p1 * p0r
    return (p0, p1, f0, Ω, Γ1 + Γ_na + Γ_cs, Γ2)
end

function gen_pa_param(p, typ, pidx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r = p
    p1 = p[9 + pidx]
    if typ == 1
        return (p1, Γ_na)
    elseif typ == 2
        return (p1, Γ_cs)
    elseif typ == 3
        return (p1, Γ_na + Γ_cs)
    elseif typ == 4
        return (p1, Γ1 + Γ_na + Γ_cs)
    end
end

function gen_ramsey_param(p, pidx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r, p0rm, df = p
    p1 = p[9 + pidx]
    return (p1, Γ2, p0rm * p1, f0 - 570.4 + df)
end

function model(xs, p)
    # p: f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r, p1s...
    function wrapper(x)
        typ = x[1]
        if typ == 1
            _, f, t, idx = x
            return model_2d(t, f, get_model_param(p, idx))
        elseif typ == 6
            _, t, _, idx = x
            return model_ramsey(t, gen_ramsey_param(p, idx))
        else
            _, t, _, idx = x
            return model_exp(t, gen_pa_param(p, typ - 1, idx))
        end
    end
    return wrapper.(xs)
end

fit = fit_survival(model, data_fit, [570.4, 2π * 1.5, 0, 2π / 0.2, 0, 0, 0.1, 0.2, 0,
                                     0.3, 0.3, 0.8, 0.3, 0.3, 0.3, 0.15],
                   plotx=false, lower=[-Inf; zeros(15)])
@show fit.uncs
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)
const param_2 = get_model_param(fit.param, 2)
const scale_1 = 1 / (param_1[1] + param_1[2])
const scale_2 = 1 / (param_2[1] + param_2[2])

const plot_freq_lo = 570.4 - 12
const plot_freq_hi = 570.4 + 12
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)
const plot_time = linspace(0, 0.79, 1000)
const plot_pa_time = linspace(0, 21, 1000)
const plot_m_time = linspace(0, 1.55, 1000)

const data_nacs_12′ = NaCsData.map_params((i, v)->v - param_1[3], data_nacs_12)
const data_nacs_25′ = NaCsData.map_params((i, v)->v - param_1[3], data_nacs_25)

ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, param_1)

figure()
errorbar([0.0], [ratio_00], [uncs_00], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi] .- param_1[3], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_12′, fmt="C1.", label="0.12 ms")
plot(plot_freq .- param_1[3], model_2d.(0.12, plot_freq, (param_1,)), "C1")
NaCsPlot.plot_survival_data(data_nacs_25′, fmt="C2.", label="0.25 ms")
plot(plot_freq .- param_1[3], model_2d.(0.25, plot_freq, (param_1,)), "C2")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("Detuning from resonance (kHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_raman_det")

figure()
errorbar([0.0], [ratio_00 * scale_1], [uncs_00 * scale_1], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi] .- param_1[3], [v00 * scale_1, v00 * scale_1], "C0-")
NaCsPlot.plot_survival_data(data_nacs_12′, scale_1, fmt="C1.", label="0.12 ms")
plot(plot_freq .- param_1[3], model_2d.(0.12, plot_freq, (param_1,)) .* scale_1, "C1")
NaCsPlot.plot_survival_data(data_nacs_25′, scale_1, fmt="C2.", label="0.25 ms")
plot(plot_freq .- param_1[3], model_2d.(0.25, plot_freq, (param_1,)) .* scale_1, "C2")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("Detuning from resonance (kHz)")
ylabel("% Atomic Ground State")
NaCsPlot.maybe_save("$(prefix)_norm_raman_det")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
fit_l12 = fit_survival(model_lorentzian, data_nacs_12, [1 / scale_1, 0.25, param_1[3], 5])
const data_nacs_12_shift = NaCsData.map_params((i, v)->v - fit_l12.param[3], data_nacs_12)

x₋ = fit_l12.param[3] - fit_l12.param[4] / 2
x₊ = fit_l12.param[3] + fit_l12.param[4] / 2
y_pm = (model_lorentzian(x₋, fit_l12.param) + model_lorentzian(x₊, fit_l12.param)) / 2

function plot_lorentz()
    NaCsPlot.plot_survival_data(data_nacs_12_shift, fmt="C0.", label="0.12 ms")
    plot(fit_l12.plotx .- fit_l12.param[3], fit_l12.ploty, "C0")
    legend(fontsize="x-small", loc="lower right")
    text(-5.83, 0.27, "\$f_{Raman}=$(770 + fit_l12.uncs[3] / 1000) \\mathrm{MHz}\$",
         color="C0", fontsize="small")
    grid()
    xlabel("Detuning from resonance (kHz)")
    ylabel("Two-body survival")
end

figure()
plot_lorentz()
NaCsPlot.maybe_save("$(prefix)_lortz_f")

figure()
plot_lorentz()
ax = gca()
ax.annotate("", (-fit_l12.param[4] / 2 - 0.1, y_pm),
            xytext=(-fit_l12.param[4] / 2 - 11 / 3, y_pm),
            arrowprops=Dict(:color=>"C3"))
ax.annotate("", (fit_l12.param[4] / 2 + 0.1, y_pm),
            xytext=(fit_l12.param[4] / 2 + 11 / 3, y_pm),
            arrowprops=Dict(:color=>"C3"))
text(-5.14, 0.137, "\$\\mathbf{\\Gamma_{FWHM}=$(fit_l12.uncs[4]) kHz}\$",
     fontsize=17, color="C3")
NaCsPlot.maybe_save("$(prefix)_lortz_f_text")

function plot_lorentz_norm()
    NaCsPlot.plot_survival_data(data_nacs_12_shift, scale_1, fmt="C0.", label="0.12 ms")
    plot(fit_l12.plotx .- fit_l12.param[3], fit_l12.ploty .* scale_1, "C0")
    legend(fontsize="x-small", loc="lower right")
    text(-5.83, 0.27 * scale_1, "\$f_{Raman}=$(770 + fit_l12.uncs[3] / 1000) \\mathrm{MHz}\$",
         color="C0", fontsize="small")
    grid()
    xlabel("Detuning from resonance (kHz)")
    ylabel("% Atomic Ground State")
end

figure()
plot_lorentz_norm()
NaCsPlot.maybe_save("$(prefix)_norm_lortz_f")

figure()
plot_lorentz_norm()
ax = gca()
ax.annotate("", (-fit_l12.param[4] / 2 - 0.1, y_pm * scale_1),
            xytext=(-fit_l12.param[4] / 2 - 11 / 3, y_pm * scale_1),
            arrowprops=Dict(:color=>"C3"))
ax.annotate("", (fit_l12.param[4] / 2 + 0.1, y_pm * scale_1),
            xytext=(fit_l12.param[4] / 2 + 11 / 3, y_pm * scale_1),
            arrowprops=Dict(:color=>"C3"))
text(-5.14, 0.137 * scale_1, "\$\\mathbf{\\Gamma_{FWHM}=$(fit_l12.uncs[4]) kHz}\$",
     fontsize=17, color="C3")
NaCsPlot.maybe_save("$(prefix)_norm_lortz_f_text")

figure()
NaCsPlot.plot_survival_data(data_nacs_t, fmt="C0.", label="770.5704 MHz")
plot(plot_time, model_2d.(plot_time, 570.4, (param_2,)), "C0")
xlim([0, 0.8])
ylim([0, ylim()[2]])
legend(fontsize="x-small", loc="upper right")
annotate("Atom", (0.01, 0.27), xytext=(0.13, 0.26), arrowprops=Dict(:color=>"C3"), color="C3")
annotate("Molecule", (0.155, 0.030),
         xytext=(0.097, 0.155), arrowprops=Dict(:color=>"C3"), color="C3")
annotate("Atom", (0.295, 0.140),
         xytext=(0.42, 0.16), arrowprops=Dict(:color=>"C3"), color="C3")
grid()
xlabel("Raman time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_t")

figure()
NaCsPlot.plot_survival_data(data_nacs_t, scale_2,
                            fmt="C0.", label="770.5704 MHz")
plot(plot_time, model_2d.(plot_time, 570.4, (param_2,)) .* scale_2, "C0")
xlim([0, 0.8])
ylim([0, ylim()[2]])
legend(fontsize="x-small", loc="upper right")
annotate("Atom", (0.01, 1), xytext=(0.13, 0.96), arrowprops=Dict(:color=>"C3"), color="C3")
annotate("Molecule", (0.155, 0.111),
         xytext=(0.097, 0.573), arrowprops=Dict(:color=>"C3"), color="C3")
annotate("Atom", (0.295, 0.517),
         xytext=(0.42, 0.591), arrowprops=Dict(:color=>"C3"), color="C3")
grid()
xlabel("Raman time (ms)")
ylabel("% Atomic Ground State")
NaCsPlot.maybe_save("$(prefix)_norm_t")

figure()
NaCsPlot.plot_survival_data(data_pa_m, fmt="C0.")
plot(plot_pa_time, model_exp(plot_pa_time, gen_pa_param(fit.param, 4, 6)), "C0-")
text(7, 0.12, "\$\\tau_{atom}=$(1 / uncs_1[5])\\ \\mathrm{ms}\$", color="C0")
xlim([0, 22])
ylim([0, 0.3])
grid()
xlabel("Atom hold time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_a_t")

const param_pa = gen_pa_param(fit.param, 4, 6)

figure()
NaCsPlot.plot_survival_data(data_pa_m, fmt="C0.")
plot(plot_pa_time, model_exp(plot_pa_time, param_pa), "C0-")
text(7, 0.12, "\$\\tau_{atom}=$(1 / uncs_1[5])\\ \\mathrm{ms}\$", color="C0")
xlim([0, 22])
ylim([0, 0.3])
grid()
xlabel("Atom hold time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_a_t")

figure()
NaCsPlot.plot_survival_data(data_pa_m, 1 / param_pa[1], fmt="C0.")
plot(plot_pa_time, model_exp(plot_pa_time, param_pa) ./ param_pa[1], "C0-")
text(7, 0.12 / param_pa[1], "\$\\tau_{atom}=$(1 / uncs_1[5])\\ \\mathrm{ms}\$", color="C0")
xlim([0, 22])
ylim([0, 1.05])
grid()
xlabel("Atom hold time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_norm_a_t")

const param_ramsey = gen_ramsey_param(fit.param, 7)
const scale_ramsey = 1 / (param_ramsey[1] + param_ramsey[3])

figure()
NaCsPlot.plot_survival_data(data_nacs_m[1:end - 1], fmt="C0s")
plot(plot_m_time, model_ramsey(plot_m_time, param_ramsey), "C0-")
text(0.218, 0.07, "\$\\tau_{molecule}=$(1 / uncs_1[6])\\ \\mathrm{ms}\$", color="C0")
xlim([0, 0.9])
grid()
xlabel("Molecule Hold Time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_m_t")

figure()
NaCsPlot.plot_survival_data(data_nacs_m[1:end - 1], scale_ramsey, fmt="C0s")
plot(plot_m_time, model_ramsey(plot_m_time, param_ramsey) .* scale_ramsey, "C0-")
text(0.218, 0.07 * scale_ramsey,
     "\$\\tau_{molecule}=$(1 / uncs_1[6])\\ \\mathrm{ms}\$", color="C0")
xlim([0, 0.9])
grid()
xlabel("Molecule Hold Time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_norm_m_t")

const img_freq = param_1[3] .+ linspace(-20, 20, 201)
const img_time = linspace(0, 0.5, 201)
const mol_2d = [model_2d(t, f, param_1, molecule=true) .* scale_1
                for f in img_freq, t in img_time]

figure()
imshow(mol_2d, aspect="auto", interpolation="none", origin="lower",
       extent=[img_time[1] - step(img_time) / 2, img_time[end] + step(img_time) / 2,
               img_freq[1] - step(img_freq) / 2 - param_1[3],
               img_freq[end] + step(img_freq) / 2 - param_1[3]])
xlabel("Raman time (ms)")
ylabel("Detuning from resonance (kHz)")
title("Molecule Population")
grid()
colorbar()
NaCsPlot.maybe_save("$(prefix)_norm_2d")

NaCsPlot.maybe_show()
