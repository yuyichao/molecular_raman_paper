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

# For anyone that see this script, the analysis done here makes no sense
# See the version as of 8dc37f814f6978f9fe76c549cc75dfe4196fb6b9
# for how this should be done correctly.

const res_freqs = [
    # 668.35 21 772.1367 0.0082
    # 668.35 21 772.1376 0.0049
    # 668.35 21 772.1172 0.0048
    # 668.35 18 771.924 0.033
    668.35 12 771.323 0.050

    560 15 770.59417 0.00013
    560 6 770.369081 0.000037
    560 3 770.287276 0.000015

    565 15 770.61656 0.00021
    565 6 770.378274 0.000080
    565 3 770.291891 0.000038

    555 15 770.57744 0.00013
    555 6 770.362165 0.000043
    555 3 770.283665 0.000017

    503 15 770.483553 0.000049
    503 12 770.43248 0.00011
    503 9 770.379505 0.000028
    503 6 770.324463 0.000019
    503 3 770.2648307 0.0000084

    492 15 770.47571 0.00010
    492 12 770.42599 0.00032
    492 9 770.375007 0.000069
    492 6 770.321355 0.000037
    492 3 770.263402 0.000028

    510 15 770.50117 0.00083
    510 12 770.44659 0.00021
    510 9 770.390205 0.000079
    510 6 770.331510 0.000057
    510 3 770.258511 0.000023

    605 15 770.76208 0.00027
    605 12 770.65624 0.00046
    605 9 770.54816 0.00058
    605 6 770.437052 0.000072
    605 3 770.321367 0.000038

    625 15 770.89410 0.00024
    625 6 770.490130 0.000059
    625 3 770.347812 0.000022
]

const rabi_freqs = [
    560 15 3.971 0.089
    560 6 1.182 0.033
    560 3 0.488 0.015

    503 15 2.073 0.051
    # 503 9 1.074 0.036
    503 6 0.680 0.022
    503 3 0.2544 0.0077

    605 15 6.39 0.15
    605 6 1.985 0.043
    605 3 0.845 0.042

    625 15 8.35 0.25
    625 6 2.513 0.036
    625 3 0.996 0.029
]

function model_sqr(x, p)
    return p[1] .+ (p[2] .+ p[3] .* x) .* x
end

function gen_pwr_model(pwr)
    return (x, p)->p[1] .* x.^pwr
end

broacast_array(f, x) = f(x)
broacast_array(f, x::AbstractArray) = f.(x)

# x: freq, power, idx
# p: fpa0, (offset, slope)...
function model_res(xs, p)
    function real_model(x)
        freq, power, idx = x
        fpa0 = p[1]
        offset = p[idx * 2]
        slope = p[idx * 2 + 1]
        return offset + slope / (freq - fpa0)
    end
    return broacast_array(real_model, xs)
end

# x: freq, power, idx
# p: (offset, strength)...
function gen_model_rabi(fpa0)
    function model(xs, p)
        function real_model(x)
            freq, power, idx = x
            offset = p[idx * 2 - 1]
            slope = p[idx * 2]
            return offset + slope / (freq - fpa0)
        end
        return broacast_array(real_model, xs)
    end
end

function gen_data(data_in)
    power_ids = Dict{Float64,Int}()
    xs = Tuple{Float64,Float64,Int}[]
    ys = Float64[]
    uncs = Float64[]
    idmax = 0
    for i in 1:size(data_in, 1)
        freq = data_in[i, 1]
        power = data_in[i, 2]
        power_id = get!(power_ids, power) do
            idmax += 1
            return idmax
        end
        push!(xs, (freq, power, power_id))
        push!(ys, data_in[i, 3])
        push!(uncs, data_in[i, 4])
    end
    return (x=xs, y=ys, unc=uncs, power_ids=power_ids)
end

function fit_power_model(model, data, _p0)
    p0 = [_p0; zeros(length(data.power_ids) * 2)]
    return fit_data(model, data.x, data.y, data.unc, p0, plotx=false)
end

get_power_idx(data, power) = data.power_ids[power]

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
        return model((x, power, data.power_ids[power]), fit.param)
    end
    plotx = get_plot_range(xs)
    return (x=xs, y=ys, unc=uncs, plotx=plotx, ploty=plot_func.(plotx))
end

function get_offsets_and_slopes(data, fit, powers=sort(collect(keys(data.power_ids))))
    offsets = Float64[]
    offsets_unc = Float64[]
    slopes = Float64[]
    slopes_unc = Float64[]
    idx_offset = length(fit.param) - length(data.power_ids) * 2
    for p in powers
        id = get_power_idx(data, p)
        push!(offsets, fit.param[idx_offset + id * 2 - 1])
        push!(offsets_unc, fit.unc[idx_offset + id * 2 - 1])
        push!(slopes, fit.param[idx_offset + id * 2])
        push!(slopes_unc, fit.unc[idx_offset + id * 2])
    end
    return (powers=powers, offsets=offsets, offsets_unc=offsets_unc,
            slopes=slopes, slopes_unc=slopes_unc)
end

function calc_omega_ma(fits_res, fits_rabi)
    @assert fits_res.powers == fits_rabi.powers
    npowers = length(fits_res.powers)
    omega_ms = Vector{Float64}(undef, npowers)
    omega_ms_unc = Vector{Float64}(undef, npowers)
    omega_as = Vector{Float64}(undef, npowers)
    omega_as_unc = Vector{Float64}(undef, npowers)
    for idx in 1:npowers
        a = fits_res.slopes[idx] * 1e15 # MHz GHz
        a_s = fits_res.slopes_unc[idx] * 1e15 # MHz GHz
        b = fits_rabi.slopes[idx] * 1e12 # kHz GHz
        b_s = fits_rabi.slopes_unc[idx] * 1e12 # kHz GHz
        # We'll just assume Omega_m >> Omega_a here since that is the case.
        # It's equivalent to takeing the approximation that a^2 >> b^2 here
        # so we can still claim that we didn't use the Omega_m >> Omega_a approximation
        # to stop ppl from complaining!
        Ωm² = Unc(-a, a_s) * 2
        Ωm = sqrt(Ωm²)
        Ωa = 2 * Unc(-b, b_s) / Ωm
        omega_ms[idx] = Ωm.a
        omega_ms_unc[idx] = Ωm.s
        omega_as[idx] = Ωa.a
        omega_as_unc[idx] = Ωa.s
    end
    return (powers = fits_res.powers,
            omega_ms=omega_ms, omega_ms_unc=omega_ms_unc,
            omega_ms_uncs=Unc.(omega_ms, omega_ms_unc),
            omega_as=omega_as, omega_as_unc=omega_as_unc,
            omega_as_uncs=Unc.(omega_as, omega_as_unc))
end

const data_res = gen_data(res_freqs)
const fit_res = fit_power_model(model_res, data_res, [705])
@show fit_res.uncs
# Hmm, ok, now we got 3 numbers (shown in the figure) to fit a quadratic curve for the offset.
# It'll of course work out great.
const fits_res = get_offsets_and_slopes(data_res, fit_res, [3, 6, 15])
@show fits_res
const fit_offset_res = fit_data(model_sqr, fits_res.powers, fits_res.offsets,
                                fits_res.offsets_unc, [770.0, 0.0, 0.0])
const fit_slope_res = fit_data(gen_pwr_model(1), fits_res.powers, fits_res.slopes,
                               fits_res.slopes_unc, [0.0])
@show fit_offset_res.uncs
@show fit_slope_res.uncs

const model_rabi = gen_model_rabi(fit_res.param[1])
const data_rabi = gen_data(rabi_freqs)
const fit_rabi = fit_power_model(model_rabi, data_rabi, Float64[])
@show fit_rabi.uncs
const fits_rabi = get_offsets_and_slopes(data_rabi, fit_rabi, [3, 6, 15])
@show fits_rabi
const fit_offset_rabi = fit_data(gen_pwr_model(1.29), fits_rabi.powers, fits_rabi.offsets,
                                 fits_rabi.offsets_unc, [0.0])
const fit_slope_rabi = fit_data(gen_pwr_model(1.29), fits_rabi.powers, fits_rabi.slopes,
                                fits_rabi.slopes_unc, [0.0])
@show fit_offset_rabi.uncs
@show fit_slope_rabi.uncs
const omegas = calc_omega_ma(fits_res, fits_rabi)
@show omegas
const fit_omega_m = fit_data(gen_pwr_model(0.5), omegas.powers, omegas.omega_ms,
                             omegas.omega_ms_unc, [0.0])
const fit_omega_a = fit_data(gen_pwr_model(0.79), omegas.powers, omegas.omega_as,
                             omegas.omega_as_unc, [0.0])
@show fit_omega_m.uncs
@show fit_omega_a.uncs

const prefix = joinpath(@__DIR__, "../imgs", "scaling")

const powers_plot = [15, 6, 3]

figure(figsize=[6.4, 4.8] * 1.11308)
for i in 1:length(powers_plot)
    power = powers_plot[i]
    pd = get_plot_data_power(data_res, fit_res, model_res, power)
    errorbar(1000 ./ (pd.x .- fit_res.param[1]), (pd.y .- fit_offset_res.param[1]) .* 1000,
             pd.unc .* 1000, fmt="C$(i - 1)o", label="$(power / 4) mW")
    plot(1000 ./ (pd.plotx .- fit_res.param[1]), (pd.ploty .- fit_offset_res.param[1]) .* 1000,
         "C$(i - 1)")
end
# xticks([288500, 288530, 288560, 288590, 288620])
# p: framan0, fpa0, offset, strength, bs...
# text(560, 1700, ("\$f_{Raman0} + a\\cdot P\$\n" *
#                  "  \$-\\dfrac{b}{f-f_{PA0}}\\cdot P\$"))
# text(480, 740, ("\$f_{Raman0}=$(fit_res.uncs[1])\$ MHz\n" *
#                 "\$f_{PA0}=$(fit_res.uncs[2] + 288000)\$ GHz\n" *
#                 "\$a=$(fit_res.uncs[3] * 1000)\$ kHz/mW\n" *
#                 "\$b=$(fit_res.uncs[4])\$ MHz\$\\cdot\$GHz/mW"), fontsize="small")
legend(loc="upper right", fontsize="small", handlelength=0.6, handletextpad=0.3)
grid()
# xticks([-200, -160, -120, -80])
yticks([0, 200, 400, 600])
ylim([0, 790])
xlabel("\$2\\pi / \\Delta~(\\mathrm{THz^{-1}})\$")
ylabel("\$\\omega_R-\\omega_{R0}\$ (\$\\mathrm{2\\pi\\times kHz}\$)")
NaCsPlot.maybe_save("$(prefix)_light_shift")

figure(figsize=[6.4, 4.8] * 1.11308)
for i in 1:length(powers_plot)
    power = powers_plot[i]
    pd = get_plot_data_power(data_rabi, fit_rabi, model_rabi, power)
    plot(1000 ./ (pd.plotx .- fit_res.param[1]), pd.ploty, "C$(i - 1)")
    errorbar(1000 ./ (pd.x .- fit_res.param[1]), pd.y, pd.unc,
             fmt="C$(i - 1)o", label="$(power / 4) mW")
end
# xticks([288500, 288530, 288560, 288590, 288620])
# p: offset, strength
# text(500, 5.0, ("\$\\left(a-\\dfrac{b}{f-705 GHz}\\right)\\cdot P^{1.29}\$"))
# text(500, 2.8, ("\$a=$(fit_rabi.uncs[1] * 1000)\$ Hz/mW\$^{1.29}\$\n" *
#               "\$b=$(fit_rabi.uncs[2])\$ kHz\$\\cdot\$GHz/mW\$^{1.29}\$"), fontsize="small")
legend(loc="upper right", fontsize="small", handlelength=0.6, handletextpad=0.3)
grid()
yticks([0, 2, 4, 6, 8])
# xticks([-200, -160, -120, -80])
ylim([0, 9.7])
xlabel("\$2\\pi / \\Delta~(\\mathrm{THz^{-1}})\$")
ylabel("\$\\Omega_{R}~(\\mathrm{2\\pi\\times kHz})\$")
NaCsPlot.maybe_save("$(prefix)_rabi")

fig = figure(figsize=[3.6, 5.4] .* 1.0516)

ax1 = fig.add_subplot(211)
errorbar(omegas.powers[1] ./ 4, omegas.omega_ms[1] ./ 1e9,
         omegas.omega_ms_unc[1] ./ 1e9, fmt="C0o")
errorbar(omegas.powers[2] ./ 4, omegas.omega_ms[2] ./ 1e9,
         omegas.omega_ms_unc[2] ./ 1e9, fmt="C1o")
errorbar(omegas.powers[3] ./ 4, omegas.omega_ms[3] ./ 1e9,
         omegas.omega_ms_unc[3] ./ 1e9, fmt="C2o")
plot(fit_omega_m.plotx ./ 4, fit_omega_m.ploty ./ 1e9, "k")
xscale("log")
yscale("log")
grid()
xticks([1, 2, 4], ["1", "2", "4"])
yticks([0.2, 0.3], ["0.2", "0.3"])
xlim([0.625, 4.25])
ylim([0.14, 0.39])
ylabel("\$\\Omega_m\$ (\$\\mathrm{2\\pi\\times GHz}\$)")
setp(ax1.get_xticklabels(), visible=false)
ax1.tick_params(axis="x", length=0)
ax1.set_xticks(0.7:0.1:4, minor=true)
ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax2 = fig.add_subplot(212)
subplots_adjust(hspace=0.0)
errorbar(omegas.powers[1] ./ 4, omegas.omega_as[1] ./ 1e6,
         omegas.omega_as_unc[1] ./ 1e6, fmt="C0o")
errorbar(omegas.powers[2] ./ 4, omegas.omega_as[2] ./ 1e6,
         omegas.omega_as_unc[2] ./ 1e6, fmt="C1o")
errorbar(omegas.powers[3] ./ 4, omegas.omega_as[3] ./ 1e6,
         omegas.omega_as_unc[3] ./ 1e6, fmt="C2o")
plot(fit_omega_a.plotx ./ 4, fit_omega_a.ploty ./ 1e6, "k")
xscale("log")
yscale("log")
grid()
xticks([1, 2, 4], ["1", "2", "4"])
yticks([2, 3, 4, 5], ["2", "3", "4", "5"])
xlim([0.625, 4.25])
ylim([1.25, 6])
xlabel("Tweezer Power (mW)")
ylabel("\$\\Omega_a\$ (\$\\mathrm{2\\pi\\times MHz}\$)")
ax2.set_xticks(0.7:0.1:4, minor=true)
ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
NaCsPlot.maybe_save("$(prefix)_slopes")

NaCsPlot.maybe_show()
