module CodeFigures

using Makie, CairoMakie
using Quantikz
using QuantumClifford: stab_to_gf2, stabilizerplot_axis
using QuantumClifford.ECC: TableDecoder, parity_checks, iscss, parity_matrix_z, parity_matrix_x, code_n, naive_encoding_circuit, naive_syndrome_circuit, shor_syndrome_circuit

using ..Helpers: logrange, instancenameof, skipredundantfix, typenameof, CircBuffer
using ..DBHelpers: dbrow, dbnarray, dbrow!

ispositivefinite(x) = x > 0 && isfinite(x)
skipzeronan(xs) = (x for x in xs if ispositivefinite(x))
plottable_log_rate(x) = ispositivefinite(x) ? x : NaN

function plottable_logical_max(x, z)
    xok = ispositivefinite(x)
    zok = ispositivefinite(z)
    xok && zok && return max(x, z)
    xok && return x
    zok && return z
    return NaN
end

function padded_loglimits(values; default=(1e-5, 1.0), lower=nothing, upper=nothing, pad_decades_fraction=0.08)
    vals = collect(skipzeronan(values))
    isempty(vals) && return default

    lo = Float64(minimum(vals))
    hi = Float64(maximum(vals))
    loglo = log10(lo)
    loghi = log10(hi)
    pad = loglo == loghi ? 0.5 : max((loghi - loglo) * pad_decades_fraction, 0.04)

    lo = 10.0^(loglo - pad)
    hi = 10.0^(loghi + pad)
    !isnothing(lower) && (lo = max(lo, lower))
    !isnothing(upper) && (hi = min(hi, upper))

    if lo >= hi
        center = (log10(lo) + log10(hi)) / 2
        lo = 10.0^(center - 0.5)
        hi = 10.0^(center + 0.5)
        !isnothing(lower) && (lo = max(lo, lower))
        !isnothing(upper) && (hi = min(hi, upper))
    end

    return lo, hi
end

function plotted_limit_values(phys_errors, results, single_error)
    physical_values = Float64[]
    logical_values = Float64[]

    if single_error
        for iᵉ in eachindex(phys_errors)
            ispositivefinite(phys_errors[iᵉ]) || continue
            for iᶜ in axes(results, 3), iᵈ in axes(results, 4), iˢ in axes(results, 5)
                logical_error = plottable_logical_max(results[iᵉ, 1, iᶜ, iᵈ, iˢ], results[iᵉ, 2, iᶜ, iᵈ, iˢ])
                ispositivefinite(logical_error) || continue
                push!(physical_values, phys_errors[iᵉ])
                push!(logical_values, logical_error)
            end
        end
    else
        for iᵉ in eachindex(phys_errors)
            ispositivefinite(phys_errors[iᵉ]) || continue
            for iˡ in axes(results, 2), iᶜ in axes(results, 3), iᵈ in axes(results, 4), iˢ in axes(results, 5)
                logical_error = results[iᵉ, iˡ, iᶜ, iᵈ, iˢ]
                ispositivefinite(logical_error) || continue
                push!(physical_values, phys_errors[iᵉ])
                push!(logical_values, logical_error)
            end
        end
    end

    return physical_values, logical_values
end

function make_decoder_figure(phys_errors, results;
    title="",
    colors=CircBuffer(Makie.wong_colors()),
    linestyles=CircBuffer([:solid, :dash, :dot, :dashdot, :dashdotdot, Linestyle([0.5, 1.0, 1.5, 2.5])]),
    markers=CircBuffer(['●', '■', '▲', '▼', '◆', '★']),
    single_error=false,
    codelabels=[], decoderlabels=[], setuplabels=[],
)
    physical_values, logical_values = plotted_limit_values(phys_errors, results, single_error)
    positive_phys_errors = collect(skipzeronan(phys_errors))
    physical_default = isempty(positive_phys_errors) ? (1e-5, 1.0) : (minimum(positive_phys_errors), min(1.0, maximum(positive_phys_errors)))
    xlimits = padded_loglimits(physical_values; default=physical_default, lower=physical_default[1], upper=1.0)
    ylimits = padded_loglimits(logical_values; upper=1.0)

    fresults = plottable_log_rate.(results)
    combined_results = map(plottable_logical_max, results[:,1,:,:,:], results[:,2,:,:,:])

    f = Figure(size=(1000,600))
    a = Axis(f[1:7,1:6],
        xscale=log10, yscale=log10,
        limits=(xlimits[1], xlimits[2], ylimits[1], ylimits[2]),
        xlabel="physical error rate",
        ylabel="logical error rate",
        title=title,
        titlesize=26,
        xlabelsize=18,
        ylabelsize=18,
        xticklabelsize=13,
        yticklabelsize=13,
        xgridcolor=(:gray, 0.18),
        ygridcolor=(:gray, 0.18),
        xminorgridvisible=true,
        yminorgridvisible=true,
        xminorgridcolor=(:gray, 0.08),
        yminorgridcolor=(:gray, 0.08),
        topspinevisible=false,
        rightspinevisible=false,
        )

    singlecode = size(results,3) == 1
    plotcolor(iᶜ,iᵈ) = singlecode ? colors[iᵈ] : colors[iᶜ]
    decoderlegendcolor(iᵈ) = singlecode ? colors[iᵈ] : :gray

    reflim = (max(xlimits[1], ylimits[1]), min(xlimits[2], ylimits[2]))
    reflim[1] < reflim[2] && lines!(a, collect(reflim), collect(reflim), color=(:black, 0.75), linewidth=1.5)
    for (iᶜ,iᵈ,iˢ) in Iterators.product(axes.((fresults,), (3,4,5))...)
        if single_error
            logical_error = combined_results[:,iᶜ,iᵈ,iˢ]
            scatter!(a, phys_errors, logical_error, marker=markers[iˢ], color=plotcolor(iᶜ,iᵈ), markersize=9)
            lines!(  a, phys_errors, logical_error, color=plotcolor(iᶜ,iᵈ), linestyle=linestyles[iᵈ], linewidth=2.6)
        else
            scatter!(a, phys_errors, fresults[:,1,iᶜ,iᵈ,iˢ], marker=:+, color=plotcolor(iᶜ,iᵈ), markersize=9)
            scatter!(a, phys_errors, fresults[:,2,iᶜ,iᵈ,iˢ], marker=:x, color=plotcolor(iᶜ,iᵈ), markersize=9)
            lines!(  a, phys_errors, fresults[:,1,iᶜ,iᵈ,iˢ], color=plotcolor(iᶜ,iᵈ), linestyle=linestyles[iᵈ], linewidth=2.6)
            lines!(  a, phys_errors, fresults[:,2,iᶜ,iᵈ,iˢ], color=plotcolor(iᶜ,iᵈ), linestyle=linestyles[iᵈ], linewidth=2.6)
        end
    end
    ca = []
    for (iᶜ,label) in enumerate(codelabels)
        push!(ca, lines!(a, [NaN], [NaN], color=plotcolor(iᶜ,1), label=label))
    end
    Legend(f[1:2,7],ca,codelabels, "Code", framevisible = false, halign=:left, titlehalign=:left, valign=:top, nbanks=2, titlesize=24, labelsize=18)
    la = []
    for (iᵈ,label) in enumerate(decoderlabels)
        push!(la, lines!(a, [NaN], [NaN], linestyle=linestyles[iᵈ], color=decoderlegendcolor(iᵈ), label=label))
    end
    Legend(f[3:6,7],la,decoderlabels, "Decoder", framevisible = false, halign=:left, titlehalign=:left, valign=:top, nbanks=1, titlesize=24, labelsize=18)
    ma = []
    if single_error
        for (iˢ,label) in enumerate(setuplabels)
            push!(ma, scatter!(a, [NaN], [NaN], marker=markers[iˢ], color=:gray, label=label))
        end
        Legend(f[7,7],ma,setuplabels, "Circuit Type", framevisible = false, halign=:left, titlehalign=:left, valign=:top, nbanks=2, titlesize=24, labelsize=18)
    else
        push!(ma, scatter!(a, [NaN], [NaN], marker=:+, color=:gray, label="X"))
        push!(ma, scatter!(a, [NaN], [NaN], marker=:x, color=:gray, label="Z"))
        Legend(f[7,7],ma,["X", "Z"], "Logical Error", framevisible = false, halign=:left, titlehalign=:left, valign=:top, nbanks=2, titlesize=24, labelsize=18)
    end
    f
end

function prep_figures(code_metadata)
    #Threads.@threads :greedy for (codetype, metadata) in code_metadata
    for (codetype, metadata) in code_metadata
        codetypename = typenameof(codetype)
        @info "Plotting figures for $(codetypename) ..."
        codes = [codetype(instance_args...) for instance_args in metadata[:family]]
        decoders = metadata[:decoders]
        setups = metadata[:setups]
        errrange = metadata[:errrange]
        e = logrange(errrange...)
        r = dbnarray(codes, decoders, setups, e)

        single_error = length(decoders)>1 || decoders==[TableDecoder]

        # Plotting benchmarking summary fig
        f = make_decoder_figure(e, r;
            title="$(codetypename)",
            codelabels=instancenameof.(codes),
            decoderlabels=skipredundantfix.(decoders),
            setuplabels=skipredundantfix.(setups),
            single_error
        )
        save("codes/$(codetypename)/totalsummary.png", f)

        # Plotting code instances
        for c in codes
            # Plotting stabilizer
            f = Figure(size=(600,300))
            sf = f[1:2,1]
            _,ax,p = stabilizerplot_axis(sf, parity_checks(c))
            ax.title = "Parity Check Tableau\n(a.k.a. Stabilizer Generators)"
            cm = Makie.cgrad([:lightgray,:black], 2, categorical = true)
            hz, hx, tz, tx = if iscss(c) === true # now 'nothing' type will be recognized as non-css
                parity_matrix_z(c)[end:-1:1,:]', parity_matrix_x(c)[end:-1:1,:]', "Z parity checks", "X parity checks"
            else
                h = stab_to_gf2(parity_checks(c))
                h[:,end÷2+1:end][end:-1:1,:]', h[:,1:end÷2][end:-1:1,:]', "Z components", "X components"
            end
            axz = Axis(f[1,2], title=tz)
            hmz = Makie.heatmap!(
                axz,
                collect(hz);
                colorrange = (0, 1),
                colormap=cm
            )
            axx = Axis(f[2,2], title=tx)
            hmx = Makie.heatmap!(
                axx,
                collect(hx);
                colorrange = (0, 1),
                colormap=cm
            )
            linkxaxes!(axx, axz)
            for ax in (axx, axz)
                Makie.hidedecorations!(ax)
                Makie.hidespines!(ax)
                ax.aspect = Makie.DataAspect()
            end
            colsize!(f.layout, 1, Relative(0.6))
            colsize!(f.layout, 2, Aspect(1, 1))
            colgap!(f.layout, 1, Relative(0.15))
            save("codes/$(codetypename)/$(instancenameof(c)).png", f)
            # Plotting circuits
            if code_n(c) <= 10
                try
                    savecircuit(naive_encoding_circuit(c), "codes/$(codetypename)/$(instancenameof(c))_encoding.png")
                catch
                    @error "$(c) failed to plot `naive_encoding_circuit`"
                end
                try
                    savecircuit(naive_syndrome_circuit(c)[1], "codes/$(codetypename)/$(instancenameof(c))_naive_syndrome.png")
                catch
                    @error "$(c) failed to plot `naive_syndrome_circuit`"
                end
                try
                    error("shor syndrome circuit plotting is problematic, fix it") #TODO
                    savecircuit(vcat(shor_syndrome_circuit(c)[1:2]...), "codes/$(codetypename)/$(instancenameof(c))_shor_syndrome.png")
                catch
                    @error "$(c) failed to plot `shor_syndrome_circuit`"
                end
            end
        end
    end
end

end
