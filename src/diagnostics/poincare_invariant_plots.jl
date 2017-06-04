
if Pkg.installed("PyPlot") != nothing
    using PyPlot
end


function plot_integral_error(t, I, filename; plot_title=L"\Delta I")
    if Pkg.installed("PyPlot") != nothing
        function power10ticks(x, pos)
            if x == 0
                return "\$ 0 \$"
            end

            # xpower      = log10(x)
            xpower      = log10(t[end])
            exponent    = @sprintf("%2d",   floor(Int64, xpower))
            coefficient = @sprintf("%2.1f", x / (10^floor(xpower)))

            return "\$ $coefficient \\times 10^\{ $exponent \} \$"
        end

        xf = matplotlib[:ticker][:FuncFormatter](power10ticks)

        yf = matplotlib[:ticker][:ScalarFormatter]()
        yf[:set_powerlimits]((-1,+1))
        yf[:set_scientific](true)
        yf[:set_useOffset](true)

        I_error = (I - I[1]) / I[1]

        fig = figure(figsize=(6,5))
        subplots_adjust(left=0.15, right=0.93, top=0.94, bottom=0.12)
        plot(t, I_error, ".", markersize=2)
        xticks(linspace(t[1], t[end], 6))
        xlim(t[1], t[end])
        xlabel(L"$t$")
        ylabel(plot_title)
        ax = gca()
        ax[:xaxis][:set_major_formatter](xf)
        ax[:yaxis][:set_major_formatter](yf)
        ax[:yaxis][:set_label_coords](-0.1,0.5)
        savefig(filename)
        close(fig)
    end
end



function plot_loop(sol, nplot, filename, dpi=100)
    if Pkg.installed("PyPlot") != nothing
        fig = figure(figsize=(10,8))
        subplots_adjust(left=0.0, right=0.7, top=1.0, bottom=0.0)
        ax  = gca(projection="3d")

        for i in 1:sol.nt+1
            if mod(i-1, div(sol.nt, nplot)) == 0
                tpower      = log10(sol.t.t[end])
                exponent    = @sprintf("%2d", tpower)
                coefficient = @sprintf("%2.2f", sol.t.t[i] / 10^tpower)
                tstr = "$coefficient \\times 10^\{ $exponent \}"
                plot3D(sol.q.d[1,i,:], sol.q.d[2,i,:], sol.q.d[3,i,:], label="\$ t = $tstr \$")
            end
        end

        xlabel(L"$x$", fontsize=18, labelpad=10)
        ylabel(L"$y$", fontsize=18, labelpad=10)
        zlabel(L"$z$", fontsize=18, labelpad=10)

        xlim(-0.6, +0.6)
        ylim(-0.6, +0.6)

        legend(loc=2, bbox_to_anchor=(1.05, 0.8))
        savefig(filename, dpi=dpi)
        close(fig)
    end
end



function plot_surface(sol, nplot, filename, dpi=100)
    if Pkg.installed("PyPlot") != nothing && nplot > 0
        fig = figure(figsize=(10,8))
        subplots_adjust(left=0.0, right=0.7, top=1.0, bottom=0.0)
        ax  = gca(projection="3d")

        nplot > sol.nt ? nplot = sol.nt : nothing

        for i in 1:sol.nt+1
            if mod(i-1, div(sol.nt, nplot)) == 0
                tpower      = log10(sol.t.t[end])
                exponent    = @sprintf("%2d", tpower)
                coefficient = @sprintf("%2.2f", sol.t.t[i] / 10^tpower)
                tstr = "$coefficient \\times 10^\{ $exponent \}"
                plot3D(sol.q.d[1,i,:], sol.q.d[2,i,:], sol.q.d[3,i,:], ".", label="\$ t = $tstr \$")
           end
        end

        xlabel(L"$x$", fontsize=18, labelpad=10)
        ylabel(L"$y$", fontsize=18, labelpad=10)
        zlabel(L"$z$", fontsize=18, labelpad=10)

        xlim(-0.15, +0.15)
        ylim(-0.15, +0.15)

        legend(loc=2, bbox_to_anchor=(1.05, 0.8))
        savefig(filename, dpi=dpi)
        close(fig)
    end
end
