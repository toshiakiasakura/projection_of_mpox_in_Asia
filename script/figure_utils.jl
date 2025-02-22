## Functions for preparing figures.

"""Visualise an epidemic curve in Japan from 2022 to 2023
stratified by travel history.
"""
function Japan_weekly_epicurve_2022_2023(df::DataFrame)
	# Obtain from 2022-07-25
	tab_7d = convert_to_weekly_case(df, med_date)

	cond = coalesce.(df[:, "Travel history"], "None") .!= "None"
	tab_7d_travel = convert_to_weekly_case(df[cond, :], med_date)

	tab_vis = leftjoin(tab_7d, tab_7d_travel, on = :timestamp, makeunique = true)
	tab_vis[:, :count_sum_1] = coalesce.(tab_vis[:, :count_sum_1], 0)
	tab_vis[:, :no_history] = tab_vis[:, :count_sum] .- tab_vis[:, :count_sum_1]
	date = tab_vis[:, :timestamp]
	cnts_no_history = tab_vis[:, :count_sum]
	cnts_history = tab_vis[:, :count_sum_1]

	xtick_label = [Date(2022, 7, 25), Date(2022, 10, 10), Date(2023, 1, 2), Date(2023, 3, 27), Date(2023, 6, 19)]
	pl = bar(date, cnts_no_history, label = "No travel history")
	bar!(pl, date, cnts_history, label = "Travel history")
	plot!(pl,
		label = :none,
		dpi = 300,
		fmt = :png,
		xtickfontsize = 10,
		ytickfontsize = 10,
		guidefontsize = 14,
		xlabel = "Medical attendance date",
		ylabel = "Weekly incidence",
		foreground_color_legend = nothing,
		background_color_legend = nothing,
		left_margin = 15Plots.pt,
		bottom_margin = 15Plots.pt,
		size = (800, 400),
		xticks = (xtick_label, xtick_label),
		ylim = [0, 20],
	)
	# vline!(pl, [Date("2023-04-15")]) # For behavioural change timing.
	savefig(pl, "../fig/epi_curve_JPN.png")
	display(pl)
end

function visualise_sexual_partner_distribution(ks_4w, Pk_4w, ks_1y, Pk_1y)
	x = 1:50
	rate_kwds = (yaxis = :log10, label = :none,
		xlim = [0.5, 51], ylim = [1e-3, 35],
		marker = (:circle, 3),
		# fillrange=1e-3,
		yticks = ([1e-2, 1e-1, 1.0, 10], ["0.01", "0.1", "1.0", "10"]),
		bar_edges = false, bar_width = 1.0, lw = 0.3,
		xlabel = "Sexual activity group",
	)
	prob_bar_kwds = (yaxis = :log10, label = :none, ylim = [1e-9, 0.5],
		fillrange = 1e-9,
		bar_width = 1.0, lw = 0.3,
		ylabel = "Proportion of \neach sexual activity group",
		xlabel = "Sexual activity group",
		xlim = [0.5, 50.5],
	)

	pl1 = plot(x, ks_4w / 365;
		ylabel = "Rate of new sexual partnership \nformulations per day",
		rate_kwds...)
	annotate!(pl1, 5, 10, text("A", :black, :left, 18))

	pl2 = bar(x, Pk_4w; #title="Natsal 4 week",
		prob_bar_kwds...)
	annotate!(pl2, 40, 8 * 1e-2, text("B", :black, :left, 18))

	pl3 = bar(x, Pk_1y; #title="Natsal 1 year",
		prob_bar_kwds...)
	annotate!(pl3, 40, 8 * 1e-2, text("C", :black, :left, 18))
	#pl4 = bar(x, ks_1y/365; rate_kwds...)
	layout = @layout [a b c]
	pl = plot(pl1, pl2, pl3, layout = layout,
		dpi = 300, size = (1000, 300),
		left_margin = 10Plots.mm,
		bottom_margin = 7Plots.mm,
	)
	savefig(pl, "../fig/sexual_distribution.png")
	display(pl)

end
