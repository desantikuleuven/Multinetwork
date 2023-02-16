using PowerPlots
using Setfield
using DataFrames
using DataStructures
using VegaLite, Vega, Query
a = Dict( 83 => [650, 797, 672, 639, 670, 735, 731, 746], 
159 => [650, 756, 668, 631, 654, 668, 659, 687],
85 => [650, 756, 650, 631, 644, 639, 659, 670],
81 => [382, 377, 406, 434, 424, 408, 356, 378],
82 =>[650, 628, 609, 588, 591, 593, 603, 625],
)
d = DataFrame()
for (branch, loading) in a
    if isempty(d)
        d[!,"Red_hours"] = loading
        d[!,"type"] = repeat([branch], length(loading))
        d[!,:repr_day_number] = [12, 15, 20, 23, 25, 28, 30, 35]
    else
        append!(d[!,"Red_hours"],  loading)
        append!(d[!,"type"] , repeat([branch], length(loading)))
        append!(d[!,:repr_day_number] , [12, 15, 20, 23, 25, 28, 30, 35])
    end
end

d |> @vlplot(
    :line, 
    height = 800, width = 800,
    y = "Red_hours:q",
    x = "repr_day_number",
    color = "type:n"
)