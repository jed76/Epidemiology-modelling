using CSV
using DataFrames
using Plots
using Distances




plantagoPops = CSV.File("Data/Plantago/abundance.csv") |> DataFrame
pops = length(plantagoPops[!,"PatchID"])

plot(plantagoPops[!,"Longitude"],plantagoPops[!,"Latitude"], seriestype=:scatter)

distanceMatrix = zeros(pops,pops)

# !!! need to calculate great circle distance for latlong data
# from: https://discourse.julialang.org/t/geodesy-how-to-calculate-the-straight-line-distance-between-two-locations-which-is-represented-by-longitude-and-latitude/19984/6?u=beroid
# using Distances
# julia> l1 = (-27.468937, 153.023628)
# julia> l2 = (-27.465933, 153.025900)
# julia> haversine(l1, l2, 6372.8)
# 0.39054922275889736
# our data appears to be great circle distance in metres from equator and date line (!)

for i in eachindex(plantagoPops[!,"PatchID"])
    for j in eachindex(plantagoPops[!,"PatchID"])
        # distanceMatrix[i,j] = haversine((plantagoPops[i, "Longitude"], plantagoPops[i, "Latitude"]),
        #                                 (plantagoPops[j, "Longitude"], plantagoPops[j, "Latitude"]),
        #                                 6371) # in metres
        distanceMatrix[i,j] = hypot((plantagoPops[i, "Longitude"] - plantagoPops[j, "Longitude"]),
                                    (plantagoPops[i, "Latitude"] - plantagoPops[j, "Latitude"]))
    end
end
# I think this gives the distances in metres

plot(plantagoPops[!,"Longitude"],plantagoPops[!,"Latitude"], seriestype=:scatter, color = plantagoPops[!,"Abundance_2012"])
histogram(distanceMatrix[:])

occupancies = CSV.File("Data/Plantago/occupancies.csv") |> DataFrame
years = unique(occupancies[!,"Year"])

occupancyByYear = zeros(length(years), 2)

for i in eachindex(years)
    year = years[i]
    occupancy_year = tryparse.(Int,occupancies[occupancies.Year .== year,:].PA)
    replace!(x -> isnothing(x) ? 0 : x, occupancy_year) # note we are converting missing data to 0s here

    total = sum(occupancy_year)
    println(year, " ", total)
    occupancyByYear[i,1] = Int(year)
    occupancyByYear[i,2] = Int(total)
end

plot(occupancyByYear[:,1], occupancyByYear[:,2])
# YES – this matches the graph in the paper. PA is clearly presence/absence.

# Putting data in an easy to use format -------------------------------------------------------------------------------------------
# CSV.write("Data/Plantago/distanceMatrix.csv", DataFrame(distanceMatrix, :auto)) 



