using CSV
using DataFrames
using Plots
using Distances
using Random
using Distributions
using StatsBase
function decision(probability)
    return rand() < probability
end

# Data from https://www.science.org/doi/full/10.1126/science.1253621

occupancies = CSV.File("Data/Plantago/occupancies.csv") |> DataFrame
pops = unique(occupancies[!,"ID"])
years = unique(occupancies[!,"Year"])

# distanceMatrix = zeros(length(pops), length(pops))
# for i in eachindex(pops)
#     pop1 = pops[i]
#     for j in eachindex(pops)
#         pop2 = pops[j]
#         distanceMatrix[i,j] = hypot((occupancies[occupancies.ID .== pop1, "Longitude"][1] - occupancies[occupancies.ID .== pop2, "Longitude"][1]),
#                                     (occupancies[occupancies.ID .== pop1, "Latitude"][1] - occupancies[occupancies.ID .== pop2, "Latitude"][1]))
#     end
#     println(i, "out of", 4367)
# end
# CSV.write("Data/Plantago/distanceMatrix.csv", DataFrame(distanceMatrix, :auto))

distanceMatrix = CSV.File("Data/Plantago/distanceMatrix.csv") |> DataFrame

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

plot(occupancyByYear[:,1], occupancyByYear[:,2], xticks = 2000:2:2012)
# YES – this matches the graph in Fig. 1B in the paper. PA is clearly presence/absence.

# Probability of extinction given occupancy ----------------------------------------------------------------------------------------------
actualOccupancies = occupancies[occupancies.PA .== "1",:] # site-years where present
extinctions = CSV.File("Data/Plantago/extinctions.csv") |> DataFrame # site-years where absent after present

PExtinction = length(extinctions[:,1]) / length(actualOccupancies[:,1])

# Epidemiological model ------------------------------------------------------------------------------------------------------------------
# take a subset of populations
popLoss = 0.0 # fraction of populations lost
popSurviving = Int(floor(length(pops)-popLoss*length(pops)))

alpha = 500 # lower alpha = lower dispersal
λ = PExtinction
hosts = popSurviving
timesteps = 100

probabilityMatrix = exp.(-(Matrix(distanceMatrix) / alpha)) # apply element-wise!!!
idx = sample(axes(probabilityMatrix,1),popSurviving, replace=false)
probabilityMatrix = probabilityMatrix[idx,idx]

A = sum(probabilityMatrix) / popSurviving
probabilityMatrix = probabilityMatrix / A
sum(probabilityMatrix) / popSurviving # should be 1

infected = zeros(hosts)
susceptible = ones(hosts)
occupancy_2001 = tryparse.(Int,occupancies[occupancies.Year .== 2001,:].PA)
replace!(x -> isnothing(x) ? 0 : x, occupancy_2001) 
infected = occupancy_2001[idx]


infectionRecord = zeros((timesteps, hosts))
totalInfected = zeros(timesteps)

# t = 1
for t in 1:timesteps
    infectionRecord[t,1:hosts] = infected
    totalInfected[t] = sum(infected)

    recovery = decision.(rand(Poisson(λ), hosts)) # disease goes extinct in a population with probability λ per unit time
    susceptible = decision.(susceptible + infected .* recovery)

    infected = decision.(probabilityMatrix*infected + infected - 100*recovery)
end

plot(1:timesteps, totalInfected)
plot!(xlabel = "time", ylabel = "Infected Populations")
savefig("Figures/Plantago model single run.png") 

# Iterate for different levels of population loss ---------------------------------------------------------------------------------------
model = function(popLoss)
    popSurviving = Int(floor(length(pops)-popLoss*length(pops)))

    alpha = 500 # lower alpha = lower dispersal
    λ = PExtinction
    hosts = popSurviving
    timesteps = 100

    probabilityMatrix = exp.(-(Matrix(distanceMatrix) / alpha)) # apply element-wise!!!
    idx = sample(axes(probabilityMatrix,1),popSurviving, replace=false)
    probabilityMatrix = probabilityMatrix[idx,idx]

    A = sum(probabilityMatrix) / popSurviving
    probabilityMatrix = probabilityMatrix / A
    sum(probabilityMatrix) / popSurviving # should be 1

    infected = zeros(hosts)
    susceptible = ones(hosts)
    occupancy_2001 = tryparse.(Int,occupancies[occupancies.Year .== 2001,:].PA)
    replace!(x -> isnothing(x) ? 0 : x, occupancy_2001) 
    infected = occupancy_2001[idx]


    infectionRecord = zeros((timesteps, hosts))
    totalInfected = zeros(timesteps)

    # t = 1
    for t in 1:timesteps
        infectionRecord[t,1:hosts] = infected
        totalInfected[t] = sum(infected)

        recovery = decision.(rand(Poisson(λ), hosts)) # disease goes extinct in a population with probability λ per unit time
        susceptible = decision.(susceptible + infected .* recovery)

        infected = decision.(probabilityMatrix*infected + infected - 100*recovery)
    end

    return Int(totalInfected[timesteps] > 0)
end

model(0.8)

samples = 10
losses = 0:0.05:0.95

extinctionSamples = zeros(length(losses), samples)

for e in eachindex(losses)
    popLoss = losses[e]
    for sample in 1:samples
        extinctionSamples[e, sample] = model(popLoss)
        println(popLoss, " ", sample)
    end
end

using CSV, Tables
CSV.write("Data/extinctionSamples2.csv", Tables.table(extinctionSamples))
