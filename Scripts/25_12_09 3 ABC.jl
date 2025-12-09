using CSV
using DataFrames
using Plots
using Random
using Distributions
using StatsBase
function decision(probability)
    return rand() < probability
end

occupancies = CSV.File("Data/Plantago/occupancies.csv") |> DataFrame
pops = unique(occupancies[!,"ID"])
years = unique(occupancies[!,"Year"])

# Probability of extinction given occupancy
actualOccupancies = occupancies[occupancies.PA .== "1",:] # site-years where present
extinctions = CSV.File("Data/Plantago/extinctions.csv") |> DataFrame # site-years where absent after present
PExtinction = length(extinctions[:,1]) / length(actualOccupancies[:,1])

distanceMatrix = CSV.File("Data/Plantago/distanceMatrix.csv") |> DataFrame

model = function(alpha)
    popLoss = 0

    popSurviving = Int(floor(length(pops)-popLoss*length(pops)))

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

    return totalInfected
end

model(100)

plot(1:100, model(500))
