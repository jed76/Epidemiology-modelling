function decision(probability)
    return rand() < probability
end


hostDensity = 1 # number of hosts per unit area
fieldsize = 10 # side length of field
timesteps = 500
λ = 0.05 # the number of infected hosts that recover per time step
λ2 = 0.02 # number of recovered hosts becoming susceptible again per time step
alpha = 0.7 # lower alpha = lower dispersal

hosts = hostDensity * fieldsize^2
infected = zeros(hosts)

hostCoords = rand(Float64, (hosts,2)) * fieldsize

distanceMatrix = zeros(hosts,hosts)

for i in 1:hosts
    for j in 1:hosts
        distanceMatrix[i,j] = hypot((hostCoords[i, 1] - hostCoords[j, 1]), (hostCoords[i, 2] - hostCoords[j,2]))
    end
end

probabilityMatrix = exp.(-(distanceMatrix / alpha)) # apply element-wise!!!

A = sum(probabilityMatrix) / hosts
probabilityMatrix = probabilityMatrix / A
sum(probabilityMatrix) / hosts # should be 1

infected = zeros(hosts)
recovered = zeros(hosts)
infected[1] = 1

infectionRecord = zeros((timesteps, hosts))
totalInfected = zeros(timesteps)
totalRecovered = zeros(timesteps)

for t in 1:timesteps
    infectionRecord[t,1:hosts] = infected
    totalInfected[t] = sum(infected)
    totalRecovered[t] = sum(recovered)

    recovery = decision.(rand(Poisson(λ), hosts)) # infected hosts recover with probability λ per unit time
    recovered = decision.(recovered + infected .* recovery)

    infected = decision.(probabilityMatrix*infected + infected - 100*recovered)

    recovered = decision.(recovered - rand(Poisson(λ2), hosts))
end

totalRecovered[timesteps] == hosts # have all hosts been infected?

using Plots
plot(1:timesteps, totalInfected)
plot!(1:timesteps, totalRecovered)
plot!(xlabel = "time", ylabel = "number in box")
savefig("Figures/basic spatial model.png") 

