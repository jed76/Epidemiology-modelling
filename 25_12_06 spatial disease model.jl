function decision(probability)
    return rand() < probability
end


hosts = 100 # number of hosts
fieldsize = 10 # side length of field
timesteps = 500

infected = zeros(hosts)

hostCoords = rand(Float64, (hosts,2)) * fieldsize

distanceMatrix = zeros(hosts,hosts)

for i in 1:hosts
    iCoords = hostCoords[i,1:2]
    for j in 1:hosts
        jCoords = hostCoords[j,1:2]

        distanceMatrix[i,j] = sqrt((iCoords[1]-jCoords[1])^2 + (iCoords[2]-jCoords[2])^2)
    end
end

alpha = 0.3

probabilityMatrix = exp.(-(distanceMatrix / alpha)) # apply element-wise!!!

A = sum(probabilityMatrix) / hosts
probabilityMatrix = probabilityMatrix / A
sum(probabilityMatrix) / hosts # should be 1

infected = zeros(hosts)
infected[1] = 1

infectionRecord = zeros((timesteps, hosts))
totalInfected = zeros(timesteps)

for t in 1:timesteps
    infectionRecord[t,1:hosts] = infected
    totalInfected[t] = sum(infected)

    infected = decision.(probabilityMatrix*infected + infected)
end

using Plots
plot(1:timesteps, totalInfected)
plot!(xlabel = "time", ylabel = "% infected")
