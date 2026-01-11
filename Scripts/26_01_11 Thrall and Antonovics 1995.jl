function pops((X_1,X_2,Y,N))
    # params from paper
    λ_1 = 2.0
    λ_2 = 1.75
    d = 0.5
    β_1 = 0.658
    β_2 = 0.202
    γ = 0.5

    b_1 = λ_1 / (γ*N + 1)
    b_2 = λ_2 / (γ*N + 1)
    X_1t = X_1*(b_1+(1-β_1*(Y/N)*(1-d)))      # healthy susceptible
    X_2t = X_2*(b_2+(1-β_2*(Y/N)*(1-d)))      # healthy resistant
    Y_t = Y*(1+β_1*(X_1/N)+β_2*(X_2/N)*(1-d)) # infected
    N_t = X_1t + X_2t + Y_t                   # total population

    return((X_1t, X_2t, Y_t, N_t))
end

# NB need to implement variable γ

#############
# RUN MODEL #
#############

# timesteps = 100

# popsizes = zeros(timesteps, 4)
# popsizes[1,:] .= (50, 50, 50, 150)

# for t in 2:timesteps
#     popsizes[t,:] .= pops(popsizes[t-1,:])
#     for i in 1:4
#         if popsizes[t,i] < 1
#             popsizes[t,i] = 0
#         end
#     end
# end

# using Plots
# plot(1:timesteps, popsizes[:,1])
# plot!(1:timesteps, popsizes[:,2])
# plot!(1:timesteps, popsizes[:,3])
# plot!(1:timesteps, popsizes[:,4])

function dispersal_probability(i)
    # params from paper
    α = 1
    θ = 1

    return exp(-(i/α)^θ) - exp(-((i+1)/α)^θ)
end

# dispersal probability matrix ---------------------
fieldSize = 10
distanceMatrix = zeros(fieldSize,fieldSize)
for i in 1:fieldSize
    for j in 1:fieldSize
        distanceMatrix[i,j] = hypot((hostCoords[i, 1] - hostCoords[j, 1]), (hostCoords[i, 2] - hostCoords[j,2]))
    end
end

probabilityMatrix = exp.(-(distanceMatrix / alpha)) # apply element-wise!!!

example = zeros(10,10)
example[1,1] = 1
example[1,3] = 4
example[3,2] = 7

x = vec(example)
y = reshape(x, (10, 10))
