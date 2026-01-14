function pops((X_1,X_2,Y,N))
    if N == 0
        return((0,0,0,0))
    else
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

# build kernel ------------------------------------------
fieldSize = 10
kernelWidth = 2 * fieldSize - 1
centre = fieldSize
m = 0.05 # of seeds disperse

K = zeros(kernelWidth, kernelWidth)
for i in 1:kernelWidth
    for j in 1:kernelWidth
        K[i,j] = dispersal_probability(hypot(i - centre, j - centre))
    end
end

K[centre, centre] = 0
K = m * K ./ sum(K)
K = m * K ./ sum(K)
sum(K) # should equal m
K[centre, centre] = 1 - m
sum(K) # should equal 1

# we want an array of the same dimensions as our field (10×10)
# where each cell contains the dispersal kernel (also 10×10)
# for its corresponding cell in the field

# cut kernel to size for a given cell
function kernelCut(i,j,K)
    width = Int((length(K[1,:]) + 1) / 2)

    return K[(width-i+1):(2*width-i), (width-j+1):(2*width-j)]
end

# Kernel array
KArray = Array{Any}(undef, fieldSize, fieldSize, fieldSize, fieldSize)
for i in 1:fieldSize
    for j in 1:fieldSize
        KArray[i,j,:,:] = kernelCut(i,j,K)
    end
end


# Run model
time = 1000


field = zeros(fieldSize, fieldSize, time, 3)
field[1,1,1,1] = 75
field[1,1,1,2] = 75
field[1,1,1,3] = 50
dispersalX1 = zeros(fieldSize, fieldSize, fieldSize, fieldSize)
dispersalX2 = zeros(fieldSize, fieldSize, fieldSize, fieldSize)
dispersalY = zeros(fieldSize, fieldSize, fieldSize, fieldSize)

for t in 1:time
    for i in 1:fieldSize
        for j in 1:fieldSize
            # i = 1
            # j = 1
            # print(i, j)
            kernel = KArray[i,j,:,:]
            if t > 1
                X1 = field[i,j,t-1,1]
                X2 = field[i,j,t-1,2]
                Y = field[i,j,t-1,3]
            else
                X1 = field[i,j,1,1]
                X2 = field[i,j,1,2]
                Y = field[i,j,1,3]
            end
            N = X1 + X2 + Y

            popst = pops((X1, X2, Y, N))
            X1 = popst[1]
            X2 = popst[2]
            Y = popst[3]

            dispersalX1[i,j,:,:] = X1 * kernel
            dispersalX2[i,j,:,:] = X2 * kernel
            dispersalY[i,j,:,:] = Y * kernel
        end
    end
    for i in 1:fieldSize
        for j in 1:fieldSize
            field[i,j,t,1] = sum(dispersalX1[:,:,i,j])
            field[i,j,t,2] = sum(dispersalX2[:,:,i,j])
            field[i,j,t,3] = sum(dispersalY[:,:,i,j])
        end
    end
end

using Plots
plot(field[1,1,:,:])

# TO DO:
# add extinction
# if population in a cell < 1, round down to 0

# KArray[1,4,:,:] # kernel for [1,4]


# for i in 1:fieldSize
#     for j in 1:fieldSize
#         distanceMatrix[i,j] = hypot((hostCoords[i, 1] - hostCoords[j, 1]), (hostCoords[i, 2] - hostCoords[j,2]))
#     end
# end

# probabilityMatrix = exp.(-(distanceMatrix / alpha)) # apply element-wise!!!

# example = zeros(10,10)
# example[1,1] = 1
# example[1,3] = 4
# example[3,2] = 7

# x = vec(example)
# y = reshape(x, (10, 10))
