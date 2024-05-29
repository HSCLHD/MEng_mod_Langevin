include("Langevin_validation.jl")
using CSV, DataFrames

using LsqFit
using PyPlot

function best_fit(x, y, p)
    model = p[1].- p[1] * exp.(x./p[2])

    fit = curve_fit(model, x, y, p)

    # best_params = fit.param

    # # Define a function to evaluate the fitted exponential
    # fitted_function(x) = best_params[1] * exp(best_params[2] * x)
    return fit.param
end



# y_fit = 0.00025 .- 0.00025* exp.(-x_data/12.0)
# size(y_data)
# limit = y_data[500]
# y_scaled = zeros(999,1)

# for i in 1:999
#     y = y_data[i]
#     if y < limit
#         y_scaled[i,1] =  1 .- y/limit
#     else
#         y_scaled[i,1] = 0.000001
#     end
# end


# log.(y_scaled)
# tau_m = -x_data[1:999]./log.(y_scaled)

# tau_m, limit = find_tau(x_data[1:999], y_data,1000)

function find_tau_opt(x, y, tau_list)
    y_limit= y[999]
    # dif = zeros(999,1)
    # for in in 1:999
    # cross_list = Float64[]
    cross_list = zeros(length(tau_list),1)
    for i in 1:length(tau_list)
        y_model = y_limit .- y_limit * exp.(-x./tau_list[i])
        cross_list[i,:] = crosscor(y[2:end],y_model[2:end],[0])
    end
    tau_opt = tau_list[argmax(cross_list)]
        # dif[i] = y_data 
    # diff = y_data - y_model
    return cross_list,tau_opt
end


y_data = Matrix(CSV.read("/Users/hanshichen/Desktop/engineering/SoftActiveCells/sims_x/positions_non_averaged_langevin/positions_maxtime=1000,dt=0.1,gamma=5.0,kBT=0.25,k=10.0,dev=0.1.csv", DataFrame))
y_data = sigmasquare(1000, 256, y_data)
y_vec = vec(y_data)
x_data = collect(1:999)
limit = y_data[999]
tau_test = 2:0.5:50
cross_list,tau = find_tau_opt(x_data,y_vec,limit,tau_test)
y_model_test = limit .- limit * exp.(-x_data./7.5)
cross_list
test_cor = crosscor(y_vec[2:end],y_model_test[2:end],[100,101,102,200,300,400,500])

test_set = zeros((2,1))
test_set[1,:] = test_cor
# fit = best_fit(x_data[1:999],y_data,[limit,20])
# y_fit = fit[1] - fit[1] * exp.(x_data/fit[2])

# PyPlot.plot(x_data[1:999], tau_m)
# y_fit = limit - limit*exp.(-x_data/tau_m)
PyPlot.plot(x_data, y_model_test)
PyPlot.plot(y_data)
PyPlot.xlabel("time lag")
PyPlot.ylabel("σ²")
PyPlot.legend()