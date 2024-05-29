using Pkg
# Pkg.add("Conda")
# using Conda
# Conda.add("nomkl")
# Conda.add("numpy")
# Conda.add("scipy")
# Pkg.add("PyCall")
# Pkg.add("PyPlot")
Pkg.add("LsqFit")
Pkg.add("CurveFit")
using LsqFit
using CurveFit
using PyPlot
using StatsBase, Distributions
using Statistics
using DelimitedFiles
using CSV, DataFrames
using Random
pygui(true)
# BAOAB(maxtime, fdt, fgamma,fkBT, fNcells, fk, fdev)

script_dir = @__FILE__
sac_dir = dirname(dirname(dirname(script_dir)))
# source_dir = abspath(joinpath(sac_dir, "validation_code_position_values"))
source_dir = sac_dir
# Read all combinations
# all_combs = Matrix(CSV.read(abspath(joinpath(source_dir,"combinations.csv")), DataFrame))
# all_combs = Matrix(CSV.read(abspath(joinpath(source_dir,"combinations_sigmasquared_check.csv")), DataFrame))


my_Ncells = 256
"""
Autocorrelationi function, hard to find the fitting curve
"""

"""
Following is the function to validate and characterize the langevine dynamics
strength B of the random noise or fluctuating force to the magni- tude of the friction or dissipation.
Need to gather more data for validating the final time. need variation of gamma, k and 
"""

# the standard deviation squared reveals how dispersed the data is. 
# The autocorrelatino and the areasquared shows how far off before the dynamics lost correlation. 
function sigmasquare(T::Int64,Ncells::Int64,pos::Matrix{Float64})
    
    std_target = zeros(Float64,(T-1, Ncells))
    
    for i in 1:T-1
        for j in 1:Ncells
            std_i = std(pos[1:i,j]) 
            std_target[i,j] = std_i
        end
    end
    std_target = mean(std_target.^2,dims=2)
    return std_target
end

function autovar(T::Int64,pos::Matrix{Float64})
    autocor_targ = mean(autocor(pos, 1:T-1),dims = 2)
    # autocor_targ = autocor(pos, 1:T-1)
    return autocor_targ
end

function areasquarevar(dt::Float64,T::Int64,pos::Matrix{Float64})
    N = Int(T/dt)
    dA_squared = zeros(Float64,(N, 256))

    for i in 1:N-1
        dA = pos[1:N-1-i,:] - pos[1+i:N-1,:]
        dA_squared[i,:] = mean(dA.^2,dims= 1)
        # dA_squared[i,:] = dA.^2
    end
    dA_squared = mean(dA_squared,dims  = 2)

    return dA_squared
end

# define the model
function best_fit(x, y, p)
    model = p[1].- p[1] * exp.(x./p[2])

    fit = curve_fit(model, x, y, p)

    # best_params = fit.param

    # # Define a function to evaluate the fitted exponential
    # fitted_function(x) = best_params[1] * exp(best_params[2] * x)
    return fit.param
end

function find_tau(x,y, maxtime::Int64)
    limit = y[Int(maxtime/2)]
    tau_m = -x./(log.(1 .- y/limit))
    return tau_m, limit
end
#TEST
# Pos = Matrix(CSV.read("/Users/hanshichen/Desktop/engineering/SoftActiveCells/sims_x/positions_non_averaged_langevin/positions_maxtime=1000,dt=0.1,gamma=0.5,kBT=0.25,k=1.0,dev=0.1.csv", DataFrame))
# PyPlot.figure()
# Test_auto = autovar(1000,Pos)
# PyPlot.plot(Test_auto, label = "non average, 1000")

# Pos_averaged = Matrix(CSV.read("/Users/hanshichen/Desktop/engineering/SoftActiveCells/sims_x/positions_from_langevin/positions_maxtime=1000,dt=0.1,gamma=0.5,kBT=0.25,k=1.0,dev=0.1.csv", DataFrame))
# test_ave_auto = autovar(1000,Pos_averaged)
# PyPlot.plot(test_ave_auto, label = "averaged, 1000")

# PyPlot.legend()
# PyPlot.xlabel("t")
# PyPlot.ylabel("Area")

# function plot_time_series(T::Int64,series::Matrix{Float64}, text::String)
#     PyPlot.plot((1:T-1),series, label = text)
#     PyPlot.xlabel("dt")
#     PyPlot.ylabel("dA squared")
# end
# cols = Dict(5.0 => "red",20.0 => "blue",40.0=>"black")

# validate gamma, k for 
# PyPlot.figure()
# for i in 1:107
#     # [maxtime, fdt, fgamma, fkBT, fNcells, fk, fdev] = all_combs[i,:]
    
#     maxtime = floor(Int,all_combs[i,1]) ;fdt =  all_combs[i,2]; fgamma=all_combs[i,3]; fkBT=all_combs[i,4]; fNcells=all_combs[i,5]; fk=all_combs[i,6]; fdev=all_combs[i,7];
#     # println(maxtime) 
    
#     if fgamma == 5 && fdev ==0.1 && fk == 5. && fNcells ==256
    
#     # if maxtime == 100. && fdev ==0.1 && fgamma == 5. && fNcells ==256
#         Pos = Matrix(CSV.read(abspath(joinpath(source_dir,"positions_maxtime=$maxtime,dt=$fdt,gamma=$fgamma,kBT=$fkBT,k=$fk,dev=$fdev.csv")), DataFrame))
#         # plot 
#         PyPlot.plot(autovar(maxtime,Pos),label = "maxtime = $maxtime")

#     end
    
# end

# for i in 1:107
#     # [maxtime, fdt, fgamma, fkBT, fNcells, fk, fdev] = all_combs[i,:]
    
#     maxtime = floor(Int,all_combs[i,1]) ;fdt =  all_combs[i,2]; fgamma=all_combs[i,3]; fkBT=all_combs[i,4]; fNcells=all_combs[i,5]; fk=all_combs[i,6]; fdev=all_combs[i,7];
#     # println(maxtime) 
    
#     # if maxtime == 100 && fdev ==0.1 && fk == 5. && fNcells ==256
    
#     if maxtime == 100. && fdev ==0.1 && fgamma == 5. && fNcells ==256
#         Pos = Matrix(CSV.read(abspath(joinpath(source_dir,"positions_maxtime=$maxtime,dt=$fdt,gamma=$fgamma,kBT=$fkBT,k=$fk,dev=$fdev.csv")), DataFrame))
#         # plot 
#         PyPlot.plot(sigmasquare(maxtime,256,Pos),label = "k = $fk")

#     end
    
# end

# for i in [1,5,10,15,19]
#     # [maxtime, fdt, fgamma, fkBT, fNcells, fk, fdev] = all_combs[i,:]
    
#     maxtime = floor(Int,all_combs[i,1]) ;fdt =  all_combs[i,2]; fgamma=all_combs[i,3]; fkBT=all_combs[i,4]; fNcells=all_combs[i,5]; fk=all_combs[i,6]; fdev=all_combs[i,7];
#     # println(maxtime) 
    
#     # if maxtime == 100 && fdev ==0.1 && fk == 5. && fNcells ==256
    
#     if maxtime == 1000. && fdev ==0.1 && fNcells ==256
#         Pos = Matrix(CSV.read(abspath(joinpath(source_dir,"positions_maxtime=$maxtime,dt=$fdt,gamma=$fgamma,kBT=$fkBT,k=$fk,dev=$fdev.csv")), DataFrame))
# #         # plot 
#         PyPlot.plot(areasquarevar(maxtime,Pos,1),label = "gamma = $fgamma")
# #         # PyPlot.plot(sigmasquare(maxtime,256,Pos),label = "k = $fk")


#     end
    
# end
# Pos = []
# PyPlot.plot(areasquarevar(1000,Pos),label = "average")
# PyPlot.legend()
# # PyPlot.title("Simulation time ")
# PyPlot.xlabel("t")
# PyPlot.ylabel("sigma ^2") 
# PyPlot.xlim(0.0,101)
# maxtime = 100;fdev =0.1; fk = 5.;fdev = 0.1;fkBT= 0.25; fgamma=5.; fdt = 0.1;
# Pos = Matrix(CSV.read(abspath(joinpath(source_dir,"positions_maxtime=$maxtime,dt=$fdt,gamma=$fgamma,kBT=$fkBT,k=$fk,dev=$fdev.csv")), DataFrame))
# # plot 
# PyPlot.plot(sigmasquare(maxtime,256,Pos),label = "gamma = $fgamma")
"""
Plot only the position wrt time
"""
function r_time(T::Int64,dt::Float64, series::Matrix{Float64}, text::String)
    PyPlot.plot(series, label = text)
    PyPlot.xlabel("t")
    PyPlot.ylabel("Area")
end

"""
Find the optimal timescale for the Simulation
"""

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
    return cross_list,tau_opt[1]
end

# Pos =  Matrix(CSV.read("/Users/hanshichen/Desktop/engineering/SoftActiveCells/boltzman2positions_maxtime=1000,dt=0.1,gamma=1.0,kBT=0.25,k=10.0,dev=0.1,validate.csv", DataFrame))
# r_time(1000,0.1,Pos,"one time series")

# xlim(1000,10000)
# ylim(1.01,1.014)


############################################################################

"""
Plot the time evolution of the radii
"""
# Pos = CSV.read("/Users/hanshichen/Desktop/engineering/SoftActiveCells/boltzman2positions_maxtime=1000,dt=0.1,gamma=1.0,kBT=0.25,k=10.0,dev=0.1,validate.csv")
# P = Matrix(P)
# J = Matrix(J)
# PyPlot.plot(mean(P,dims=2), label = "Python")
# PyPlot.plot(mean(J,dims = 2), label = "Julia")
# PyPlot.xlabel("tau")
# PyPlot.ylabel("autocorrelation")
# PyPlot.legend()
# PyPlot.figure()
# PyPlot.plot(P[:,1:3], label = ["py first cell", "py second cell", "py third cell"])

# PyPlot.plot(J[:,1:3], label = ["Ju first cell", "ju second cell", "ju third cell"])

# PyPlot.plot(P[:,1], label = ["py first cell"])

# PyPlot.plot(J[:,1], label = ["Ju first cell"])

# axhline(y = 0)
# PyPlot.xlabel("tau")
# PyPlot.ylabel("autocorrelation")

# PyPlot.legend()

# K = zeros(Float64,(5,997, 256))
# for i in 1:5
#     K[i,:,:] = crosscor(J[:,i],P,1:997)
# end

"""
Crosscorrelation of the radii update, gives high negative values
"""
# K_pos = zeros(Float64,(5,997, 256))
# for i in 1:5
#     K_pos[i,:,:] = crosscor(pos_j[:,i],pos_p,1:997)
# end
# PyPlot.figure()
# PyPlot.plot(K_pos[1,:,1:3], label = ["first cell", "second cell", "third cell"])
# axhline(y = 0)
# PyPlot.xlabel("tau")
# PyPlot.ylabel("crosscorrelation")

# PyPlot.legend()

"""
The following is more correct
sample the (delta A) ^2 

The following is incorrect
Diffusion coefficient searching, find std of t = 50, 100, 250, 500, 999

For the current radii distribution, it is essentially a 1D random walk
D = sigma^2 / (2*t)

"""
# Ncells = 256
# dt = 2
# # timestep
# T = [10,50, 100,250,500,999]
# Dlist_py = zeros(Float64,(length(T), Ncells))
# Dlist_ju = zeros(Float64,(length(T), Ncells))
# stdlist_py = zeros(Float64,(length(T), Ncells))
# stdlist_ju = zeros(Float64,(length(T), Ncells))
# stdlist_py = zeros(Float64,(999, Ncells))
# stdlist_ju = zeros(Float64,(999, Ncells))
# dt = 0.1
# for j in 1:999
#     for i in 1:Ncells
        
# #         t = T[j]
#         std_py = std(pos_p[1:j,i]) 
#         std_ju = std(pos_j[1:j,i])
#         stdlist_ju[j,i] = std_ju
#         stdlist_py[j,i] = std_py
#         # println(D_py)
#         # println(D_ju)
#         # D_ju = std_ju^2/(2*t*dt)
#         # D_py = std_py^2/(2*t*dt)
#         # Dlist_ju[j,i] = D_ju
#         # Dlist_py[j,i] = D_py
#         # push!(Dlist_ju[j,i],D_ju)
#         # push!(Dlist_py[j,i],D_py)
#     end
# end
# take the average across all cells given that they are currently independent to one another
# D_ave_py = mean(Dlist_py,dims=2)
# D_ave_ju = mean(Dlist_ju,dims=2)

# really low mean square error across the whole dataset at given time interval, of the order e-5

# MSE =sum(((Dlist_ju .- Dlist_py).^2), dims=2)/Ncells
# X = collect(-10:10)
# 0.000850081
# PyPlot.plot(Normal(1.0,stdlist_ju[1,1]))
# d = Normal(0, 1)
# n=rand(d,1000)
# PyPlot.plot(pdf.(d, X))
# PyPlot.scatter(X,pdf.(d, X))
# PyPlot.plot((1:999),mean(stdlist_ju.^2,dims=2), label = "Julia")
# PyPlot.plot((1:999),mean(stdlist_py.^2,dims=2), label = "Python")
# PyPlot.legend()
# PyPlot.xlabel("t")
# PyPlot.ylabel("sigma ^2")

"""
Plot sum of (dA)^2 against different dt?
"""
dA_squared_j = zeros(Float64,(998, Ncells))
dA_squared_p = zeros(Float64,(998, Ncells))

# for i in 1:998
#     for j in 1:256

# test = pos_j[1:998,:] - pos_j[2:999,:]

# for dt in 1:998
#     # for i in 1:256
#     K_j = pos_j 
#     K_p = pos_p
#     dA_j = pos_j[1:999-dt,:] - pos_j[1+dt:999,:]
#     dA_p = pos_p[1:999-dt,:] - pos_p[1+dt:999,:]
#     # K_j[1+dt:999,:] = dA_j
#     # K_p[1+dt:999,:] = dA_p
#     dA_squared_j[dt,:] = mean(dA_j.^2,dims= 1)
#     dA_squared_p[dt,:] = mean(dA_p.^2,dims= 1)
# end
        

# dA_squared_j = mean(dA_squared_j,dims  = 2)
# dA_squared_p = mean(dA_squared_p,dims  = 2)

# PyPlot.plot((1:998),dA_squared_j, label = "julia")
# PyPlot.plot((1:998),dA_squared_p, label = "python")
# PyPlot.xlabel("dt")
# PyPlot.ylabel("dA squared")

# PyPlot.legend()
# axhline(y = 0)
