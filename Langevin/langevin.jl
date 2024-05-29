# a = [1,2,3]
# b = [2,4,6]
# a + b
# https://hockygroup.hosting.nyu.edu/exercise/langevin-dynamics.html
"""
initial area(radii) being uniform (all ones)

The noise is inserted by R = rand(Normal(0.1,dev),Ncells)
where dev is set to 0.1 and 0.2 as usual 
require fairly heavy damping to start with for a more chaotic (unpredictable) pattern of radii distribution both in time and across cells
need to tune it so that the mean across both time and cells would be close to 1.0, the current mean is actually equal to 1.0281540565111256
would provide a std ~0.10375612298080972, which aligns with what we have in the static case, wherre we set the deviation to 0.1 to observe the impact of static noise(heterogeneity).

when dev of the noise = 0.8, at uniform case, the std is 0.2012 (which is pretty good for the sake of the project), mean = 1.028
dev = 0.45, std is 0.100858. mean is 1.027. so in the similar range of the analysis

I will keep gamma the same value as possible, 10 seems to be a good number. the underdamped case doesn't yield good results however as it is not powerful enough to overcome the noise itself possibly?
"""
# using Pkg
# import Pkg; Pkg.add("PlotlyJS")

using PyCall
using StatsBase
using Random
using Plots
using Distributions, LinearAlgebra
using PyPlot
# using Plots
using StatsBase
using Statistics
# Pkg.add("DelimitedFiles")
using CSV, DataFrames
using DelimitedFiles
using PlotlyJS
using Combinatorics
# using 
pygui(true)
Random.seed!(16)
Ncells = 256
Means= ones(Ncells)
initial_position = rand(Normal(1.0,0.1),Ncells)
initial_velocity = Means*0.2
# For harmonic force equilibrium, not actual force
script_dir = @__FILE__
sac_dir = dirname(dirname(dirname(script_dir)))
save_dir = abspath(joinpath(sac_dir, "sims_x"))

function harmonic_energy(x::Union{Float64,Vector{Float64}}, x0::Union{Float64,Vector{Float64}}, k::Float64)
    # if typeof(x)== Vector{Float64}
    dim_x = size(x)
    dim_x0 = size(x0)
    if dim_x != dim_x0
        println("Error: Vectors must have the same dimension.")
        return
    end
        energy = 0.5*k*(x - x0).^2
        force = -k*(x - x0)
    # else
    #     energy = 0.5*k*(x - x0)^2
    #     force = -k*(x - x0)
    # end
return energy,force
end

# function plot_enenrgy_force(spacing::Float64,k::Float64,r0::Float64, rmin::Float64, rmax::Float64)
#     PyPlot.figure()
#     r_points = collect(rmin:spacing:rmax)
#     energies, forces = harmonic_energy(r_points,r0,k)
#     println(forces)
#     println(size(forces)[1])
#     reshape(forces,(size(forces)[1],1))
#     # label = "U(r)"
#     # label = label + ", k=%s"%(k)
#     p = PyPlot.plot(r_points, energies, label = "U(r), k = $k")

#     # PyPlot.plot(r_points, forces,label="",color=p[1].get_color(),linestyle='o')
#     PyPlot.plot(r_points, forces, color = "steelblue", linestyle=":")
#     PyPlot.legend()
# end
# plot_enenrgy_force(0.001,1.,1.,0.8,1.2)

# # plot_energy_force(spacing = 0.001, k=2, r0 = 1)
# BAOAB method, differentiate half a timestep, uses the symmetry of the system
# B
function position_update(x::Vector{Float64}, xdot::Vector{Float64}, dt::Float64)
    # r the vector made of the radii of all cells
    # rdot the small update for each step 
    x_new = x + xdot*dt/2
    return x_new
end
# A
function velocity_update(xdot::Vector{Float64}, F::Vector{Float64}, dt::Float64)
    xdot_new = xdot + F*dt/2
    return xdot_new
end
# O
function random_velocity_update(Ncells::Int64,xdot::Vector{Float64}, gamma::Float64, kBT::Float64,dt::Float64,dev::Float64)
    
    # introduce noise/fluctuation into the temporal variation
    R = rand(Normal(0,dev),Ncells)
    c1 = exp(-gamma*dt)
    # try
    c2 = sqrt((1 - c1^2)*kBT)
    xdot_new = c1*xdot + c2 * R
    # catch y
    #     if isa(y, DomainError)
    #         println("negative product")
    #     end
    # end

    return xdot_new
end

function BAOAB(max_time::Int64, dt::Float64, gamma::Float64, kBT::Float64, Ncells::Int64, k::Float64,dev::Float64)
    x = initial_position
    x0 = initial_position
    v = initial_velocity
    t = 0
    step_number = 1
    # positions = Vector{Vector{Float64}}()
    total_step = ceil(Int,max_time/dt)
    positions = zeros(Float64,(total_step, Ncells))
    # velocities = 
    # velocities = Vector{Vector{Float64}}()
    # total_energies = Vector{Vector{Float64}}()
    velocities = zeros(Float64,(total_step, Ncells))
    total_energies = zeros(Float64,(total_step, Ncells))
    save_times = zeros(Float64,total_step)
    
    while(t<max_time)
            
            # B
            potential_energy, force = harmonic_energy(x,x0,k)
            v = velocity_update(v,force,dt)
            
            #A
            x = position_update(x,v,dt)

            #O
            v = random_velocity_update(Ncells,v,gamma,kBT,dt,dev)
            
            #A
            x = position_update(x,v,dt)
            
            # B
            potential_energy, force = harmonic_energy(x,x0,k)
            v = velocity_update(v,force,dt)
            
                # if step_number%save_frequency == 0 && step_number>0
            e_total = .5*v.^2 .+ potential_energy
                    # push!(positions,x)
                    # push!(velocities,v)
                    # push!(total_energies,e_total)
                    # push!(save_times,t)
                # end

            
            
        if step_number <= total_step    
            
            save_times[step_number] = t
            positions[step_number, :] = x
            velocities[step_number, :] = v
            total_energies[step_number, :] = e_total
        end
        t = t+dt
        step_number = step_number + 1
    end
    # return save_times, positions, velocities, total_energies,step_number   
    return positions,save_times, velocities, total_energies
end

my_k = 5.
my_max_time = 100.
my_dev = 1.0
Ncells = 256
# # simlen
# initial_position = Base.vect(UInt8(1), 0.8, 1//2)
# initial_position= 0.8*ones(Ncells)

# initial_position = Base.vect(1.0)
# initial_velocity = Base.vect(0.1)
my_gamma=10.
# random walk case
# my_kBT=0.25
my_kBT = 2.07*10^(-21)
my_dt = 0.1
# my_dt=0.8

# for i in 1:100
# calculate the radii matrices for the model
# positions, times, velocities, total_energies= BAOAB(my_max_time, my_dt, my_gamma,my_kBT,initial_position,initial_velocity,Ncells,my_k,my_dev)
# writedlm( "positions_julia_1.csv",  positions, ',')
# println(size(positions))

# println(step_number)

# Plot the mean and the one cell autocorrelation
# PyPlot.figure()
# PyPlot.plot(times,positions[:,1].-1.037,label="cell 1")
# across time dims = 1
# across cells dims = 2
# median_cell = median(positions,dims=1)
# mean_cell = mean(positions,dims=1)
# autocor(x, [lags]; demean=true)

# min_cell = min(positions, 2)
# max_cell = max(positions,2)
# median1 = median_cell[1]
# mean1 = mean_cell[1]

# println(std(positions))
# println() 
# dif = abs(mean1)
# println(mean_cell[1])
# , mean = $mean_cell[1]
# PyPlot.figure()

# PyPlot.plot(times,positions[:,1],label="cell 1")
# sigma = std(positions[:,1])
# axhline(y = 1 - std(positions[:,1]))
# axhline(y = 1 + std(positions[:,1]))

# PyPlot.title("dev = $my_dev, gamma = $my_gamma, dt = $my_dt, median = $median1, mean = $mean1")

# axhline(y = mean1)
# PyPlot.plot(times,positions[:,2],label="cell 2")
# PyPlot.plot(times,positions[:,3],label="cell 3")
# xlim(20, 100)
# ylim(1.2,1.5)
# PyPlot.legend
# display()

# PyPlot.plot(times,velocities)
# PyPlot.plot(times,positions)

# PyPlot.xlabel("time")
# plt[:show]()
# show()
# PyPlot.legend(loc="upper center")
# PyPlot.legend()
# PyPlot.show(block=false)
# PyPlot.figure()
# PyPlot.plot(times,total_energies[:,1],label="Simulated E")
# PyPlot.xlabel("time")
# PyPlot.ylabel("Total Energy")
# PyPlot.legend()
# histogram(positions[1,:])
# PyPlot.show()

# L = length(arr) 
# 1:L
# A = autocor(positions, 0:length(positions))
# f = open("/Users/hanshichen/Desktop/engineering/SoftActiveCells/scripts_x/positions.txt", "r")
# py_auto = readdlm("/Users/hanshichen/Desktop/engineering/SoftActiveCells/scripts_x/positions.txt")
# py_pos = readdlm("/Users/hanshichen/Desktop/engineering/SoftActiveCells/scripts_x/positions.txt", Float64)
# # println(py_auto)
# # # few file operations 
# # close(f)
# py_auto = autocor(py_pos, 1:999)
# ju_auto = autocor(positions, 1:999)
# writedlm( "autocor_py_1.csv",  py_auto, ',')
# PyPlot.figure()


# plot(heatmap(z = py_auto))
# PyPlot.imshow(py_auto, cmap="hot", interpolation="nearest")
# PyPlot.plot(py_auto[:,3:5])
# PyPlot.plot(ju_auto[:,1:2])


# ave_ju = mean(ju_auto, dims=2)
# ave_py = mean(py_auto, dims=2)
# PyPlot.imshow(ave_ju, cmap="hot", interpolation="nearest")
# PyPlot.imshow(ave_py, cmap="hot", interpolation="nearest")
# PyPlot.plot(ave_ju, label="Julia")
# PyPlot.plot(ave_py, label="Python")
# PyPlot.plot(py_auto[:,3:5])
# PyPlot.legend()

# PyPlot.show()

# K = 1.:10.; Gamma = 0.5:0.5:10. ;Max_time = 50:50:1000; Dev = 0.05:0.2; my_dt = 0.1; 
# K = [10.]; Gamma = [1.,5.,10.] ;Max_time = [50,100,1000]; Dev = 0.05:0.05:0.2; my_dt = 0.1; 
# K = [1.,2.,5.,10.]; Gamma = 0.5:0.5:10 ;Max_time = [50,200,500,1000]; Dev = 0.1; my_dt = 0.1; 
K = 1.:0.5:10; Gamma = 5.0 ;Max_time = [1000]; Dev = 0.1; my_dt = 0.1; 

my_k = 5.
my_max_time = 100.
my_dev = 1.0
# Ncells = 256
# # simlen
# initial_position = Base.vect(UInt8(1), 0.8, 1//2)
# initial_position= 0.8*ones(Ncells)
Means= ones(Ncells)
initial_position = rand(Normal(1.0,0.1),Ncells)
initial_velocity = Means*0.2
# initial_position = Base.vect(1.0)
# initial_velocity = Base.vect(0.1)
my_gamma=10.
# random walk case
my_kBT=0.25
my_dt = 0.1
# Ncells = 256
all_com = collect(Iterators.product(Max_time,my_dt,Gamma,my_kBT,256,K,Dev))
# f = open("/Users/hanshichen/Desktop/engineering/SoftActiveCells/scripts_x/positions.txt", "r")

writedlm( "k_update_combinations.csv",  x, ',')
# (my_max_time, my_dt, my_gamma,my_kBT,initial_position,initial_velocity,Ncells,my_k,my_dev)

# A1 = all_com[1]
# for set in all_com
#     positions, times, velocities, total_energies= BAOAB(my_max_time, my_dt, my_gamma,my_kBT,initial_position,initial_velocity,Ncells,my_k,my_dev)
# end

# BAOAB(max_time::Int64, dt::Float64, gamma::Float64, kBT::Float64, Ncells::Int64, k::Float64,dev::Float64)
# nosim = 10
# full_pos = zeros(Float64,(nosim,10000,256))
# full_pos = zeros(Float64,(10000,256))
for (maxtime, fdt, fgamma, fkBT, fNcells, fk, fdev) in all_com
    # i = 
    # full_pos = zeros(Float64,(10000,256))

    # for i in 1:nosim
        # positions, times, velocities, total_energies= BAOAB(1000, 0.1, 1.0, 0.25, 256, 10.0, 0.1)
        positions, times, velocities, total_energies= BAOAB(maxtime, fdt, fgamma,fkBT, fNcells, fk, fdev)
        # full_pos[:,:] = positions
        # push!(full_pos,positions)

        
    # end
    try
        writedlm(abspath(joinpath(save_dir, "positions_from_langevin/positions_maxtime=$maxtime,dt=$fdt,gamma=$fgamma,kBT=$fkBT,k=$fk,dev=$fdev.csv")),  mean(positions,dims=2), ',')
        writedlm(abspath(joinpath(save_dir, "positions_non_averaged_langevin/positions_maxtime=$maxtime,dt=$fdt,gamma=$fgamma,kBT=$fkBT,k=$fk,dev=$fdev.csv")), positions,  ',')
      
    catch
    end

end
# mean()
# writedlm( "boltzman_positions_maxtime=$maxtime,dt=$fdt,gamma=$fgamma,kBT=$fkBT,k=$fk,dev=$fdev.csv", mean(full_pos,dims=1), ',')
# writedlm( "boltzman2positions_maxtime=1000,dt=0.1,gamma=1.0,kBT=0.25,k=10.0,dev=0.1,validate.csv", mean(full_pos,dims=(1,3)), ',')