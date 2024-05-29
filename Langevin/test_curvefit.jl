using PyPlot
using CurveFit
using CSV, DataFrames

x = 0.0:0.02:2.0
y0 = @. 1 + x + x*x + randn()/10
fit = CurveFit.curve_fit(Polynomial, x, y0, 2)
y0b = fit.(x) 
plot(x, y0, "o", x, y0b, "r-", linewidth=3)


x = [1,2,3,4,5,6,7,8,9]; y0 = [1,2,3,4,5,6,7,8,9]
fit = CurveFit.curve_fit(LinearFit,x,y0)
y0b = fit.(x)

plot(x, y0, "o", x, y0b, "r-", linewidth=3)


K = Matrix(CSV.read(abspath("/Users/hanshichen/Desktop/engineering/SoftActiveCells/sims_x/positions_non_averaged_langevin/positions_maxtime=1000,dt=0.1,gamma=9.0,kBT=0.25,k=1.0,dev=0.1.csv"), DataFrame))
first_cell = K[1,:]