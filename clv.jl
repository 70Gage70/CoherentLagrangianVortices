using ThreadsX
using MAT
using StaticArrays
using AxisArrays
using Dates
using Plots
using NetCDF, NCDatasets
using Interpolations
using ScatteredInterpolation

using CoherentStructures, OrdinaryDiffEq, Tensors

# define metric
function metric(x, y)
   a = 6387.0 # [km]
   sqrt2 = a^2 - x^2 - y^2 # (x,y) meaningful inside circle of radius a
   sqrt2 < 0 && return one(Tensor{2, 2}) # ensure entries are positive outside circle
   g11 = 1 + x^2/sqrt2
   g12 = x*y/sqrt2
   g22 = 1 + y^2/sqrt2
   return Tensors.SymmetricTensor{2, 2, Float64}((g11, g12, g22)) 
end

# define metric representation of Cauchy-Green tensor
function cauchygreen(odefun, xy0, tspan, δ; kwargs...)
   lf = linearized_flow(odefun, xy0, [tspan[1], tspan[end]], δ; kwargs...)
   xyt, dF = lf[1][end], lf[2][end]
   G0 = metric(xy0...)
   G = metric(xyt...)
   C = G0^(-1/2) * (dF' ⋅ G ⋅ dF) * G0^(-1/2)
   return Tensors.SymmetricTensor{2, 2, Float64}((C[1, 1], C[1, 2], C[2, 2])) 
end

# fix issue with interpolateVF in original CoherentStructures.jl
function interpolateVF(
   X::AbstractRange{S1},
   Y::AbstractRange{S1},
   T::AbstractRange{S1},
   U::AbstractArray{S2,3},
   V::AbstractArray{S2,3},
   itp_type = BSpline(Cubic(Free(OnGrid()))) ) where {S1 <: Real, S2 <: Real}

   UV = map(SVector{2,S2}, U, V)::Array{SVector{2,S2},3}
   return extrapolate(Interpolations.scale(Interpolations.interpolate(UV, itp_type), X, Y, T), SVector{2}([0.0,0.0]))
end

# load velocity data in coordinates
fileIn = "ERA5-600K-May.nc"
x = ncread(fileIn, "x")    # [km]
y = ncread(fileIn, "y")
t = ncread(fileIn, "t")    # [days]
vx = ncread(fileIn, "u")   # [km/days]
vy = ncread(fileIn, "v")
replace!(vx, NaN => 0)
replace!(vy, NaN => 0)

# velocity interpolant
xmin, xmax = extrema(x)
ymin, ymax = extrema(y)
tmin, tmax = extrema(t)
xrange = range(xmin, stop = xmax, length = length(x))
yrange = range(ymin, stop = ymax, length = length(y))
trange = range(tmin, stop = tmax, length = length(t))
Fv = interpolateVF(xrange, yrange, trange, vx, vy)

# deal with coasts (land) if present...
landxy = vx[:,:,1] .== 0
land = map(Iterators.product(1:size(landxy, 1)-1, 1:size(landxy, 2)-1)) do (i, j)
   return landxy[i,j] & landxy[i+1,j] & landxy[i,j+1] & landxy[i+1,j+1]
end
dx, dy = step.((xrange, yrange))
is_land = let land=land, dx=dx, dy=dy, xmin=xmin, ymin=ymin
x -> land[Int((x[1] - xmin)÷dx + 1), Int((x[2] - ymin)÷dy + 1)]
end

# x0,y0
xrangemin, xrangemax = extrema(xrange)
yrangemin, yrangemax = extrema(yrange)
x0min, x0max, y0min, y0max = xrangemin+5, xrangemax-5, yrangemin+5, yrangemax-5
Nx = 128
Ny = floor(Int, (y0max - y0min) / (x0max - x0min) * Nx)
x0 = range(x0min, stop=x0max, length=Nx)
y0 = range(y0min, stop=y0max, length=Ny)
P = AxisArray(SVector{2}.(x0, y0'), x0, y0)

# t0
t0 = (Date(2002, 5, 10) - Date(0, 1, 1)).value + 1

# vortex parameters
#known_centers = [Singularity((0,0), 1)]
p = LCSParameters(
   boxradius = 5000, 
   n_seeds = 500,	
   pmin = .8,
   pmax = 2,
   rdist = 1, 
   merge_heuristics = []
   )
T = 20

# compute C
function vF(du, u, p, t)
   x, y = u
   du[1:2] .= Fv(x, y, t)
end
tspan = range(t0, stop = t0+T, length = 2)
CG = u -> cauchygreen(vF, u, tspan, 1; tolerance = 1e-6, solver= Tsit5())
C = ThreadsX.map(CG, P) |> c -> AxisArray(c, row = P.axes[1].val, col = P.axes[2].val)
@info "CG computation::Done!"

# extract CLVs
vortices, singularities = ellipticLCS(C, p; 
	outermost = true,
   verbose = true,
   #suggested_centers = known_centers
)
@info "CLV extraction(s)::Done!"

# plot
plot_singularities!(singularities)
plot_barrier!(vortices[1].barriers[1],lc=:black)
png("output")

### BEGIN: save to MATLAB .mat
# fileOut = "output.mat"
# rm(fileOut, force = true)
# file = matopen(fileOut, "w")
# N = size(vortices,1)
# for i in 1:N
#    xy = vortices[i].barriers[1].curve
#    G0inv = metric(xy[1]...)^(-1/2)
#    xy_out = map(x -> G0inv * x, xy) |> stack
#    xout, yout = xy_out[1,:], xy_out[2,:]
#    pout  = vortices[i].barriers[1].p
#    write(file, "x$i", xout)
#    write(file, "y$i", yout)
#    write(file, "p$i", pout)
# end
# close(file)
# @info "Saving to $(fileOut)::Done!"
### END: save to MATLAB .mat

### BEGIN: save to NetCDF .nc
fileOut = "output.nc"
rm(fileOut, force = true)
ds = NCDataset(fileOut, "c")
N = size(vortices,1)
xnames = ["x$i" for i = 1:N]
ynames = ["y$i" for i = 1:N]
pnames = ["p$i" for i = 1:N]
for i = 1:N
   xy = vortices[i].barriers[1].curve
   G0inv = metric(xy[1]...)^(-1/2)
   xy_out = map(x -> G0inv * x, xy) |> stack
   xout, yout = xy_out[1,:], xy_out[2,:]
   pout  = vortices[i].barriers[1].p  
   defVar(ds, xnames[i], xout, ("$(xnames[i])_dim",))
   defVar(ds, ynames[i], yout, ("$(ynames[i])_dim",))
   defVar(ds, pnames[i], pout, ("$(pnames[i])_dim",))
end
close(ds)
@info "Saving to $(fileOut)::Done!"
### END: save to NetCDF .nc