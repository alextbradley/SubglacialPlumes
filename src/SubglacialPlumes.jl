module SubglacialPlumes

using Parameters, DifferentialEquations

export Params

#Abstract types
abstract type AbstractModel{T <: Real, N<: Integer} end

#Struct to hold model parameters.
#Format: fieldname::Type = default_value.
#T & N are type parameters, usually real numbers (e.g. Float64) and integers
#(e.g. Int64) respectively. Variable names and default values from Jenkins 2011.
@with_kw struct Params{T <: Real}
E0::T = 3.6e-2 #entrainment coefficient
Cd::T = 2.5e-3 #drag coefficient
St::T = 5.9e-4 #Stanton number
λ1::T = -5.73e-2 #Seawater freezing point slope
λ2::T = 8.32e-2 #Seawater freezing point offset
λ3::T = 7.61e-4 #depth dependence freezing point
L::T  = 3.35e5 #Latent heat of fusion for ice
ci::T  = 2.009e3 #Specific heat capacity for ice
cw::T = 3.974e3 #Specific heat capacity for water
βs::T = 7.86e-4 #haline contraction coefficient
βt::T = 3.87e-5 #thermal expansion coefficnet
g::T  = 9.81    #gravitational acceleration
StT::T= 1.1e-3  #Thermal Stanton number
StS::T= 3.1e-5  #Haline Stanton number
secs_per_year:: T = 365.25*24*60*60
Si::T = 0.0       #salinity of ice
Ti::T = 0.0     #tempatature of ice
ρ0::T = 1000.0    #reference density of water
initial_geometry::Array{T,2} #initial geoemetry of the configuration
initial_Sa::Array{T,1} = 34.6 * ones(size(initial_geometry)[2],); @assert size(initial_Sa) == size(initial_geoemtry)[2]#default initial ambient is 34.6 everywhere
initial_Ta::Array{T,1}  = 0.5  * ones(size(initial_geometry)[2],); @assert size(initial_Ta) == size(initial_geoemtry)[2]#default initial ambient is 0.5C everywhere
end

#grid for storing variables on the grid specified by the geometry
@with_kw struct Grid{R <: Real, N <: Integer}
n::N #number of grid points
zgl::R #depth of the grounding line
geometry_length::R #arc length of the geometry
x::Array{R,1}; @assert size(x) == (n,) #x-coordinates of geometry
z::Array{R,1}; @assert size(z) == (n,) #z-coordinates of geometry
dz::Array{R,1}; @assert size(dz) == (n,) #slope of geometry (dZb/dx)
s::Array{R,1}; @assert size(s) == (n,) #arc length parameter at grid points
Ta::Array{R,1}; @assert size(Ta) == (n,)
Sa::Array{R,1}; @assert size(Sa) == (n,)
U::Array{R,1} = zeros(n,); @assert size(U) == (n,)
D::Array{R,1} = zeros(n,); @assert size(D) == (n,)
S::Array{R,1} = zeros(n,); @assert size(S) == (n,)
T::Array{R,1} = zeros(n,); @assert size(T) == (n,)
Δρ::Array{R,1} = zeros(n,); @assert size(Δρ) == (n,)
ΔT::Array{R,1} = zeros(n,); @assert size(ΔT) == (n,)
end

#structure to hold information
@with_kw struct State{T <: Real, N<: Integer} <: AbstractModel{T,N}
params::Params{T}
grid::Grid{T,N}
end

#useful functions

#local freezing point of water
temp_freezing(S,Z,params) = params.λ1 * S + params.λ2 + params.λ3 * Z

#ambient freezing point
temp_af(Sa, Zb,params) = temp_freezing(Sa, Zb,params)

#meltwater freezing point
temp_if(Zb,params) = temp_freezing(params.Si, Zb, params)

#effective meltwater temperature
temp_ief(Z, params) = temp_freezing(params.Si,Z,params) - params.L/params.cw - params.ci*(temp_freezing(params.Si,Z,params) - params.Ti)/params.cw

#density contrast (linear equation of state)
density_contrast(S,Sa, T, Ta, params) = params.ρ0 * (params.βs * (Sa - S) - params.βt*(Ta - T))

#effective meltwater density contrast
density_contrast_effective(Sa, Ta, Z,params) = density_contrast(params.Si, Sa,  temp_ief(Z,params), Ta, params)

#melt rate factor
m_naught(S,Z,params) = params.St / (temp_freezing(S,Z,params) - temp_ief(S,Z,params))

#return the length of the geometry using linear interpolation
get_geometry_length(geometry) = sum( sqrt.(diff(geometry[1,:]).^2 .+ diff(geometry[2,:]).^2))


#create plume from input parameters (parameters and geometry?)
function start(params) where {T <: Real}
    #accepts input Z = Z_b(X) in form geometry = [Xb; Zb(Xb)].
    # Xb(1) specifies grounding line position, Zb(1) specifies grounding line position
    grid = Grid(
        n = size(params.initial_geometry)[2],
        zgl = params.initial_geometry[2,1],
        geometry_length = get_geometry_length(params.initial_geometry),
        x = params.initial_geometry[1,:],
        z = params.initial_geometry[2,:],
        s = get_arclength(params.initial_geometry),
        dz = get_slope(params.initial_geometry),
        Sa = params.initial_Sa,
        Ta = params.initial_Ta)
    plume=State(params, grid)
    return plume
end

function get_slope(geometry)
    #returns the slope Z_b'(X) at grid points X
    dx = diff(geometry[1,:])[1] #we use a regular grid
    slope = zeros(1,size(geometry)[2])
    slope[2:end-1] = (-geometry[2,1:end-2] + geometry[2,3:end])/2 /dx

    #one sided fd for slope at both ends using 23 point stencil
    slope[1] = (-3/2 * geometry[2,1] + 2*geometry[2,2] -1/2 *geometry[2,3])/dx
    slope[end] =(1/2* geometry[2,end-2] - 2*geometry[2,end-1] + 3/2* geometry[2,end])/dx
    return slope[1,:]
end

function get_arclength(geometry)
    #returns the arc length at grid points
    s = zeros(1,size(geometry)[2])
    section_lengths = sqrt.(diff(geometry[1,:]).^2 .+ diff(geometry[2,:]).^2)
    s[2:end] =  [sum(section_lengths[1:i]) for i in 1:length(s)-1] #take sum over all elements with index below in section_lengths
    return s[1,:]
end


function update_geometry!(plume)
    return nothing

end

function get_local_variables(s, grid)
    #return the slope of the base, Z co-ordinate, ambient salinity and
    #temperature associated with arc length parameter s
    α = get_local_var(s, grid.s, grid.dz)
    Z = get_local_var(s, grid.s, grid.z)
    Sa = get_local_var(s, grid.s, grid.Sa)
    Ta = get_local_var(s, grid.s, grid.Ta)
    return α, Z, Sa, Ta
end

function get_local_var(s_local, s, var)
    #returns the value of variable Var (defined on the same grid as S) at the
    #local arclength position s
    mxval, mxindx= findmax(-broadcast(abs, s .-s_local); dims = 1); #returns index of grid point with haf closest to zeros
    mxindx = mxindx[] #remove array packaging
    if mxindx > length(s) - 1 #we may encounter this case close to the front
        local_var = var[end]
        return local_var
    end

    if s[mxindx] > s_local #s_local is smaller than nearest grid point, so s_local ∈ (s[mxindx -1], s[mxindx])
        local_var = ((s[mxindx] - s_local)*var[mxindx-1] + (s_local - s[mxindx .- 1])*var[mxindx])/(s[mxindx] - s[mxindx - 1])
    else #s_local is largest than s at nearest grid point so s_local ∈ (s[mxindx], s[mxindx +1])
        local_var = ((s[mxindx+1] - s_local)*var[mxindx] + (s_local - s[mxindx])*var[mxindx+1])/(s[mxindx+1] - s[mxindx])
    end
    return local_var
end





################### ODE Functions #############################################
function get_rhs(du, v, plume, s)
    #unpack parameters
    D = v[1]
    U = v[2]
    Δρ = v[3]
    ΔT = v[4]
    @unpack params, grid = plume

    #recover salinity and temperature
    #T, S = get_temperature_and_salinity(Δ)
    #M0 = ..
    α, Z, Sa, Ta = get_local_variables(s, grid)

    M0 = params.cw * params.St / params.L #update to include salinity and temperautr

    du[1] = params.E0*U*α + M0*U*ΔT
    du[2] = params.g*D*Δρ*α/params.ρ0 - params.Cd*U^2
    du[3] = M0*U*ΔT*density_contrast_effective(Sa, Ta, Z,params) #+ effective density contrast gradient term
    du[4] = ((Ta - temp_af(Sa, Z, params))*params.E0*α*U +
            (temp_ief(Z, params) - temp_freezing(params.Si, Z, params))*M0*U*ΔT -
            params.λ3*α*D*U)
end

function update_func(A, v, plume, s)
    #unpack parameters
    D = v[1]
    U = v[2]
    Δρ = v[3]
    ΔT = v[4]

    A[1,1] = U
    A[1,2] = D
    A[1,3] = 0
    A[1,4] = 0
    A[2,1] = U^2
    A[2,2] = 2*D*U
    A[2,3] = 0
    A[2,4] = 0
    A[3,1] = U*Δρ
    A[3,2] = D*Δρ
    A[3,3] = D*U
    A[3,4] = 0
    A[4,1] = U*ΔT
    A[4,2] = D*ΔT
    A[4,3] = 0
    A[4,4] = D*U
end

function get_ic(plume)
    #wrapper function to return the initial condition depending on
    #specification.Currently only no discharge initial condition implemented
    return  no_discharge_ic(plume)
end


function no_discharge_ic(plume; xc = 100)
    #returns the initial condition associated with no subglacial discharge
    @unpack params, grid = plume
    αlocal = grid.dz[1]
    E  = params.E0 * αlocal
    Ltilde = params.L + params.ci*(temp_freezing(params.Si,grid.zgl,params) - params.Ti)
    tau = grid.Ta[1] - temp_freezing(grid.Sa[1], grid.zgl, params)

    v0 = zeros(4,)
    v0[4] = E/(E + params.St) * tau;
    v0[3] = params.St * params.cw  * v0[4] * density_contrast_effective(grid.Sa[1], grid.Ta[1], grid.zgl,params) / E /Ltilde
    v0[2] = (E * v0[3] * αlocal * params.g / params.ρ0 /(2*E + (3/2*params.Cd)))^(1/2) * xc^(1/2)
    v0[1] = 2/3 * E * xc
    return v0
end

function solve_plume(plume; sspan = (0,plume.grid.geometry_length))
    #get initial condition
    v0 = get_ic(plume)

    #set up differential equation
    M = DiffEqArrayOperator(ones(4,4), update_func = update_func)
    prob = ODEProblem(ODEFunction(get_rhs, mass_matrix = M), v0, sspan, plume)

    #set condition to trigger integration termination when velocity small
    condition(v,t, integrator) = v[2] - 1e-3
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    sol = solve(prob, Rodas4(), callback = cb)
    return sol
end



end # module
