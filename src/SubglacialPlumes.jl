module SubglacialPlumes

using Parameters, DifferentialEquations

export Params, solve!, parametrize!

#Abstract types
abstract type AbstractModel{T <: Real, N<: Integer} end

"""
Struct to hold model parameters.
#Format: fieldname::Type = default_value.
#T & N are type parameters, usually real numbers (e.g. Float64) and integers
#(e.g. Int64) respectively. Variable names and default values from Jenkins 2011.
Also store the specified ambient temperature and salinity and the initial geometry here.
"""
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
geometry::Array{T,2}  = linear_geometry(0.01, -1000) #default geometry is one with slope o.01 and grounding line depth -1000
ambient_grid::Array{T,1} = Array(range(geometry[2,1], stop = 0, length = 100)); @assert ~any(diff(diff(ambient_grid)) .> 1e-5) #default ambient grid is linearly space between grounding line depth and surface
#assert regular grid in the ambient salinity
Sa::Array{T,1} = 34.6 * ones(size(ambient_grid)); @assert size(Sa) == size(ambient_grid)#default initial ambient is 34.6 everywhere
Ta::Array{T,1}  = 0.5  *ones(size(ambient_grid)); @assert size(Ta) == size(ambient_grid)#default initial ambient is 0.5C everywhere
end

"""
Structure for holding variables defined on the (user-specified) grid and grid properties
"""
@with_kw struct Grid{R <: Real, N <: Integer}
n::N #number of grid points
geometry_length::R #arc length of the geometry
x::Array{R,1}; @assert size(x) == (n,) #x-coordinates of geometry
z::Array{R,1}; @assert size(z) == (n,) #z-coordinates of geometry
dz::Array{R,1}; @assert size(dz) == (n,) #slope of geometry (dZb/dx)
s::Array{R,1}; @assert size(s) == (n,) #arc length parameter at grid points
U::Array{R,1} = zeros(n,); @assert size(U) == (n,)
D::Array{R,1} = zeros(n,); @assert size(D) == (n,)
S::Array{R,1} = zeros(n,); @assert size(S) == (n,)
T::Array{R,1} = zeros(n,); @assert size(T) == (n,)
Δρ::Array{R,1} = zeros(n,); @assert size(Δρ) == (n,)
ΔT::Array{R,1} = zeros(n,); @assert size(ΔT) == (n,)
m::Array{R,1} = zeros(n,); @assert size(m) == (n,)
mparam::Array{R,1} = zeros(n,); @assert size(m) == (n,) #store a parametrization of the melt rate
end

"""
structure for storing useful quantities computed from the input parameters
"""
@with_kw struct Store{R <: Real}
    zgl::R #depth of the grounding line
    v0::Array{R,1} = zeros(4,) #store initial conditions
    tau::R #thermal driving at the grounding line
    l0::R #z -length scale associated with freezing point dependence
    dρa_dz::Array{R,1}; #ambient density gradient, defined on ambient grid
end

""" 
Structure to hold each of the above structures when plume is instantiated
"""
@with_kw struct State{T <: Real, N<: Integer} <: AbstractModel{T,N}
params::Params{T}
grid::Grid{T,N}
store::Store{T}
end

#useful functions

#define for testing
f(x,y) = 2x + 3y

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
m_naught(S,Z,params) = params.St / (temp_freezing(S,Z,params) - temp_ief(Z,params))
m_naught(params) =  params.cw * params.St / params.L #approximation using L/c >> Tf - Tfi

#return the length of the geometry using linear interpolation
get_geometry_length(geometry) = sum( sqrt.(diff(geometry[1,:]).^2 .+ diff(geometry[2,:]).^2))

linear_geometry(α, zgl; n = 1000) = vcat(Array(range(0, stop = abs(zgl)/α, length = n))',Array( zgl .+ α* range(0, stop = abs(zgl)/α, length = n))') #initial geoemetry of the configuration. Default is linear with GL depth -1000 and slope 0.01

get_tau(Ta0,Sa0, zgl, params) = Ta0 - temp_freezing(Sa0, zgl, params)

#create plume state from input parameters
function start(params) 
    #accepts input Z = Z_b(X) in form geometry = [Xb; Zb(Xb)].
    # Xb(1) specifies grounding line position, Zb(1) specifies grounding line position
    grid = Grid(
        n = size(params.geometry)[2],
        geometry_length = get_geometry_length(params.geometry),
        x = params.geometry[1,:],
        z = params.geometry[2,:],
        s = get_arclength(params.geometry),
        dz = get_slope(params.geometry))

    store = Store(
        zgl = params.geometry[2,1],
        tau = get_tau(params.Ta[1], params.Sa[1], params.geometry[2,1], params),
        l0  = get_tau(params.Ta[1], params.Sa[1], params.geometry[2,1], params)/params.λ3,  #z lengthscale associated with freezing point dependence. Uses slope at first grid point as angle scale
        dρa_dz = get_ambient_density_gradient(params.Sa, params.Ta, params)
    )
    plume=State(params, grid, store)
    return plume
end




"""
Returns the gradient in ambient density at grid points. Computed using centered finite differences, except for the endpoint, which use one sided FD.
"""
function get_ambient_density_gradient(Sa, Ta, params)
    #compute the associated ambient density using linear equation of state
    dz = diff(params.ambient_grid)[1] #we use a regular grid
    ρa = zeros(size(Sa))
    dρa_dz = zeros(size(Sa))
    @. ρa = params.ρ0 * (params.βs * Sa - params.βt * Ta)
    dρa_dz[2:end-1] = (-ρa[1:end-2] + ρa[3:end])/2/dz
    dρa_dz[1] = (-3/2 * ρa[1] + 2*ρa[2] -1/2 *ρa[3])/dz
    dρa_dz[end] = (1/2* ρa[end-2] - 2*ρa[end-1] + 3/2* ρa[end])/dz
    return dρa_dz
end

"""
returns the slope Z_b'(X) at grid points X. Computed using centered finite differences, except for the end points which use one-sided differences
"""
function get_slope(geometry)
    dx = diff(geometry[1,:])[1] #we use a regular grid
    slope = zeros(1,size(geometry)[2])
    slope[2:end-1] = (-geometry[2,1:end-2] + geometry[2,3:end])/2 /dx

    #one sided for all grid points before final
    slope[1:end-1] = (geometry[2,2:end] - geometry[2,1:end-1])/ dx

    #one sided fd for slope at both ends using 23 point stencil
    #slope[1] = (-3/2 * geometry[2,1] + 2*geometry[2,2] -1/2 *geometry[2,3])/dx
    slope[end] =(1/2* geometry[2,end-2] - 2*geometry[2,end-1] + 3/2* geometry[2,end])/dx
    return slope[1,:]
end

"""
Returns the arc length of the curve specified by the geometry 
"""
function get_arclength(geometry)
    #returns the arc length at grid points
    s = zeros(1,size(geometry)[2])
    section_lengths = sqrt.(diff(geometry[1,:]).^2 .+ diff(geometry[2,:]).^2)
    s[2:end] =  [sum(section_lengths[1:i]) for i in 1:length(s)-1] #take sum over all elements with index below in section_lengths
    return s[1,:]
end


"""
Returns the value of variable Var (defined on the same grid as S) at the local position s
"""
function get_local_var(s_local, s, var)

    mxval, mxindx= findmax(-broadcast(abs, s .-s_local); dims = 1); #returns index of grid point with haf closest to zeros
    mxindx = mxindx[] #remove array packaging
    if mxindx > length(s) - 1 #we may encounter this case close to the front
        local_var = var[end]
        return local_var
    end

    if (mxindx == 1) && s_local < s[mxindx] #if local point is smaller than all grid points in s
        local_var = var[1]
    end

    if s[mxindx] > s_local #s_local is smaller than nearest grid point, so s_local ∈ (s[mxindx -1], s[mxindx])
        local_var = ((s[mxindx] - s_local)*var[mxindx-1] + (s_local - s[mxindx .- 1])*var[mxindx])/(s[mxindx] - s[mxindx - 1])
    else #s_local is largest than s at nearest grid point so s_local ∈ (s[mxindx], s[mxindx +1])
        local_var = ((s[mxindx+1] - s_local)*var[mxindx] + (s_local - s[mxindx])*var[mxindx+1])/(s[mxindx+1] - s[mxindx])
    end
    return local_var
end

"""
Converts a value of Δρ, ΔT, Z, Sa(Z), Ta(Z) to the corresponding values of salinity S and temperature T
"""
function convert_to_st(Δρ, ΔT, Z, Sa, Ta, params)
    #convert a value of Δρ and ΔT to T, S
    S = params.βs * Sa - Δρ./params.ρ0 - params.βt * (Ta - ΔT - params.λ2 - params.λ3 * Z) 
    S = S/(params.βs - params.λ1 * params.βt)
    T = ΔT + temp_freezing(S, Z, params)
    return T, S
end


################### ODE Functions #############################################
"""
Returns the right hand side of the plume equations (equations 14-17 in Jenkins 2011)
"""
function get_rhs(du, v, plume, s)
    #unpack parameters
    D = v[1]
    U = v[2]
    Δρ = v[3]
    ΔT = v[4]
    @unpack params, grid, store = plume

    #recover salinity and temperature
    #M0 = .. 
    α = get_local_var(s, grid.s, grid.dz)
    Z = get_local_var(s, grid.s, grid.z)
    Sa = get_local_var(Z, params.ambient_grid, params.Sa)
    Ta = get_local_var(Z, params.ambient_grid, params.Ta)
    dρa_dz = get_local_var(Z, params.ambient_grid, store.dρa_dz)
    M0 = m_naught(params) #approximation using L/c >> Tf - Tfi

    du[1] = params.E0*U*α + M0*U*ΔT
    du[2] = params.g*D*Δρ*α/params.ρ0 - params.Cd*U^2
    du[3] = (M0*U*ΔT*density_contrast_effective(Sa, Ta, Z,params)  + 
            α * dρa_dz * D * U); 
    du[4] = ((Ta - temp_af(Sa, Z, params))*params.E0*α*U +
            (temp_ief(Z, params) - temp_freezing(params.Si, Z, params))*M0*U*ΔT -
            params.λ3*α*D*U)
end

"""
Specifies the mass matrix A, where plume equations are specified by A * dv/dt = f (v = D, U, Δρ, ΔT)
"""
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


""" 
    Update the initial conditions stored in grid. 
    Overload this method to use initial conditions other than zero flux initial conditions.
"""
function update_ic!(plume::AbstractModel; xc = 10) 
    @unpack params, grid, store = plume
    αlocal = grid.dz[1]
    E  = params.E0 * αlocal
    Ltilde = params.L + params.ci*(temp_freezing(params.Si,store.zgl,params) - params.Ti)
    
    #get the ambient temperature and salnity at the grounding line
    Sa = get_local_var(store.zgl, params.ambient_grid, params.Sa)
    Ta = get_local_var(store.zgl, params.ambient_grid, params.Ta)
    
    v0 = zeros(4,) #explicitly define to ensure size compatibility
    v0[4] = E/(E + params.St) * store.tau;
    v0[3] = params.St * params.cw  * v0[4] * density_contrast_effective(Sa, Ta, store.zgl,params) / E /Ltilde
    v0[2] = (E * v0[3] * αlocal * params.g / params.ρ0 /(2*E + (3/2*params.Cd)))^(1/2) * xc^(1/2)
    v0[1] = 2/3 * E * xc
    store.v0 .= v0
    return plume
end

"""
    Solve the plume equations. Variables stored on grid are modified. solution span set span sspan: sovles on the whole s domain if there are no points at which dz goes negative, else solves until the first point where dz goes negative. 
"""
function solve!(plume; sspan = (0, (any(plume.grid.dz .< 1e-5) ? plume.grid.s[findfirst(plume.grid.dz .< 1e-5)] : plume.grid.s[end])))
    @unpack params, grid, store = plume

    #endow plume with an initial condition (overload to change from no-discharge ic)
    update_ic!(plume)

    #set up differential equation
    M = DiffEqArrayOperator(ones(4,4), update_func = update_func)
    prob = ODEProblem(ODEFunction(get_rhs, mass_matrix = M), store.v0, sspan, plume)

    #set condition to trigger integration termination when velocity small
    condition(v,t, integrator) = v[2] - eps()
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    #note that the below does not account for the case where the plume terminates before edge of ice shelf
    sol = solve(prob, Rodas4(), callback = cb, saveat = plume.grid.s)
    nout = length(sol[1,:]) #number of cells solution fills
    grid.D[1:nout] .= sol[1,:]
    grid.U[1:nout] .= sol[2,:]
    grid.Δρ[1:nout] .= sol[3,:]
    grid.ΔT[1:nout] .= sol[4,:]

    #find corresponding salinity and temperature (no vectorization at params not vectorized)
    #for i = 1:length(grid.n)
    #    grid.T[i], grid.S[i] = SubglacialPlumes.convert_to_st(plume.grid.Δρ[i], plume.grid.ΔT[i], plume.grid.z[i], plume.grid.Sa[i], plume.grid.Ta[i], params)
    #end
    
    update_melt_rates!(plume)
    return plume
end

"""
Updates the melt rates stored in the grid structure to correspond to the values of U, ΔT stored in same structure
"""
function update_melt_rates!(plume)
    @unpack params, grid = plume
    M0 = m_naught(params) #approximation using L/c >> Tf - Tfi
    grid.m .= params.secs_per_year * M0.* grid.U .* grid.ΔT
    return plume
end

"""
    Update the mparam field to include a parametrization of the melt rate. 
    Second argument specifies the method of paramtrization (default is Lazeroms w/out ad hoc geometric dependence)
"""
function parametrize!(plume; method = "Lazeroms")
    if method == "Lazeroms"
        parametrization_lazeroms!(plume)
    elseif method == "Lazeroms_AHG"
        parametrization_lazeroms_AHG!(plume)
    else
        error("Parametrization method must be one of following: Lazeroms, Lazeroms_AHG")
    end

end

"""
    Modify the melt rate parametrization field to the standard Lazeroms parametrization
"""
function parametrization_lazeroms!(plume)
    @unpack params, grid, store = plume
    ΔTscale   = params.E0 *grid.dz[1] * store.tau/params.St; 
    Uscale = sqrt(params.βs * grid.Sa[1] * params.g * store.l0 * store.tau * params.E0 * grid.dz[1]/(params.L/params.cw) / params.Cd);
    snd = grid.s./(store.l0 / grid.dz[1]) #dimensionless version of arclength parameter s
    
    #need to add control for snd > 1
    M0 = m_naught(params) #approximation using L/c >> Tf - Tf
    @. grid.mparam = params.secs_per_year * M0 * Uscale *  ΔTscale * (3*(1 - snd)^(4/3) - 1).*(1 - (1 - snd)^(4/3))^(1/2) /2 /sqrt(2)
    return plume
end

"""
    Modify the melt rate parametrization field to the Lazeroms ad-hoc geometry parametrization
"""
function parametrization_lazeroms_AHG!(plume)
    @unpack params, grid, store = plume
    ΔTscale = zeros(grid.n,)
    Uscale  = zeros(grid.n,)
    @. ΔTscale = params.E0 *grid.dz * store.tau/params.St; #vector: different angle in scale at each grid point
    @. Uscale = sqrt(params.βs * grid.Sa[1] * params.g * store.l0 * store.tau * params.E0 * grid.dz/(params.L/params.cw) / params.Cd);
    snd = grid.s./(store.l0 / grid.dz[1]) #dimensionless version of arclength parameter s
    #need to add control for snd > 1
    M0 = m_naught(params) #approximation using L/c >> Tf - Tf
    @. grid.mparam = params.secs_per_year * M0 * Uscale *  ΔTscale * (3*(1 - snd)^(4/3) - 1).*(1 - (1 - snd)^(4/3))^(1/2) /2 /sqrt(2)
    return plume
end


end # module
