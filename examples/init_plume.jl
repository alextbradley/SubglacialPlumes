#initializing wavi example
α = 0.01
geometry = zeros(2,1000)
geometry[1,:] = range(0, stop = 1000/α, length = 1000)
geometry[2,:] = -1000 .+ α*geometry[1,:]

params = Params(initial_geometry = geometry)
plume = SubglacialPlumes.start(params)
