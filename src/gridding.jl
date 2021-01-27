# Generate double resolution grid
function get_2x_axis(lx, ly, lz, dx, dy, dz)
    nx, ny, nz = cld.([lx,ly,lz], [dx,dy,dz]) .|> Int
    nx, ny, nz = max.([nx,ny,nz], 1)
    x2s = range(-1/2, stop=1/2, length=2nx+1)[1:end-1] * nx*dx
    y2s = range(-1/2, stop=1/2, length=2ny+1)[1:end-1] * ny*dy
    z2s = range(-1/2, stop=1/2, length=2nz+1)[1:end-1] * nz*dz
    return x2s, y2s, z2s
end

# Functions for downsampling ep and mu
function downsample_2x_grid(ep2, mu2)
    @assert all(iseven.(size(ep2)))
    @assert size(ep2) == size(mu2)
    nx, ny, nz = div.(size(ep2), 2)
    ep = ones(ComplexF64, nx, ny, nz, 3)
    ep[:,:,:,1] = ep2[2:2:end, 1:2:end, 1:2:end]
    ep[:,:,:,2] = ep2[1:2:end, 2:2:end, 1:2:end]
    ep[:,:,:,3] = ep2[1:2:end, 1:2:end, 2:2:end]
    mu = ones(ComplexF64, nx, ny, nz, 3)
    mu[:,:,:,1] = mu2[1:2:end, 2:2:end, 2:2:end]
    mu[:,:,:,2] = mu2[2:2:end, 1:2:end, 2:2:end]
    mu[:,:,:,3] = mu2[2:2:end, 2:2:end, 1:2:end]
    return ep, mu
end

# functions for PMLs?