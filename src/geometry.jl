struct Geometry
    ep
    mu
    dx
    dy
    dz
end

function Geometry(lx, ly, lz, dx, dy, dz, ep::Function, mu::Function)
    x2s, y2s, z2s = get_2x_axis(lx, ly, lz, dx, dy, dz)
    ep2 = [ep(x,y,z) for x in x2s, y in y2s, z in z2s]
    mu2 = [mu(x,y,z) for x in x2s, y in y2s, z in z2s]
    ep, mu = downsample_2x_grid(ep2, mu2)
    return Geometry(ep, mu, dx, dy, dz)
end