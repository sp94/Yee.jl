using SparseArrays

function getyeeops(nx::Int, ny::Int, nz::Int, dx::Real, dy::Real, dz::Real; kxlx=0, kyly=0, kzlz=0)
    n = nx*ny*nz
    Dx = [ix!=nx for ix in 1:nx, iy in 1:ny, iz in 1:nz]
    Dy = [iy!=ny for ix in 1:nx, iy in 1:ny, iz in 1:nz]
    Dz = [iz!=nz for ix in 1:nx, iy in 1:ny, iz in 1:nz]
    DxE = spdiagm(
        0 => -ones(n),
        1 => Dx[1:n-1],
        1-nx => exp(1im*kxlx)*(1 .- Dx)[nx:n]
    ) / dx
    DyE = spdiagm(
        0 => -ones(n),
        nx => Dy[1:n-nx],
        nx-nx*ny => exp(1im*kyly)*(1 .- Dy)[1+nx*ny-nx:n]
    ) / dy
    DzE = spdiagm(
        0 => -ones(n),
        nx*ny => Dz[1:n-nx*ny],
        nx*ny-nx*ny*nz => exp(1im*kzlz)*(1 .- Dz)[1+nx*ny*nz-nx*ny:n]
    ) / dz
    DxH, DyH, DzH = -DxE', -DyE', -DzE'
    return DxE, DyE, DzE, DxH, DyH, DzH, spzeros(n,n)
end