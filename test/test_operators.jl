nx, ny, nz = 2, 3, 4
dx, dy, dz = 0.1, 0.2, 0.3

kx, ky, kz = rand(3)
F = zeros(ComplexF64, nx, ny, nz)
for ix in 1:nx, iy in 1:ny, iz in 1:nz
    F[ix,iy,iz] = exp(1im*kx*ix*dx + 1im*ky*iy*dy + 1im*kz*iz*dz)
end
F = F[:]

DxE, DyE, DzE, DxH, DyH, DzH, Z = getyeeops(nx, ny, nz, dx, dy, dz, kxlx=kx*nx*dx, kyly=ky*ny*dy, kzlz=kz*nz*dz)
@test DxE*F ≈ (exp(+1im*kx*dx)*F - F)/dx
@test DyE*F ≈ (exp(+1im*ky*dy)*F - F)/dy
@test DzE*F ≈ (exp(+1im*kz*dz)*F - F)/dz
@test DxH*F ≈ (F - exp(-1im*kx*dx)*F)/dx
@test DyH*F ≈ (F - exp(-1im*ky*dy)*F)/dy
@test DzH*F ≈ (F - exp(-1im*kz*dz)*F)/dz