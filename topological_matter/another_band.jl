
# TODO: find reference 
# course: topological matter 

using LinearAlgebra
using Plots


sig_x=[0. 1;1 0.]
sig_y=[0. -im;im 0.]
sig_z=[1 0.;0. -1]
iden =[1. 0.; 0. 1.]
# Define Hamiltonian
function H(gam, k0, kx, ky, kz, m, tx, t)
    return eigvals(gam*(cos(kx)-cos(k0))*iden-(m*(2-cos(ky)-cos(kz))+2*tx*(cos(kx)-cos(k0)))*sig_x-2*t*sin(ky)*sig_y-2*t*sin(kz)*sig_z)
end


H_max=2
plot(xlims=(-H_max,H_max),ylims=(-H_max,H_max),camera=(10, 10))
I=60
vals_kx = collect(-H_max:H_max/I:H_max)
vals_kz = collect(-H_max:H_max/I:H_max)
z = zeros(length(vals_kx), length(vals_kz), 2)
#parametros para el calculo del Hamiltoniano
gamma=3
k0=pi/2
ky=0
t=1
tx=1 #tx=t
m=2 #m=2*t
for i in 1:length(vals_kx)
    for j in 1:length(vals_kz)
        z[i,j,:]=H(gamma, k0, vals_kx[i], ky, vals_kz[j], m, tx, t)
        display(scatter!([vals_kx[i]], [vals_kz[j]], [z[i,j,1]], markercolor = "red", markerstrokecolor = "blue", markersize = 1, color=:inferno, legend=false))
        display(scatter!([vals_kx[i]], [vals_kz[j]], [z[i,j,2]], markercolor = "red", markerstrokecolor = "red", markersize = 1,legend=false))
    end
end
#surface(x, y, z[:,:,1])
#surface!(x, y, z[:,:,2])