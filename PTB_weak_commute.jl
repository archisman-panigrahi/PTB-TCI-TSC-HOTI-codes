using Plots
using LinearAlgebra
using LaTeXStrings
using StatsBase
using CSV
using DataFrames

mu=0.01;
t = 1;
t0 = 1;

z1_array = range(0.1,sqrt(2)-0.1,5);
z2_array = range(-sqrt(2)+0.1,-0.1,5);

for i_vals = 1:length(z1_array)
    for j_vals = 1:length(z2_array)
        z1 = z1_array[i_vals];
        z2 = z2_array[j_vals];

m_0 = (z1+z2)/sqrt(2);
Delta = (z1-z2)/sqrt(2);

# m_0 = -0.2;
# Delta = 0.5;

x_periodic = 0;
y_periodic = 0;

L = 80;

Lx = L;
Ly = L;

exportdata=1;


sigma_x = [0 1;
           1 0];
sigma_y = [0 -im;
           im 0];
sigma_z = [1 0;
           0 -1];
eye2 = [1 0;
        0 1];
gamma_01 = kron(eye2, sigma_x);
gamma_02 = kron(eye2, sigma_y);
gamma_03 = kron(eye2, sigma_z);
gamma_13 = kron(sigma_x,sigma_z);
gamma_30 = kron(sigma_z, eye2);
# gamma_31 = kron(sigma_z, sigma_x);
# gamma_32 = kron(sigma_z, sigma_y);
# gamma_33 = kron(sigma_z, sigma_z);
gamma_12 = kron(sigma_x, sigma_y);
gamma_23 = kron(sigma_y, sigma_z)

#m = (-1+sqrt(5))/2;
m=2/3
c1 = 6
c2 = 18

function line_down(x::Float64)
    return c1 + m*x
end

function line_up(x::Float64)
    return c2 + m*x
end

function line_middle(x::Float64)
    return (line_up(x) + line_down(x))/2
end

line_x = [1.0,L]
line_down_y = line_down.(line_x)
line_up_y = line_up.(line_x)

points_x_array = range(1, L,L)
points_y_array = range(1, L,L)

points2 = zeros(2*L^2)
for ii = 1:L
    for jj = 1:L
        points2[2*((jj-1)*L+ii)-1] = points_x_array[ii]
        points2[2*((jj-1)*L+ii)] = points_y_array[jj]
    end
end


### These arrays have N^2 elements, and hold the x and y coordinates of nth site
points_x_array = points2[1:2:end]
points_y_array =  points2[2:2:end]

site_index = range(1,L^2, L^2)

### PTB_index will contain ordered indices of the sites in PTB

PTB_index_dummy = collect(range(1,L^2,L^2))
for ii = 1:L^2
    if (points_y_array[ii] >= line_down(points_x_array[ii])) && (points_y_array[ii] <= line_up(points_x_array[ii]))
        ## Do nothing
    else
        PTB_index_dummy[ii] = 0.0
    end
end

## Here we store non-zero indices
PTB_index = findall(!iszero, PTB_index_dummy)

N_PTB = size(PTB_index)[1]

points_PTB_array_x = zeros(N_PTB)
points_PTB_array_y = zeros(N_PTB)

for ii = 1:N_PTB
    points_PTB_array_x[ii] = points_x_array[PTB_index[ii]]
    points_PTB_array_y[ii] = points_y_array[PTB_index[ii]]
end

### Draw lines separating PTB
line_x = [1.0,L]
line_down_y = line_down.(line_x)
line_up_y = line_up.(line_x)

println("Number of sites in PTB = ",N_PTB)
println("Amount of sites in PTB = ",100* N_PTB/L^2, " %")

println(PTB_index)

projected_points_PTB_array_x = (m*points_PTB_array_y + points_PTB_array_x - m*ones(N_PTB)*c1)/(m^2 + 1);
projected_points_PTB_array_y = m*projected_points_PTB_array_x + c1 * ones(N_PTB);

## This array contains distance along the line
## This should be x coordinates on the line * sec(theta), where the slope is tan(theta) = m
distance_along_projected_line = projected_points_PTB_array_x * sqrt(m^2 +1)


include("generate_matrices2D.jl")
include("angmom.jl")

println("t = ",t)
println("t0 = ", t0)
println("m_0 = ", m_0)
println("Delta = ", Delta)


println("Lx = ", Lx)
println("Ly = ", Ly)

if(x_periodic^2 != x_periodic || y_periodic^2 != y_periodic)
    println("x_periodic and y_periodic must be 0 or 1")
    println("You entered x_periodic = ", x_periodic, " y_periodic = ", y_periodic)
    return
end


const_2D, CX2D, SX2D, CY2D, SY2D, CXCY2D, SXCY2D, CXSY2D = generate_matrices_2D(Lx,Ly,x_periodic,y_periodic)

# h_SSH = t*kron(SX2D,gamma_01) + 0*kron(SY2D,gamma_02) + kron(m_0*const_2D- t0*(CX2D + const_2D),gamma_03) + Delta*kron(const_2D,gamma_13) -mu*kron(const_2D,gamma_30)

h_SSH_without_local = t*kron(SX2D,gamma_01) + 0*kron(SY2D,gamma_02) - kron(t0*(CX2D),gamma_03)
# h_SSH_without_local = copy(h_SSH)

if(Hermitian_Check(h_SSH_without_local) == false)
    println("Error! Hamiltonian is not Hermitian")
    return
end

const_2D= Nothing
CX2D= Nothing
SX2D= Nothing
CY2D= Nothing
SY2D= Nothing
CXCY2D= Nothing
SXCY2D= Nothing
CXSY2D= Nothing

PTB_orbital_index = ones(4*N_PTB);
for ii = 1:N_PTB
    PTB_orbital_index[4*ii-3] = 4*PTB_index[ii]-3
    PTB_orbital_index[4*ii-2] = 4*PTB_index[ii]-2
    PTB_orbital_index[4*ii-1] = 4*PTB_index[ii]-1
    PTB_orbital_index[4*ii] = 4*PTB_index[ii]
end
println(PTB_orbital_index)

outside_orbital_index = setdiff(1:4*Lx*Ly,PTB_orbital_index);
println(outside_orbital_index)

NOrbitalsOutside = 4*Lx*Ly - 4*N_PTB

NOrbitalsInside = 4*N_PTB


H_PTB_renor = zeros(ComplexF64, NOrbitalsInside, NOrbitalsInside);

H_11 = zeros(ComplexF64, NOrbitalsInside, NOrbitalsInside);

H_22 = zeros(ComplexF64, NOrbitalsOutside, NOrbitalsOutside);

H_21 = zeros(ComplexF64, NOrbitalsOutside, NOrbitalsInside);

H_12 = zeros(ComplexF64, NOrbitalsInside, NOrbitalsOutside);

println("Matrix dimensions:")
println("H_11: $(size(H_11))")
println("H_22: $(size(H_22))")
println("H_21: $(size(H_21))")
println("H_12: $(size(H_12))")

for ii = 1:NOrbitalsInside
    for jj = 1:NOrbitalsInside
        H_11[ii,jj] = h_SSH_without_local[Int(PTB_orbital_index[ii]), Int(PTB_orbital_index[jj])]
    end
end

Threads.@threads for ii = 1:NOrbitalsOutside
    for jj = 1:NOrbitalsOutside
        H_22[ii,jj] = h_SSH_without_local[Int(outside_orbital_index[ii]), Int(outside_orbital_index[jj])]
    end
end

Threads.@threads for ii = 1:NOrbitalsInside
    for jj = 1:NOrbitalsOutside
        H_12[ii,jj] = h_SSH_without_local[Int(PTB_orbital_index[ii]), Int(outside_orbital_index[jj])]
    end
end

Threads.@threads for ii = 1:NOrbitalsOutside
    for jj = 1:NOrbitalsInside
        H_21[ii,jj] = h_SSH_without_local[Int(outside_orbital_index[ii]), Int(PTB_orbital_index[jj])]
    end
end

h_SSH= Nothing
h_SSH_without_local = Nothing

H_PTB_renor_before_adding_local = H_11 - Hermitian(H_12*inv(Hermitian(H_22 + 10^(-8)*Matrix(1.0I, NOrbitalsOutside, NOrbitalsOutside)))*H_21);

### generate identity Matrix in PTB space
identity_PTB = zeros(N_PTB,N_PTB);
for ii = 1:N_PTB
    identity_PTB[ii,ii] = 1.0;
end

# h_SSH = t*kron(SX2D,gamma_01) + 0*kron(SY2D,gamma_02) + kron(m_0*const_2D- t0*(CX2D + const_2D),gamma_03) + Delta*kron(const_2D,gamma_13) -mu*kron(const_2D,gamma_30)
H_PTB_renor = H_PTB_renor_before_adding_local + m_0 * kron(identity_PTB, gamma_03) - t0*kron(identity_PTB, gamma_03) + Delta*kron(identity_PTB, gamma_13) - mu*kron(identity_PTB, gamma_30);

H_11= Nothing
H_12= Nothing
H_21= Nothing
H_22= Nothing

if(Hermitian_Check(H_PTB_renor) == false)
    println("Error! Hamiltonian is not Hermitian")
    return
end

## Store the x and y coordinates of the sites in PTB
YList_PTB = ceil.(PTB_index/Lx);
XList_PTB = PTB_index - (YList_PTB - ones(N_PTB))*Lx;

## Kronecker product to get the orbitals
XListKron_PTB = kron(XList_PTB,[1,1,1,1]);
YListKron_PTB = kron(YList_PTB,[1,1,1,1]);

X1_PTB = diagm(XListKron_PTB);
Y1_PTB = diagm(YListKron_PTB);

(energy_eigenvalues_PTB, eigenstates_PTB) = eigen(H_PTB_renor);
#println(size(eigenstates))
filled_eigenstates_PTB = eigenstates_PTB[:,1:2*N_PTB]

## We create a projector P, which projects to the space of filled eigenstates (half-filled)
P_PTB = conj(filled_eigenstates_PTB) * transpose(filled_eigenstates_PTB)
## Q projects to the empty eigenstates
Q_PTB = kron(diagm(ones(N_PTB)),kron(eye2,eye2)) - P_PTB
W_PTB = kron(diagm(ones(N_PTB)),gamma_12)
## We define the local Chern operators, whose diagonal elements in the Wannier basis are the local Chern numbers per orbital
local_pol_PTB = W_PTB*(Q_PTB*X1_PTB*P_PTB + P_PTB*X1_PTB*Q_PTB);

local_pol_PTBSiteWiseList_PTB = zeros(N_PTB)*im;

## Here we add the two chern numbers for the two orbitals
for ii = 1:N_PTB
    local_pol_PTBSiteWiseList_PTB[ii] = local_pol_PTB[4*ii - 3,4*ii - 3] + local_pol_PTB[4*ii - 2,4*ii - 2] + local_pol_PTB[4*ii - 1,4*ii - 1] + local_pol_PTB[4*ii,4*ii]
end
local_pol_PTBSiteWiseList_PTB = real(local_pol_PTBSiteWiseList_PTB)

Gap_PTB = 2*minimum(abs.(energy_eigenvalues_PTB))

function state_to_real_space_LDoS(v::Array{ComplexF64})
    n_sites = Int(size(v)[1]/4)
    prob_dist = zeros(n_sites)
    for ii = 1:n_sites
        prob_dist[ii] = abs(v[4*ii-3])^2 + abs(v[4*ii-2])^2 +  abs(v[4*ii-1])^2 +  abs(v[4*ii])^2
    end
    return prob_dist
end 

#boundary_state_nearest_zero = eigenstates_PTB[:,N_PTB];
probability_boundary_state = state_to_real_space_LDoS(eigenstates_PTB[:,2*N_PTB]) + state_to_real_space_LDoS(eigenstates_PTB[:,2*N_PTB+1]);

color_map = cgrad([RGB(1,1,1), RGB(0,0,1), RGB(1,0,0)])

# Calculate alphas based on y values
function alpha_function(p, pmax)
    if p < pmax/10
        return 0
    else
        return p/pmax
    end
end
alphas = alpha_function.(probability_boundary_state, maximum(probability_boundary_state))

#plt_PTB_BD_states = scatter(distance_along_projected_line, probability_boundary_state)
plt_PTB_BD_states = scatter(distance_along_projected_line, ones(N_PTB), 
                            zcolor=probability_boundary_state, legend=false, xlabel="x", colorbar_title="Value",
                            c=color_map, ms=15, seriesalpha=alphas, markerstrokewidth=0, grid=false,
                            yaxis=false)

display(plt_PTB_BD_states)

plt_PTB_local_marker_sitewise = scatter(1:N_PTB, local_pol_PTBSiteWiseList_PTB,ylims=(-3,3),yticks=range(-3,3,7))

### I generate an array which are the points closest to the middle line
middle_line_x_coordinates = 1:Lx;
middle_line_y_coordinates = round.(line_middle.(float(middle_line_x_coordinates)));

Indices_of_PTB_middle_points = zeros(Lx);
local_marker_PTB_list_middle_points = zeros(Lx);

for ii = 1:Lx
    for jj = 1:N_PTB
        if (XList_PTB[jj] == middle_line_x_coordinates[ii]) && (YList_PTB[jj] == middle_line_y_coordinates[ii])
            Indices_of_PTB_middle_points[ii] = PTB_index[jj];
            local_marker_PTB_list_middle_points[ii] = local_pol_PTBSiteWiseList_PTB[jj];
        end
    end
end


mode_chern = modes(round.(local_marker_PTB_list_middle_points))[1]

### We will plot the local Chern marker along the line y = Ly/2
plt_local_PTB = scatter(1:Lx, local_marker_PTB_list_middle_points, ylims=(-3,3), legend=:none,
                xlabel="x", ylabel="PTB Local Z2 topo marker at middle line", 
                title=string("Lx = ", string(Lx), ", t = ", string(t), ", t0 = ", string(t0), ", m_0 = ", string(m_0), ", Delta = ", string(Delta), ", mu = ", string(mu)),titlefontsize=10)
plt_local_PTB = plot!(1:Lx, mode_chern*ones(Lx), linestyle=:dash, thickness=2, linewidth=2)
display(plt_local_PTB)

# Define the folder path (relative or absolute)
folder_path = "data/topo_super/"  # Change to your desired folder
# Ensure the folder exists
isdir(folder_path) || mkdir(folder_path)

filename = "t0=$(t0)_t=$(t)_m_0=$(m_0)_mu=$(mu)_Delta=$(Delta)_x_periodic=$(x_periodic)_y_periodic=$(y_periodic)_L=$(L)_m=$(m)_c1=$(c1)_c2=$(c2).csv"
filename_without_extension = replace(filename, ".csv" => "")
isdir(string(folder_path,"weak_invariant/")) || mkdir(string(folder_path,"weak_invariant/"))

if exportdata==1
    CSV.write(string(folder_path,"weak_invariant/",filename), (; local_marker_PTB_list_middle_points),writeheader=false)
end

if exportdata == 1
    savefig(plt_local_PTB, string(folder_path, "weak_invariant/plots/central", filename_without_extension, ".png"))
    savefig(plt_PTB_local_marker_sitewise, string(folder_path, "weak_invariant/plots/all", filename_without_extension, ".png"))
end

    end
end