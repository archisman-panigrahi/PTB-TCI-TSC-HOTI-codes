using LinearAlgebra
using DelimitedFiles
using Plots;
gr()
using LaTeXStrings

#Plots.default(show = true)
#pyplot()

include("generate_matrices1D.jl")
include("generate_matrices2D.jl")
include("angmom.jl")

function energy_eigenvalues_hscti(Lx::Int64, Ly::Int64, t::Float64, t_prime::Float64, m_0::Float64, x_periodic::Int64, y_periodic::Int64)
    #close("all")
    t0 = 1.0;

    println("t = ",t)
    println("t_prime = ", t_prime)
    println("t_0 = ", t0)
    println("m_0 = ", m_0)

    if(x_periodic^2 != x_periodic || y_periodic^2 != y_periodic)
        println("x_periodic and y_periodic must be 0 or 1")
        println("You entered x_periodic = ", x_periodic, " y_periodic = ", y_periodic)
        return
    end
    t_start = time()

    sigma_x, sigma_y, sigma_z = generate_spin_matrices(1/2)
    const_2D, CX2D, SX2D, CY2D, SY2D, CXCY2D, SXCY2D, CXSY2D = generate_matrices_2D(Lx,Ly,x_periodic,y_periodic)

    h_scti = t*kron((SX2D + CXSY2D),sigma_x) + t*kron(SY2D+SXCY2D,sigma_y) + kron((m_0-2*t_prime)*const_2D + t0*(CX2D + CY2D) + 2*t_prime*CXCY2D,sigma_z)
    if(Hermitian_Check(h_scti) == false)
        println("Error! Hamiltonian is not Hermitian")
        return
    end
    
    energies = eigvals(h_scti)

    t_done = time()
    println("time taken = ", t_done - t_start, "seconds")

    plt1 = scatter(1:2*Lx*Ly, energies)
    display(plt1)
end

function energy_eigenvalues_hscti_ky_good_quantum_number(Lx::Int64, t::Float64, t_prime::Float64, m_0::Float64, x_periodic::Int64, nvals_ky::Int64)
    #close("all")
    #gr()
    t0 = 1.0;
    println("t = ",t)
    println("t_prime = ", t_prime)
    println("t_0 = ", t0)
    println("m_0 = ", m_0)

    if(x_periodic^2 != x_periodic)
        println("x_periodic must be 0 or 1")
        println("You entered x_periodic = ", x_periodic)
        return
    end
    t_start = time()

    sigma_x, sigma_y, sigma_z = generate_spin_matrices(1/2)
    const_1D, CX1D, SX1D = generate_matrices_1D(Lx,x_periodic)

    ky_array = range(-pi,pi,nvals_ky)
    Ey_array = zeros(2*Lx, nvals_ky)

    h_scti = zeros(ComplexF64,4*Lx,4*Lx)
    
    E_array_long = zeros(4*Lx*nvals_ky, 1)
    ky_array_long = zeros(4*Lx*nvals_ky, 1)

    for jj = 1:nvals_ky
        h_scti = t*kron((SX1D + CX1D*sin(ky_array[jj])),sigma_x) + t*kron(sin(ky_array[jj])*const_1D+SX1D*cos(ky_array[jj]),sigma_y) + kron((m_0-2*t_prime)*const_1D + t0*(CX1D + cos(ky_array[jj])*const_1D) + 2*t_prime*CX1D*cos(ky_array[jj]),sigma_z)
        if(Hermitian_Check(h_scti) == false)
            println("Error! Hamiltonian is not Hermitian")
            return
        end
        Ey_array[:,jj] = eigvals(h_scti)

        ## For plotting, I make two long vectors, of length nvals_ky*(4*Lx)
        for ii = 1:(2*Lx)
            ky_array_long[(jj-1)*2*Lx+ii] = ky_array[jj]
            E_array_long[(jj-1)*2*Lx+ii] = Ey_array[ii,jj]
        end
    end

    t_done = time()
    println("time taken = ", t_done - t_start, "seconds")

    plt1 = scatter(ky_array_long,E_array_long,xlabel="ky", ylabel="E", title = string("t = ", string(t), ", t_prime = ", string(t_prime), ", m0 = ", string(m_0)))
    xlims!(minimum(ky_array_long), maximum(ky_array_long))
    ylims!(minimum(E_array_long), maximum(E_array_long))
    display(plt1)
end

function LocalChernMarker(Lx::Int64,Ly::Int64,t::Float64,t_prime::Float64,m_0::Float64)
    t0 = 1.0;

    println("t = ",t)
    println("t_prime = ", t_prime)
    println("t_0 = ", t0)
    #println("m_0 = ", m_0)

    t_start = time()

    sigma_x, sigma_y, sigma_z = generate_spin_matrices(1/2)
    const_2D, CX2D, SX2D, CY2D, SY2D, CXCY2D, SXCY2D, CXSY2D = generate_matrices_2D(Lx,Ly,1,1)

    eye2 = diagm([1,1])

    XList = kron(range(1,1,Ly),range(1,Lx,Lx))
    YList = kron(range(1,Ly,Ly),range(1,1,Lx))
    #println(XList)
    #println(YList)

    X1 = kron(diagm(XList),eye2)
    Y1 = kron(diagm(YList),eye2)

    h_scti = t*kron((SX2D + CXSY2D),sigma_x) + t*kron(SY2D+SXCY2D,sigma_y) + kron((m_0-2*t_prime)*const_2D + t0*(CX2D + CY2D) + 2*t_prime*CXCY2D,sigma_z)

    (energy_eigenvalues, eigenstates) = eigen(h_scti);
    #println(size(eigenstates))
    filled_eigenstates = eigenstates[:,1:Lx*Ly]

    ## We create a projector P, which projects to the space of filled eigenstates (half-filled)
    P = conj(filled_eigenstates) * transpose(filled_eigenstates)
    ## Q projects to the empty eigenstates
    Q = kron(const_2D,eye2) - P
    ## We define the local Chern operators, whose diagonal elements in the Wannier basis are the local Chern numbers per orbital
    ChernMatrix = -4*pi*imag(P * X1 * Q * Y1 * P);
    ChernMatrixDiagonalList = diag(ChernMatrix)
    ChernMatrixSiteWiseList = zeros(Lx*Ly)
    ## Here we add the two chern numbers for the two orbitals
    for ii = 1:Lx*Ly
        ChernMatrixSiteWiseList[ii] = ChernMatrixDiagonalList[2*(ii-1) + 1] + ChernMatrixDiagonalList[2*(ii-1) + 2]
    end

    println(ChernMatrixSiteWiseList)

    println("time taken for calculation = ", time() - t_start, " seconds")

    ### We will plot the local Chern marker along the line y = Ly/2
    plt1 = scatter(1:Lx, ChernMatrixSiteWiseList[(Int(round(Ly/2)) - 1)*Lx+1:(Int(round(Ly/2)) - 1)*Lx+Lx], ylims=(-3,3), xlabel="x", ylabel="Local Chern marker at y=L/2", title=string("Lx = ", string(Lx), ", t = ", string(t), ", t_prime = ", string(t_prime), ", m0 = ", string(m_0)))
    display(plt1)
end

function BottIndexHalfFilled_hscti(Lx::Int64, Ly::Int64, t::Float64, t_prime::Float64)
    t0 = 1.0;

    no_of_points = 21;

    println("t = ",t)
    println("t_prime = ", t_prime)
    println("t_0 = ", t0)
    #println("m_0 = ", m_0) h_scti = t*kron((SX2D + CXSY2D),gamma_1) + t*kron(SY2D+SXCY2D,gamma_2) + kron((m_array[a]-2*t_prime)*const_2D + t0*(CX2D + CY2D) + 2*t_prime*CXCY2D,gamma_3)


    t_start = time()

    sigma_x, sigma_y, sigma_z = generate_spin_matrices(1/2)
    const_2D, CX2D, SX2D, CY2D, SY2D, CXCY2D, SXCY2D, CXSY2D = generate_matrices_2D(Lx,Ly,1,1)

    eye2 = diagm([1,1])


    xarray2D = kron(kron(range(1,1,Ly),range(1,Lx,Lx)), [1,1])
    yarray2D = kron(kron(range(1,Ly,Ly), range(1,1,Lx)), [1,1])

    U_x = diagm(exp.(2 * pi * im * xarray2D/Lx));
    V_y = diagm(exp.(2 * pi * im * yarray2D/Ly));

    bott_index_store = zeros(no_of_points)
    m_array = range(-3,3,no_of_points)

    for a=1:no_of_points
        h_scti = t*kron((SX2D + CXSY2D),sigma_x) + t*kron(SY2D+SXCY2D,sigma_y) + kron((m_array[a]-2*t_prime)*const_2D + t0*(CX2D + CY2D) + 2*t_prime*CXCY2D,sigma_z)

        (energy_eigenvalues, eigenstates) = eigen(h_scti);
        #println(size(eigenstates))
        #display(scatter(1:2*Lx*Ly,energy_eigenvalues))
        filled_eigenstates = eigenstates[:,1:Lx*Ly]

        P = conj(filled_eigenstates) * transpose(filled_eigenstates)


        U = P * U_x * P + (kron(const_2D,eye2) - P);
        V = P * V_y * P + (kron(const_2D,eye2) - P);
        bott = U * V * U' * V'
        #println(size(bott))
        #bott_index = real((-im/(2*pi)) * sum(log.(eigvals(bott))))

        bott_index_store[a] = real((-im/(2*pi)) * sum(log.(eigvals(bott))))
        println(m_array[a],bott_index_store[a])
    end

    #println(bott_index)
    println("time taken for calculation = ", time() - t_start, " seconds")
    scatter(m_array, bott_index_store)
end