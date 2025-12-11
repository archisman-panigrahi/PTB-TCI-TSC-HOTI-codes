include("Hermitian_Check.jl")

function generate_matrices_double_dislocation2D(Lx::Int64, Ly::Int64, l1::Int64, l2::Int64, m_dislocation_lattice::Int64, periodic_x::Int64, periodic_y::Int64)
    # if(periodic_x != 0 || periodic_y != 0)
    #     println("Periodic dislocation matrices not implemented yet")
    #     return
    # end
    nSites = Lx*Ly + l2
    CX2D = zeros(nSites, nSites)
    CY2D = zeros(nSites, nSites)
    SX2D = zeros(ComplexF64,nSites, nSites)
    SY2D = zeros(ComplexF64,nSites, nSites)
    CXCY2D = zeros(nSites, nSites)
    SXCY2D = zeros(ComplexF64,nSites, nSites)
    CXSY2D = zeros(ComplexF64,nSites, nSites)

    ### Constant matrix
    const_2D = Diagonal(ones(nSites))

    # l1 denotes the row in which dislocation is placed.
    # l1 belongs to the set (1,2,... Ly -1)

######### kx matrices ###########
#kx matrices, lower block
for a = 1:(Lx - 1)
	for b = 0:(l1-1)
		SX2D[a + b*Lx, a + b*Lx + 1] = im/2;
		SX2D[a + b*Lx + 1, a + b*Lx] = -im/2;
		CX2D[a + b*Lx, a + b*Lx + 1] = 1/2;
		CX2D[a + b*Lx + 1, a + b*Lx] = 1/2;
	end
end

#kx matrices, middle block
for a = 1:Lx
	for b = 0:(l2 - 1)
		SX2D[a + l1*Lx + b * (Lx + 1), a + l1*Lx + b * (Lx + 1) + 1] = im/2;
		SX2D[a + l1*Lx + b * (Lx + 1) + 1, a + l1*Lx + b * (Lx + 1)] = -im/2;
		CX2D[a + l1*Lx + b * (Lx + 1), a + l1*Lx + b * (Lx + 1) + 1] = 1/2;
		CX2D[a + l1*Lx + b * (Lx + 1) + 1, a + l1*Lx + b * (Lx + 1)] = 1/2;
	end
end

#kx matrices, upper block
for a = 1:(Lx-1)
	for b = 0:(Ly - (l1 + l2) - 1)
		SX2D[a + (l1 + l2)*Lx + l2 + b * Lx, a + (l1 + l2)*Lx + l2 + b * Lx + 1] = im/2;
		SX2D[a + (l1 + l2)*Lx + l2 + b * Lx + 1, a + (l1 + l2)*Lx + l2 + b * Lx] = -im/2;
		CX2D[a + (l1 + l2)*Lx + l2 + b * Lx, a + (l1 + l2)*Lx + l2 + b * Lx + 1] = 1/2;
		CX2D[a + (l1 + l2)*Lx + l2 + b * Lx + 1, a + (l1 + l2)*Lx + l2 + b * Lx] = 1/2;
	end
end

### Periodic Boundary Condition ###
if(periodic_x == 1)
#lower block
    for b = 0:(l1 - 1)
    	SX2D[(b+1)*Lx, b*Lx + 1] = im/2;
    	SX2D[b*Lx + 1, (b+1)*Lx] = -im/2;
    	CX2D[(b+1)*Lx, b*Lx + 1] = 1/2;
    	CX2D[b*Lx + 1, (b+1)*Lx] = 1/2;
    end

#middle block
    for b = 0:(l2 - 1)
    	SX2D[l1*Lx + (b+1)*(Lx+1), l1*Lx + b*(Lx+1) + 1] = im/2;
    	SX2D[l1*Lx + b*(Lx+1) + 1, l1*Lx + (b+1)*(Lx+1)] = -im/2;
    	CX2D[l1*Lx + (b+1)*(Lx+1), l1*Lx + b*(Lx+1) + 1] = 1/2;
    	CX2D[l1*Lx + b*(Lx+1) + 1, l1*Lx + (b+1)*(Lx+1)] = 1/2;
    end

#upper block
    for b = 0:(Ly - (l1 + l2) -1)
    	SX2D[(b+1)*Lx + Lx*(l1+l2) + l2, b*Lx + 1 + Lx*(l1+l2) + l2] = im/2;
    	SX2D[b*Lx + 1 + Lx*(l1+l2) + l2, (b+1)*Lx + Lx*(l1+l2) + l2] = -im/2;
    	CX2D[(b+1)*Lx + Lx*(l1+l2) + l2, b*Lx + 1 + Lx*(l1+l2) + l2] = 1/2;
    	CX2D[b*Lx + 1 + Lx*(l1+l2) + l2, (b+1)*Lx + Lx*(l1+l2) + l2] = 1/2;
    end

end


########## ky matrices ###########
#ky matrices, lower block (without connectors)
for a = 1:Lx
	for b = 0:(l1-2)
		SY2D[a + b*Lx, a + (b+1)*Lx] = im/2;
		SY2D[a + (b+1)*Lx, a + b*Lx] = -im/2;
		CY2D[a + b*Lx, a + (b+1)*Lx] = 1/2;
		CY2D[a + (b+1)*Lx, a + b*Lx] = 1/2;
	end
end

#middle block, without connectors
for a = 1:(Lx+1)
	for b = 0:(l2-2)
		SY2D[a + b*(Lx+1) + l1*Lx, a + (b+1)*(Lx+1) + l1*Lx] = im/2;
		SY2D[a + (b+1)*(Lx+1) + l1*Lx, a + b*(Lx+1) + l1*Lx] = -im/2;
		CY2D[a + b*(Lx+1) + l1*Lx, a + (b+1)*(Lx+1) + l1*Lx] = 1/2;
		CY2D[a + (b+1)*(Lx+1) + l1*Lx, a + b*(Lx+1) + l1*Lx] = 1/2;
	end
end

#ky matrices, upper block (without connectors)
for a = 1:Lx
	for b = 0:(Ly - (l1 + l2) - 2)
		SY2D[a + b*Lx + (l1+l2)*Lx + l2, a + (b+1)*Lx + (l1+l2)*Lx + l2] = im/2;
		SY2D[a + (b+1)*Lx + (l1+l2)*Lx + l2, a + b*Lx + (l1+l2)*Lx + l2] = -im/2;
		CY2D[a + b*Lx + (l1+l2)*Lx + l2, a + (b+1)*Lx + (l1+l2)*Lx + l2] = 1/2;
		CY2D[a + (b+1)*Lx + (l1+l2)*Lx + l2, a + b*Lx + (l1+l2)*Lx + l2] = 1/2;
	end
end

#ky matrices, lower connectors
for a = 1:m_dislocation_lattice
	SY2D[a + (l1-1)*Lx, a + l1*Lx] = im/2;
	SY2D[a + l1*Lx, a + (l1-1)*Lx] = -im/2;
	CY2D[a + (l1-1)*Lx, a + l1*Lx] = 1/2;
	CY2D[a + l1*Lx, a + (l1-1)*Lx] = 1/2;
end

for a = 1:(Lx - m_dislocation_lattice)
	SY2D[a + m_dislocation_lattice + (l1-1)*Lx, a + m_dislocation_lattice + 1 + l1*Lx] = im/2;
	SY2D[a + m_dislocation_lattice + 1 + l1*Lx, a + m_dislocation_lattice + (l1-1)*Lx] = -im/2;
	CY2D[a + m_dislocation_lattice + (l1-1)*Lx, a + m_dislocation_lattice + 1 + l1*Lx] = 1/2;
	CY2D[a + m_dislocation_lattice + 1 + l1*Lx, a + m_dislocation_lattice + (l1-1)*Lx] = 1/2;
end
#ky matrices, upper connectors
for a = 1:m_dislocation_lattice
	SY2D[a + (l1 + l2)*Lx + l2 - (Lx + 1), a + (l1 + l2)*Lx + l2] = im/2;
	SY2D[a + (l1 + l2)*Lx + l2, a + (l1 + l2)*Lx + l2 - (Lx + 1)] = -im/2;
	CY2D[a + (l1 + l2)*Lx + l2 - (Lx + 1), a + (l1 + l2)*Lx + l2] = 1/2;
	CY2D[a + (l1 + l2)*Lx + l2, a + (l1 + l2)*Lx + l2 - (Lx + 1)] = 1/2;
end

for a = 1:(Lx - m_dislocation_lattice)
	SY2D[a + (m_dislocation_lattice+1) + (l1 + l2)*Lx + l2 - (Lx + 1), a + m_dislocation_lattice + (l1 + l2)*Lx + l2] = im/2;
	SY2D[a + m_dislocation_lattice + (l1 + l2)*Lx + l2, a + (m_dislocation_lattice+1) + (l1 + l2)*Lx + l2 - (Lx + 1)] = -im/2;
	CY2D[a + (m_dislocation_lattice+1) + (l1 + l2)*Lx + l2 - (Lx + 1), a + m_dislocation_lattice + (l1 + l2)*Lx + l2] = 1/2;
	CY2D[a + m_dislocation_lattice + (l1 + l2)*Lx + l2, a + (m_dislocation_lattice+1) + (l1 + l2)*Lx + l2 - (Lx + 1)] = 1/2;
end

####Periodic Boundary condition#####
if(periodic_y == 1)

    for a = 1:Lx
	    SY2D[Lx*(Ly-1) + l2 + a, a] = im/2;
	    SY2D[a, Lx*(Ly-1) + l2 + a] = -im/2;
	    CY2D[Lx*(Ly-1) + l2 + a, a] = 1/2;
	    CY2D[a, Lx*(Ly-1) + l2 + a] = 1/2;
    end

end

### CXCY Matrices
for ii = 1:Lx-1
	for jj = 0:l1-2
		CXCY2D[ii + jj* Lx, (ii+1) + (jj+1)*Lx] = 1/4;
		CXCY2D[(ii+1) + (jj+1)*Lx, ii + jj* Lx] = 1/4;
	end
end

for ii = 1:Lx-1
	for jj = 0:l1-2
		CXCY2D[ii + 1 + jj*Lx, ii + (jj+1)*Lx] = 1/4;
		CXCY2D[ii + (jj+1)*Lx, ii + 1 + jj*Lx] = 1/4;
	end
end

for ii = 1:Lx
	for jj = 0:l2-2
		CXCY2D[ii + jj*(Lx+1) + l1*Lx, (ii+1) + (jj+1)*(Lx+1) + l1*Lx] = 1/4
		CXCY2D[(ii+1) + (jj+1)*(Lx+1) + l1*Lx, ii + jj*(Lx+1) + l1*Lx] = 1/4
	end
end

for ii = 1:Lx
	for jj = 0:l2-2
		CXCY2D[ii + 1 + jj*(Lx+1) + l1*Lx, ii + (jj+1)*(Lx+1) + l1*Lx] = 1/4
		CXCY2D[ii + (jj+1)*(Lx+1) + l1*Lx, ii + 1 + jj*(Lx+1) + l1*Lx] = 1/4
	end
end

for ii = 1:Lx-1
	for jj = 0:Ly-(l1+l2)-2
		CXCY2D[ii + jj*Lx + (l1 + l2)*Lx + l2, (ii+1) + (jj+1)*Lx + (l1 + l2)*Lx +l2] = 1/4;
		CXCY2D[(ii+1) + (jj+1)*Lx + (l1 + l2)*Lx + l2, ii + jj*Lx + (l1 + l2)*Lx + l2] = 1/4;
	end
end

for ii = 1:Lx-1
	for jj = 0:Ly-(l1+l2)-2
		CXCY2D[ii+1 + jj*Lx + (l1 + l2)*Lx + l2, ii + (jj+1)*Lx + (l1+l2)*Lx + l2] = 1/4
		CXCY2D[ii + (jj+1)*Lx + (l1+l2)*Lx + l2, ii+1 + jj*Lx + (l1 + l2)*Lx + l2] = 1/4
	end
end

for ii = 1:m_dislocation_lattice
	CXCY2D[ii + (l1-1)*Lx, ii+1+l1*Lx] = 1/4;
	CXCY2D[ii+1+l1*Lx, ii + (l1-1)*Lx] = 1/4;
end

for ii = 1:Lx-m_dislocation_lattice-1
	CXCY2D[ii+m_dislocation_lattice + (l1-1)*Lx, ii+2 + m_dislocation_lattice + l1*Lx] = 1/4;
	CXCY2D[ii+2 + m_dislocation_lattice + l1*Lx, ii+m_dislocation_lattice + (l1-1)*Lx] = 1/4;
end

for ii = 1:m_dislocation_lattice-1
	CXCY2D[ii + 1 + (l1-1)*Lx, ii + l1*Lx] = 1/4;
	CXCY2D[ii + l1*Lx, ii + 1 + (l1-1)*Lx] = 1/4;
end

for ii = 1:Lx-m_dislocation_lattice
	CXCY2D[ii + m_dislocation_lattice + (l1-1)*Lx, ii + m_dislocation_lattice + l1*Lx] = 1/4
	CXCY2D[ii + m_dislocation_lattice +  l1*Lx, ii + m_dislocation_lattice + (l1-1)*Lx] = 1/4
end

for ii = 1:m_dislocation_lattice -1
	CXCY2D[ii + Lx*(l1+l2) + l2 - (Lx+1), ii + 1 + Lx*(l1+l2) + l2] = 1/4;
	CXCY2D[ii + 1 + Lx*(l1+l2) + l2, ii + Lx*(l1+l2) + l2 - (Lx+1)] = 1/4;
end

for ii = 1:(Lx+1)-m_dislocation_lattice-1
	CXCY2D[ii + m_dislocation_lattice + Lx*(l1+l2) + l2 - (Lx+1), ii + m_dislocation_lattice + Lx*(l1+l2) +l2] = 1/4
	CXCY2D[ii + m_dislocation_lattice + Lx*(l1+l2) + l2, ii + m_dislocation_lattice + Lx*(l1+l2) + l2 - (Lx+1)] = 1/4
end

for ii = 1:m_dislocation_lattice
	CXCY2D[ii+1 + Lx*(l1+l2) + l2 - (Lx+1), ii + Lx*(l1+l2) + l2] = 1/4;
	CXCY2D[ii + Lx*(l1+l2) + l2, ii+1 + Lx*(l1+l2) + l2 - (Lx+1)] = 1/4;
end

for ii = 1:Lx-m_dislocation_lattice-1
	CXCY2D[ii + m_dislocation_lattice + 2 + Lx*(l1+l2) + l2 - (Lx+1), ii + m_dislocation_lattice + Lx*(l1+l2) + l2] = 1/4
	CXCY2D[ii + m_dislocation_lattice + Lx*(l1+l2) + l2, ii + m_dislocation_lattice + 2 + Lx*(l1+l2) + l2 - (Lx+1)] = 1/4
end

if (periodic_x == 1 && periodic_y == 1)
	for jj = 1:l1
		CXCY2D[jj*Lx, jj*Lx + 1] = 1/4
		CXCY2D[jj*Lx + 1, jj*Lx] = 1/4
	end

	for jj = 1:l2
		CXCY2D[jj*(Lx+1) + l1*Lx, jj*(Lx+1) + l1*Lx + 1] = 1/4;
		CXCY2D[jj*(Lx+1) + l1*Lx + 1, jj*(Lx+1) + l1*Lx] = 1/4;
	end

	for jj = 1:Ly-(l1+l2) -1
		CXCY2D[jj*Lx + (l1+l2)*Lx + l2, jj*Lx + (l1+l2)*Lx + l2 + 1] = 1/4;
		CXCY2D[jj*Lx + (l1+l2)*Lx + l2 + 1, jj*Lx + (l1+l2)*Lx + l2] = 1/4;
	end

	for jj = 1:l1-1
		CXCY2D[(jj+1)*Lx, (jj-1)*Lx + 1] = 1/4
		CXCY2D[(jj-1)*Lx + 1, (jj+1)*Lx] = 1/4
	end

	CXCY2D[(l1+1)*Lx + 1, (l1-1)*Lx + 1] = 1/4
	CXCY2D[(l1-1)*Lx + 1, (l1+1)*Lx + 1] = 1/4

	for jj = 1:l2-1
		CXCY2D[(jj+1)*(Lx+1) + l1*Lx, (Lx+1)*(jj-1) + l1*Lx + 1] = 1/4
		CXCY2D[(Lx+1)*(jj-1) + l1*Lx + 1, (jj+1)*(Lx+1) + l1*Lx] = 1/4
	end

	CXCY2D[(l1+l2+1)*Lx + l2, (l1+l2)*Lx + l2 - (Lx+1)] = 1/4
	CXCY2D[(l1+l2)*Lx + l2 - (Lx+1), (l1+l2+1)*Lx + l2] = 1/4

	for jj=1:Ly-(l1+l2)-1
		CXCY2D[(jj+1)*Lx + (l1+l2)*Lx + l2, (jj-1)*Lx + (l1+l2)*Lx + l2 + 1] = 1/4
		CXCY2D[(jj-1)*Lx + (l1+l2)*Lx + l2 + 1, (jj+1)*Lx + (l1+l2)*Lx + l2] = 1/4
	end

	for ii = 1:Lx-1
		CXCY2D[Lx*Ly + l2 - Lx + ii, ii + 1] = 1/4
		CXCY2D[ii + 1, Lx*Ly + l2 - Lx + ii] = 1/4

		CXCY2D[(Ly-1)*Lx + l2 + ii + 1, ii] = 1/4
		CXCY2D[ii, (Ly-1)*Lx + l2 + ii + 1] = 1/4
	end

end


    if (Hermitian_Check(CX2D)==1 && Hermitian_Check(CY2D)==1 && Hermitian_Check(SX2D)==1 && Hermitian_Check(SY2D)==1 && Hermitian_Check(CXCY2D)==1 && Hermitian_Check(SXCY2D)==1 && Hermitian_Check(CXSY2D)==1 && Hermitian_Check(const_2D)==1)
        println("Verified: Building Block Matrices are Hermitian")
        return const_2D, CX2D, SX2D, CY2D, SY2D, CXCY2D
    else
        println("Error!!! Building Block Matrices are not Hermitian!")
    end
end