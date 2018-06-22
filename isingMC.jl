#Monte Carlo algorithm for the Ising model on an LxL square lattice
#Example 8.1 from Giordano/Nakanishi Computational Physics
#Emanuel Casiano-Diaz


using PhysConsts: k
using Plots

kB = k

#"Ising Model Hamiltonian"
#"H = J \sum_{<ij>} s_i s_j -  \mu H \sum_{i} s_i"

function createIsingLattice(L)
	#'''Create a random square Ising lattice of dimensions LxL
	s = rand((L,L))
	s[s .>= 0.5] = 1
	s[s .!= 1] = -1
	s
	end
 
function isingEnergy(s,J)
	#'''Count the total energy of the given 2D Ising Lattice'''
	E = 0 #Energy
	#Count horizontally
	for j in 1:size(s)[1]
		for i in 1:size(s)[2]
			if i == size(s)[2]
				E += J * s[j,i]*s[j,1] - mu * H * s[j,i] 	
			else
				E += J * s[j,i]*s[j,i+1] - mu * H * s[j,i] 	
			end
		end
	end

	#Count vertically (field contribution ignored; previously counted)
	for i in 1:size(s)[2]
		for j in 1:size(s)[1]
			if j == size(s)[1]
				E += -J * s[j,i]*s[1,i] 	
			else
				E += -J * s[j,i]*s[j+1,i] 	
			end
		end
	end
	E
end

function isingMagnetization(s)
	#'''Count the total energy of the given 2D Ising Lattice'''
	M = 0 #Magnetization
	#Count horizontally
	for j in 1:size(s)[1]
		for i in 1:size(s)[2]
			M += s[j,i]
		end
	end
	M
end

function isingMC(s,J,T,MC_steps)
	#Sweep through the lattice, count the nearest neighbor interaction
	#contribution of each spin to the energy, update spins accordingly
	
	println("MC Step: 0")
	display(s)
	println("")
	println("---------------------------------")
	h = size(s)[1] #height of the Ising lattice
	w = size(s)[2] #width
	for k in 1:MC_steps
		for j in 1:h
			for i in 1:w
				if j == 1 && i == 1      #Top left corner
					Eflip = -J*(s[1,1]*s[1,w] + s[1,1]*s[1,2] + s[1,1]*s[h,1] + s[1,1]*s[2,1])
				elseif j == 1 && i == w  #Top right corner
					Eflip = -J*(s[1,w]*s[1,w-1] + s[1,w]*s[1,1] + s[1,w]*s[h,w] + s[1,w]*s[2,w])
				elseif(j == h && i == 1)  #Bottom left corner
					Eflip = -J*(s[h,1]*s[h,w] + s[h,1]*s[h,2] + s[h,1]*s[h-1,1] + s[h,1]*s[1,1])
				elseif j == h && i == w  #Bottom right corner
					Eflip = -J*(s[h,w]*s[h,w-1] + s[h,w]*s[h,1] + s[h,w]*s[h-1,w] + s[h,w]*s[1,w])
				elseif j == 1           #Top row
					Eflip = -J*(s[1,i]*s[1,i-1] + s[1,i]*s[1,i+1] + s[1,i]*s[h,i] + s[1,i]*s[2,i])
				elseif j == h           #Bottom row
					Eflip = -J*(s[h,i]*s[h,i-1] + s[h,i]*s[h,i+1] + s[h,i]*s[h-1,i] + s[h,i]*s[1,i])
				elseif i == 1           #Leftmost column
					Eflip = -J*(s[j,1]*s[j,w] + s[j,1]*s[j,2] + s[j,1]*s[j-1,1] + s[j,1]*s[j+1,1])
				elseif i == w           #Rightmost column
					Eflip = -J*(s[j,w]*s[j,w-1] + s[j,w]*s[j,1] + s[j,w]*s[j-1,w] + s[j,w]*s[1,w])
				else                    #Inner spins
					Eflip = -J*(s[j,i]*s[j,i-1] + s[j,i]*s[j,i+1] + s[j,i]*s[j-1,i] + s[j,i]*s[j+1,i])
				end
	
				#Flip spin if energy required to flip it is negative (this reduces total energy)
				if Eflip <= 0
					s[j,i] = (-1)*s[j,i]
				else
					#Generate uniformly distributed random number for comparison with Eflip
					r = rand()
					if r <= exp(-Eflip/(kB*T))
						s[j,i] = (-1)*s[j,i]
					else
						s[j,i] = s[j,i]
					end
				end
			end
		end
	println("MC Step: ",@sprintf("%d", k))
	display(s)
	println("")
	println("---------------------------------")
	end
s
end

#Parameters	
T = 2.25
J = 1
L = 6
MC_steps = 10

s = createIsingLattice(L)

s = isingMC(s,J,T,MC_steps)
M = isingMagnetization(s)

#Create heatmap representing Square Ising Lattice
println("")
print("Magnetic Moment: ")
println(M)
