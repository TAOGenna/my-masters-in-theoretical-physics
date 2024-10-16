#=
This code consider the case of a Kitaev chain and a
quantum dot attached at one end. There is a repultion considered
between the quantum dot and the first site of the chain. 
This repultion is thought as a Hartree-Fock such that there 
will be shiftings in the occupation of site d and 1, as well
as shifting to the cooper pairs, etc. and all of this changes
are considered in the function ini!(...). The iteration
is because we are solving this problem with a self-consistent
approach. Finally, we plot the results to see how the majorana
fermions are affected. 
THE PLOT:
The plot involves seeing how the 
eigenenergies of the entire system change with respect to the value of e_d, 
the occupation energy of the quantum dot. 
=#
using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings
using DataFrames
using IterTools
using Printf
using PrettyTables
gr()
global N=40+2 # always an even number +2 (+2 because of the quantum dot)
#-------------Values for the Kitaev Chain-----------
# delta=0.2 | t=1 | mu=0 
delta=0.2
t=1
mu=0
#-------------Values for the quantum dot------
t_prime=0.2
e_d=-0.01
H_max=1.2
#----------------------------------------------
dt(x,y)= ==(x,y) # like a delta dirac
string_array2 = [" | d , c >", " | d , a >"," | 1 , c >", " | 1 , a >", " | 2 , c >", " | 2 , a >", " | 3 , c >", " | 3 , a >", " | 4 , c >", " | 4 , a >", " | 5 , c >", " | 5 , a >", " | 6 , c >", " | 6 , a >", " | 7 , c >", " | 7 , a >", " | 8 , c >", " | 8 , a >", " | 9 , c >", " | 9 , a >", " | 10 , c >", " | 10 , a >", " | 11 , c >", " | 11 , a >", " | 12 , c >", " | 12 , a >", " | 13 , c >", " | 13 , a >", " | 14 , c >", " | 14 , a >", " | 15 , c >", " | 15 , a >", " | 16 , c >", " | 16 , a >", " | 17 , c >", " | 17 , a >", " | 18 , c >", " | 18 , a >", " | 19 , c >", " | 19 , a >", " | 20 , c >", " | 20 , a >"]
#---------------Plot Settings-----------------
#myplot=plot(xlims=(-H_max,H_max))
#title!(L"$\Delta=0.2,t=1,t'=0.2,\mu=0, N=40$")
#xlabel!(L"$e_d$")
#ylabel!(L"$\langle c^\dagger_d c^\dagger_1\rangle$")
#---------------------------------------------

function ini!(H::Matrix{Float64},val,hopping_shift,ocd,oc1,vv)
    #display(vv)
    H[4,1] += -vv*val
    H[1,4] += -vv*val
    H[2,3] += vv*val
    H[3,2] += vv*val
    # <1,d|H|0,c>=-delta
    # <1,c|H|0,d>=delta
    #----------------------------------------------------
    H[1,3] += -hopping_shift*vv
    H[3,1] += -hopping_shift*vv
    H[2,4] +=  hopping_shift*vv
    H[4,2] +=  hopping_shift*vv
    #----------------------------------------------------
    
    H[1,1] +=  oc1*vv
    H[2,2] += -oc1*vv
    H[3,3] +=  ocd*vv
    H[4,4] += -ocd*vv
    
end
function generate_array(length::Int,H_max)
    return LinRange(-H_max, H_max, length)
end
function separator()
    display("-------------------------------------------------------------------------")
    display("-------------------------------------------------------------------------")
end

H=zeros(N,N)
legend_added = false
for i in 3:N#     j'
   for j in 3:N# j
        i1=ceil(i/2)
        j1=ceil(j/2)
        if !iseven(i) && !iseven(j)#     <i',c|H|j,c>
             H[i,j]=dt(j1-1,i1)*(-t)+dt(j1+1,i1)*(-t)
        elseif iseven(i) && !iseven(j)#  <i',a|H|j,c>
             H[i,j]=dt(i1,j1-1)*(delta)-dt(i1,j1+1)*(delta)
        elseif !iseven(i) && iseven(j)#  <i',c|H|j,a>
             H[i,j]=dt(i1,j1+1)*(delta)-dt(i1,j1-1)*(delta)
        else#                            <i',a|H|j,a>
             H[i,j]=dt(i1,j1+1)*(t)+dt(i1,j1-1)*(t)
        end
    end
end
through=generate_array(400,H_max)
#display(through)
H[3,1]+=-t_prime             # <1,c|H|0,c>=t_prime
H[4,2]+=t_prime              # <1,d|H|0,d>=-t_prime
H[1,3]+=-t_prime
H[2,4]+=t_prime
#H[4,1]=H[1,4]=-delta        # <1,d|H|0,c>=-delta
#H[2,3]=H[3,2]=delta            # <1,c|H|0,d>=delta
epsilon_expec=1e-4
num_elements = 2000
iterations=50
#V=1
V=0
initial_ev = [0.133222,0.315101,0.531797,0.489767]
legend_added=false
#---------------PLOT SETTINGS---------------------------#
min_energy=1# aumentar esto para el caso de N=40
myplot=plot(xlims=(-H_max,H_max),ylims=(-min_energy,min_energy))
#title!(L"$e_d=0.0,\quad 0.5,\quad 1.0$")
title!(L"$N=100+2,\Delta=0.2,t=1,t'=0.2, V=0$")
xlabel!(L"$e_d$")
ylabel!(L"$Eigenenergies$")
#-------------------------------------------------------#
for values_ed in through
    H[1,1]=values_ed                  # <0,c|H|0,c>
    H[2,2]=-values_ed                 # <0,d|H|0,d>
    matrix_plot=[]
    for x in 1:iterations
       ini!(H,initial_ev[1],initial_ev[2],(initial_ev[3]-0.5),(initial_ev[4]-0.5),V)
       espec=eigvals(H)
       espec_vec=eigvecs(H)
       #display(H)
       if x==iterations
            matrix_plot=deepcopy(H)
       end
       ini!(H,-initial_ev[1],-initial_ev[2],-(initial_ev[3]-0.5),-(initial_ev[4]-0.5),V)
       ev=[0.0,0.0,0.0,0.0]
       for i in 1:N
           if espec[i]<0.
              if abs(espec[i])<10^(-8) #we're using element-wise multiplication and addition
                  ev[1] += conj(espec_vec[1,i]) * espec_vec[4,i]*(1/2) # < c+d c+1 >
                  ev[2] += conj(espec_vec[1,i]) * espec_vec[3,i]*(1/2) # < c+d c_1 >
                  ev[3] += conj(espec_vec[1,i]) * espec_vec[1,i]*(1/2) # < c+d c_d >
                  ev[4] += conj(espec_vec[3,i]) * espec_vec[3,i]*(1/2) # < c+1 c_1 >
              else 
                  ev[1] += conj(espec_vec[1,i]) * espec_vec[4,i]
                  ev[2] += conj(espec_vec[1,i]) * espec_vec[3,i]
                  ev[3] += conj(espec_vec[1,i]) * espec_vec[1,i] # = < c+d c_d >
                  ev[4] += conj(espec_vec[3,i]) * espec_vec[3,i]
              end
           end
       end #at the end of this loop we get the new expectation values, they are supposed to match the initial ones
       global initial_ev = ev
    end
    #-----------------plotting eigenenergies---------------------------------------------
    espec_plot=eigvals(matrix_plot)
    for tmr in espec_plot
        if tmr>-min_energy && tmr<min_energy
            if !legend_added
                display(scatter!([values_ed], [tmr], markercolor = "blue", markerstrokecolor = "blue", markersize = 0.9,label="",markeralpha = 0.5))
                global legend_added=true
            else 
                display(scatter!([values_ed], [tmr], markercolor = "blue", markerstrokecolor = "blue", markersize = 0.9,label="",markeralpha = 0.5))
            end
        end
    end  
    #-------------------------------------------------------------------------------------   
end
#------------------------------------------------------------------------------------------------
#---------------------------DISPLAYING BASIC INFORMATION-----------------------------------------
#=
display("Paramenters value -> delta=$delta | t=$t | V=$V | N=$N")
string_array = ["< c+_d c+_1 >", "< c+_d c_1 >", "< c+_d c_d >", "< c+_1 c_1 >"]
initial_ev=reverse(initial_ev)
string_array=reverse(string_array)
df = DataFrame(strings = string_array, similar_values = initial_ev)
# Display the DataFrame
#display("the bigger values are at the end, this is the last site")
pretty_table(df, header = ["equation", "self-consisten value"], tf = tf_markdown)
=#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
