using Plots
using LinearAlgebra

#This animation is an approach to describe the motion of a double pendulum in vacume using lagrangian mechanics. It was a project for the first semester of physics at Universidad EIA. 
#I used Runge-Kutta method for ordinary, second order differential equations.
#The initial conditions are established within the function main, and can be altered arbitrarely to produce different behaviors on the pendulums. 
#The masses of each pendulum and the length of the rods can also be changed.
#the cartesian plane has been inverted so it is easier to work with 
#if you want to remove the line for the mass1 all you have to do is change the plot command for the graph of each mass, thus having plot! fro mass1 and plot for mass 2


#global constants
g = 9.81

r1, r2, m1, m2= 1, 1, 1, 1 #length of rod1, length of rod 2, mass of body 1, mass of body 2

#Function to relate each initial condition to it's derivative
function derivatives(t, u)

#initial conditions: angle1, angle2, angular vel1, angular vel2
    ang1, ang2, omg1, omg2 = u[1:4]

#some useful statements
    d_a = ang1 - ang2
    z = m1 + m2

#Factors for angular acceleration of rod1
    num1_1 = g * m2 * sin(ang1) * cos(d_a)
    num1_2 = r1 * m2 * omg1^2 * sin(d_a) * cos(d_a)
    num1_3 = z * g * sin(ang1)
    num1_4 = m2 * r2 * omg2^2 * sin(d_a)
    denom1 = r1 * (z - m2 * cos(d_a)^2)

#factors for angular acceleration of rod2
    num2_1 = m2 * r2 * omg2^2 * sin(d_a) * cos(d_a);
    num2_2 = z * g * sin(ang1) * cos(d_a);
    num2_3 = z * r1 * omg1^2 * sin(d_a);
    num2_4 = z * g * sin(ang2);
    denom2 = r2 * (z + cos(d_a)^2);

    acc1 = (num1_1 - num1_2 - num1_3 - num1_4) / denom1 #angular acceleration for body 1

    acc2 = (num2_1 + num2_2 + num2_3 - num2_4) / denom2 #angular acceleration for body 2

    du = Vector([omg1, omg2, acc1, acc2]) #derivatives vector

    return du

end

#numerical approximation to the solution using Runge_Kutta
function Runge_Kutta(u, t, dt)

    k1 = derivatives(t, u);
    k2 = derivatives(t + dt/2, (u[1] + dt/2*k1[1], u[2] + dt/2*k1[2], u[3] + dt/2*k1[3], u[4] + dt/2*k1[4]));
    k3 = derivatives(t + dt/2, (u[1] + dt/2*k2[1], u[2] + dt/2*k2[2], u[3] + dt/2*k2[3], u[4] + dt/2*k2[4]));
    k4 = derivatives(t + dt, (u[1] + dt*k3[1], u[2] + dt*k3[2], u[3] + dt*k3[3], u[4] + dt*k3[4]));

    for i in 1:4
        u[i] += dt/6 * (k1[i]+2*k2[i]+2*k3[i]+k4[i])
    end
end

#graphing and initial conditions 
function main()

#intiial conditions

  #Angles for each rod
    ang1 = pi/2 
    ang2 = 0    
  #Angular velocities for each rod
    om1 = 0     
    om2 = 0     

    u0 = Vector([ang1, ang2, om1, om2]) #vector for initial conditions

    #time parameters
    t0 = 0
    dt = 0.05
    t_max = 25

    #making sure the steps of the iteration are integers to avoid limiting dt
    steps = convert(Int64, t_max/dt)

    #vectors for all position from t0 to t_max
    pos1 = Vector{Array{Float64}}(undef, steps)
    pos2 = Vector{Array{Float64}}(undef, steps)

    #Iteration for the position of each pendulum starting at t0
    for i in 1:steps 

        #Cartesian coordinates for each mass
        Cart1 = Vector([r1*sin(u0[1]), r1*cos(u0[1])])
        Cart2 = Cart1 + Vector([r2*sin(u0[2]), r2*cos(u0[2])])

        #Store each position per iteration
        pos1[i] = copy(Cart1) 
        pos2[i] = copy(Cart2)

        Runge_Kutta(u0, t0, dt)   #Update u0

        t0 += dt  #update t0

    end 

@gif for i in 1:steps
        
        plot(map(p -> p[1], pos1[1:i]), map(p -> p[2], -pos1[1:i]), xlims=(-2.5, 2.5), ylims=(-2.5, 0.5), label="m1", linecolor=:red) #Draw the first pendulum
        plot!(map(p -> p[1], pos2[1:i]), map(p -> p[2], -pos2[1:i]),  xlims=(-2.5, 2.5), ylims=(-2.5, 0.5), label="m2", linecolor=:blue) #Draw the second pendulum with traceline
        plot!([pos1[i][1], pos2[i][1]],[-pos1[i][2], -pos2[i][2]], label = false, lw = 1, color=:black) #Draw line from mass1 to mass2
        plot!([0, pos1[i][1]], [0, -pos1[i][2]], label = false, lw = 1, color=:black) #draw line from the origin to mass1
        scatter!([pos1[i][1]], [-pos1[i][2]], label="", markercolor=:red) #Draw mass1
        scatter!([pos2[i][1]], [-pos2[i][2]], label="", markercolor=:blue) #Draw mass2
        scatter!([0], [0], label = false, color=:black)

    end fps = 15 #fps to play with the speed of the animation

end

anim = main() #run animation
