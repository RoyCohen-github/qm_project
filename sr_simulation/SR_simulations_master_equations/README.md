## [SR_system](./SR_system.py) contains the objects that implemets the atomic states and their evolutoin.
The basic object is 'SR_system' which can be instanciate by:

    new_system = SR_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)
    
You have to pass the initial S and M quantum number of the system, which are the  first and second argument.
You can also specify the total number of atoms in the system, but  is will default to '20'

#### time evolution:
To start the computatoin of the system according to the master equatoin you can envoke:

    # sim_duratoin: is the time you want to let the system evolve
    new_system.creat_time_evolution(sim_duratoin)
    
calling the time evolution function multyple times will continue the evolutoin from the last knowen result.
