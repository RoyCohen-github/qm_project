## [SR_system](./SR_system.py) contains the objects that implemets the atomic states and their evolutoin.
The basic object is 'SR_system' which can be instanciate by:

    from SR_system import *
    new_system = SR_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)
    
You have to pass the initial S and M quantum number of the system, which are the  first and second argument.
You can also specify the total number of atoms in the system, but  is will default to '20'

### Time evolution:
To start the computatoin of the system according to the master equatoin you can envoke:

    # sim_duratoin: is the time you want to let the system evolve
    new_system.create_time_evolution(sim_duratoin)
    
calling the time evolution function multyple times will continue the evolutoin from the last knowen result.
It is best to pass time in the characteristic time unit of the system, this can be done like the following:

    sim_duratoin_unitless = 10
    new_system.create_time_evolution(sim_duratoin_unitless * new_system.natural_time_unit)
    
### Additional functions:
There are a lot of other usefull functions to help you visuallize the results fast, the best resource is the [generate_sims](./final_sims/generate.py) file.
The above file is an example of  how to use the different visualization fucntions

#### here is  a list of the all the functoins:
        create_time_evolution
        get_f_expectation_value
        get_mandel_Q_param
        plot_m_expectation_value
        plot_scatter_p_evolution
        plot_Gamma
        plot_Q
        stem_Pm
        plot_Pn
        plot_Pm
        stem_Pn_vs_poiss_equivalent
        animate_Pm
   

## In addition, this file contains an objects that implemets the evolutoin for equivalent classical system.
The usage is similar to the SR evolutoin, usage example can be found in [generate_sims](./final_sims/generate.py) file.
