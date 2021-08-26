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
        create_time_evolution:
        get_f_expectation_value(f):
            """ given functoin f(m) return its expectation value over time"""
        get_mandel_Q_param():
            """ returns mandel_Q_param over time """
        plot_m_expectation_value():
            """ plot m expectation value over time """
        plot_scatter_p_evolution():
            """ plot the system probability evolution over time, together with a line to indicate the expectation value over time """
            # plot example below
        plot_Gamma():
            """ plot emission rate value over time """
        plot_Q():
            """ plot mandel_Q_param over time """
        stem_Pm(t):
            """ stem the system m state probabilities at a given time """
        plot_Pn(t):
            """ plot the system n state (photonic) probabilities at a given time """
        plot_Pm(t):
            """ plot the system m state (atomic) probabilities at a given time """
        stem_Pn_vs_poiss_equivalent(t):
            """ plot the system n state (photonic) probabilities at a given time vs equivalent pois distribution with the same expectation value """
        animate_Pm(dest):
            """ generate a .gif file  to animate the porbability evolution of the system at 'dest' location """
   
scatter_p_plot_example:
[generate_sims](./final_sims/generate.py) file.

## In addition, this file contains an objects that implemets the evolutoin for equivalent classical system.
The usage is similar to the SR evolutoin, usage example can be found in [generate_sims](./final_sims/generate.py) file.
