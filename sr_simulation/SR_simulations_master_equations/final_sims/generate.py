from SR_simulations_master_equations.SR_system import *

N = 40

SR_m_max = SR_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)
SR_m_0 = SR_system(N * SR_system.spin, m0=0, N=N)
CL_m_max = CL_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)
CL_m_0 = CL_system(N * SR_system.spin, m0=0, N=N)

ntu = 30
SR_m_max.create_time_evolution(ntu * SR_m_max.natural_time_unit)
SR_m_0.create_time_evolution(ntu * SR_m_0.natural_time_unit)
CL_m_max.create_time_evolution(20 * ntu * CL_m_max.natural_time_unit)
CL_m_0.create_time_evolution(20 * ntu * CL_m_0.natural_time_unit)


##########################################################################################
##########################################################################################


print('gen SR_probability_evolution_m0_max')
SR_m_max.plot_scatter_p_evolution()
plt.savefig('SR_probability_evolution_m0_max.svg')
plt.close()

print('gen SR_probability_evolution_m0_0')
SR_m_0.plot_scatter_p_evolution()
plt.savefig('SR_probability_evolution_m0_0.svg')
plt.close()

print('gen CL_probability_evolution_m0_max')
CL_m_max.plot_scatter_p_evolution()
plt.savefig('CL_probability_evolution_m0_max.svg')
plt.close()

print('gen CL_probability_evolution_m0_0')
CL_m_0.plot_scatter_p_evolution()
plt.savefig('CL_probability_evolution_m0_0.svg')
plt.close()

##########################################################################################
##########################################################################################

print('gen m_sr_evolution_m0_max_vs_m0_0')
SR_m_max.plot_m_expectation_value(normlize_m_range=False, label='m_0=max')
SR_m_0.plot_m_expectation_value(normlize_m_range=False, label='m_0=0')
plt.legend()
plt.savefig('m_sr_evolution_m0_max_vs_m0_0.svg')
plt.close()

print('gen m_evolution_m0_max_sr_vs_cl')
SR_m_max.plot_m_expectation_value(normlize_m_range=False, label='m_0=max (with SR)')
CL_m_max.plot_m_expectation_value(normlize_m_range=False, label='m_0=max (no SR)')
plt.legend()
plt.savefig('m_evolution_m0_max_sr_vs_cl.svg')
plt.close()


##########################################################################################
##########################################################################################

print('gen gamma_per_m_for_m0_max_vs_m0_0_vs_cl')
SR_m_max.plot_Gamma(as_functoin_of_m=True, label='m=max')
SR_m_0.plot_Gamma(as_functoin_of_m=True, label='m=0')
CL_m_max.plot_Gamma(as_functoin_of_m=True, label='(without SR)')
plt.legend()
plt.savefig('gamma_per_m_for_m0_max_vs_m0_0_vs_cl.svg')
plt.close()


##########################################################################################
##########################################################################################


print('gen Q_param_m0_max')
SR_m_max.plot_Q(label='with sr', as_functoin_of_m=True)
CL_m_max.plot_Q(label='without sr', as_functoin_of_m=True)
plt.legend()
plt.title('$m_0 = max$')
plt.savefig('Q_param_m0_max.svg')
plt.close()


print('gen Q_param_m0_0')
SR_m_0.plot_Q(label='with sr', as_functoin_of_m=True)
CL_m_0.plot_Q(label='without sr', as_functoin_of_m=True)
plt.legend()
plt.title('$m_0 = 0$')
plt.savefig('Q_param_m0_0.svg')
plt.close()


##########################################################################################
##########################################################################################


selected_system = SR_m_max

t_vals = [selected_system.Pm.t[selected_system.get_mandel_Q_param() == selected_system.get_mandel_Q_param().max()].iloc[0],
          selected_system.Pm.t[selected_system.get_m_expectation_value().abs().min() == selected_system.get_m_expectation_value().abs()].iloc[0],
          selected_system.Pm.t[(selected_system.get_mandel_Q_param()-(-.5)).abs() == (selected_system.get_mandel_Q_param()-(-.5)).abs().min()].iloc[0],
          ]
t_vals = [float(val) for val in t_vals]
t_proposes = ['Q_max',
             'm_is_0',
             'neqative_Q'
            ]

for t_val, t_propose in zip(t_vals, t_proposes):
    selected_system.stem_Pn_vs_poiss_equivalent(t_val)
    Q = np.interp(t_val, selected_system.Pm.t, selected_system.get_mandel_Q_param()).round(2)
    m = np.interp(t_val, selected_system.Pm.t, selected_system.get_m_expectation_value()).round(2)
    plt.title(f'$m_0 = {selected_system.m0}, Q = {Q}, m = {m}$')
    print(f'gen Pn_vs_poiss_equivalent_{t_propose}')
    plt.savefig(f'Pn_vs_poiss_equivalent_{t_propose}.svg')
    plt.close()



