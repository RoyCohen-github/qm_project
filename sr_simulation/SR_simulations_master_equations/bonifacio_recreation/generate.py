from SR_simmulations.SR_system import *

N = 200

SR_m_max = SR_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)

ntu = 30
SR_m_max.create_time_evolution(ntu * SR_m_max.natural_time_unit)

# t is in natural time units
t_in_ntu = [2.65, 3.97, 5.3, 6.62, 7.95]
for t in t_in_ntu:
    SR_m_max.plot_Pm(t * SR_m_max.natural_time_unit)
plt.savefig(f'bonifacio_Fig_1_recreation.svg')
plt.close()

a = 1 + 1
a = a + 1
