from SR_simmulations.SR_system import *

N = 40

# SR_m_max = SR_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)
# SR_m_0 = SR_system(N * SR_system.spin, m0=0, N=N)
# CL_m_max = CL_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)
CL_m_0 = CL_system(N * SR_system.spin, m0=0, N=N)

ntu = 30
# SR_m_max.creat_time_evolution(ntu * SR_m_max.natural_time_unit)
# SR_m_0.creat_time_evolution(ntu * SR_m_0.natural_time_unit)
# CL_m_max.creat_time_evolution(20 * ntu * SR_m_0.natural_time_unit)
CL_m_0.creat_time_evolution(20 * ntu * CL_m_0.natural_time_unit)