from brian2 import *

def visualise(S):
    Ns = len(S.source)
    Nt = len(S.target)
    figure(figsize=(10, 4), dpi= 80, facecolor='w', edgecolor='k')
    subplot(121)
    plot(zeros(Ns), arange(Ns), 'ok', ms=10)
    plot(ones(Nt), arange(Nt), 'ok', ms=10)
    for i, j in zip(S.i, S.j):
        plot([0, 1], [i, j], '-k')
    xticks([0, 1], ['Source', 'Target'])
    ylabel('Neuron index')
    xlim(-0.1, 1.1)
    ylim(-1, max(Ns, Nt))
    subplot(122)
    plot(S.i, S.j, 'ok')
    xlim(-1, Ns)
    ylim(-1, Nt)
    xlabel('Source neuron index')
    ylabel('Target neuron index')

def rand_params(Parameter,Unit,N_Cells,Step):
    Nn = [int(N_Cells/2), N_Cells-int(N_Cells/2)] 
    shuffle(Nn)
    Base = int(1/Step)
    Start = int(Base*Parameter)
    Begin = Start - Nn[0]
    End = Start + Nn[1]
    Param_vector = [x / float(Base) for x in range(Begin, End, 1)]*Unit
    shuffle(Param_vector)
    return Param_vector

def IO_equations():
  eqs_IO_V = '''
  dVs/dt = (-(I_ds + I_ls + I_Na + I_Ca_l + I_K_dr + I_h + I_as) + Iapp_s)/Cm : volt
  dVd/dt = (-(I_sd + I_ld + I_Ca_h + I_K_Ca + I_c) + Iapp_d)/Cm : volt
  dVa/dt = (-(I_K_a + I_sa + I_la + I_Na_a))/Cm : volt
  I_c : metre**-2*amp
  Iapp_s : metre**-2*amp
  Iapp_d : metre**-2*amp
  '''
  eqs_IO_Ca = '''
  dCa/dt = (-3*I_Ca_h*((uamp / cm**2)**-1)*mM - 0.075*Ca)/ms : mM
  '''
  eqs_IO_Isom = '''
  I_as    = (g_int/(1-p2))*(Vs-Va)     : metre**-2*amp
  I_ls    = g_ls*(Vs-V_l)              : metre**-2*amp
  I_ds    = (g_int/p)*(Vs-Vd)          : metre**-2*amp
  I_Na    = g_Na*m_inf**3*h*(Vs-V_Na)  : metre**-2*amp
  I_Ca_l  = g_Ca_l*k*k*k*l*(Vs-V_Ca)   : metre**-2*amp
  I_K_dr  = g_Kdr*n*n*n*n*(Vs-V_K)     : metre**-2*amp
  I_h     = g_h*q*(Vs-V_h)             : metre**-2*amp
  I_K_s   = g_K_s*(x_s**4)*(Vs-V_K)    : metre**-2*amp
  '''
  eqs_IO_Iden = '''
  I_sd    = (g_int/(1-p))*(Vd-Vs)      : metre**-2*amp
  I_ld    = g_ld*(Vd-V_l)              : metre**-2*amp
  I_Ca_h  = g_Ca_h*r*r*(Vd-V_Ca)       : metre**-2*amp
  I_K_Ca  = g_K_Ca*s*(Vd-V_K)          : metre**-2*amp
  '''
  eqs_IO_Iax = '''
  I_K_a  = g_K_a *x_a**4*(Va-V_K)      : metre**-2*amp
  I_sa   = (g_int/p2)*(Va-Vs)          : metre**-2*amp
  I_la   = g_la*(Va-V_l)               : metre**-2*amp
  I_Na_a = g_Na_a*m_a**3*h_a*(Va-V_Na) : metre**-2*amp
  '''
  eqs_IO_activation = '''
  dh/dt = (h_inf - h)/tau_h : 1
  dk/dt = (k_inf - k)/tau_k : 1
  dl/dt = (l_inf - l)/tau_l : 1
  dn/dt = (n_inf - n)/tau_n : 1
  dq/dt = (q_inf - q)/tau_q : 1
  dr/dt = (r_inf - r)/tau_r : 1
  ds/dt = (s_inf - s)/tau_s : 1
  m_a = m_inf_a : 1
  dh_a/dt = (h_inf_a - h_a)/tau_h_a : 1
  dx_a/dt = (x_inf_a - x_a)/tau_x_a : 1
  dx_s/dt = (x_inf_s - x_s)/tau_x_s : 1
  '''
  eqs_IO_inf = '''
  m_inf   = alpha_m /(alpha_m+beta_m)        : 1
  h_inf   = alpha_h/(alpha_h+beta_h)         : 1
  k_inf   = 1/(1+e**(-(Vs/mvolt+61)/4.2))    : 1
  l_inf   = 1/(1+e**((Vs/mvolt+85.5)/8.5))   : 1
  n_inf   = alpha_n/(alpha_n+beta_n)         : 1
  q_inf   = 1/(1+e**((Vs/mvolt+75)/(5.5)))   : 1
  r_inf   = alpha_r/(alpha_r + beta_r)       : 1
  s_inf   = alpha_s/(alpha_s+beta_s)         : 1
  m_inf_a = 1/(1+(e**((-30-Va/mvolt)/ 5.5))) : 1
  h_inf_a = 1/(1+(e**((-60-Va/mvolt)/-5.8))) : 1
  x_inf_a = alpha_x_a/(alpha_x_a+beta_x_a)   : 1
  x_inf_s = alpha_x_s/(alpha_x_s + beta_x_s) : 1
  '''
  eqs_IO_tau = '''
  tau_h   = 170*msecond/(alpha_h+beta_h)                                          : second
  tau_k   = 5*msecond                                                             : second
  tau_l   = 1*msecond*(35+(20*e**((Vs/mvolt+160)/30/(1+e**((Vs/mvolt+84)/7.3))))) : second
  tau_n   = 5*msecond/(alpha_n+beta_n)                                            : second
  tau_q   = 1*msecond/(e**((-0.086*Vs/mvolt-14.6))+e**((0.07*Vs/mvolt-1.87)))     : second
  tau_r   = 5*msecond/(alpha_r + beta_r)                                          : second
  tau_s   = 1*msecond/(alpha_s + beta_s)                                          : second
  tau_h_a = 1.5*msecond*e**((-40-Va/mvolt)/33)                                    : second
  tau_x_a = 1*msecond/(alpha_x_a + beta_x_a)                                      : second
  tau_x_s = 1*msecond/(alpha_x_s + beta_x_s)                                      : second
  '''
  eqs_IO_alpha = '''
  alpha_m   = (0.1*(Vs/mvolt + 41))/(1-e**(-(Vs/mvolt+41)/10)) : 1
  alpha_h   = 5.0*e**(-(Vs/mvolt+60)/15) : 1
  alpha_n   = (Vs/mvolt + 41)/(1-e**(-(Vs/mvolt+41)/10)) : 1
  alpha_r   = 1.7/(1+e**(-(Vd/mvolt - 5)/13.9)) : 1
  alpha_s   = ((0.00002*Ca/mM)*int((0.00002*Ca/mM)<0.01) + 0.01*int((0.00002*Ca/mM)>=0.01)) : 1
  alpha_x_a = 0.13*(Va/mvolt + 25)/(1-e**(-(Va/mvolt+25)/10)) : 1
  alpha_x_s = 0.13*(Vs/mvolt + 25)/(1-e**(-(Vs/mvolt+25)/10)) : 1
  '''

  eqs_IO_beta = '''
  beta_m = 9.0*e**(-(Vs/mvolt+60)/20)                        : 1
  beta_h = (Vs/mvolt+50)/(1-e**(-(Vs/mvolt+50)/10))          : 1
  beta_n = 12.5*e**(-(Vs/mvolt+51)/80)                       : 1
  beta_r = 0.02*(Vd/mvolt + 8.5)/(e**((Vd/mvolt + 8.5)/5)-1) : 1
  beta_s = 0.015                                             : 1
  beta_x_a  = 1.69*e**(-0.0125*(Va/mvolt + 35))              : 1
  beta_x_s  = 1.69*e**(-0.0125*(Vs/mvolt+ 35))               : 1
  '''

  eqs_vector = '''
  V_Na : volt
  V_K  : volt
  V_Ca : volt
  V_l  : volt
  V_h  : volt
  Cm : farad*meter**-2
  g_Na   : siemens/meter**2
  g_Kdr  : siemens/meter**2
  g_Ca_l : siemens/meter**2
  g_h    : siemens/meter**2
  g_Ca_h : siemens/meter**2
  g_K_Ca : siemens/meter**2
  g_ls : siemens/meter**2
  g_ld : siemens/meter**2
  g_int  : siemens/meter**2
  g_Na_a   : siemens/meter**2
  g_K_a   : siemens/meter**2
  g_la   : siemens/meter**2
  g_K_s   : siemens/meter**2
  p : 1
  p2 : 1
  '''

  eqs_IO = eqs_IO_beta
  eqs_IO += eqs_IO_alpha
  eqs_IO += eqs_IO_tau
  eqs_IO += eqs_IO_inf
  eqs_IO += eqs_IO_activation
  eqs_IO += eqs_IO_Iax
  eqs_IO += eqs_IO_Iden
  eqs_IO += eqs_IO_Isom
  eqs_IO += eqs_IO_Ca
  eqs_IO += eqs_IO_V
  eqs_IO += eqs_vector
  return eqs_IO

def cell_values(N_Cells_IO, IO_response):
  IO_V_Na = rand_params(55,mvolt ,N_Cells_IO,(1.0/N_Cells_IO))  #55*mvolt
  IO_V_K = rand_params(-75,mvolt ,N_Cells_IO,(1.0/N_Cells_IO))  #-75*mvolt
  IO_V_Ca = rand_params(120,mvolt ,N_Cells_IO,(1.0/N_Cells_IO))  #120*mvolt
  IO_V_l = rand_params(10,mvolt ,N_Cells_IO,(1.0/N_Cells_IO))  #10*mvolt 
  IO_V_h = rand_params(-43,mvolt ,N_Cells_IO,(1.0/N_Cells_IO))  #-43*mvolt 
  IO_Cm = rand_params(1,uF*cm**-2 ,N_Cells_IO,(0.1/N_Cells_IO))  #1*uF*cm**-2 
  IO_g_Na = rand_params(150,mS/cm**2,N_Cells_IO,(1.0/N_Cells_IO))  #150*mS/cm**2
  IO_g_Kdr = rand_params(9.0,mS/cm**2,N_Cells_IO,(0.1/N_Cells_IO))  #9.0*mS/cm**2
  IO_g_K_s = rand_params(5.0,mS/cm**2,N_Cells_IO,(0.1/N_Cells_IO))  #5.0*mS/cm**2
  IO_g_h = rand_params(0.12,mS/cm**2,N_Cells_IO,(0.01/N_Cells_IO))  #0.12*mS/cm**2
  IO_g_Ca_h = rand_params(4.5,mS/cm**2,N_Cells_IO,(0.1/N_Cells_IO))  #4.5*mS/cm**2
  IO_g_K_Ca = rand_params(35,mS/cm**2,N_Cells_IO,(0.5/N_Cells_IO))  #35*mS/cm**2
  IO_g_Na_a = rand_params(240,mS/cm**2,N_Cells_IO,(1.0/N_Cells_IO))  #240*mS/cm**2
  IO_g_K_a = rand_params(20,mS/cm**2,N_Cells_IO,(0.5/N_Cells_IO))  #20*mS/cm**2
  IO_g_ls = rand_params(0.016,mS/cm**2,N_Cells_IO,(0.001/N_Cells_IO))  #0.016*mS/cm**2
  IO_g_ld = rand_params(0.016,mS/cm**2,N_Cells_IO,(0.001/N_Cells_IO))  #0.016*mS/cm**2
  IO_g_la = rand_params(0.016,mS/cm**2,N_Cells_IO,(0.001/N_Cells_IO))  #0.016*mS/cm**2
  IO_g_int = rand_params(0.13,mS/cm**2,N_Cells_IO,(0.001/N_Cells_IO))  #0.13*mS/cm**2
  IO_p = rand_params(0.25,1,N_Cells_IO,(0.01/N_Cells_IO))  #0.25
  IO_p2 = rand_params(0.15,1,N_Cells_IO,(0.01/N_Cells_IO))   #0.15
  if IO_response=='oscillatory':    
      IO_g_Ca_l =  [.5*mS/cm**2]*N_Cells_IO
  elif IO_response=='non':    
      IO_g_Ca_l =  [.1*mS/cm**2]*N_Cells_IO
  elif IO_response=='spiking':    
      IO_g_Ca_l =  [1.1*mS/cm**2]*N_Cells_IO
  elif IO_response=='both':    
      IO_g_Ca_l =  [.75*mS/cm**2]*N_Cells_IO

  return IO_V_Na, IO_V_K, IO_V_Ca, IO_V_l, IO_V_h, IO_Cm, IO_g_Na, IO_g_Kdr, IO_g_K_s, IO_g_h, IO_g_Ca_h, IO_g_K_Ca, IO_g_Na_a, IO_g_K_a, IO_g_ls, IO_g_ld, IO_g_la, IO_g_int, IO_p, IO_p2, IO_g_Ca_l

def neurons(N_Cells_IO,IO_response,dt,dt_rec):
  eqs_IO = IO_equations()
  IO_V_Na, IO_V_K, IO_V_Ca, IO_V_l, IO_V_h, IO_Cm, IO_g_Na, IO_g_Kdr, IO_g_K_s, IO_g_h, IO_g_Ca_h, IO_g_K_Ca, IO_g_Na_a, IO_g_K_a, IO_g_ls, IO_g_ld, IO_g_la, IO_g_int, IO_p, IO_p2, IO_g_Ca_l = cell_values(N_Cells_IO,IO_response)

  ####### Coupled group
  IO_Coupled = NeuronGroup(N_Cells_IO, model = eqs_IO, threshold='Vs>20*mV' , method = 'euler',name = 'SchweighoferOlive_Coupled',dt=dt)
  IO_Statemon_Coupled = StateMonitor(IO_Coupled, variables = ['Vs','Vd','I_c'], record = True, dt=dt_rec)
  IO_Spikemon_Coupled = SpikeMonitor(IO_Coupled)
  IO_rate_Coupled = PopulationRateMonitor(IO_Coupled)

  ####### Uncoupled group
  IO_Uncoupled = NeuronGroup(N_Cells_IO, model = eqs_IO, threshold='Vs>20*mV' , method = 'euler',name = 'SchweighoferOlive',dt=dt)
  IO_Statemon_Uncoupled = StateMonitor(IO_Uncoupled, variables = ['Vs','Vd','I_c'], record = True, dt=dt_rec)
  IO_Spikemon_Uncoupled = SpikeMonitor(IO_Uncoupled)
  IO_rate_Uncoupled = PopulationRateMonitor(IO_Uncoupled)

  for ii in range(0, N_Cells_IO, 1):
    IO_Coupled.V_Na[ii] = IO_V_Na[ii]
    IO_Coupled.V_K[ii] = IO_V_K[ii]
    IO_Coupled.V_Ca[ii] = IO_V_Ca[ii]
    IO_Coupled.V_l[ii] = IO_V_l[ii]
    IO_Coupled.V_h[ii] = IO_V_h[ii]
    IO_Coupled.Cm[ii] = IO_Cm [ii]
    IO_Coupled.g_Na[ii] = IO_g_Na[ii]
    IO_Coupled.g_Kdr[ii] = IO_g_Kdr[ii]
    IO_Coupled.g_K_s[ii] = IO_g_K_s[ii]
    IO_Coupled.g_h[ii] = IO_g_h[ii]
    IO_Coupled.g_Ca_h[ii] = IO_g_Ca_h[ii]
    IO_Coupled.g_K_Ca[ii] = IO_g_K_Ca[ii]
    IO_Coupled.g_Na_a[ii] = IO_g_Na_a[ii]
    IO_Coupled.g_K_a[ii] = IO_g_K_a[ii]
    IO_Coupled.g_ls[ii] = IO_g_ls[ii]
    IO_Coupled.g_ld[ii] = IO_g_ld[ii]
    IO_Coupled.g_la[ii] = IO_g_la[ii]
    IO_Coupled.g_int[ii] = IO_g_int[ii]
    IO_Coupled.p[ii] = IO_p[ii]
    IO_Coupled.p2[ii] = IO_p2[ii]
    IO_Coupled.g_Ca_l[ii] =  IO_g_Ca_l[ii]
    IO_Uncoupled.V_Na[ii] = IO_V_Na[ii]
    IO_Uncoupled.V_K[ii] = IO_V_K[ii]
    IO_Uncoupled.V_Ca[ii] = IO_V_Ca[ii]
    IO_Uncoupled.V_l[ii] = IO_V_l[ii]
    IO_Uncoupled.V_h[ii] = IO_V_h[ii]
    IO_Uncoupled.Cm[ii] = IO_Cm [ii]
    IO_Uncoupled.g_Na[ii] = IO_g_Na[ii]
    IO_Uncoupled.g_Kdr[ii] = IO_g_Kdr[ii]
    IO_Uncoupled.g_K_s[ii] = IO_g_K_s[ii]
    IO_Uncoupled.g_h[ii] = IO_g_h[ii]
    IO_Uncoupled.g_Ca_h[ii] = IO_g_Ca_h[ii]
    IO_Uncoupled.g_K_Ca[ii] = IO_g_K_Ca[ii]
    IO_Uncoupled.g_Na_a[ii] = IO_g_Na_a[ii]
    IO_Uncoupled.g_K_a[ii] = IO_g_K_a[ii]
    IO_Uncoupled.g_ls[ii] = IO_g_ls[ii]
    IO_Uncoupled.g_ld[ii] = IO_g_ld[ii]
    IO_Uncoupled.g_la[ii] = IO_g_la[ii]
    IO_Uncoupled.g_int[ii] = IO_g_int[ii]
    IO_Uncoupled.p[ii] = IO_p[ii]
    IO_Uncoupled.p2[ii] = IO_p2[ii]
    IO_Uncoupled.g_Ca_l[ii] =  IO_g_Ca_l[ii]    
  
  return IO_Coupled, IO_Statemon_Coupled, IO_Uncoupled, IO_Statemon_Uncoupled

def syn(g_c_coupled,g_c_uncoupled,IO_Coupled,IO_Uncoupled):
  eqs_IO_syn_Coupled = ''' I_c_pre = (g_c_coupled)*(0.6*e**(-((Vd_pre/mvolt-Vd_post/mvolt)/50)**2) + 0.4)*(Vd_pre-Vd_post) : metre**-2*amp (summed)'''
  IO_synapse_Coupled = Synapses(IO_Coupled, IO_Coupled, eqs_IO_syn_Coupled, name = 'IO_Synapse_Coupled')

  eqs_IO_syn_Uncoupled = '''I_c_pre = (g_c_uncoupled)*(0.6*e**(-((Vd_pre/mvolt-Vd_post/mvolt)/50)**2) + 0.4)*(Vd_pre-Vd_post) : metre**-2*amp (summed)'''
  IO_synapse_Uncoupled = Synapses(IO_Uncoupled, IO_Uncoupled, eqs_IO_syn_Uncoupled, name = 'IO_Synapse_Uncoupled')
  
  return IO_synapse_Coupled, IO_synapse_Uncoupled

def plot_neurons(N_Cells_IO,IO_Statemon_Coupled,IO_Statemon_Uncoupled):
  plt.figure(figsize=(14, 10), dpi= 80, facecolor='w', edgecolor='k')
  plt.subplot(2, 1, 1)
  for ii in range(0,N_Cells_IO):
      plot(IO_Statemon_Coupled.t/ms,IO_Statemon_Coupled.Vs[ii]/mV, label='IO_Cell'+str(ii+1))
  title('Membrane Potential Coupled')
  ylabel('V [mV]')
  plt.subplot(2, 1, 2)
  for ii in range(0,N_Cells_IO):
      plot(IO_Statemon_Uncoupled.t/ms,IO_Statemon_Uncoupled.Vs[ii]/mV, label='IO_Cell'+str(ii+1))
  title('Membrane Potential Uncoupled')
  xlabel('Time [ms]')
  ylabel('V [mV]')
  show()

def plot_coupling_connectivity(N_Cells_IO,IO_Statemon_Coupled,IO_Statemon_Uncoupled):
  plt.figure(figsize=(14, 10), dpi= 80, facecolor='w', edgecolor='k')
  plt.subplot(2, 1, 1)
  for ii in range(0,N_Cells_IO):
      plot(IO_Statemon_Coupled.t/ms,IO_Statemon_Coupled.I_c[ii]/(mS/meter**2), label='I_c_IO_Cell'+str(ii+1))
  title('Coupling Connectivity Coupled')
  ylabel('I_c [mS/m**2]')
  plt.subplot(2, 1, 2)
  for ii in range(0,N_Cells_IO):
      plot(IO_Statemon_Uncoupled.t/ms,IO_Statemon_Uncoupled.I_c[ii]/(mS/meter**2), label='I_c_IO_Cell'+str(ii+1))
  title('Coupling Connectivity Uncoupled')
  xlabel('Time [ms]')
  ylabel('I_c [mS/m**2]')
  show()