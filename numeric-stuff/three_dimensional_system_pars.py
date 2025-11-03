# a complete (ish) store of all parameter values used to develop the plots
# noting that the calculation of particular points, fold lines and separatrix
# values mean that some parts of the codebase may need to be changed for more accurate results
kv = 100

three_spikes = {
    'I': 82.35/(20*kv),
    'epsilon': 0.00015,
    'N': 50/kv,
    'a': kv/18,
    'm': -1.2/18,
    'K': -100/kv,
    'a_k': 20/20,
    'sub': 8/20,
    'S': -100/kv,
    'a_l': 2.0/20,
    'L': -70/kv,
    'gamma_y': kv/10,
    'beta_y': -10/10,
    'gamma_z': kv/15,
    'beta_z': -21/15,
}

two_spikes = {
    'I': 80/(20*kv),
    'epsilon': 0.00015,
    'N': 50/kv,
    'a': kv/18,
    'm': -1.2/18,
    'K': -100/kv,
    'a_k': 20/20,
    'sub': 8/20,
    'S': -100/kv,
    'a_l': 2.0/20,
    'L': -70/kv,
    'gamma_y': kv/10,
    'beta_y': -10/10,
    'gamma_z': kv/15,
    'beta_z': -21/15,
}

one_spike_high_sub = {
    'I': 70/(20*kv),
    'epsilon': 0.00015,
    'N': 50/kv,
    'a': kv/18,
    'm': -1.2/18,
    'K': -100/kv,
    'a_k': 20/20,
    'sub': 8/20,
    'S': -100/kv,
    'a_l': 2.0/20,
    'L': -70/kv,
    'gamma_y': kv/10,
    'beta_y': -10/10,
    'gamma_z': kv/15,
    'beta_z': -21/15,
}

type_two = {
    'I': 0/(20*kv),
    'epsilon': 0.00015,
    'N': 50/kv,
    'a': kv/18,
    'm': -1.2/18,
    'K': -100/kv,
    'a_k': 20/20,
    'sub': 2/20,
    'S': -100/kv,
    'a_l': 2.0/20,
    'L': -70/kv,
    'gamma_y': kv/10,
    'beta_y': -10/10,
    'gamma_z': kv/15,
    'beta_z': -21/15,
}

type_two_one_spike = {
    'I': 45/(20*kv),
    'epsilon': 0.00015,
    'N': 50/kv,
    'a': kv/18,
    'm': -1.2/18,
    'K': -100/kv,
    'a_k': 20/20,
    'sub': 2/20,
    'S': -100/kv,
    'a_l': 2.0/20,
    'L': -70/kv,
    'gamma_y': kv/10,
    'beta_y': -10/10,
    'gamma_z': kv/15,
    'beta_z': -21/15,
}

type_one = {
    'I': 0/(20*kv),
    'epsilon': 0.00015,
    'N': 50/kv,
    'a': kv/18,
    'm': -1.2/18,
    'K': -100/kv,
    'a_k': 20/20,
    'sub': 3/20,
    'S': 50/kv,
    'a_l': 2.0/20,
    'L': -70/kv,
    'gamma_y': kv/10,
    'beta_y': -10/10,
    'gamma_z': kv/15,
    'beta_z': -21/15,
}