import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

import numpy as np
from pathlib import Path
import os
import probeinterface as pi
from probeinterface.plotting import plot_probe
import probeinterface.plotting as piplot


def load_probe(probe_file: Path, force: bool = False):
    '''
    Load a probe file

    Parameters
    ----------
    probe_file : Path
        Path to the probe .json file.
    force : bool
        When True, will recreate a rewrite over the given probe_file
        When False, will load the probe_file if it exists, and create it if it doesn't

    Returns
    -------
    probes : si.ProbeGroup object
        Gives channel geometry and device_ids

    '''
    if probe_file.exists() and not force:
        probes = pi.read_probeinterface(probe_file)
        return probes
    # if probe.json doesn't exist yet, generate probes from scratch
    # linear probe neuronexus A1x32-6mm-50-177
    linear = pi.generate_linear_probe(num_elec=32, ypitch=50, contact_shapes='circle', contact_shape_params={'radius': 7.5})
    # set the device channel indices, based on wiring diagram
    l_ch_ind = np.array([92, 70, 84, 72, 85, 71, 93, 74, 86, 73, 83, 75, 87, 67, 82, 76,
                         88, 66, 81, 78, 89, 77, 80, 79, 90, 69, 91, 68, 94, 65, 95, 64])
    linear.set_device_channel_indices(l_ch_ind)

    # multishank probe neuronexus A4x16-Poly2-5mm-23s-200-177
    model_shank = pi.generate_multi_columns_probe(num_columns=2, num_contact_per_column=[8, 8],
                                            xpitch=30, ypitch=46, y_shift_per_column=[23, 0],
                                            contact_shapes='circle', contact_shape_params={'radius': 7.5})
    ch_ind = [[41, 34, 32, 45, 53, 47, 33, 43, 51, 42, 35, 39, 49, 38, 36, 37],
              [56, 62, 58, 54, 40, 48, 46, 55, 44, 50, 59, 57, 61, 52, 63, 60],
              [0, 6, 8, 4, 14, 22, 9, 16, 12, 18, 7, 5, 10, 3, 2, 1],
              [28, 23, 19, 30, 17, 11, 21, 31, 20, 13, 25, 29, 24, 15, 27, 26]]
    # rearrange channels slightly:
    for i, ch in enumerate(ch_ind):
        ch_ind[i] = ch[1::2] + ch[0::2]
    # flatten
    dev_ind = [i for j in ch_ind for i in j]

    shanks = []
    # add rest of shanks
    for s in np.arange(4):
        c_shank = model_shank.copy()
        c_shank.move([200*s, 0])
        shanks.append(c_shank)

    multi_shank = pi.combine_probes(shanks)
    multi_shank.set_device_channel_indices(dev_ind)
    multi_shank.move([385, 0])

    probes = pi.ProbeGroup()
    # need to be added in order, according to which port on the intan box they're plugged into
    probes.add_probe(multi_shank)
    probes.add_probe(linear)

    pi.write_probeinterface(probe_file, probes)

    return probes
    #return [multi_shank, linear]


