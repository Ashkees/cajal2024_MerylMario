import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()


import spikeinterface.full as si
import numpy as np
from pathlib import Path
import os
import probeinterface as pi
from probeinterface.plotting import plot_probe
import probeinterface.plotting as piplot
import probe
from spikeinterface.postprocessing import compute_spike_amplitudes, compute_principal_components
from spikeinterface.exporters import export_to_phy

def load_sort_save(raw_file, probe_file, force=True):
    # set some params
    base_folder = raw_file.parent
    fs = 20000
    n_cpus = os.cpu_count()
    n_jobs = n_cpus - 4
    job_kwargs = dict(n_jobs=n_jobs, chunk_duration="1s", progress_bar=True)
    # load from raw or saved
    if (base_folder / "preprocessed").is_dir() and not force:
        recording = si.load_extractor(base_folder / "preprocessed")
    else:
        full_raw_rec = si.read_binary(raw_file, fs, 'int16', num_channels=96, gain_to_uV=0.195)
        probes = probe.load_probe(probe_file, force=False)
        raw_rec = full_raw_rec.set_probegroup(probes)
        raw_rec.annotate(is_filtered=True)
        recording = raw_rec.save(folder=base_folder / "preprocessed", overwrite=True, **job_kwargs)
    # sort using spykingcircus2 and save
    sorting_SC2 = si.run_sorter_by_property('spykingcircus2', recording, grouping_property='group',
                                            working_folder=base_folder / 'results_SC2',
                                            verbose=True, job_kwargs=job_kwargs)
    sorting_SC2.save(folder=base_folder / "agg_SC2")
    # extract waveforms and export to Phy
    wav_ex = si.extract_waveforms(recording, sorting_SC2, folder=base_folder / "waveforms_sparse",
                                  sparse=True, overwrite=True, **job_kwargs)
    compute_spike_amplitudes(waveform_extractor=wav_ex)
    si.export_to_phy(wav_ex, output_folder=base_folder / 'phy_SC2', remove_if_exists=True,
                     compute_amplitudes=True, compute_pc_features=True, copy_binary=False,
                     **job_kwargs)
    return sorting_SC2, wav_ex





go_on = False

if __name__ == '__main__':

    raw_f = Path(r'E:\Mulle\Ashley_data\Clustering\Anael\Anael_analysis\CM_200902_151604_Anael\amplifier.dat')
    base_folder = raw_f.parent
    probe_f = Path(r'E:\Mulle\Ashley_data\Clustering\probe_files\SI_4shank_linear.json')
    sorting_SC2, we = load_sort_save(raw_f, probe_f)

if go_on:
    probes = probe.load_probe(pf, force=True)
    probe_df = probes.to_dataframe(complete=True)
    # plot the probe layout
    #piplot.plot_probe_group(probes, same_axes=False, with_device_index=True)
    #plot_probe(probes[0], with_device_index=True)
    #plot_probe(probes[1], with_device_index=True)


    raw_rec = full_raw_rec.set_probegroup(probes)

    fig, ax = plt.subplots()
    si.plot_probe_map(raw_rec, with_channel_ids=True, ax=ax)

    # could filter and apply CMR, but Anael's data are already filtered
    recording_f = si.bandpass_filter(raw_rec, freq_min=300, freq_max=9000)
    recording_cmr = si.common_reference(recording_f, reference='global', operator='median')

    w = si.plot_traces(lin_rec, channel_ids=np.arange(64, 96), order_channel_by_depth=True,
                       show_channel_ids=True, time_range=[10, 13])

    # take a section of 5 min for example
    raw_rec_sub = raw_rec.frame_slice(start_frame=0 * fs, end_frame=300 * fs)

    # separate out the two probes


    # to save a si object
    n_cpus = os.cpu_count()
    n_jobs = n_cpus - 4
    job_kwargs = dict(n_jobs=n_jobs, chunk_duration="1s", progress_bar=True)
    if (base_folder / "preprocessed").is_dir():
        recording_saved = si.load_extractor(base_folder / "preprocessed")
    else:
        recording_saved = raw_rec_sub.save(folder=base_folder / "preprocessed", **job_kwargs)
    # to load a si object
    recording_loaded = si.load_extractor(base_folder / "preprocessed")
    # add an annotation to say that recording is already filtered, since it was done before by spykingcircus
    recording_loaded.annotate(is_filtered=True)

    c_shank = '1'
    channel_ids = probe_df['device_channel_indices'][probe_df['shank_ids'] == c_shank].values
    w = si.plot_traces(recording_loaded, channel_ids=channel_ids, order_channel_by_depth=True,
                       show_channel_ids=True, time_range=[10, 13])


    # sorting_SC2 = si.run_sorter('spykingcircus2', recording_saved,
    #                             output_folder=base_folder / 'results_SC2',
    #                             verbose=True, job_kwargs=job_kwargs)

    # to do the sorting
    sorting_SC2 = si.run_sorter_by_property('spykingcircus2', recording_loaded, grouping_property='group',
                                working_folder=base_folder / 'aggregated_sorting',
                                verbose=True, job_kwargs=job_kwargs)
    sorting_SC2 = sorting_SC2.save(folder=base_folder /"agg_sorting")
    # to load from saved
    sorting_SC2_loaded = si.read_sorter_folder(base_folder / "results_SC2" / '0')  # multishank
    sorting_SC2_loaded_1 = si.read_sorter_folder(base_folder / "results_SC2" / '1')  # linear
    sorting_SC2_loaded_aggregate = si.aggregate_units([sorting_SC2_loaded, sorting_SC2_loaded_1])

    sorting_SC2_loaded_agg = si.load_extractor(base_folder / "agg_sorting")

    we = si.extract_waveforms(recording_loaded, sorting_SC2_loaded_agg, folder=base_folder / "waveforms_sparse",
                              sparse=True, overwrite=True, **job_kwargs)
    from spikeinterface.postprocessing import compute_spike_amplitudes, compute_principal_components
    from spikeinterface.exporters import export_to_phy
    # supposed to compute amplitudes and PCs before extraction, but couldn't get it to extract with the phy stuff
    compute_spike_amplitudes(waveform_extractor=we)
    compute_principal_components(waveform_extractor=we, n_components=3, mode='by_channel_global')

    si.export_to_phy(we, output_folder=base_folder / 'phy_SC2', remove_if_exists=True,
                     compute_amplitudes=True, compute_pc_features=True, copy_binary=False,
                     **job_kwargs)

    sorting_phy_curated = si.read_phy(base_folder / 'phy_SC2/', exclude_cluster_groups=['noise', 'mua'])
    sorting_phy_Anael = si.read_phy(base_folder / 'amplifier/amplifier.GUI/', exclude_cluster_groups=['noise', 'mua'])
    sorting_phy_Cath = si.read_phy(Path(r'E:\Mulle\Ashley_data\Clustering\Anael\Catherine_analysis\CM_200902_151604_Catherine\amplifier\amplifier.GUI'), exclude_cluster_groups=['noise', 'mua'])


    w_rs = si.plot_rasters(sorting_phy_curated, time_range=(0, 60), backend="matplotlib")

    comp_Ash_Cath = si.compare_two_sorters(sorting_phy_curated, sorting_phy_Cath, 'Ash', 'Cath')
    si.plot_agreement_matrix(comp_Ash_Cath, unit_ticks=False)







