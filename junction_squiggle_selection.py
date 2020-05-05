'''
Functions to select junction squiggles without the requirement of transcript reference

Input:
    BED line
    Cigar string
    fast5 filename

output:
    Normalised junction squiggles
'''

import h5py
import re

def tombo_squiggle_to_basecalls(read_fast5_fn):
    '''
    script from tombo (https://gist.github.com/marcus1487/5cf083644143aafc770dceb213783346)
    Input:
        read_fast5: fast5 filename
    Returns:
        rsqgl_results <class>: Resquiggle result
           # useful rsqgl_results attributes:
                #   - genome_seq: clipped basecalls to model-able positions (can't get k-mers and thus expected levels for first and last few bases)
                #   - raw_signal: median normalized raw signal values (only for portion of raw signal mapped to model-able basecalls)
                #   - read_start_rel_to_raw: start of rsqgl_results.raw_signal and rsqgl_results.segs within original all_raw_signal
                #   - segs: start position for each model-able base (should be one longer than length of genome_seq to include last end point) 
    '''

    from tombo import tombo_helper, tombo_stats, resquiggle
    # extract read info
    fast5_data = h5py.File(read_fast5_fn)
    seq_samp_type = tombo_helper.get_seq_sample_type(fast5_data)
    seq_data = resquiggle.get_read_seq(fast5_data, 'Basecall_1D_000', 'BaseCalled_template', seq_samp_type, 0)

    # prep tombo objects
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    start_clip = std_ref.central_pos
    end_clip = std_ref.kmer_width - std_ref.central_pos - 1

    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)

    # extract raw signal
    all_raw_signal = tombo_helper.get_raw_read_slot(fast5_data)['Signal'][:]
    # if this is a direct RNA read, flip raw signal to process from 5' to 3'
    if seq_samp_type.rev_sig:
        all_raw_signal = all_raw_signal[::-1]

    # spoof mapping read to its own basecalls (with pseudo-clipped bases in genome_location for model-able postions)
    map_results = tombo_helper.resquiggleResults(
        align_info=tombo_helper.alignInfo(ID=seq_data.id, Subgroup='BaseCalled_template', ClipStart=0, ClipEnd=0, Insertions=0,
                                        Deletions=0, Matches=len(seq_data.seq), Mismatches=0),
        genome_loc=tombo_helper.genomeLocation(Start=std_ref.central_pos, Strand='+', Chrom='NA'),
        genome_seq=seq_data.seq, raw_signal=all_raw_signal, mean_q_score=seq_data.mean_q_score)

    # align raw signal to basecalls
    try:
        rsqgl_results = resquiggle.resquiggle_read(
        map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal)
    except:
        rsqgl_params = rsqgl_params ._replace(bandwidth=2000)
        rsqgl_results = resquiggle.resquiggle_read(
        map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal)

    return rsqgl_results, start_clip, end_clip


def genome_to_read_pos_conversion(cigar):
    '''
    Input:
        g_pos: 0-base position within genome
        cigar: CIGAR string from BAM/SAM file
    
    Returns:
        r_pos: 0-base position within a read
    '''

    cigar_long = []
    for count, type in re.findall('(\d+)([A-Za-z])', cigar):
        cigar_long += int(count) * [type]

    r_index = -1
    g_r_mapping = []
    for i in cigar_long:
        if i in "MIS":
            r_index += 1
        if i in "MD":
            g_r_mapping.append(r_index)
        if i == "N":
            g_r_mapping.append(-1)

    return g_r_mapping

