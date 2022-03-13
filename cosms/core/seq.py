#!/usr/bin/env python3
# coding: utf-8
"""seq
sequencing of COS based on given fragment ion peaks
"""
import copy
from collections import defaultdict as ddict
# verbose = True
verbose = False

def determine_sequences(precursor = None,
                        fragments = None,
                        prec_seq = None,
                        ignored_frag = None,
                        ignored_seq = None):
    """ Sequencing of COS based on fragment ions
    """
    if verbose:
        print('\n-- Begin of function --', precursor)
        print('precursor.nonred_maxonly', precursor.nonred_maxonly)
    ignore_ion_if_conflict_threshold = 0.05
    ignore_ion_threshold = 0.01
    ignore_unresolved_seqs_threshold = 0.06 # 0.06 is standard
    ignore_single_seq_threshold = 0.00
    prec_seq_copy = copy.deepcopy(prec_seq) # work on copy
    dp_minus1 = precursor.dp - 1
    if precursor.nonred_maxonly:
        return_value = False
        for nonred_type in precursor.nonredend_types:
            if (nonred_type, dp_minus1) in fragments:
            # do not sequence if no largest b-ion present TODO what happens then?
                return_value = True
                break
        if not return_value:
            return False
    if ignored_frag == None: # fragment ions that are not considered
        ignored_frag = set()
    if ignored_seq == None: # low abundant sequences that are considered as absent
        ignored_seq = set()
    
    allowed_ion_types = set() # considered types of fragment ions
    bc_seq = ddict(set) # sequences based on b- or c-ions
    if precursor.nonred_maxonly:
        for nonred_type in precursor.nonredend_types:
            allowed_ion_types.add((nonred_type, dp_minus1))
            for f_key, (intensity, rel_intensity) in fragments[nonred_type, 
                                                               dp_minus1].items():
                if rel_intensity < ignore_ion_threshold:
                    ignored_frag.add(f_key)
                if f_key not in ignored_frag:
                    # set union without ignored sequences
                    bc_seq[nonred_type, dp_minus1] |= prec_seq_copy[f_key] - ignored_seq 
        
    yz_seq = ddict(set) # sequences based on different y-ions
    for (frag_type, dp), frag2intensity in fragments.items():
        if frag_type in precursor.redend_types and dp > 1: 
            # y1 are not considered (and usually outside measured m/z range)
            for f_key, (intensity, rel_intensity) in frag2intensity.items():
                if rel_intensity < ignore_ion_threshold:
                    ignored_frag.add(f_key)
                if (f_key not in ignored_frag and 
                    len(prec_seq_copy[f_key] - ignored_seq) > 0):
                    allowed_ion_types.add((frag_type, dp))
                    yz_seq[frag_type, dp] |= prec_seq_copy[f_key] - ignored_seq
        elif frag_type in precursor.nonredend_types and dp > 1 and not precursor.nonred_maxonly:
            for f_key, (intensity, rel_intensity) in frag2intensity.items():
                if rel_intensity < ignore_ion_threshold:
                    ignored_frag.add(f_key)
                if (f_key not in ignored_frag and 
                    len(prec_seq_copy[f_key] - ignored_seq) > 0):
                    allowed_ion_types.add((frag_type, dp))
                    bc_seq[frag_type, dp] |= prec_seq_copy[f_key] - ignored_seq
                
                
    allowed_seq = set().union(*[s for s in bc_seq.values()]
                              ).intersection(*[s for s in yz_seq.values()])
    if verbose:
        print('Considered sequences:', 
              ', '.join(sorted([''.join(_) for _ in allowed_seq])))
        print('Check fragment ions:')
    """ only sequences that are supported by one of largest b-ions and a 
    complete y-ion series (except y1 and ignored y ions) are considered 
    """
    seq2i = ddict(list) 
    """
    seq2i = {sequences1: [(abs_intensity1, abundance1), 
                          (abs_intensity2, abundance2), ...], 
             sequences2: ...}
    """ 
    for (frag_type, dp), frag2intensity in sorted(fragments.items(), 
                                                  key = lambda x: (x[0][0], 
                                                                   -x[0][1]), 
                                                  reverse = True):
        # y-ions first, then b-ions; sort by increasing DP
        if not (frag_type, dp) in allowed_ion_types:
            continue
        print(frag_type, dp)
        print(frag2intensity)
        raw_i = [intensity for f_key, (intensity, rel_intensity)
                 in frag2intensity.items() if f_key not in ignored_frag]
        total_i = sum(raw_i) # used to weight contribution of fragment types
        max_i = max(raw_i)
        for f_key, (intensity, rel_intensity) in sorted(frag2intensity.items(),
                                                        key = lambda x: x[1][1]
                                                        ): 
            # evaluate fragments from low to high abundant
            if f_key in ignored_frag:
                continue
            prec_seq_copy[f_key] &= allowed_seq
            if verbose:
                print('\t{0}{1}'.format(frag_type, dp), f_key)
            if prec_seq_copy[f_key]:
                seq2i[tuple(sorted(prec_seq_copy[f_key]))
                      ].append((total_i, 
                                intensity / total_i)) # adjusted relative intensity
                if verbose:
                    print('\t   Intersection of potential precursor sequences and considered sequences:')
                    print('\t  ', ' + '.join(sorted([''.join(_) for _ in prec_seq_copy[f_key]])), '=', '{0:.2f}'.format(intensity / total_i))
                    for w, i in seq2i[tuple(sorted(prec_seq_copy[f_key]))]:
                        print('\t\t{1:.2f}\t{0:>10.0f}'.format(w,i))
            else: 
                """ fragment ion corresponds only to sequences that are not
                supported by a complete series of fragment ions
                """
#                 print('Fragment ion creates conflict')
                rel_to_max = intensity / max_i
                if verbose:
                    print('\t\tIntersection with considered sequences is empty -> conflict')
                    print('\t\tHeight is {0:.0%} of height of largest ion of this type'.format(rel_to_max))
                if (rel_to_max > ignore_ion_if_conflict_threshold and
                    (frag_type, dp) not in (('b', dp_minus1), 
                                            ('c', dp_minus1))):
                    """ conflicting fragment ion is not that small;
                    complete (y) fragment type will be ignored
                    """
                    if verbose:
                        print('\t\t-> all fragments of this type will be ignored')
                    ignored_frag |= set(frag2intensity.keys()) 
                    return determine_sequences(precursor = precursor,
                                               fragments = fragments,
                                               prec_seq = prec_seq,
                                               ignored_frag = ignored_frag)
                else:
                    """ conflicting fragment ion is small and will be ignored
                    """
                    if verbose:
                        print('\t\t-> this fragment will be ignored')
                    ignored_frag.add(f_key)
                    return determine_sequences(precursor = precursor,
                                               fragments = fragments,
                                               prec_seq = prec_seq,
                                               ignored_frag = ignored_frag)
                    
    """ actual sequencing starts here
    """
    if verbose:
        print('Start sequencing:')
    last_seq2i = []
    while last_seq2i != sorted(list(seq2i)):
        if verbose:
            print('\t..........')
        last_seq2i = sorted(list(seq2i))
        seq2del = set()
        new_seq2i = ddict(list)
        for seq1, i1 in seq2i.items():
            # seq1 and seq2 can be the same only once (dict keys)
            new_seq2i[seq1] += i1 
            for seq2, i2 in seq2i.items():
                if set(seq2) < set(seq1): 
                    """ seq2 is subset of seq1 and seq2 != seq1
                    -> seq2 can be subtracted from seq1 and intensity of 
                    remaining sequences in seq1 can be calculated;
                    original seq1 can be removed from equation system
                    """
                    seq2del.add(seq1)
                    set_difference = tuple(sorted(set(seq1) - set(seq2)))
                    if verbose:
                        print('\t', seq2, 'is subset of', seq1)
                        print('\t\t Remaining sequence(s):', set_difference)
                    for total_i1, rel_intensity1 in i1:
                        for total_i2, rel_intensity2 in i2:
                            new_seq2i[set_difference].append((
                                # new weighting -> sum of both equations
                                (total_i1 + total_i2),# / 2, 
                                # new abundance -> subtract intensity of subset
                                # from intensity of larger set
                                rel_intensity1 - rel_intensity2 
                                ))
                            if verbose:
                                print('\t\t {1:.2f}\t{0:>10.0f}'.format(total_i1 + total_i2,
                                                                       rel_intensity1 - rel_intensity2))
                            
        seq2i = ddict(list)
        for seqs, i in new_seq2i.items():
            if seqs in seq2del:
                continue
            # values have to averaged here; otherwise list can get too long
            seq2i[seqs].append((
                # sum weightings
                sum([total_i for total_i, rel_intensity in i]),# / len(i),
                # weighted average
                sum([total_i * rel_intensity for total_i, rel_intensity in i]) / 
                sum([total_i for total_i, rel_intensity in i])
                ))
            if verbose:
                print('\t  ', ' + '.join([''.join(_) for _ in seqs]), '=', '{0:.2f}'.format(sum([total_i * rel_intensity for total_i, rel_intensity in i]) / 
                sum([total_i for total_i, rel_intensity in i])))
        """ this loop is repeated until there are no more subsets of other sets
        found within the equation system 
        """
        
    solved_seq2i = {seqs: i[0][1] for seqs, i in seq2i.items()} 
    # all sequences in solved equation system (even multiples)
    all_seq = [s for seqs in solved_seq2i for s in seqs] 
    out_rows = []
    if verbose:
        print('\tNo more subsets found within equation system')
        print('Check results:')
    for seqs, rel_intensity in sorted(solved_seq2i.items(), key = lambda x: x[1]):
        # sort by increasing intensity
        if (rel_intensity < ignore_single_seq_threshold or 
            (len(seqs) > 1 and rel_intensity < ignore_unresolved_seqs_threshold)):
            """ if the intensity is very low the corresponding sequence is 
            neglected; especially when the corresponding sequence is not 
            completely resolved
            """
            ignored_seq |= set(seqs)
            if verbose:
                print('\t Unresolved sequences -> start sequencing again, neglecting sequences:')
                print('\t', sorted(seqs), round(rel_intensity,4))
            return determine_sequences(precursor = precursor,
                                       fragments = fragments,
                                       prec_seq = prec_seq,
                                       ignored_frag = ignored_frag,
                                       ignored_seq = ignored_seq)
            
        notes = set()
        for s in seqs:
            if all_seq.count(s) > 1:
                notes.add("not completely resolved")
        if len(seqs) > 1:
            # more than one possible sequence
            pos2monomers = ddict(set)
            for s in seqs:
                for pos, monomer in enumerate(s):
                    pos2monomers[pos].add(monomer)
            merged_sequence = ''
            for pos, monomers in sorted(pos2monomers.items()):
                if len(monomers) > 1:
                    merged_sequence += 'x'
                else:
                    merged_sequence += list(monomers)[0]
        else: # single sequence possible
            merged_sequence = ''.join(seqs[0])
        if len(ignored_frag) > 0:
            notes.add('neglected {0}'.format(
                ', '.join(['{0:.2f} ({1})'.format(mz, name)
                                                for mz, name in sorted(ignored_frag)])))
        if len(ignored_seq) > 0:
            notes.add('neglected {0}'.format(', '.join(
                list(sorted([''.join(_) for _ in ignored_seq])))))
        
        out_rows.append({
            'Relative intensity (w.r.t. oligo)': rel_intensity,
            'Notes': '; '.join(sorted(notes)),
#             'All possible sequences': ', '.join(seqs),
            'All possible sequences': ', '.join([_ for _ in 
                                                 [''.join(__) for __ in 
                                                  sorted(seqs)]]),
            'Sequence': merged_sequence
            })
        
    return out_rows
                    
            