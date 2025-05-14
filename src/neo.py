import numpy as np
import os, cooler
#def cool2mar(clr, chunksize=int(1e7), exclude=['M', 'Y', 'MT', 'EBV'] ):
#    map_ = map
#    nnz = clr.info['nnz']
#    n_bins = clr.info['nbins']
#    edges = np.arange(0, nnz+chunksize, chunksize)
#    spans = list(zip(edges[:-1], edges[1:]))
#    marg = (
#        balance.split(clr, spans=spans, map=map_, use_lock=False)
#            .prepare(balance._init)
#            .pipe([])
#            .pipe(balance._marginalize)
#            .reduce(balance.add, np.zeros(n_bins))
#    )
#    table = clr.bins()[:][['chrom', 'start', 'end']]
#    table['Coverage'] = marg.astype(int)
#    pool = []
#    chroms = [c for c in clr.chromnames if ((not c.lstrip('chr') in exclude) and (not '_' in c))]
#    for chrom in chroms:
#        pool.append(table[table['chrom']==chrom])
#    table = pd.concat(pool)
#    return table, clr.binsize

import cooler
import numpy as np
import pandas as pd

def cool2mar(clr, chunksize=int(1e7), exclude=['M', 'Y', 'MT', 'EBV']):
    nnz = clr.info['nnz']
    n_bins = clr.info['nbins']

    # Initialize marginals
    marg = np.zeros(n_bins)

    # Chunk-based marginal aggregation
    for start in range(0, nnz, chunksize):
        stop = min(start + chunksize, nnz)
        chunk = clr.pixels()[start:stop]
        marg += np.bincount(chunk['bin1_id'], weights=chunk['count'], minlength=n_bins)
        marg += np.bincount(chunk['bin2_id'], weights=chunk['count'], minlength=n_bins)

    # Create bin table
    table = clr.bins()[:][['chrom', 'start', 'end']].copy()
    table['Coverage'] = marg.astype(int)

    # Filter chromosomes
    chroms = [c for c in clr.chromnames if c not in exclude and '_' not in c.lstrip('chr')]
    table = table[table['chrom'].isin(chroms)]

    return table.reset_index(drop=True), clr.binsize


def cool2marginals(clr_path, chunksize=int(1e7), exclude=['M', 'Y', 'MT', 'EBV']):
    clr = cooler.Cooler(clr_path)
    bins = clr.bins()[:]

    # Filter bins by excluded chromosomes
    valid_bins = bins[~bins['chrom'].isin(exclude)].index.values
    marginals = np.zeros(clr.info['nbins'])

    # Read interactions in chunks and compute marginals
    for start in range(0, clr.info['nnz'], chunksize):
        stop = min(start + chunksize, clr.info['nnz'])
        chunk = clr.pixels()[start:stop]
        marginals_chunk = np.bincount(chunk['bin1_id'], weights=chunk['count'], minlength=clr.info['nbins'])
        marginals += marginals_chunk
        # Since the matrix is symmetric
        marginals_chunk = np.bincount(chunk['bin2_id'], weights=chunk['count'], minlength=clr.info['nbins'])
        marginals += marginals_chunk

    # Mask excluded chromosomes
    marginals[~np.isin(np.arange(len(marginals)), valid_bins)] = np.nan

    return marginals
