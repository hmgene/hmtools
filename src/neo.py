
def cool2mar(clr, chunksize=int(1e7), exclude=['M', 'Y', 'MT', 'EBV'] ):
    map_ = map
    nnz = clr.info['nnz']
    n_bins = clr.info['nbins']
    edges = np.arange(0, nnz+chunksize, chunksize)
    spans = list(zip(edges[:-1], edges[1:]))
    marg = (
        balance.split(clr, spans=spans, map=map_, use_lock=False)
            .prepare(balance._init)
            .pipe([])
            .pipe(balance._marginalize)
            .reduce(balance.add, np.zeros(n_bins))
    )
    table = clr.bins()[:][['chrom', 'start', 'end']]
    table['Coverage'] = marg.astype(int)
    pool = []
    chroms = [c for c in clr.chromnames if ((not c.lstrip('chr') in exclude) and (not '_' in c))]
    for chrom in chroms:
        pool.append(table[table['chrom']==chrom])
    table = pd.concat(pool)
    return table, clr.binsize
