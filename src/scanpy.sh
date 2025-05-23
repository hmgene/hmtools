
scanpy-mtxs2h5ad(){
usage="$FUNCNAME <input.csv: sample_name,mtx_file> <output.h5ad>"
if [ $# -lt 1 ];then echo "$usage";return;fi

cat $1 | python <( echo '
import pandas as pd
import scanpy as sc
import sys
inp=pd.read_csv(sys.stdin,header=None)
out="'${2%.h5ad}.h5ad'";

tt=None;

for i,row in inp.iterrows():
    s=row[0]
    f=row[1]
    x=sc.read_mtx(f).T
    b=pd.read_csv(f.replace("matrix.mtx","barcodes.tsv"),sep="\t",header=None)
    g=pd.read_csv(f.replace("matrix.mtx","features.tsv"),sep="\t",header=None)
    x.obs.index = b[0].values
    x.var.index = g[1].values   ## Gene
    x.var_names_make_unique()
    x.obs["sample"] = s; #re.search(r"(AXH[O\d]+)", f).group(1)
    if tt is None:
        tt=x
    else:
        tt=tt.concatenate(x)
    tt.write_h5ad(out)
')
}

scanpy-mtxs2h5ad-test(){
echo "
AXHO16,bigdata/local_sc/AXHO16_COL_1/filtered_feature_bc_matrix/matrix.mtx.gz
AXHO16,bigdata/local_sc/AXHO16_COL_2/filtered_feature_bc_matrix/matrix.mtx.gz
AXHO16,bigdata/local_sc/AXHO16_COL_3/filtered_feature_bc_matrix/matrix.mtx.gz
AXHO16,bigdata/local_sc/AXHO16_COL_4/filtered_feature_bc_matrix/matrix.mtx.gz
" | grep -v "^$" | scanpy-mtxs2h5ad  - AXHO16.h5ad
}

scanpy-prepro(){
usage="$FUNCNAME <h5ad> [<o.h5ad> [mt=25,ngenes=200,ncells=3]"];
if [ $# -lt 1 ];then echo "$usage";return;fi
i=${1%.h5ad}.h5ad
o=${2:-$i} ## override

python <( echo '
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
inp="'$i'"
out="'$o'"
sc.settings.figdir = inp.replace(".h5ad","") + "_figs";

tt=sc.read_h5ad(inp)
tt.var["mt"] = tt.var_names.str.startswith("MT-")
tt.obs["pct_counts_mt"] = np.sum(tt[:, tt.var["mt"]].X, axis=1).A1 / np.sum(tt.X, axis=1).A1 * 100
tt.var["ribo"] = tt.var_names.str.startswith(("RPS", "RPL"))
tt.var["hb"] = tt.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics( tt, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
sc.pp.filter_cells(tt, min_genes=200)
sc.pp.filter_genes(tt, min_cells=3)
tt= tt[tt.obs["pct_counts_mt"] < 25, :]
tt.layers["counts"] = tt.X.copy()
sc.pl.violin(tt, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True,save="_filter_summary.png",show=False)

tt.obs["n_genes"] = (tt.X > 0).sum(axis=1).A1
sc.pp.filter_cells(tt, min_genes=200)
sc.pp.filter_genes(tt, min_cells=3)
tt.layers["counts"] = tt.X.copy()
sc.pp.normalize_total(tt)
sc.pp.log1p(tt)
sc.pp.highly_variable_genes(tt, n_top_genes=2000, batch_key="sample")
#sc.pl.highly_variable_genes(tt)
sc.tl.pca(tt)
#sc.pl.pca_variance_ratio(tt, n_pcs=50, log=True)
tt.write_h5ad(out)
')
}

scanpy-cluster(){
usage="$FUNCNAME <adata.h5ad> [res=0.2] [out.h5ad] ";if [ $# -lt 1 ];then echo "$usage";return; fi
i=${1%.h5ad}.h5ad; r=${2:-0.2}; o=${3:-${i%.h5ad}_res$r.h5ad}

python <( echo '
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
inp="'$i'";
res='$r';
out="'$o'";
sc.settings.figdir = out.replace(".h5ad","") + "_figs";

tt=sc.read_h5ad(inp)
sc.pp.neighbors(tt)
sc.tl.umap(tt)
sc.tl.leiden(tt,resolution=res, flavor="igraph", n_iterations=2)
sc.tl.dendrogram(tt, groupby="leiden")
sc.pl.umap( tt, color="leiden", save="_leiden.png",show=False)
sc.pl.umap( tt, color="sample", save="_sample.png",show=False)

sc.tl.rank_genes_groups(tt, groupby="leiden", method="t-test")
sc.pl.rank_genes_groups_dotplot(tt, groupby="leiden", standard_scale="var", n_genes=5,save="_leiden_5top.png",show=False)
#sc.pl.umap(tt, color="CFTR",save="_CFTR.png")

result = tt.uns["rank_genes_groups"]
groups = result["names"].dtype.names  # list of cluster/group names
top_n = 10
records = []
for group in groups:
    for i in range(top_n):
        records.append({
            "group": group,
            "gene": result["names"][group][i],
            "logfoldchange": result["logfoldchanges"][group][i],
            "pval_adj": result["pvals_adj"][group][i],
            "score": result["scores"][group][i],
        })

df = pd.DataFrame(records)
df.to_csv( out.replace(".h5ad","")+f"_top{top_n}.tsv",sep="\t")
tt.write_h5ad(out)
')

}

