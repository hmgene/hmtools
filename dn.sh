#wget https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar; mv picard.jar bigdata/
#curl -o bigdata/cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1747707394&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=ZuPMlQkS~R4f04JEvQPDHi2dYjJsz8AasjSZzlRroH~o8K4O2-PPoApa79gCy5FyTpUFMlrOWlbksZJBKczyY8wh24kgiiY21ckD8ovl6euzOi2uU3N-kurWpgY4gYhQtTqJiDi1ApNrKCROHsITBTpO1~W6MuKhXNezeBTQyJhAmFO5rBcXDDZqIzfOTbbNT3z58yOf26v~v8-f9CgkRvPkElOxB1jKCpPw4X-Cbrda~dv0PajxDz6CisECrPmAkM~vv3IlQJe-fraJjG72UC~Fi3LHjghp9f39InpNWM~6rOIDqd0f7C7aejweVYOj5QXJDqdFyOm8WEgf8tdgpg__"
#curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" -o bigdata/refdata-gex-GRCh38-2024-A.tar.gz

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes; mv hg19.chrom.sizes data/
