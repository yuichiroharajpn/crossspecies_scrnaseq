library(Seurat)
#library(stringr)
require("Seurat.utils")
options(Seurat.object.assay.version = "v3")

args <- commandArgs(trailingOnly = T)

dir <- args[1]
file <- paste(dir, "out", "filtered_feature_bc_matrix", sep="/")
prefix <- args[3]
out_prim_rds <- paste(prefix, "prim.rds", sep=".")
out_rds <- paste(prefix, "ortho1to1.rds", sep=".")

objs_in <- Read10X(file, unique.features = TRUE)
objs_org <- CreateSeuratObject(counts=objs_in, min.cells=0, min.features = 200, project="P2")
objs_org <- NormalizeData(objs_org, verbose = FALSE)

saveRDS(objs_org, out_prim_rds)

df_ortho <- read.delim(args[2])
if(length(args)>4){
	ref <- c(paste(gsub(args[5], ' ', '.'), "gene.name", sep="."))
	new <- as.list(df_ortho[,ref])
}else{
	new <- as.list(df_ortho[,2])
}
	
tgt <- c(paste(gsub('[[:blank:]]', '.', args[4]), "gene.name", sep="."))
names(new)<-as.list(df_ortho[,tgt])

objs_org_s <- subset(objs_org, features=df_ortho[,tgt])
common <- as.character(unlist(new[rownames(objs_org_s@assays$RNA@counts)]))

objs<-RenameGenesSeurat(obj = objs_org_s, newnames = common)
rownames(objs@assays$RNA@meta.features) <-as.character(unlist(new[rownames(objs@assays$RNA@meta.features)]))

saveRDS(objs, out_rds)

