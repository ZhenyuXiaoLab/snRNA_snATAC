##use liger with all Granja's cicero gene activity matrix



library(liger)# 0.5.0
##library(Seurat,lib.loc = "/home/mjwang/.conda/envs/myenv/lib/R/library_opt/") #3.2.2
library('Seurat') #3.2.3
library('ggplot2')

#library(harmony)

library('magrittr')

library('RANN')

library('patchwork')

library(hexbin)
library(RColorBrewer)

#library("ggpointdensityplot")
library(viridis)


library(grid)
library(gridExtra)


library('Peacock.test')



###look up color sets in ArchR###
library('ArchR')
color_list <- ArchR::ArchRPalettes

options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_list)){ 
  len = length(color_list[[name]])
  color = color_list[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}

saveRDS(color_list,'ArchR.color_list.rds')

##


##python yellowbrick color palettes##
library(rlist)
color_set_yellowbrick <- readRDS('color_set_yellowbrick.rds')
color_set_yellowbrick.flat <- list.flatten(color_set_yellowbrick)
options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_set_yellowbrick.flat)){ 
  len = length(color_set_yellowbrick.flat[[name]])
  color = color_set_yellowbrick.flat[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}

saveRDS(color_set_yellowbrick.flat,'color_set_yellowbrick.flat.rds')

############Buen colors########
library(BuenColors)
color_set0 <- jdb_color_maps #17 different colors
names(color_set0) <- NULL
#plot(1:17,1:17,pch = 19, cex = 5,col=jdb_color_maps)

#discrete colors
color_set1 <- jdb_palette("solar_extra") #9 discrete but gradient colors
color_set2 <- jdb_palette("brewer_spectra") #9 discrete but gradient colors
color_set3 <- jdb_palette("flame_light") #9 discrete but gradient colors, good!

color_set3_ext12 <- colorRampPalette(colors = as.character(color_set3))(12)
color_set3_ext17 <- colorRampPalette(colors = as.character(color_set3))(17)

#############ArchR colors############
#hmcols <- colorRamps::blue2green2red(length(bks) ) #colors
color_peak <- ArchR::paletteContinuous(set = 'solarExtra',n=256,reverse=FALSE)  
color_tfdev = ArchR::paletteContinuous(set = 'blueYellow',n=257,reverse=FALSE)                       
#color_ga <- paletteContinuous(set='solarExtra',n=257,reverse=FALSE) 
#color_ga <- paletteContinuous(set='horizon',n=257,reverse=FALSE)                      
#color_ga <- paletteContinuous(set='horizonExtra',n=257,reverse=FALSE)  #good        
color_rna <- ArchR::paletteContinuous(set='greenBlue',n=256,reverse=FALSE)
#color_ga <- paletteContinuous(set='blueYellow',n=257,reverse=FALSE)
color_ga <- ArchR::paletteContinuous(set='greyMagma',n=257,reverse=FALSE)

color_rna <- colorRampPalette(c('grey','red'))(10) 

#########customized colors########
color_snap = c('1'='grey','2'='#E31A1C','3'='#FFD700','4'='#771122','5'='#777711','6'='#1F78B4','7'='#68228B','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')
#names(color_snap) <- NULL

#modified for CTB with dark red colors
color_snap_mod1 = c('1'='#777711','2'='#E31A1C','3'='#68228B','4'='#771122','5'='grey','6'='#1F78B4','7'='#FFD700','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')

color_signac = c(
'0'='#E6D55E','1'='#792B8A','2'='#DA703D','3'='#9DC8E5','4'='#BA273C','5'='#C2C184','6'='#7F8084','7'='#65AB53','8'='#D082AF','9'='#496EAB','10'='#DE896D','11'='#491F8B','12'='#E1AD49','13'='#8E1B85','14'='#E7EE77','15'='#7D1A1D','16'='#96B355')
names(color_signac) <- NULL



color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" )

color_good <- c("#E7D654", "#6F1482" ,"#DC7035", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,"#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,
 "#CA362E" ,"#2B3918","#1E1E1E" )


color_good <- c("#E7D654", "#6F1482" ,"#DC7035", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", "#63AC4E", "#D181B0" ,
                "#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,"#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,
                "#CA362E" ,"#2B3918","#1E1E1E" )


color_set <- c('0'   =    '#2679b4',
'1'   =    '#fd7f28',
'2'   =    '#2f9d6a',
'3'    =   '#d42a2f',
'4'    =   '#a94af8',
'5'    =   '#8b564c',
'6'    =   '#e17ac1',
'7'   =    '#b5bc66',
'8'   =    '#29bece',
'9'    =   '#b1cae9',
'10'  =    '#fdba7d',
'11'   =   '#9ade8d',
'12'   =   '#fd9998',
'13'  =    'grey',
'14'  =    'navy',
'15'  =  'brown',
'16'   =   'black',
'17'   =   'yellow'
              )



color_gradient_my <- c(
    rgb(5,48,97,maxColorValue = 255),
    rgb(42,113,178,maxColorValue = 255),
    rgb(147,198,222,maxColorValue = 255),
    rgb(239,243,245,maxColorValue = 255),
    rgb(253,219,199,maxColorValue = 255),
    rgb(214,96,77,maxColorValue = 255),
    rgb(121,5,34,maxColorValue = 255)

)

color_cellranger <-c('#820610','#C50F1E','#F42428','#F86D30','#FBB33D','#FCFB4E','#C0FB61','#87FB8F','#41FEF9','#2BAED7','#155CB1','#08238D')

###select one global color set###
color <- color_good



sample <- 'placenta 8 weeks RNA integrate with ATAC'


######read in predefined cluster
#rna
#seurat.obj <- readRDS('../../02.seurat_harmony/PLA-8w-RNA.final.rds')


# seurat.obj  <- RenameIdents(seurat.obj,'0'='1', '1'='2','2'='3','3'='4','4'='5',
#            '5'='6','6'='7','7'='8','8'='9','9'='10',
#            '11'='11','12'='12','13'='13','14'='14','15'='15','16'='16' )


#cluster.df.add <- readRDS('../../02.seurat_harmony/cluster.df.add.rds')
##length(Idents(seurat.obj) ) #11206
#all.equal(cluster.df.add$cluster,Idents(seurat.obj),check.attributes = FALSE) #TRUE

##rna_cluster <- Idents(seurat.obj)


rna_cluster.tab <- readRDS("../../02.seurat_harmony/cluster.df.add.rds") #11206 full
rna_cluster <- rna_cluster.tab$cluster
names(rna_cluster) <- rownames(rna_cluster.tab)
table(rna_cluster)

   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74 

#rna te only
rna_cluster.te.tab <- readRDS("../../02.seurat_harmony/cluster.df.add.final.rds") #10198 te only
rna_cluster.te <- rna_cluster.te.tab$cluster
names(rna_cluster.te) <- rownames(rna_cluster.te.tab)

table(rna_cluster.te)
   1    2    3    4    5    6    7    8    9   10 
1720 1428 1225 1093 1083  962  892  852  655  288 


table(rownames(rna_cluster.te.tab) %in% rownames(rna_cluster.tab) )

 TRUE 
10198

# table(rna_cluster)
# rna_cluster
#    0    1    2    3    4    5    6    7    8    9   11   12   13   14   15   16 
# 1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74 



# rna_cluster <- cluster.df.add$cluster
# names(rna_cluster) = rownames(cluster.df.add) 

# all.equal(Idents(seurat.obj) ,rna_cluster) #TRUE


# table(rna_cluster)
# rna_cluster
#    0    1    2    3    4    5    6    7    8    9   11   12   13   14   15   16 
# 1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74 

# #modify levels to remove 0
# cluster.new <- as.numeric(as.character(rna_cluster)) + 1
# cluster.new <- factor(cluster.new,levels=sort(unique(cluster.new)) )
# names(cluster.new) = rownames(cluster.df.add) 
# rna_cluster <- cluster.new

# table(rna_cluster)
# rna_cluster
#    1    2    3    4    5    6    7    8    9   10   12   13   14   15   16   17 
# 1720 1428 1225 1093 1083  963  893  853  656  365  288  189  173  123   80   74

# #rename cluster id

# mapid <- c('1'='1', '2'='2','3'='3','4'='4','5'='5',
#            '6'='6','7'='7','8'='8','9'='9','10'='10',
#            '12'='11','13'='12','14'='13','15'='14','16'='15','17'='16'   )


# rna_cluster <- mapid[rna_cluster]



# names(rna_cluster) <- rownames(cluster.df.add)

# table(rna_cluster)



#add prefix 
names(rna_cluster) <- paste0('placenta.rna_',names(rna_cluster))

##saveRDS(rna_cluster,'rna_cluster.rds')

names(rna_cluster.te) <- paste0('placenta.rna_',names(rna_cluster.te))
saveRDS(rna_cluster.te,'rna_cluster.te.rds')


# rna_cluster.tab = read.table("scanpy_seurat_cluster_all_df.txt",header=T,stringsAsFactors = TRUE,row.names = 1)
# rna_cluster = factor(rna_cluster.tab$cluster, levels=sort(unique(rna_cluster.tab$cluster))  )
# names(rna_cluster) = rownames(rna_cluster.tab) 



#atac 
cluster.atac.obj <- readRDS("../../02.snapATAC_harmony/cluster.df.add.final.rds") #12526 full version
atac_cluster <- cluster.atac.obj$cluster
names(atac_cluster) = rownames(cluster.atac.obj)

##rename to -1, -2
cellid <- names(atac_cluster) 
idy <- which(grepl(pattern = "^placenta_donor2",x = cellid    ))
cellid[idy] <- gsub( pattern = '-1$',replacement = '-2' ,x =  cellid [idy]  )
cellid  <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = cellid )

names(atac_cluster) <- cellid 

#add prefix
names(atac_cluster) <- paste0('placenta.atac_',names(atac_cluster))

table(atac_cluster)

atac_cluster
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
1832 1783 1090 1546 1475 1509 1019  444  395  352  290  253  235  171  132

##saveRDS(atac_cluster,'atac_cluster.rds')

##atac te only
cluster.atac.te.obj <- readRDS("../../02.snapATAC_harmony/cluster.df.add.final.final.rds") #11093 te only
atac_cluster.te <- cluster.atac.te.obj$cluster
names(atac_cluster.te) = rownames(cluster.atac.te.obj)

##rename to -1, -2
cellid <- names(atac_cluster.te) 
idy <- which(grepl(pattern = "^placenta_donor2",x = cellid    ))
cellid[idy] <- gsub( pattern = '-1$',replacement = '-2' ,x =  cellid [idy]  )
cellid  <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = cellid )

names(atac_cluster.te) <- cellid 

#add prefix
names(atac_cluster.te) <- paste0('placenta.atac_',names(atac_cluster.te))

table(atac_cluster.te)

atac_cluster.te
   1    2    3    4    5    6    7    8    9 
1832 1783 1090 1546 1475 1509 1019  444  395 

saveRDS(atac_cluster.te,'atac_cluster.te.rds')


table(rownames(cluster.atac.te.obj) %in% rownames(cluster.atac.obj) )
TRUE 
11093


# #atac_cluster.tab = read.table("seuratObject.PLA-8w-ATAC-1.umap.cl.8029.txt",header=T,stringsAsFactors = TRUE)
# atac_cluster.tab = read.table("snapATAC.PLA-8w-ATAC-1.umap.bin50k.filtered.txt",header=T,stringsAsFactors = TRUE)

# atac_cluster = factor(atac_cluster.tab$cluster, levels = 0:15 ) #remove c16, but not removed, become NA
# names(atac_cluster) = rownames(atac_cluster.tab)

# #atac_cluster.tab[is.na(atac_cluster.tab),] ##get NA records

# #atac_cluster.new <- droplevels(atac_cluster)

# atac_cluster.tab$cluster = atac_cluster



# #cisTopic
# #atac_cluster.tab.ct = read.table("cisTopicObject.PLA-8w-ATAC-1.umap.louvaincl.txt",header=T,stringsAsFactors = TRUE)
# #colnames(atac_cluster.tab.ct) = c('cluster',"UMAP_1","UMAP_2")
# atac_cluster.tab.ct = read.table("cisTopicObject.PLA-8w-ATAC-1.tsne.louvaincl.txt",header=T,stringsAsFactors = TRUE)
# colnames(atac_cluster.tab.ct) = c('cluster',"tSNE_1","tSNE_2")

# #atac_cluster.ct = factor(atac_cluster.tab.ct$cluster-1, levels = 0:17 )
# atac_cluster.ct = factor(atac_cluster.tab.ct$cluster-1, levels = 0:20 )
# names(atac_cluster.ct) = rownames(atac_cluster.tab.ct)
# atac_cluster.tab.ct$cluster = atac_cluster.ct

# #cellranger-atac
# atac_cluster.tab.cr = read.table("cellranger-atac.cl.txt",header=T,stringsAsFactors = TRUE,row.names = 1)
# colnames(atac_cluster.tab.cr) = c('cluster',"tSNE_1","tSNE_2")

# atac_cluster.cr = factor(atac_cluster.tab.cr$cluster-1, levels = 0:11 )
# names(atac_cluster.cr) = rownames(atac_cluster.tab.cr)
# atac_cluster.tab.cr$cluster =atac_cluster.cr

# #episcanpy
# atac_cluster.tab.sp = read.table("adatabin400reg.umap.tsne.cl.txt",header=T,stringsAsFactors = TRUE,row.names = 1)
# colnames(atac_cluster.tab.sp) = c('cluster',"tSNE_1","tSNE_2","UMAP_1","UMAP_2")

# atac_cluster.sp = factor(atac_cluster.tab.sp$cluster, levels = 0:12 )
# names(atac_cluster.sp) = rownames(atac_cluster.tab.sp)
# atac_cluster.tab.sp$cluster =atac_cluster.sp



#####gene activity matrix from fragments counts#######

genepromoter.bc <- read.table(file = "atac_activity_new/out.gr.findoverlap.aggre.useR.txt", sep = "\t", as.is = c(4,7), header = FALSE) #use all genes, new findOverlap method instead of bedmap(may have some unknow problem with DNMT1 count..)

#genepromoter.bc <- read.table(file = "atac_activity/genebodyandpromoter.GRCh38.all.counts", sep = "\t", as.is = c(4,7), header = FALSE) #use all genes
#genepromoter.bc <- read.table(file = "atac_activity_norm/genebodyandpromoter.GRCh38.all.counts", sep = "\t", as.is = c(4,7), header = FALSE) #use all genes
barcodes = cellid
#head(barcodes)

#placenta.atac.new <- liger:::makeFeatureMatrix(genepromoter.bc.new, barcodes)

placenta.atac <- liger:::makeFeatureMatrix(genepromoter.bc, barcodes) #57840 x 12526 #matrix 63970 x 12526 #matrix: 63970 x 8029 , ##63970 x 7416 #22285 x 7416
max(placenta.atac) #336 #176 #336
#max(placenta.atac.new) #336

##add prefix to rowname
all.equal(colnames(placenta.atac),cellid) #TRUE
colnames(placenta.atac) <- names(atac_cluster)

#all.equal(colnames(placenta.atac.new),cellid) #TRUE
#colnames(placenta.atac.new) <- names(atac_cluster)


#saveRDS(placenta.atac,"liger.gene_activity.allgene.rds")
#placenta.atac = readRDS("liger.gene_activity.allgene.rds")

saveRDS(placenta.atac,"placenta.atac.rds")
#saveRDS(placenta.atac,"placenta.atac.new.rds")


#placenta.atac <- placenta.atac.new
#57840 x 12526

###read in 10x rna raw count table###
#rna.dir = '../../01.aggregation/PLA-8w-RNA-aggre/filtered_feature_bc_matrix/'
rna.dir = '../../01.aggregation/PLA-8w-RNA-aggre_nonorm/filtered_feature_bc_matrix/'
placenta.rna = liger::read10X(sample.dirs = list(rna.dir), sample.names = list('placenta.rna') ) 
#dgCMatrix 33538 x 15407
#dgCMatrix 33538 x 9496



# ####gene activity from SnapATAC gmat and gene expression mat from Seurat data
# gmat <- readRDS('gmat.rds') #magic normalized SnapATAC gene body activity score
# exprMat <- seurat.obj@assays$RNA@counts




##create liger object
#check colnames atac
all.equal (colnames(placenta.atac), names(atac_cluster)) #TRUE
#check colnames rna
table(names(rna_cluster) %in% colnames(placenta.rna) )
 TRUE 
11206 

table(colnames(placenta.rna) %in% names(rna_cluster)  )
FALSE  TRUE 
 4201 11206 

saveRDS(placenta.rna,'placenta.rna.rds')


##a.placenta = createLiger (list(atac=placenta.atac[,names(atac_cluster)], rna=placenta.rna[,names(rna_cluster)]  ) ) # atac: 30448 x 12526 rna: 26386 x 11206   #atac: 30475 x 8029, rna: 25401 9461; #atac: 30472 x 7416, rna: 25401 x 9461


atac_cluster_d1 <-  names(atac_cluster)[which( grepl(pattern = "-1$",x = names(atac_cluster)) )]
atac_cluster_d2 <-  names(atac_cluster)[which( grepl(pattern = "-2$",x = names(atac_cluster)) )]
length(atac_cluster_d1) + length(atac_cluster_d2) == length(atac_cluster) #TRUE, 12526


rna_cluster_d1 <-  names(rna_cluster)[which( grepl(pattern = "-1$",x = names(rna_cluster)) )]
rna_cluster_d2 <-  names(rna_cluster)[which( grepl(pattern = "-2$",x = names(rna_cluster)) )]
length(rna_cluster_d1) + length(rna_cluster_d2) == length(rna_cluster) #TRUE, 11206


a.placenta = createLiger (
    list(atac1=placenta.atac[,atac_cluster_d1],  #57840 x 5042 #63970 x 5042
         atac2=placenta.atac[,atac_cluster_d2],  #57840 x 7484 #63970 x 7484
         rna1=placenta.rna[,rna_cluster_d1],  #33538 x 6658
         rna2=placenta.rna[,rna_cluster_d2]   #33538 x 4548
         
        ) 
) 


##real dim after filter:
#atac1 57518 x 5042
#atac2 57204 x 7484
#rna1 24993 x 6658
#rna2 25066 x 4548

#former liger raw.data
#30373 5042
#30087 7484
#24993 6658
#25066 4548

a.placenta <- normalize(a.placenta ) 
#a.placenta <- selectGenes(a.placenta , datasets.use = 2,do.plot = T) #only use var genes of rna
a.placenta <- selectGenes(a.placenta , datasets.use = 3,do.plot = T,var.thresh = 0.2) #use RNA from d1, 1584 genes with var.thresh=0.2,  3286 default 

#or use given variable genes?
#var.genes <- readRDS("var.genes.old.rds") #1857 of the old liger var.gene (use problemic gene activity mat)
#a.placenta@var.genes <- var.genes #1857 var genes

a.placenta <- scaleNotCenter(a.placenta )

rm(genepromoter.bc)#,placenta.atac,placenta.rna)
gc()

# running suggestK on multiple cores can greatly decrease the runtime
#k.suggest <- suggestK(a.placenta, num.cores = 15, gen.new = T, plot.log2 = F, ##unused argument (return.results = T)
#   nrep = 5)
a.placenta <- optimizeALS(a.placenta, k=35) #first use k = 20 for a try

#a.placenta <- optimizeALS(a.placenta, k=30, lambda = 2) #default lambda = 5, tune down if datasets difference are expected to be relatively small

#a.placenta <- optimizeALS(a.placenta, k=35, lambda = 1) #default lambda = 5, tune down if datasets difference are expected to be relatively small



##before quantile norm for each lib
a.placenta <- runTSNE(a.placenta, use.raw = T)
p1 <- plotByDatasetAndCluster(a.placenta, return.plots = T)
print(p1[[1]])



#a.placenta <- quantileAlignSNF(a.placenta)
a.placenta <- quantile_norm(a.placenta) #resolution = 0.4, small.clust.thresh = 20 #including assign cluster steps
#Warning message in regularize.values(x, y, ties, missing(ties)):
#“collapsing to unique 'x' values”


# ###########harmony to NMF matrix? no, do not over normalize############
# V <- a.placenta@H.norm #not a.placenta@V
# 23732 x 20
# #11560 x 50
# #15107 x 50

# meta_data <- cluster.all.coord #see below code
# 23732 x 6
# #11560 x 6
# #15107 x 6

# #meta_data$sample = ifelse(grepl(pattern = "-1$", x=rownames(meta_data) ), 'D1', 'D2'  )
# #meta_data$sample <- factor(meta_data$sample,levels=c('D1','D2'))
# #table(meta_data$sample)
# #  D1   D2 
# #6996 4564 

# #  D1   D2 
# #9461 5646 


# all.equal(rownames(V),rownames(meta_data)) #TRUE
# V_harmony <- HarmonyMatrix(V,meta_data, 'anno', do_pca = FALSE)
# #Harmony converged after 3 iterations
# #23732 x 20

# #Harmony converged after 6 iterations
# #11560 x 50

# #do not converged after 10 iterations
# #15107 x 50

# #a.placenta@H.norm <- V_harmony #overwrite the NMF

# ####compare dimred before and after harmony
# # dimred.beforeharmony <- uwot::umap(X = a.placenta@H.norm ,n_components = 2,n_threads = 4)
# # dimred.afterharmony <- uwot::umap(X = V_harmony,n_components = 2,n_threads = 4)

# # plot(dimred.beforeharmony,pch=19,cex=0.1)
# # plot(dimred.afterharmony,pch=19,cex=0.1)

# dimred.beforeharmony <- uwot::umap(X = a.placenta@H.norm ,n_components = 2,n_threads = 4,min_dist=0.01, spread=1)#,  min.prob.lower=1e-3)
# dimred.afterharmony <- uwot::umap(X = V_harmony,n_components = 2,n_threads = 4,min_dist=0.01, spread=1)

# options(repr.plot.width=5,repr.plot.height=5)
# plot(dimred.beforeharmony,pch=19,cex=0.1)
# plot(dimred.afterharmony,pch=19,cex=0.1)
# ##############



##do clustering
a.placenta <- louvainCluster(a.placenta, resolution = 0.8) #will rewrite the cluster slot
names(a.placenta@clusters) <- rownames(a.placenta@H.norm)

table(cluster.df.add.bk$cluster) #1
  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
3628 3157 2805 2500 2426 2339 2226 1335  992  529  405  368  321  288  176  139 
  16   17 
  71   27 

table(a.placenta@clusters) #0.8, use this
   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
4649 3267 2844 2672 2472 2453 2210  993  518  405  368  321  288  175   70   27

table(a.placenta@clusters) #0.6
   0    1    2    3    4    5    6    7    8    9   10   11   12   13 
4763 4690 3889 2817 2295 2140 1063  493  405  365  321  288  176   27


######do dimension reduction and clustering by liger function####
#tSNE
# a.placenta <- runTSNE(a.placenta)
# plots <- plotByDatasetAndCluster(a.placenta,pt.size = .1,text.size = 8,axis.labels=c('tSNE1','tSNE2'),return.plots = T) 

#UMAP with tuning
#a.placenta <- runUMAP(a.placenta)
##a.placenta <- runUMAP(a.placenta,distance = 'cosine', n_neighbors = 30, min_dist = 0.2,rand.seed = 10)


plots <- plotByDatasetAndCluster(a.placenta,pt.size = .1,text.size = 8,axis.labels=c('UMAP1','UMAP2'),return.plots = T) 

##mixability of ATAC two donor and RNA two donor
options(repr.plot.width=6.5,repr.plot.height=5)
plots[[1]] +   
  scale_color_manual(values = c('atac1'='#94C6DD','atac2'='#1273AE','rna1'='#F3A585','rna2'='#C80927') ) + #'#94C6DD', '#1273AE'  vs '#F3A585','#C80927'
  theme(legend.position = 'right',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines"),
        legend.text = element_text(size  = 15),
        legend.title = element_text(size  = 15)
       )

##cluster id on umap
options(repr.plot.width=6.5,repr.plot.height=5)
plots[[2]] +   
  scale_color_manual(values = color_good) + 
  theme(legend.position = 'right',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines"),
        legend.text = element_text(size  = 15),
        legend.title = element_text(size  = 15)
       )




######customized dimension reduction with uwot::umap and RANN build graph then clustering with leidenlg#####

#all.equal(a.placenta@H.norm,seurat_obj@reductions$inmf@cell.embeddings,check.attributes = FALSE)#TRUE
#matINMF <- seurat_obj@reductions$inmf@cell.embeddings

matINMF <- a.placenta@H.norm

##tuning umap

set.seed(10)
dimred.umap <- uwot::umap(X = a.placenta@H.norm ,n_components = 2,n_threads = 2,min_dist=0.2, spread=0.5)

options(repr.plot.width=15,repr.plot.height=15)
plot(dimred.umap,pch=19,cex = 0.2)




# ####tuning cluster
# getCluster <- function(data.mat = , mat.use = matINMF ){ 
#     ##input a cell x dim matrix (the first round of dr, as pca, lsi, diffusion map,NMF...)
#     ##get knn graph (getGraph) then clustering with tuning (leidenlg)
#     ##quick plot the cluster.df result and output
    

#     ##get snn graph and adjacency matrix
#     #https://github.com/TomKellyGenetics/leiden
#     snn <- RANN::nn2(mat.use,k=30)$nn.idx
#     adjacency_matrix <- matrix(0,nrow(mat.use),nrow(mat.use) )
#     rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- rownames(mat.use)
#     for(ii in 1:nrow(mat.use)) {
#         adjacency_matrix[ii,rownames(data_mat)[snn[ii,]]] <- 1L
#     }
#     #check that rows add to k
#     sum(adjacency_matrix[1,]) == 30
#     table(apply(adjacency_matrix, 1, sum))

#     #do leiden
#     partition <- leiden(adjacency_matrix)
# }
# #############


# ###modify a.placenta tsne.coord and cluster slots ?




########save obj####
##saveRDS(a.placenta,'a.placenta.rds')
#a.placenta <- readRDS('a.placenta.rds')
saveRDS(a.placenta@raw.data,'a.placenta.raw.data.rds')



############get the  cluster.df.add meta table############
coords.int = a.placenta@tsne.coords #23732 #16877

clusters.int = a.placenta@clusters #23732 #16877
table(a.placenta@clusters)

   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
4649 3267 2844 2672 2472 2453 2210  993  518  405  368  321  288  175   70   27 

#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 3628 3157 2805 2500 2426 2339 2226 1335  992  529  405  368  321  288  176  139 
#   16   17 
#   71   27

#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 3428 2876 2707 2627 2078 2015 1763 1582 1065  916  768  470  402  332  285  231 
#   16   17 
#  158   29 


#  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 4247 3924 3756 2620 2399 2146  909  726  528  523  462  378  289  259  163  162 
#   16   17   18   19 
#  123  116    1    1 

#     0     1     2     3     4     5     6     7     8     9    10    11 
# 15735  4644   908   765   436   349   310   259   162   162     1     1 

levels(clusters.int) = as.numeric(gtools::mixedsort(levels(clusters.int)) ) + 1
table(clusters.int)
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
4649 3267 2844 2672 2472 2453 2210  993  518  405  368  321  288  175   70   27 

#for ( i in seq_len(length(clusters.int)) ){clusters.int[i]<- clusters.int[i]+1)}

all.equal(rownames(a.placenta@H.norm),rownames(a.placenta@tsne.coords)) #TRUE
all.equal (rownames(a.placenta@tsne.coords),names(a.placenta@clusters) )#TRUE
all.equal(names(a.placenta@clusters)  ,  c(names(atac_cluster),names(rna_cluster)) ) #TRUE


# alignment.clusters = a.placenta@alignment.clusters #23732 #16877
# #levels(alignment.clusters) = 0:14
# levels(alignment.clusters) = gtools::mixedsort(levels(alignment.clusters ))
# all.equal (a.placenta@clusters, a.placenta@alignment.clusters ) ##TRUE, the same

#combine predefined clusters
cluster.all = rbind (data.frame(cluster = atac_cluster, type = 'atac'  ), data.frame(cluster = rna_cluster, type = 'rna'  ) )
#cluster.all = rbind (data.frame(cluster = rna_cluster, type = 'rna'  ),data.frame(cluster = atac_cluster, type = 'atac'  ) )
all.equal(rownames(cluster.all),rownames(coords.int) ) #TRUE
all.equal(rownames(cluster.all),names(clusters.int) ) #TRUE
#all.equal(rownames(cluster.all),names(alignment.clusters) ) #TRUE

#cluster.all.coord = cbind(cluster.all,coords.int,"clusters_int" = clusters.int, "alignment_clusters" = alignment.clusters)
cluster.all.coord = cbind(cluster.all,coords.int,"clusters_int" = clusters.int)

head(cluster.all.coord)
#colnames(cluster.all.coord) = c("cluster","type",'UMAP1','UMAP2',"clusters_int","alignment_clusters")
colnames(cluster.all.coord) = c("cluster","type",'UMAP1','UMAP2',"clusters_int")
#all.equal (cluster.all.coord$clusters_int, cluster.all.coord$alignment_clusters ) #TRUE, can omit alignment_clusters

table(is.na(cluster.all.coord$cluster)) #FALSE


#####add donor info####
cluster.all.coord$anno <- 'none'

idx1 <- which(cluster.all.coord$type == 'rna' & grepl(pattern='-1$',x=rownames(cluster.all.coord)) )
idx2 <- which( cluster.all.coord$type == 'rna' & grepl(pattern='-2$',x=rownames(cluster.all.coord)) )
idx3 <- which( cluster.all.coord$type == 'atac' & grepl(pattern='-1$',x=rownames(cluster.all.coord)) )
idx4 <- which(cluster.all.coord$type == 'atac' & grepl(pattern='-2$',x=rownames(cluster.all.coord)) )

intersect(idx1,idx2) #0
intersect(idx3,idx4) #0
intersect(idx1,idx3) #0

length(idx1) + length(idx2) + length(idx3) + length(idx4) == nrow(cluster.all.coord) #TRUE


cluster.all.coord[idx1,'anno'] <- 'rna_D1'
cluster.all.coord[idx2,'anno'] <- 'rna_D2'
cluster.all.coord[idx3,'anno'] <- 'atac_D1'             
cluster.all.coord[idx4,'anno'] <- 'atac_D2'

cluster.all.coord$anno <- factor(cluster.all.coord$anno,levels=c('rna_D1','rna_D2','atac_D1','atac_D2' ) )

###save cluster.df
cluster.df.add <- cluster.all.coord[,c('clusters_int','UMAP1','UMAP2','cluster','type','anno')]
colnames(cluster.df.add) <- c('cluster','UMAP_1','UMAP_2','cluster_lib','type','anno')


#######filter cluster.df.add with te only
cluster.all.te = rbind (data.frame(cluster = atac_cluster.te, type = 'atac'  ), data.frame(cluster = rna_cluster.te, type = 'rna'  ) ) #21291

table(row.names(cluster.all.te) %in% row.names(cluster.all) )
TRUE 
21291 

all.equal( cluster.all[row.names(cluster.all.te),2],cluster.all.te[,2], check.attributes = FALSE )
all.equal( droplevels(cluster.all[row.names(cluster.all.te),1]),cluster.all.te[,1], check.attributes = FALSE )
table(cluster.all[row.names(cluster.all.te),1])
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
3552 3211 2315 2639 2558 2471 1911 1296 1050    0  288    0    0    0    0    0 

table(cluster.all.te[,1])
 1    2    3    4    5    6    7    8    9   10 
3552 3211 2315 2639 2558 2471 1911 1296 1050  288 

table(rownames(cluster.all.te) %in% rownames(cluster.df.add) )
 TRUE 
21291 

cluster.df.add.te <- cluster.df.add[rownames(cluster.all.te),] #21291

table(cluster.df.add.te$cluster)
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
4631 3252 2829 2651 2443 2429 2168    2  516    8  349    6    2    5    0    0 

table(cluster.df.add.te$cluster_lib)
  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
3552 3211 2315 2639 2558 2471 1911 1296 1050    0  288    0    0    0    0    0 


##rename and drop level for cluster_lib

cluster_id <- cluster.df.add.te$cluster_lib
#cluster_id <- droplevels(cluster_id)
table(cluster_id)
1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
3552 3211 2315 2639 2558 2471 1911 1296 1050    0  288    0    0    0    0    0 

for(i in seq_len(length(cluster_id)) ){ if(cluster_id[i] == 11){cluster_id[i] = 10}   }
table(cluster_id)
1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
3552 3211 2315 2639 2558 2471 1911 1296 1050  288    0    0    0    0    0    0 

cluster_id <- droplevels(cluster_id)
table(cluster_id)
   1    2    3    4    5    6    7    8    9   10 
3552 3211 2315 2639 2558 2471 1911 1296 1050  288 

cluster.df.add.te$cluster_lib <- cluster_id

##saveRDS(cluster.all.coord,"cluster.all.coord.rds")
saveRDS(cluster.df.add,'cluster.df.add.rds') #use this
saveRDS(cluster.df.add.te,'cluster.df.add.te.rds') #use this

#plot(readRDS('cluster.all.coord.rds.bk')[,3:4])
#plot(readRDS('old_rds/cluster.all.coord.rds')[,3:4])

############

cluster.df.add.bk <- cluster.df.add
cluster.df.add <- cluster.df.add.te



##subset a.placenta liger object

a.placenta.te <- subsetLiger(a.placenta,cells.use = rownames(cluster.all.te) ) #21291
table(cluster.df.add.te$anno)
 rna_D1  rna_D2 atac_D1 atac_D2 
   5982    4216    4557    6536 
##checked ok

saveRDS(a.placenta.te , 'a.placenta.te.rds')


############start to customized plot################

UMAP = cluster.df.add[,c('UMAP_1','UMAP_2')]

left <- 1.3*min(UMAP[,1])
right <- 0.7*max(UMAP[,1])
top <- 1*max(UMAP[,2])
bottom <- 1*min(UMAP[,2])


#######visualize integrated library#########
options(repr.plot.width=7.5,repr.plot.height=7.5)
#ggplot(data=cluster.all.coord,aes(x=UMAP1,y=UMAP2,col=type) ) +
ggplot(data=cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=anno) ) +
#ggplot(data=cluster.all.coord,aes(x=UMAP1,y=UMAP2,col=anno) ) +
  geom_point(size=0.2,alpha=0.9) +
  scale_color_manual(values = c('rna_D1'='#94C6DD', 'rna_D2'='#1273AE'  , 'atac_D1'='pink','atac_D2'='#C80927'),name='8w library') +
  xlim(left,right) + ylim(bottom,top) +
  #theme_classic() +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 18, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines"),
        legend.text = element_text(size  = 15),
        legend.title = element_text(size  = 15)
       ) +
  guides(colour=guide_legend(override.aes=list(size=8)))
ggsave(filename = "pdfs_te/PLA-8w-RNA-ATAC-liger.source_of_donor.pdf",height=7.5,width=7.5,useDingbats=FALSE)

# ggplot(data=cluster.all.coord,aes(x=dim1,y=dim2,col=clusters_int) ) +
#   geom_point(size=0.5) +
#   scale_color_manual(values = color_good) +
#   theme_classic() +
#   guides(colour=guide_legend(override.aes=list(size=6)))


# ggplot(data=cluster.all.coord,aes(x=dim1,y=dim2,col=clusters_int) ) +
#   geom_point(size=0.5) +
#   scale_color_manual(values = color_good) +
#   theme_classic() +
#   guides(colour=guide_legend(override.aes=list(size=6)))


######plot the integrated clusters
options(repr.plot.width=7.5,repr.plot.height=7.5)
ggplot(data=cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster) ) +
#ggplot(data=cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster_lib) ) +
#ggplot(data=cluster.all.coord,aes(x=UMAP1,y=UMAP2,col=clusters_int) ) +
  geom_point(size=0.2) +
  scale_color_manual(values = color_good) +
  xlim(left,right) + ylim(bottom,top) +
  #theme_classic() +
  theme(legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 25, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines"),
        legend.text = element_text(size  = 15),
        legend.title = element_text(size  = 15)
       ) +
  guides(colour=guide_legend(override.aes=list(size=8)))
ggsave(filename = "pdfs_te/PLA-8w-RNA-ATAC-liger.umap.nolabel.pdf",height=7.5,width=7.5,useDingbats=FALSE)




############customized way to plot umap-cluster with text halo###########
centers <- cluster.df.add %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))

centers_shift = data.frame()
##add little shift for text x and y, to plot text halo, borrow from snapATAC
theta= seq(0, 2*pi, length.out=50)
r=0.1
strwidth = 0.5 
strheight = 0.5
xo <- r*strwidth # r*strwidth('A')
yo <- r*strheight #r*strheight('A')
for (i in seq_len(nrow(centers))){
  for (j in theta) {
        centers_shift = rbind(centers_shift,
                              data.frame(
                                  cluster=as.character(unlist(centers[i,'cluster'])),
                                  x=centers[i,'x'] + cos(j)*xo, 
                                  y=centers[i,'y'] + sin(j)*yo
                                 )
                       )
      }
}

####the UMAP plot with annotation
#label right
options(repr.plot.height=5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
  #scale_colour_manual(values = color_snap_mod1)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  xlim(left,right) + ylim(bottom,top) +
  theme(
        legend.position = 'right',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "\ntotal cells:",nrow(cluster.df.add),  sep=" ") ) +
#   geom_text(data = centers_shift, #the halo
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "white", 
#             size = 6.5) +
#   geom_text(data = centers, 
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "black", 
#             size = 6) +
  guides(col = guide_legend(override.aes = list(size = 3))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs_te/PLA-8w-RNA-ATAC-liger.UMAP.labelright.pdf",height=5,width=6,useDingbats=FALSE)

##label on cluster

cluster.df.add.filter <- cluster.df.add

table(cluster.df.add.filter$cluster)

   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
4631 3252 2829 2651 2443 2429 2168    2  516    8  349    6    2    5    0    0 

cluster.df.add.filter <- subset(cluster.df.add,cluster %in% c('1','2','3','4','5','6','7','9','11'))

nrow(cluster.df.add.filter)
#21268

table(cluster.df.add.filter$cluster)
 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
4631 3252 2829 2651 2443 2429 2168    0  516    0  349    0    0    0    0    0 



#rename them
cluster.int <- as.character(droplevels(cluster.df.add.filter$cluster))
table(cluster.int)
#  1    2    3    4    5    6    7    9   11 
# 4631 3252 2829 2651 2443 2429 2168  516  349
  1   11    2    3    4    5    6    7    9 
4631  349 3252 2829 2651 2443 2429 2168  516

map_cellname <- list(
    '1' = '1',
    '2' = '2',
    '3' = '3',
    '4' = '4',
    '5' = '5',
    '6' = '6',
    '7' = '7',
    '9' = '8',
    #'10' = NA,
    '11' = '9'
    #'12' = NA,
    #'13' = NA,
    #'14' = NA


)


for (i in seq_len(length(cluster.int)) ){   
    cluster.int[i] <- map_cellname[[cluster.int[i]]]
}

table(cluster.int)
1    2    3    4    5    6    7    8    9 
4631 3252 2829 2651 2443 2429 2168  516  349 


cluster.df.add.filter$cluster <- cluster.int #21268
table(cluster.df.add.filter$cluster)





##subset a.placenta liger object again

a.placenta.te <- subsetLiger(a.placenta.te,cells.use = rownames(cluster.df.add.filter) ) #21268
# Warning message in FUN(X[[i]], ...):
# “Number of subsetted cells too small (less than 25), please check cells.use!”
# Warning message in FUN(X[[i]], ...):
# “Number of subsetted cells too small (less than 25), please check cells.use!”
# Warning message in FUN(X[[i]], ...):
# “Number of subsetted cells too small (less than 25), please check cells.use!”
# Warning message in FUN(X[[i]], ...):
# “Number of subsetted cells too small (less than 25), please check cells.use!”
# [1] "Removing 2 genes not expressing in atac2."
# [1] "AC131056.6" "DEFB103A" 

table(cluster.df.add.filter$anno)

 rna_D1  rna_D2 atac_D1 atac_D2 
   5978    4216    4549    6525

#cluster.df.add.te.bk <- cluster.df.add.te
cluster.df.add.te <- cluster.df.add.filter

cluster.df.add.te <- cluster.df.add.te[,-8]



####rename int cluster id########
cluster.int.te <- as.character(droplevels(a.placenta.te@clusters))
names(cluster.int.te) <- names(a.placenta.te@clusters)
table(cluster.int.te)
 0    1   10    2    3    4    5    6    8 
4631 3252  349 2829 2651 2443 2429 2168  516 

map_cellname_te <- list(
    '0' = '1',
    '1' = '2',
    '2' = '3',
    '3' = '4',
    '4' = '5',
    '5' = '6',
    '6' = '7',
    '8' = '8',
    '10' = '9'
    #'11' = '9'
    #'12' = NA,
    #'13' = NA,
    #'14' = NA


)


all.equal(names(cluster.int.te), rownames(cluster.df.add.te), check.attributes = FALSE ) #TRUE
table(cluster.df.add.te$cluster)
  1    2    3    4    5    6    7    8    9 
4631 3252 2829 2651 2443 2429 2168  516  349 


for (i in seq_len(length(cluster.int.te)) ){   
    cluster.int.te[i] <- map_cellname_te[[cluster.int.te[i]]]
}

table(cluster.int.te)
  1    2    3    4    5    6    7    8    9 
4631 3252 2829 2651 2443 2429 2168  516  349 

all.equal(cluster.int.te,cluster.df.add.te$cluster,check.attributes = FALSE)#TRUE


cluster.int.te <- factor(cluster.int.te)


a.placenta.te@clusters <- cluster.int.te


saveRDS(a.placenta.te , 'a.placenta.te.rds')
saveRDS(cluster.df.add.te , 'cluster.df.add.te.rds')


###centers
centers <- cluster.df.add.te %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))

centers_shift = data.frame()
##add little shift for text x and y, to plot text halo, borrow from snapATAC
theta= seq(0, 2*pi, length.out=50)
r=0.1
strwidth = 0.5 
strheight = 0.5
xo <- r*strwidth # r*strwidth('A')
yo <- r*strheight #r*strheight('A')
for (i in seq_len(nrow(centers))){
  for (j in theta) {
        centers_shift = rbind(centers_shift,
                              data.frame(
                                  cluster=as.character(unlist(centers[i,'cluster'])),
                                  x=centers[i,'x'] + cos(j)*xo, 
                                  y=centers[i,'y'] + sin(j)*yo
                                 )
                       )
      }
}

###

options(repr.plot.height=5,repr.plot.width=5)
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
ggplot(cluster.df.add.te,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = 0.2,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
  xlim(left,right) + ylim(bottom,top) +
  #scale_colour_manual(values = color_snap_mod1)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  theme(
        legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "\ntotal cells:",nrow(cluster.df.add.te),  sep=" ") ) +
  geom_text(data = centers_shift, #the halo
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "white", 
            size = 6.5) +
  geom_text(data = centers, 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            size = 6.5) +
  ##guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs_te/PLA-8w-RNA-ATAC-liger.UMAP.labelon.pdf",height=5,width=5,useDingbats=FALSE)
ggsave(filename = "pdfs_te/PLA-8w-RNA-ATAC-liger.UMAP.labelon.final.pdf",height=5,width=5,useDingbats=FALSE)
###############



##########density 2d with contour plot UMAP#########
ggplot(cluster.df.add.te, aes(UMAP_1, UMAP_2)) + 
  geom_point(color = "lightgray") +
  geom_density_2d(color='red') +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_gradientn(colors = c("#FFEDA0", "#FEB24C", "#F03B20")) +
  #scale_fill_distiller(palette = "Spectral", direction = -1)+
  theme(
        legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       )
#########



#########pixel style hexbin density scatterplot (grid)#######
##https://cran.r-project.org/web/packages/hexbin/vignettes/hexagon_binning.pdf
pdf( "pdfs_te/pixel_style.cluster.density.pdf",height=7.5,width=7.5,useDingbats = FALSE)

options(repr.plot.height=7.5,repr.plot.width=7.5)

bin<-hexbin(cluster.df.add.te$UMAP_1, cluster.df.add.te$UMAP_2, xbins=40,xbnds = range(cluster.df.add.te$UMAP_1)*1.3,ybnds = range(cluster.df.add.te$UMAP_2)*1.3 )

my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))

hb <- plot(bin,main="Nuclei pileup density" , colramp=my_colors , legend=F,xlab='UMAP1',ylab='UMAP2')

dev.off()

#bin@xbnds <- c(left,right)
#bin@ybnds <- c(bottom,top)

#pdf( "pdfs_te/pixel_style.cluster.density.pdf",height=7.5,width=7.5,useDingbats = FALSE)
# options(repr.plot.height=7.5,repr.plot.width=7.5)

# grid.newpage()
# grid.show.viewport( viewport(x=0.6, y=0.6,
#                              w=unit(1, "inches"), 
#                              h=unit(1, "inches"),
#                              xscale=c(left,right)+.05*c(-1,1),
#                              yscale=c(bottom,top)+.05*c(-1,1),
#                             )  
#                   )



# #pushViewport(viewport(layout = grid.layout(1, 1)))
# #pushViewport(viewport(layout.pos.col = 1,layout.pos.row = 1))

# hb <- plot(bin,main="hexagon pileup density" , colramp=my_colors , legend=F,xlab='UMAP1',ylab='UMAP2')#,xlim=c(left,right),ylim=c(bottom,top) ) 
# ##grid.points(c(left,right), c(bottom,top) ,pch = 16,gp = gpar(cex = .4))
# ##hb$plot.vp@xscale <- c(-6,15)
# ##hb$plot.vp@yscale <- c(-10,10)

# pushHexport(hb$plot.vp,xscale=c(-6,15),yscale=c(-10,10))
# popViewport(hb$plot.vp,xscale=c(-6,15),yscale=c(-10,10))


# #dev.off()

# pushViewport(dataViewport(c(-6,15), c(-10,10), 
#                layout.pos.row=1, layout.pos.col=1, name="11"))


#######



##################################post process#####################################
##https://macoskolab.github.io/liger/liger-vignette.html

##feature (meta table)
plotFeature(a.placenta.te,'nUMI')

##work cloud
word_clouds <- plotWordClouds(a.placenta.te, num.genes = 10, do.spec.plot = F, return.plots = T)
print(word_clouds[[11]])

########explore iNMF factorization gene loading (gene signature )########
##metagenes–sets of co-expressed genes###


gene_loadings_atac12 <- plotGeneLoadings(a.placenta.te, num.genes = 15, do.spec.plot = F, return.plots = T,dataset1 = 'atac1',dataset2='atac2' )
#35 factors (k=35)
#20 factors (k = 20)

gene_loadings_atac1_rna1 <- plotGeneLoadings(a.placenta.te, num.genes = 15, do.spec.plot = F,return.plots = T, dataset1 = 'atac1',dataset2='rna1' )

gene_loadings_atac2_rna2 <- plotGeneLoadings(a.placenta.te, num.genes = 15, do.spec.plot = F, return.plots = T,dataset1 = 'atac2',dataset2='rna2' )

gene_loadings_rna12 <- plotGeneLoadings(a.placenta.te, num.genes = 15, do.spec.plot = F, return.plots = T,dataset1 = 'rna1',dataset2='rna2' )


# print(gene_loadings[[4]]) #EVT
# print(gene_loadings[[5]])
# print(gene_loadings[[9]])
# print(gene_loadings[[10]]) #?
# print(gene_loadings[[12]])
# print(gene_loadings[[13]])
# print(gene_loadings[[17]])
# print(gene_loadings[[18]])

# print(gene_loadings[[1]])
# print(gene_loadings[[3]])

options(repr.plot.height=7.5,repr.plot.width=7.5)
#for(i in 1:10){
#for(i in 11:20){ 
for(i in 1:20){ 
#for(i in 21:30){ 
#for(i in 31:35){ 
#for(i in 11:length(gene_loadings)){
    #print(gene_loadings_atac12[[i]]) 
    print(gene_loadings_rna12[[i]])
    #print(gene_loadings_atac1_rna1[[i]]) 
    #print(gene_loadings_atac2_rna2[[i]])  
}

##arrange plot by patchwork
options(repr.plot.height=15,repr.plot.width=15)
gene_loadings_atac12[[2]] + gene_loadings_rna12[[2]] + gene_loadings_atac1_rna1[[2]] + gene_loadings_atac2_rna2[[2]] +
plot_layout(ncol=2,nrow=3)

options(repr.plot.height=15,repr.plot.width=15)
gene_loadings_atac12[[2]] + gene_loadings_rna12[[2]]  +
plot_layout(ncol=2,nrow=2)

 

#######visualization for iNMF factor marker gene#####
markers_atac12 = getFactorMarkers(a.placenta.te,num.genes = 15,dataset1 = 'atac1',dataset2='atac2')
markers_rna12 = getFactorMarkers(a.placenta.te,num.genes = 15,dataset1 = 'rna1',dataset2='rna2')
markers_atac1_rna1 = getFactorMarkers(a.placenta.te,num.genes = 15,dataset1 = 'atac1',dataset2='rna1')
markers_atac2_rna2 = getFactorMarkers(a.placenta.te,num.genes = 15,dataset1 = 'atac2',dataset2='rna2')

#table(markers_atac12[['shared']][,'factor_num']) #15 genes for each 35 factors

marker_gene = c('HLA-G','DNMT1','CDH1','FLT1','PAPPA','ERVFRD-1','LAIR2','VIM')
geneid = marker_gene[8]

##plot gene enrichment heatmap-umap
gene_plot = plotGene(a.placenta.te,gene=geneid,return.plots=T)

options(repr.plot.height=15,repr.plot.width=15)
gene_plot[['atac1']] + gene_plot[['atac2']]+gene_plot[['rna1']]+gene_plot[['rna2']] +
plot_layout(ncol=2,nrow=2) #patch_work


##plot gene violin among clusters
gene_plot = plotGeneViolin(a.placenta.te, gene = geneid,return.plots=T)
names(gene_plot) <- c('atac1','atac2','rna1','rna2')

#options(repr.plot.height=15,repr.plot.width=15)
gene_plot[['atac1']] + gene_plot[['atac2']]+gene_plot[['rna1']]+gene_plot[['rna2']] +
plot_layout(ncol=2,nrow=2) #patch_work





####identification of joint-cluster specific genes
dge.wilcoxon <- runWilcoxon(a.placenta.te, data.use = 'all', compare.method = 'clusters') #quick

saveRDS(object = dge.wilcoxon,file='dge.wilcoxon.te.rds')
   
##visualize for  wilcoxon differential genes with heatmap?





#######exploring factors and clusters####
#factors
plotFactors(a.placenta.te) #20 factors, 11,13,14,18,20 ##do not plot directly

#cluster proportions
plotClusterProportions(a.placenta.te)


#cluster factors
plotClusterFactors(a.placenta.te,use.aligned=T)


##comparing different cluster assignment
## calculate adjusted Rand Index between liger cluster assignments and another assignment
calcARI(a.placenta.te, c(rna_cluster.te, atac_cluster.te)) #0.21662
# calculate purity between liger cluster assignments and another assignment for just dataset 1 
calcPurity(a.placenta.te, rna_cluster.te) #0.605
calcPurity(a.placenta.te, atac_cluster.te)  #0.498

calcAlignment(a.placenta.te) #0.82169
calcAgreement(a.placenta.te) #0.0812259498741353
# see if certain clusters are more integrated than others
calcAlignmentPerCluster(a.placenta.te) 
#
                             
                                 


############liger to seurat obj and do annotation and visualization with seurat functions########
# Create Seurat object from liger object, keeping liger highly variable genes
seurat_obj = ligerToSeurat(a.placenta.te, use.liger.genes = T)
saveRDS(seurat_obj,'seurat_obj_fromliger.te.rds')

# DefaultAssay(seurat_obj) <- "RNA"


# #plot(seurat_obj@reductions$tsne@cell.embeddings,pch=19,cex=0.1)


# DimPlot(seurat_obj, reduction = "tsne", #UMAP in fact
#               #group.by = "Eday", 
#               #cols = color_set3_ext17,
#               label.size = 10,
#               label = TRUE,
#               #repel = TRUE, 
#               pt.size = 0.1) +

#               guides(colour = guide_legend(override.aes = list(size=8)))
        


# gene_list1 = c("FLT1","GCM1","CGA","PAGE4","DNMT1",'CST6') 
# gene_list2 = c("VIM",'LAIR2','HLAG','MMP9','PLAC8',"ERVFRD-1")
# gene_list3 = c("PECAM1","CD68","HLA-A",'KRT7','MKI67','PCNA')


# genes_stb1 = c('INSL4', 'SDC1', 'INHA')
# genes_stb2 = c("FLT1","GCM1","CGA",'KISS1', 'CSHL1')


# marker.genes = c(
#     "DNMT1", "CDH1", "MKI67",
#     "FLT1", "CSHL1", "PSG8", 
#     "ERVFRD-1", "LAIR2", "PLAC8",
#     'VIM','PECAM1','CD14'
#   );

# marker.genes = c( #try to distinguish STB terminals 
#     "FLT1", 'LEP','INSIG2', #FLT1-enriched
#     "CSHL1",'CSH1','PAPPA',  #PAPPA-enriched
#     'PSG1','CGA','GCM1' #general

#   );

# ##STR markers
# marker.genes = c(
#     "VIM", "PECAM1", "THY1",
#     "NR2F1", "DLK1", "HIVEP3", 
#     "CD68", "CD14", "HLA-A",
#     'HLA-DPA1','HLA-DPB1','MKI67'
#   );



# ##check gene name 

# marker.genes %in% rownames(seurat_obj) #exact match
# grep(pattern = 'CGB',x = rownames(seurat_obj),value = TRUE)


# seurat_obj <- ScaleData(seurat_obj, verbose = FALSE,features = rownames(seurat_obj@assays$RNA@data) )

# options(repr.plot.width=12,repr.plot.height=18)
# FeaturePlot(object = seurat_obj, features = marker.genes, cols = c("grey", "blue"), 
#     reduction = "tsne",slot='data',ncol=2)                         


# ##################






########compare cicero aggragated gene peak-coaccessibility activity score ,######
########gmat and promoter_genebody gene count activity score#####


##1 cicero gene activity score of aggre peaks (cellranger-atac gtf, gencode v28, ensembl v92)####
ciceroGA <- readRDS("../../02.snapATAC_harmony/cicero_Granja/ciceroGA.rds")
#19529 x 12526


##rename to -1, -2
cellid <- colnames(ciceroGA) 
idy <- which(grepl(pattern = "^placenta_donor2",x = cellid    ))
cellid[idy] <- gsub( pattern = '-1$',replacement = '-2' ,x =  cellid [idy]  )
cellid  <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = cellid )


##check core barcode
check.df <- data.frame(colid = colnames(ciceroGA), cellid = cellid   )
check.df$colid <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = check.df$colid )
check.df$colid <- gsub(pattern = '-1$',replacement = "", x = check.df$colid )
check.df$cellid <- gsub(pattern = '-1|-2$',replacement = "", x = check.df$cellid )
all.equal(check.df$colid,check.df$cellid) #TRUE, the core barcode 

#add prefix
colnames(ciceroGA)  <- paste0('placenta.atac_',cellid)



##2 SnapATAC processed gmat (cellranger-atac gtf, gencode v28, ensembl v92)#######
gmat <- t(readRDS('gmat.rds')) #SnapATAC, with RPM depth normalize and smooth by magic 
#58344 x 12526
#gmat_head100 <- gmat[1:100,1:100]

#all.equal(gmat[1:100,1:100],gmat_head100,check.attributes = FALSE) #TRUE


##rename to -1, -2
cellid <- colnames(gmat) 
idy <- which(grepl(pattern = "^placenta_donor2",x = cellid    ))
cellid[idy] <- gsub( pattern = '-1$',replacement = '-2' ,x =  cellid [idy]  )
cellid  <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = cellid )

##check core barcode
check.df <- data.frame(colid = colnames(gmat), cellid = cellid   )
check.df$colid <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = check.df$colid )
check.df$colid <- gsub(pattern = '-1$',replacement = "", x = check.df$colid )
check.df$cellid <- gsub(pattern = '-1|-2$',replacement = "", x = check.df$cellid )
all.equal(check.df$colid,check.df$cellid) #TRUE, the core barcode 

#add prefix
colnames(gmat)  <- paste0('placenta.atac_',cellid)


##3 liger bedmap count gene + promoter(upstream gene 2k) with EnsDb.Hsapiens.v86######

placenta.atac <- readRDS('placenta.atac.rds')


##3 seurat RNA scale.data

exprMat.z <- readRDS('../../02.seurat_harmony/exprMat.scaledata.rds') #24307 x 10198
#add prefix
colnames(exprMat.z)  <- paste0('placenta.rna_',colnames(exprMat.z))


table(colnames(exprMat.z) %in% colnames(placenta.rna) )
 TRUE 
10198 
table(colnames(placenta.rna) %in%  colnames(exprMat.z) )
 TRUE 
10198 


#######compare liger genePromoter activity score, SnapATAC gmat, cicero aggre_peaks gene activity score

##check cell id
all.equal(colnames(ciceroGA),colnames(gmat)) #TRUE
all.equal(colnames(gmat),colnames(placenta.atac)) #TRUE

all.equal (colnames(ciceroGA),c(colnames(a.placenta@raw.data[['atac1']]), colnames(a.placenta@raw.data[['atac2']])) ) #TRUE
all.equal (colnames(placenta.atac),c(colnames(a.placenta@raw.data[['atac1']]), colnames(a.placenta@raw.data[['atac2']])) ) #TRUE

##all cell id equal


##check gene id duplication and get shared gene set among three gene activity score mat
rowid_ga <- rownames(ciceroGA) #19529
rowid_gmat <- rownames(gmat) #58344
rowid_liger <- rownames(placenta.atac) #57840 #63970

sum(duplicated(rowid_ga)) #0
sum(duplicated(rowid_gmat)) #1549
sum(duplicated(rowid_liger)) #2332 #7327


#1 SnapATAC(no promoter region with rpm and magic smooth) vs liger (gene+2k promoter raw count)
rowid_gmat_liger <- intersect(rowid_gmat,rowid_liger) #47205 #47747

gmat_liger.df <- data.frame(gmat_rowave=Matrix::rowMeans(gmat[rowid_gmat_liger,]),
                            liger_rowave=Matrix::rowMeans(placenta.atac[rowid_gmat_liger,])
                           )
res.cor <- cor(gmat_liger.df[,1],gmat_liger.df[,2])#0.981 #0.66
gene_n <- nrow(gmat_liger.df) #47205 #47747

pdf( "pdfs_te/Correlation_of_SnapATAC_gmat-vs-liger_method.pdf",height=7.5,width=7.5,useDingbats = FALSE)

options(repr.plot.width=7.5,repr.plot.height=7.5)
plot(gmat_liger.df[,1],
     gmat_liger.df[,2],
     pch = 19,
     cex = 0.1,
     xlim =c(0,750),
     ylim = c(0,15),
     xlab = 'Gene activity (SnapATAC method)',
     ylab = 'Gene activity (liger method)',
     main = 'Correlation of different \ngene activity score method',
     cex.lab = 1.2,
     cex.main = 1.2
     
    )
text(x = 1,y=15,labels = paste("r = ",round(res.cor,digits = 2),sep=''),pos = 4 ,cex.lab=2)
text(x = 3,y=14.5,labels = paste('gene number: ',gene_n,sep=''),pos = 4, cex.lab=2 )

dev.off()

#2 ciceroGA (aggregated high coaccessible peaks to gene ) vs liger
rowid_ga_liger <- intersect(rowid_ga,rowid_liger) #14477

ga_liger.df <- data.frame(
       ga_rowave = Matrix::rowMeans(ciceroGA[rowid_ga_liger,]),
       liger_rowave = Matrix::rowMeans(placenta.atac[rowid_ga_liger,])
    )

res.cor <- cor(ga_liger.df[,1],ga_liger.df[,2]) #0.268 #0.12
gene_n <- nrow(ga_liger.df) #14472 #14477


pdf( "pdfs_te/Correlation_of_ciceroGA-vs-liger_method.pdf",height=7.5,width=7.5,useDingbats = FALSE)

options(repr.plot.width=7.5,repr.plot.height=7.5)
plot(ga_liger.df[,1],
     ga_liger.df[,2],
     pch = 19,
     cex = 0.1,
     xlim =c(0,4e-4),
     ylim = c(0,5),
     xlab = 'Gene activity (ciceroGA method)',
     ylab = 'Gene activity (liger method)',
     main = 'Correlation of different \ngene activity score method',
     cex.lab = 1.2,
     cex.main = 1.2
     
    )
text(x = 4e-4,y=5,labels = paste("r = ",round(res.cor,digits = 2),sep=''),pos = 2 )
text(x = 4e-4,y=4.8,labels = paste('gene number: ',gene_n,sep=''),pos = 2)

dev.off()



######



#########atac raw countd and rna raw count reproducibility between donors##########
#########get pseudobulk RNA correlation scatter plot and pseudobulk ATAC correlation scatter plot from raw.data of liger object### 

rowid_atac1 <- rownames(a.placenta.te@raw.data[['atac1']]) #57518 #30373
rowid_atac2 <- rownames(a.placenta.te@raw.data[['atac2']]) #57204 #30087
rowid_rna1 <- rownames(a.placenta.te@raw.data[['rna1']]) #24993
rowid_rna2 <- rownames(a.placenta.te@raw.data[['rna2']]) #25066



sum(duplicated(rowid_atac1)) #2314 #1218
sum(duplicated(rowid_atac2)) #2305 #1212

rowid_atac1 <- rowid_atac1[!duplicated(rowid_atac1)] #55204
rowid_atac2 <- rowid_atac2[!duplicated(rowid_atac2)] #54899

sum(duplicated(rowid_rna1)) #0
sum(duplicated(rowid_rna2)) #0


#atac
rowid_atac12 <- intersect(rowid_atac1 ,rowid_atac2 ) #54755 #shared gene 28801
combine_atac12 <-  union(rowid_atac1 ,rowid_atac2 ) #55348 #29229


atac12.df <- data.frame(atac1_rowave=Matrix::rowMeans(a.placenta.te@raw.data[['atac1']][rowid_atac12,]),
                         atac2_rowave=Matrix::rowMeans(a.placenta.te@raw.data[['atac2']][rowid_atac12,])
                        )
cor(atac12.df[,1],atac12.df[,2]) #0.985 #0.98
nrow(atac12.df) #54755 #28801

pdf( "pdfs_te/Reproducibility_of_snATAC_donor_library.pdf",height=7.5,width=7.5,useDingbats = FALSE)
options(repr.plot.width=7.5,repr.plot.height=7.5)
plot(atac12.df[,1],
     atac12.df[,2],
     pch = 19,
     cex = 0.5,
     xlim =c(0,20),
     ylim = c(0,8),
     xlab = "pseudobulk snATAC gene activity score (donor1)",
     ylab = "pseudobulk snATAC gene activity score (donor2)",
     main = 'Reproducibility of snATAC library',
     cex.main = 2,
     cex.lab = 1.2
    )
text(x = 1,y=8,labels = paste("r = ",round(cor(atac12.df[,1],atac12.df[,2]),digits = 2),sep='') )
text(x = 3,y=7.5,labels = paste('gene number: ',nrow(atac12.df),sep='') )
dev.off()

# atac12.df <- data.frame(atac1_rowave=Matrix::rowMeans(a.placenta.te@raw.data[['atac1']][rowid_atac12,]),
#                         atac2_rowave=Matrix::rowMeans(a.placenta.te@raw.data[['atac2']][rowid_atac12,])
#                        )
# ggplot(data = atac12.df, mapping = aes(x = atac1_rowave, y = atac2_rowave)) +
#   geom_point(size = 0.2) +
#   #geom_pointdensity() +
#   #stat_density2d(aes(fill = ..density..), contour = F, geom = 'polygon') +
#   #geom_bin2d(bins = 100, color ="white")+
#   #geom_hex(bins = 40, color = "white")+
#   #geom_density_2d() +
#   stat_density_2d(aes(fill = ..level..), geom = "polygon") +
#   #scale_color_viridis() +
#   #scale_fill_gradient(low =  "#00AFBB", high = "#FC4E07")+
#   theme_minimal()


# ggplot(diamonds, aes(carat, price)) +
#   geom_bin2d(bins = 20, color ="white")+
#   scale_fill_gradient(low =  "#00AFBB", high = "#FC4E07")+
#   theme_minimal()

#rna
rowid_rna12 <- intersect(rowid_rna1 ,rowid_rna2 ) #23673
rna12.df <- data.frame(rna1_rowave=Matrix::rowMeans(a.placenta.te@raw.data[['rna1']][rowid_rna12,]),
                         rna2_rowave=Matrix::rowMeans(a.placenta.te@raw.data[['rna2']][rowid_rna12,])
                        )
cor(rna12.df[,1],rna12.df[,2])#0.93
nrow(rna12.df) #23673

pdf( "pdfs_te/Reproducibility_of_snRNA_donor_library.pdf",height=7.5,width=7.5,useDingbats = FALSE)
options(repr.plot.width=7.5,repr.plot.height=7.5)
plot(rna12.df[,1],
     rna12.df[,2],
     pch = 19,
     cex = 0.5,
     xlim =c(0,25),
     ylim = c(0,50),
     xlab = "pseudobulk snRNA gene activity score (donor1)",
     ylab = "pseudobulk snRNA gene activity score (donor2)",
     main = 'Reproducibility of snRNA library',
     cex.main = 2,
     cex.lab = 1.2
    )
text(x = 1,y=50,labels = paste("r = ",round(cor(rna12.df[,1],rna12.df[,2]),digits = 2),sep='') )
text(x = 3,y=47,labels = paste('gene number: ',nrow(rna12.df),sep='') )
dev.off()


########

cellDecorate = function(cluster.df=cluster.df.add,data=gmat,gene='HLA-G',q.low=0,q.high=0.9,title='',color.use = color_ga){
    #prepare data df
    stopifnot(gene %in% rownames(data) )
    #sum(colnames(data) %in% rownames(cluster.df))
    
    feature.value <- data[gene,]
    cutoff <- quantile(feature.value,probs = c(q.low,q.high)) 
    cutoff.low <- cutoff[1]
    cutoff.high <- cutoff[2]
    feature.value[feature.value < cutoff.low] <- cutoff.low
    feature.value[feature.value > cutoff.high] <- cutoff.high
    
    cluster.df.data <- cbind(cluster.df[colnames(data),], value=feature.value )
    
    #plot
    p=ggplot(data=cluster.df.data,aes(x=UMAP_1,y=UMAP_2,col=value ) ) +
        geom_point(size=0.5,alpha=0.8) +
        #scale_color_manual(values = color_good) +
        #scale_color_gradientn(colours = rainbow(n = 7,rev = TRUE)) +
        scale_color_gradientn(colours = color.use) + 
        xlim(left,right) + ylim(bottom,top) +
        labs(title = title, subtitle="",caption = '',tag="",label="") +
        xlab ('UMAP1') +
        ylab ('UMAP2') +
        #scale_x_continuous(limits = c(-10,10)) +
        #scale_y_continuous(limits = c(-10,10)) +
        ggtitle(paste(gene,title,sep=' ') ) +
        #theme_classic() #+
        theme(
            legend.position = 'none',
            axis.text=element_blank(), 
            axis.title = element_text(size = 25, face = "bold"),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(color="black", fill = NA,size=1),
            #         panel.background = element_rect(fill = "white", colour = "white", 
            #                 size = rel(1)),
            #panel.border = element_blank(),
            plot.title = element_text(size = 25, face = "bold"),
            #complete = TRUE
            plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
        ) 
    #guides( col = guide_legend(override.aes = list(size = 5)) )
    p
}

#####compare three gene activity method

#marker_gene = c('DNMT1','ERVFRD-1','HLA-G', 'FLT1','PAPPA','PSG8')

#marker_gene = c('DNMT1','ERVFRD-1','MKI67', 'FLT1','PAPPA','PSG8')

marker_gene = c('LEP', 'SH3TC2', 'FLT1', 'PAPPA','PSG8', 'ERVFRD-1')

# gene_test = marker_gene[1]

# p1 <- cellDecorate(cluster.df=cluster.df.add,data=gmat,gene=gene_test,q.low=0.2,q.high=1,title='SnapATAC gmat')
# p2 <- cellDecorate(cluster.df=cluster.df.add,data=ciceroGA,gene=gene_test,q.low=0,q.high=1,title='ciceroGA')
# p3 <- cellDecorate(cluster.df=cluster.df.add,data=placenta.atac,gene=gene_test,q.low=0.3,q.high=1,title='liger method')

# options(repr.plot.width=20,repr.plot.height=7.5)
# p1 + p2 + p3 + plot_layout(ncol=3,nrow=1) #patch_work
# ####

######marker gene grid plot for liger integration UMAP
gene_test1 = marker_gene[1]
gene_test2 = marker_gene[2]
gene_test3 = marker_gene[3]

gene_test1 = marker_gene[4]
gene_test2 = marker_gene[5]
gene_test3 = marker_gene[6]

color_rna <- color_set_yellowbrick.flat[['Purples.5']]
color_rna <- color_list[['solarExtra']]

#'LEP'  atac 0.2-1, rna 0.2-1 
#'SH3TC2'  atac 0.2-1, rna 0.5-0.9 
#'FLT1',  atac 0.2-1, rna 0.2-1 

#'PAPPA', atac 0.2-1, rna 0.2-1 
#'PSG8',   atac 0.2-1, rna 0.2-1 
#'EFVFRD1'  atac 0.2-1, rna 0-1 


p1 <- cellDecorate(cluster.df=cluster.df.add,data=gmat,gene=gene_test1,q.low=0.2,q.high=1,title='Gene Activity Score')
p2 <- cellDecorate(cluster.df=cluster.df.add,data=gmat,gene=gene_test2,q.low=0.2,q.high=1,title='Gene Activity Score')
p3 <- cellDecorate(cluster.df=cluster.df.add,data=gmat,gene=gene_test3,q.low=0.2,q.high=1,title='Gene Activity Score')
#p4 <- cellDecorate(cluster.df=cluster.df.add,data=placenta.rna,gene=gene_test1,q.low=0,q.high=1,title='Gene Expression',color.use = color_rna)
#p5 <- cellDecorate(cluster.df=cluster.df.add,data=placenta.rna,gene=gene_test2,q.low=0,q.high=1,title='Gene Expression',color.use = color_rna)
#p6 <- cellDecorate(cluster.df=cluster.df.add,data=placenta.rna,gene=gene_test3,q.low=0,q.high=1,title='Gene Expression',color.use = color_rna)
p4 <- cellDecorate(cluster.df=cluster.df.add,data=exprMat.z,gene=gene_test1,q.low=0.2,q.high=1,title='Gene Expression',color.use = color_rna)
p5 <- cellDecorate(cluster.df=cluster.df.add,data=exprMat.z,gene=gene_test2,q.low=0.2,q.high=1,title='Gene Expression',color.use = color_rna)
p6 <- cellDecorate(cluster.df=cluster.df.add,data=exprMat.z,gene=gene_test3,q.low=0,q.high=1,title='Gene Expression',color.use = color_rna)

res.marker.list <- list()
res.marker.list[[paste('atac-',gene_test1,sep='')]] <- p1
res.marker.list[[paste('atac-',gene_test2,sep='')]] <- p2
res.marker.list[[paste('atac-',gene_test3,sep='')]] <- p3
res.marker.list[[paste('rna-',gene_test1,sep='')]] <- p4
res.marker.list[[paste('rna-',gene_test2,sep='')]] <- p5
res.marker.list[[paste('rna-',gene_test3,sep='')]] <- p6

for(i in names(res.marker.list)){
 ggsave(filename=paste("pdfs_te/integrated.marker.accessibility-vs-expression.",i,".pdf",sep=''),plot=res.marker.list[[i]],height=6,width=6,useDingbats = FALSE)
}



##combine to pdf
options(repr.plot.width=20,repr.plot.height=15)
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol=3,nrow=2) #patch_work
#ggsave("pdfs_te/integrated.marker.accessibility-vs-expression.pdf",height=20,width=15,useDingbats = FALSE)

#ggsave("pdfs_te/integrated.marker.accessibility-vs-expression.part1.pdf",height=20,width=15,useDingbats = FALSE)
ggsave("pdfs_te/integrated.marker.accessibility-vs-expression.part2.pdf",height=20,width=15,useDingbats = FALSE)

###


#############compariable riverplot (liger method to link atac cluster to rna cluster)##################
all.equal(names(a.placenta@clusters), c(names(atac_cluster),names(rna_cluster) )  ) #TRUE
all.equal(names(a.placenta.te@clusters), c(names(atac_cluster.te),names(rna_cluster.te) )  ) #FALSE

table(names(atac_cluster.te) %in% names(a.placenta.te@clusters))
FALSE  TRUE 
   19 11074 

flag_atac <- names(atac_cluster.te) %in% names(a.placenta.te@clusters)
atac_cluster.te <- atac_cluster.te[flag_atac]

table(names(rna_cluster.te) %in% names(a.placenta.te@clusters))
FALSE  TRUE 
    4 10194 

flag_rna <- names(rna_cluster.te) %in% names(a.placenta.te@clusters)
rna_cluster.te <- rna_cluster.te[flag_rna]

all.equal(names(a.placenta.te@clusters), c(names(atac_cluster.te),names(rna_cluster.te) )  ) #TRUE


table(a.placenta.te@clusters)

   1    2    3    4    5    6    7    8    9 
4631 3252 2829 2651 2443 2429 2168  516  349

###


#source("makeRiverPlot.r") #add 'rna- atac-'
#environment(makeRiverplot) <- asNamespace('liger')

options(repr.plot.width=12,repr.plot.height=15)
# makeRiverplot(a.placenta.te, rna_cluster.te, atac_cluster.te, min.frac = 0.2,river.yscale = 5, #node.order = set_node_order,
#               river.usr = c(0, 1, -0.6, 1.6))

# makeRiverplot(a.placenta.te, rna_cluster.te, atac_cluster.te, min.frac = 0.1,river.yscale = 5, #node.order = set_node_order,
#               river.usr = c(0, 1, -0.6, 1.6))

makeRiverplot(a.placenta.te, atac_cluster.te, rna_cluster.te,  min.frac = 0.1,river.yscale = 5, 
              #node.order = set_node_order,
river.usr = c(0, 1, -0.6, 1.6))

makeRiverplot(a.placenta.te, atac_cluster.te, rna_cluster.te,  min.frac = 0.22,river.yscale = 5, 
              #node.order = set_node_order,
river.usr = c(0, 1, -0.6, 1.6),river.lty = 0,label.cex = 2) #use min.frac 0.3 for clear trend



##matching rna - atac cluster id (riverplot roughly matched)

#atac    rna
6        9
3        5

9        10

5        6


7        4 #?
1        7


2        3

8,4      2 #?


###########do knn ratio analysis (quantitative method to linke atac cluster to rna cluster)################

########cluster distribution plot ####
dotDistri = function (cluster = NULL, id = NULL, type = NULL){
    cluster.sel = cluster[(cluster$type == type & cluster$cluster == id),]
    n_sel = nrow(cluster.sel)
    #cat ('select for cluster ',id,' n = ',n_sel," \n")
    color = ''
    #color = ifelse(type == 'atac','red','navy')
    
    if(type == 'atac'){
      color = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel)),'black' ,'red'   )
    }else if (type == 'rna'){
        color = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel)),'lightblue' ,'navy'   ) 
    }
    
    plot(cluster$UMAP1,cluster$UMAP2,pch = 16, type='p',col='grey',cex=0.5,xlab='UMAP1',ylab='UMAP2',main=paste(type," cells cluster ",id,"\nn = ",n_sel,sep=''),cex.main=2,cex.axis=2,cex.lab=2 ) 
    points(cluster.sel$UMAP1,cluster.sel$UMAP2,pch = 16, cex=0.5,col=color)
    return(paste("cluster ",id," ok",sep='') )
}

options(repr.plot.height=15,repr.plot.width=15)
par(mfrow=c(3,3))
dotDistri(cluster = cluster.all.coord, id = '6', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '3', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '9', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '5', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '7', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '1', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '2', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '4', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '12', type = 'atac')

par(mfrow=c(3,3))
dotDistri(cluster = cluster.all.coord, id = '8', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '10', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '11', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '13', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '14', type = 'atac')
dotDistri(cluster = cluster.all.coord, id = '15', type = 'atac')




par(mfrow=c(3,3))
dotDistri(cluster = cluster.all.coord, id = '8', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '9', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '5', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '11', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '6', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '4', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '7', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '1', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '3', type = 'rna')


par(mfrow=c(3,3))
dotDistri(cluster = cluster.all.coord, id = '2', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '12', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '10', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '13', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '14', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '15', type = 'rna')
dotDistri(cluster = cluster.all.coord, id = '16', type = 'rna')






######
dotDistri_both = function (cluster = NULL, id1 = NULL, id2 = NULL){
    #id1 for atac, id2 for rna
    cluster.sel1 = cluster[(cluster$type == 'atac' & cluster$cluster == id1),]
    n_sel1 = nrow(cluster.sel1)
    cluster.sel2 = cluster[(cluster$type == 'rna' & cluster$cluster == id2),]
    n_sel2 = nrow(cluster.sel2)
    
    ##stopifnot(sum(is.na(cluster.sel1)) == 0)
    ##stopifnot(sum(is.na(cluster.sel2)) == 0)
    
    ##write output for kdeplot
    ##fileout <- paste( "atac_c",id1,'_rna_c',id2,'.txt',sep=''   )
    ##write.table(x = rbind(cluster.sel1,cluster.sel2),file = fileout,sep='\t',row.names = TRUE,col.names = TRUE,quote = FALSE)
     
    color_atac = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel1)),'pink' ,'red'   )
    color_rna = ifelse(test=grepl(pattern='-1$',x=rownames(cluster.sel2)),'lightblue' ,'navy'   ) 
   
    stopifnot(length(color_atac) == nrow(cluster.sel1))
    stopifnot(length(color_rna) == nrow(cluster.sel2))
    
    plot(cluster$UMAP1,cluster$UMAP2,pch = 16, type='n',col='grey',cex=0.5,xlab='UMAP1',ylab='UMAP2',main=paste("atac cells cluster ",id1,", n = ",n_sel1,"\nrna cells cluster ",id2,", n = ",n_sel2,sep=''),cex.main=2,cex.axis=2,cex.lab=2 ) 
    #points(cluster.sel2$UMAP1,cluster.sel2$UMAP2,pch = 16, cex=0.2,col='navy') #rna
    #points(cluster.sel1$UMAP1,cluster.sel1$UMAP2,pch = 16, cex=0.2,col='red') #atac

    points(cluster.sel1$UMAP1,cluster.sel1$UMAP2,pch = 16, cex=0.5,col=color_atac) #atac
    points(cluster.sel2$UMAP1,cluster.sel2$UMAP2,pch = 16, cex=0.5,col=color_rna ) #rna
    
    return(paste("cluster ",id1,", ",id2," ok. ",sep='') )
    #return(paste("cluster ",id1,", ",id2," ok. coords output to ",fileout,sep='') )
}



#stb
options(repr.plot.width=15,repr.plot.height=15)
par(mfrow=c(2,2))
dotDistri_both(cluster = cluster.all.coord, id1 = '5',id2='6') #naive STB
dotDistri_both(cluster = cluster.all.coord, id1 = '2',id2='3') #PAPPA-path?
dotDistri_both(cluster = cluster.all.coord, id1 = '1',id2='7') #FLT1-path?
dotDistri_both(cluster = cluster.all.coord, id1 = '4',id2='2') #intermediate? ageing?


#ctb & fusion & evt
options(repr.plot.width=15,repr.plot.height=15)
par(mfrow=c(2,2))
dotDistri_both(cluster = cluster.all.coord, id1 = '3',id2='5') # TS-like CTB??
dotDistri_both(cluster = cluster.all.coord, id1 = '6',id2='9') #naive CTB
dotDistri_both(cluster = cluster.all.coord, id1 = '9',id2='11') #fusion
dotDistri_both(cluster = cluster.all.coord, id1 = '12',id2='12') #EVT

#str (four types, 5 clusters)
options(repr.plot.width=15,repr.plot.height=15)
par(mfrow=c(2,2))
dotDistri_both(cluster = cluster.all.coord, id1 = '13',id2='13')
dotDistri_both(cluster = cluster.all.coord, id1 = '10',id2='10')
dotDistri_both(cluster = cluster.all.coord, id1 = '14',id2='14')
dotDistri_both(cluster = cluster.all.coord, id1 = '15',id2='15')
dotDistri_both(cluster = cluster.all.coord, id1 = '11',id2='16')

#ctb proliferation rna only
options(repr.plot.width=15,repr.plot.height=15)
par(mfrow=c(2,2))
dotDistri_both(cluster = cluster.all.coord, id1 = 'none',id2='8')



# #detect convex then compare?
# see kde_contour.py

# #customized riverplot? no link atac dots and rna dots first by knn?



#####################evaluate cluster overlap percent by knn#################
library(FNN)
library(igraph)
##get_knnx method use k = 5




# #######################a little test###############
# set.seed(123)
# x1 = cbind(runif(10),runif(10)) ##reference xy coordinates, can extend to n dims?
# #x1.bk <- x1
# #row.names(x1) = paste0('x1_',1:nrow(x1)) #no use to add rownames
# len1 = nrow(x1)
# x2 = cbind(runif(12),runif(12)) ##query xy coordinates, can extend to n dims?
# #x2.bk <- x2
# #x2 = cbind(runif(2),runif(2))
# len2 = nrow(x2)
# #row.names(x2) = paste0('x2_',1:nrow(x2))

# ###plot in xy coordinate to view dots nearest directly####
# par(mfrow=c(1,3))
# options(repr.plot.width=15,repr.plot.height=5)
# plot(x1,cex=5,pch  = 19,col='black',main='reference',cex.main=3)
# text(x1,labels = 1:len1,col='grey',cex = 2 )
# plot(x2,cex=5,pch  = 19,col='grey',main='query',cex.main=3)
# text(x2,labels = 1:len2,col='red',cex = 2 )
# xall <- rbind(x1,x2)
# lenall <- nrow(xall)
# plot(xall,cex=5,pch  = 19,col=c(rep('black',len1),rep('grey',len2)),main='merge',cex.main=3 )#  c('black','grey'),each=10 ) )
# text(xall,labels = c(1:len1,1:len2),col=c(rep('grey',len1),rep('red',len2)),cex = 2 )#rep(c('grey','red'),each=10 ) )


# #######search and plot with igraph########

# ###plot x1 graph and x2 graph
# ##graph1 (reference)
# par(mfrow=c(1,1))
# options(repr.plot.width=5,repr.plot.height=5)
# nn1 = get.knn(x1,k)
# test.df1 = data.frame(from = rep(1:nrow(nn1$nn.index), k), 
#                     to = as.vector(nn1$nn.index), 
#                     weight = 1/(1 + as.vector(nn1$nn.dist))
#                    )
# test.df.nw1 = igraph::graph_from_data_frame(test.df1, directed = TRUE)
# plot(test.df.nw1,main='graph of reference',cex.main=3)


# ##graph2 (query)
# par(mfrow=c(1,1))
# options(repr.plot.width=5,repr.plot.height=5)
# nn2 = get.knn(x2,k)
# test.df2 = data.frame(from = rep(1:nrow(nn2$nn.index), k), 
#                     to = as.vector(nn2$nn.index), 
#                     weight = 1/(1 + as.vector(nn2$nn.dist))
#                    )
# test.df.nw2 = igraph::graph_from_data_frame(test.df2, directed = TRUE)
# plot(test.df.nw2,main='graph of query',cex.main=3)


# ##query x2 from x1
# #set.seed(100) #graph varies even set seed, but topology is similar
# k = 3
# par(mfrow=c(1,1))
# options(repr.plot.width=6,repr.plot.height=6)
# nn = get.knnx(x1,x2,k)
# #nn = get.knnx(x1,x2,k) #for one point query
# test.df = data.frame(from = paste0('q',rep(1:nrow(nn$nn.index), k)), 
#                     to = as.vector(nn$nn.index),
#                     #weight =  as.vector(nn$nn.dist)
#                     weight = 1/(1 + as.vector(nn$nn.dist))
#                    )
# test.df.nw = igraph::graph_from_data_frame(test.df, directed = TRUE)
# plot(test.df.nw,
#      main='graph of qeury to reference',
#      cex.main=3
#      )
#      #vertex.size=15,
#      #edge.width=E(test.df.nw)$weight)

# ##modify the network attributes
# get.data.frame(test.df.nw,what='vertices')

# #see vertex order
# V(test.df.nw)

# #set vertex color
# V(test.df.nw)$color <- c(rep('red',12),rep('orange',10))

# #plot the new graph
# plot(test.df.nw,
#      main='graph of qeury to reference',
#      cex.main=3
#      )

# # length(edge.betweenness(test.df.nw) )
# # length(E(test.df.nw)$weight )

# # #edge.betweenness(test.df.nw) = E(test.df.nw)$weight * 5 #not settable
# # edge_attr(test.df.nw) #see edge attribes

# # test.df.nw <- test.df.nw %>% set_edge_attr('color',value='red') 
# # plot(test.df.nw,main='graph of qeury to reference',cex.main=3)



# #get all nearest from reference
# unique(as.vector(nn$nn.index) )

# #filter eucliean dist?
# boxplot(as.vector(nn$nn.dist))
# filter <- quantile(as.vector(nn$nn.dist),probs = c(0.1,0.25,0.5,0.75,0.9,1) )[3]
# unique(as.vector(nn$nn.index) [as.vector(nn$nn.dist <= filter)])

# #####test end########






# ###start the real drive##
# id1='2'
# id2='3'
# cluster = cluster.all.coord
# cluster.sel1 = cluster[(cluster$type == 'atac' & cluster$cluster == id1),] #atac
# n_sel1 = nrow(cluster.sel1)
# cluster.sel2 = cluster[(cluster$type == 'rna' & cluster$cluster == id2),] #rna
# n_sel2 = nrow(cluster.sel2)

# stopifnot(sum(is.na(cluster.sel1)) == 0)
# stopifnot(sum(is.na(cluster.sel2)) == 0)


# x1 <- cluster.sel1[,c('UMAP1','UMAP2')] #atac, reference
# #x1 <- cluster[cluster$type == 'atac',c('dim1','dim2')] #all atac, reference, not so good
# len1 = nrow(x1) #1783 #605 #7394
# x2 <- cluster.sel2[,c('UMAP1','UMAP2')] #rna, query one cluster per time
# len2=nrow(x2) #1225 #1009

# ##plot in xy coordinate to view dots nearest directly####
# par(mfrow=c(1,3))
# options(repr.plot.width=15,repr.plot.height=5)
# plot(x1,cex=0.5,pch  = 19,col='red',main=paste('reference-atac:cluster',id1,sep = ''),cex.main=3)
# #text(x1,labels = 1:len1,col='grey',cex = 2 )
# plot(x2,cex=0.5,pch  = 19,col='navy',main=paste('query-rna:cluster',id2,sep = ''),cex.main=3)
# #text(x2,labels = 1:len2,col='red',cex = 2 )
# xall <- rbind(x1,x2)
# lenall <- nrow(xall)
# plot(xall,cex=0.5,pch  = 19,col=c(rep('red',len1),rep('navy',len2)),main='merge',cex.main=3 )#  c('black','grey'),each=10 ) )
# #text(xall,labels = c(1:len1,1:len2),col=c(rep('grey',len1),rep('red',len2)),cex = 2 )#rep(c('grey','red'),each=10 ) )


# ##query x2 from x1
# #set.seed(100) #graph varies even set seed, but topology is similar
# k = 5
# #par(mfrow=c(1,1))
# #options(repr.plot.width=5,repr.plot.height=5)
# nn = get.knnx(x1,x2,k)
# #nn = get.knnx(x1,x2,k) #for one point query
# # nn.df = data.frame(from = rep(1:nrow(nn$nn.index), k), 
# #                     to = as.vector(nn$nn.index), 
# #                     weight =  as.vector(nn$nn.dist)
# #                     #weight = 1/(1 + as.vector(nn$nn.dist))
# #                    )
# #test.df.nw = igraph::graph_from_data_frame(test.df, directed = TRUE)
# #plot(test.df.nw,main='graph of qeury to reference',cex.main=3)


# #all nearest from reference
# NN_index <- unique(as.vector(nn$nn.index) ) #1530 #545 #1702
# #545/605 = 0.90

# #see hits dot cluster distribution in atac
# par(mfrow=c(1,1))
# options(repr.plot.width=8,repr.plot.height=8)
# plot(x1,cex=0.5,pch  = 19,col='grey',main='query-rna-in-reference-atac',cex.main=3)
# points(x1[NN_index,],cex=0.5,pch  = 19,col='red')
# points(x2,cex=0.5,pch  = 19,col='navy')
# legend("topright",legend = c('scATAC','scATAC(outside)','scRNA'),fill = c('red','grey','navy'))


#cluster.hits <- cluster[rownames(x1[NN_index,]),c('cluster','type','dim1','dim2')] 

#table(cluster.hits$type)
#table(cluster.hits$cluster)




# #########read in snapATAC cluster and plot######################

# cluster.snap.filter <- read.table("snapATAC.PLA-8w-ATAC-1.umap.bin50k.filtered.txt",
#                                   sep = '\t',
#                                   header=TRUE,
#                                   quote="",
#                                   row.names=1
#                                  )
# #table(cluster.snap.filter$cluster_snap)
# cluster.snap.filter$cluster_snap <- factor(cluster.snap.filter$cluster_snap, sort(unique(cluster.snap.filter$cluster_snap)) )
# ##6539 


# #library('BuenColors')
# library('cowplot')
# ggplot(data=cluster.snap.filter,aes(x = UMAP1_snap,y = UMAP2_snap, col = cluster_snap) ) +
#   geom_point(size = 1) +
#   scale_color_manual(values=color_snap) +
#   theme_cowplot()
#   #pretty_plot()



# ##substitute snapATAC cluster 
# ids.atac <- rownames(cluster.all.coord[cluster.all.coord$type == 'atac',] ) #7394

# ids.share <- intersect(ids.atac,rownames(cluster.snap.filter)) #6491

# cluster.temp1 <- cluster.all.coord[ids.share,] 
# cluster.temp2 <- cluster.snap.filter[ids.share,]



# cluster.temp1$cluster <-  cluster.temp2$cluster_snap


# cluster.all.coord.substitute_snap <- rbind(cluster.temp1,cluster.all.coord[cluster.all.coord$type == 'rna',])


# ##plot
# par(mfrow=c(3,3))
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '2', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '10', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '3', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '1', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '4', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '5', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '8', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '6', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '13', type = 'atac')

# par(mfrow=c(3,3))
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '9', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '15', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '12', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '14', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '11', type = 'atac')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '7', type = 'atac')

# par(mfrow=c(3,3))
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '7', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '8', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '6', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '10', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '2', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '5', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '1', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '4', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '3', type = 'rna')


# par(mfrow=c(3,3))
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '9', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '13', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '11', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '12', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '14', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '15', type = 'rna')
# dotDistri(cluster = cluster.all.coord.substitute_snap, id = '16', type = 'rna')
# #dotDistri(cluster = cluster.all.coord.substitute_snap, id = '16', type = 'rna')







# options(repr.plot.width=12,repr.plot.height=5)
# par(mfrow=c(1,3))
# dotDistri_both(cluster = cluster.all.coord.substitute_snap, id1 = '4',id2='1')
# dotDistri_both(cluster = cluster.all.coord.substitute_snap, id1 = '8',id2='4')
# dotDistri_both(cluster = cluster.all.coord.substitute_snap, id1 = '5',id2='3')


doKNN <- function(cluster = cluster.all.coord.substitute_snap,id1='5',id2='3',k=5,leg.pos = 'bottomright' ){
    ###start the real drive##
    #id1='5' #'8' #'4'
    #id2='3'#'4' #'1'
    #cluster = cluster.all.coord.substitute_snap
    cluster.sel1 = cluster[(cluster$type == 'atac' & cluster$cluster == id1),]
    #n_sel1 = nrow(cluster.sel1)
    cluster.sel2 = cluster[(cluster$type == 'rna' & cluster$cluster == id2),]
    #n_sel2 = nrow(cluster.sel2)

    stopifnot(sum(is.na(cluster.sel1)) == 0)
    stopifnot(sum(is.na(cluster.sel2)) == 0)

    x1 <- cluster.sel1[,c('UMAP1','UMAP2')] #atac, reference
    #x1 <- cluster[cluster$type == 'atac',c('dim1','dim2')] #all atac, reference, not so good
    len1 = nrow(x1) #504 #424 #765 #605 #7394
    x2 <- cluster.sel2[,c('UMAP1','UMAP2')] #rna, query one cluster per time
    len2=nrow(x2) #1009

    
#     #####do some statistic test ?###
#     x <- x1
#     y <-  x2
#     ks <- peacock2(x, y)
#     cat ('two-dimensional Kolmogorov-Smirnov test,ks=:',ks,'\n')

#     par(mfrow=c(1,1))
#     options(repr.plot.width=2.5,repr.plot.height=2.5)
#     plot(x,pch=19,cex = 0.8, col = 'red', main = paste('two-dimensional Kolmogorov-Smirnov test,ks=:',ks,sep="") )
#     points(y,pch=19,cex = 0.8, col = 'navy')

    ##plot in xy coordinate to view dots nearest directly####
    par(mfrow=c(2,2))
   # cat('Before setting options:width=',options()$repr.plot.width,' height=',options()$repr.plot.height,'\n')
    options(repr.plot.width=10,repr.plot.height=10)
    #cat('After setting options:width=',options()$repr.plot.width,' height=',options()$repr.plot.height,'\n')
    
    plot(x1,cex=0.5,pch  = 19,col='red',main=paste('reference-atac','(c',id1,")"),cex.main=2,cex.lab=2)
    #text(x1,labels = 1:len1,col='grey',cex = 2 )
    plot(x2,cex=0.5,pch  = 19,col='navy',main=paste('query-rna','(c',id2,")"),cex.main=2)
    #text(x2,labels = 1:len2,col='red',cex = 2 )
    xall <- rbind(x1,x2)
    lenall <- nrow(xall)
    plot(xall,cex=0.5,pch  = 19,col=c(rep('red',len1),rep('navy',len2)),main='merge',cex.main=2 )#  c('black','grey'),each=10 ) )
    #text(xall,labels = c(1:len1,1:len2),col=c(rep('grey',len1),rep('red',len2)),cex = 2 )#rep(c('grey','red'),each=10 ) )


    ##query x2 from x1
    #set.seed(100) #graph varies even set seed, but topology is similar
    ###k = 5
    #par(mfrow=c(1,1))
    #options(repr.plot.width=5,repr.plot.height=5)
    nn = get.knnx(x1,x2,k)
    #nn = get.knnx(x1,x2,k) #for one point query
    # nn.df = data.frame(from = rep(1:nrow(nn$nn.index), k), 
    #                     to = as.vector(nn$nn.index), 
    #                     weight =  as.vector(nn$nn.dist)
    #                     #weight = 1/(1 + as.vector(nn$nn.dist))
    #                    )
    #test.df.nw = igraph::graph_from_data_frame(test.df, directed = TRUE)
    #plot(test.df.nw,main='graph of qeury to reference',cex.main=3)


    #all nearest from reference
    NN_index <- unique(as.vector(nn$nn.index) ) #644 #545 #1702
    #344/504 = 68% #396/424 = 93% 644/765 = 0.84
    perc_hit <- 100*(length(NN_index)/len1 )
    perc_hit <- round(perc_hit)
    #see hits dot cluster distribution in atac
    #par(mfrow=c(1,1))
    #options(repr.plot.width=8,repr.plot.height=8)
    plot(x1,cex=0.5,pch  = 19,col='grey',main=paste('query rna in \nreference atac. (',perc_hit,"%)",sep = ''),cex.main=2)
    points(x1[NN_index,],cex=0.5,pch  = 19,col='red')
    points(x2,cex=0.5,pch  = 19,col='navy')
    #legend("topright",legend = c('scATAC','scATAC(outside)','scRNA'),
    legend(leg.pos,legend = c('scATAC','scATAC(outside)','scRNA'),
           fill = c('red','grey','navy'),
           x.intersp = 0.8, #x space
           y.intersp = 0.8, #y space
           cex = 1,
           text.width = 3.2 #overall text width
          )
    

    return(paste('atac:c',id1,' ','rna:c',id2,' ok. ','knn hits: ',perc_hit,"%",sep='') )
}



##for STB subtype, atac vs rna
pdf('pdfs_te/nearest-neighbor-ratio.c2atac-vs-c3rna.pdf',width=10,height=10,useDingbats = FALSE)
  doKNN(cluster = cluster.all.coord,id1='2',id2='3',k=5,leg.pos = 'bottomright' ) #86% #84% 
dev.off()

pdf('pdfs_te/nearest-neighbor-ratio.c1atac-vs-c7rna.pdf',width=10,height=10,useDingbats = FALSE)
 doKNN(cluster = cluster.all.coord,id1='1',id2='7',k=5 ) #76% #93%
dev.off()

pdf('pdfs_te/nearest-neighbor-ratio.c5atac-vs-c6rna.pdf',width=10,height=10,useDingbats = FALSE)
 doKNN(cluster = cluster.all.coord,id1='5',id2='6',k=5 ) #73% #68%
dev.off()

##for CTB,fusion,  matching
pdf('pdfs_te/nearest-neighbor-ratio.c3atac-vs-c5rna.pdf',width=10,height=10,useDingbats = FALSE)
 doKNN(cluster = cluster.all.coord,id1='3',id2='5',k=5,leg.pos = 'topright'  ) #83% CTB  
dev.off()

pdf('pdfs_te/nearest-neighbor-ratio.c6atac-vs-c9rna.pdf',width=10,height=10,useDingbats = FALSE)
 doKNN(cluster = cluster.all.coord,id1='6',id2='9',k=5 ,leg.pos = 'topright') #80% CTB   
dev.off()

pdf('pdfs_te/nearest-neighbor-ratio.c9atac-vs-c11rna.pdf',width=10,height=10,useDingbats = FALSE)
  doKNN(cluster = cluster.all.coord,id1='9',id2='11',k=5,leg.pos = 'topright' ) #83% fusion
dev.off()

##for EVT matching
doKNN(cluster = cluster.all.coord,id1='12',id2='12',k=5 ) #92% EVT

##STR matching
doKNN(cluster = cluster.all.coord,id1='13',id2='13',k=5 )


###summary as table

"cluster atac id"  "cluster rna id"  "nearest-neighbor-overlap-ratio"  "cluster_integration id"
c2(stb-PAPPA) c3(stb-PAPPA)    87%                           c2,c3(part),c4(part)
c1(stb-FLT1)          c7(stb-FLT1)       68%                           c1(part), c4(part)
c5(stb-naive)         c6(stb-naive)      64%                           c5
c9(ctb-fusion compenent)  c11(ctb-fusion competent)  84%               c9 
c3(ctb-2)             c5(ctb-2)            79%                         c4(?)
c6(ctb-1)             c9(ctb-1)           75%                          c7

##to result.nearest-neighbor-ratio.txt


##corrected version
"cluster atac id"  "cluster rna id"  "nearest-neighbor-overlap-ratio"  "cluster_integration id"
c2(stb-PAPPA) c3(stb-PAPPA)    87%                           c2,c3(part)
c1(stb-FLT1)          c7(stb-FLT1)       68%                           c4(part),c1(part)
c5(stb-naive)         c6(stb-naive)      64%                           c5
c9(ctb-fusion compenent)  c10(ctb-fusion competent)  84%               c8 
c3(ctb-2)             c5(ctb-2)            79%                         c6
c6(ctb-1)             c9(ctb-1)           75%                          c7

c4(stb-PAPPA-aged?)     c2(stb-PAPPA-aged?) 




####collect and plot table output####


source('dogridTable.r')

res.tab <- read.table('result.nearest-neighbor-ratio.txt',header = TRUE, sep='\t')

res.g <- doGridTable(data=res.tab, w=15,h=3,basesize = 20,plot=TRUE) 


#####do some statistic test ?###
library('Peacock.test')
x <- matrix(rnorm(120, 10, 1), ncol=2)
y <-  matrix(rnorm(160, 10, 1), ncol=2)
ks <- peacock2(x, y)
ks

options(repr.plot.width=5,repr.plot.height=5)
plot(x,pch=19,cex = 0.8, col = 'red',xlim=c(0,15),ylim=c(0,15))
points(y,pch=19,cex = 0.8, col = 'navy')

ks.test(x,y)











######################start data imputation##########################

###########add binarized peak mat to liger obj#######################

a.placenta.ds <- a.placenta.te #21268 total cells (rna + atac)
table(a.placenta.te@cell.data$dataset)

atac1 atac2  rna1  rna2 
 4549  6525  5978  4216


pmat <- readRDS('pmat.processed.rds') #11093 x 299537
#table(rownames(pmat) %in% rownames(a.placenta.ds@tsne.coords ) )#no

##rename pmat rowname to -1, -2
cellid <- rownames(pmat)
idy <- which(grepl(pattern = "^placenta_donor2",x = cellid    ))
cellid[idy] <- gsub( pattern = '-1$',replacement = '-2' ,x =  cellid [idy]  )
cellid  <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = cellid )
cellid <- paste0('placenta.atac_',cellid)
rownames(pmat) <- cellid 

table(rownames(pmat) %in% rownames(a.placenta.ds@tsne.coords ) )
FALSE  TRUE 
   19 11074

table(colnames(a.placenta.ds@raw.data[['atac1']] ) %in% rownames(pmat) )
TRUE 
4549

table(colnames(a.placenta.ds@raw.data[['atac2']] ) %in% rownames(pmat) )
TRUE 
6525

##add to obj (with substitution)
#a.placenta.ds@raw.data[['atac']] <- t(pmat) #299537 x 11093

a.placenta.ds@raw.data[['atac1']] <- t(pmat)[,colnames(a.placenta.ds@raw.data[['atac1']] )]
a.placenta.ds@raw.data[['atac2']] <- t(pmat)[,colnames(a.placenta.ds@raw.data[['atac2']] )]


# #merge two rna slot?

# rna1.mat <- a.placenta.ds@raw.data[['rna1']] #24355 x 5978
# rna2.mat <- a.placenta.ds@raw.data[['rna2']] #24578 x 4216

# sharedid <- intersect(rownames(rna1.mat),rownames(rna2.mat) )
# #23032

# rna.mat <-  cbind(rna1.mat[sharedid,],rna2.mat[sharedid,])
# #23032 x 10194

# a.placenta.ds@raw.data[['rna']] <- rna.mat #merge two rna donor
# #a.placenta.ds@raw.data[['rna1']] <- NULL
# #a.placenta.ds@raw.data[['rna2']] <- NULL
# #a.placenta.ds@raw.data[['atac1']] <- NULL
# #a.placenta.ds@raw.data[['atac2']] <- NULL

a.placenta.ds <- liger::normalize(a.placenta.ds)

saveRDS(a.placenta.ds,'peak_gene_link_liger/a.placenta.ds.rds')

##try to imputate data, from rna to atac

#a.placenta.ds@H.norm

#a.placenta.ds_torna1 <- imputeKNN(a.placenta.ds, reference = 'atac1',queries='rna1',knn_k = 50)
#a.placenta.ds_torna2 <- imputeKNN(a.placenta.ds, reference = 'atac2',queries='rna2',knn_k = 50)

a.placenta.ds_torna <- imputeKNN(a.placenta.ds, reference = 'atac1',queries='rna1',knn_k = 50) #will modify rna1 slot
a.placenta.ds_torna <- imputeKNN(a.placenta.ds_torna, reference = 'atac2',queries='rna2',knn_k = 50)

#a.placenta.ds_toatac1 <- imputeKNN(a.placenta.ds, reference = 'rna1',queries='atac1',knn_k = 50)
#a.placenta.ds_toatac2 <- imputeKNN(a.placenta.ds, reference = 'rna2',queries='atac2',knn_k = 50)

a.placenta.ds_toatac <- imputeKNN(a.placenta.ds, reference = 'rna1',queries='atac1',knn_k = 50) #will modify atac slot
a.placenta.ds_toatac <- imputeKNN(a.placenta.ds_toatac, reference = 'rna2',queries='atac2',knn_k = 50)



##correlate peak with gene


gmap <- a.placenta.te@norm.data[['rna1']]

pmat <- a.placenta.ds_torna@norm.data[['rna1']] 


regnet <- linkGenesAndPeaks(gene_counts = gmap, peak_counts = pmat, dist = 'spearman',alpha = 0.05, path_to_coords = 'atac_activity/genebodyandpromoter.GRCh38.all.bed') #very slow

saveRDS(regnet ,'peak_gene_link_ligerregnet.rds')
#299537 x 18190 (peak x gene)


####output certain gene
PAPPA <- regnet[, 'PAPPA']
PAPPA <- PAPPA[abs(PAPPA) > 0]  #130 with positive and negative
#View(PAPPA[order(abs(PAPPA), decreasing = TRUE)])

table(PAPPA > 0)
FALSE  TRUE 
    9   121 

options(repr.plot.width=5,repr.plot.height=5)
hist(PAPPA,breaks = 50) #most >=0.2
abline(v=0.2,col='red')

#output PAPPA p2g links
table(PAPPA > 0.25)

FALSE  TRUE 
   64    66

outbed <- PAPPA[PAPPA > 0.25]

res.df <- data.frame()
for(i in seq_along(outbed)){
    v = outbed[i]
    names(v) <- NULL
    id = names(outbed[i])
    id.split <- unlist(strsplit(id,split = ':|-') )
    res.df <- rbind.data.frame(res.df ,data.frame('chr'=id.split[1],'start'=id.split[2],'end'=id.split[3],'id'=id,'score'=v,'strand'='+'))
    
}

write.table(res.df,file='peak_gene_link_liger/PAPPA.p2g.bed',sep='\t',quote = FALSE,row.names = FALSE,col.names =FALSE)



##output all gene p2g bed
cnt <- 0
ngene <- ncol(regnet)
for( gene in colnames(regnet) ) { #18190
    cnt <- cnt + 1
    #if(cnt %% 100 == 0){cat('#')}
    #if(cnt == ngene){cat('\n')}
    
    p2g <- regnet[, gene]
    ##p2g <- p2g[abs(p2g) > 0]  #130 with positive and negative
    #View(PAPPA[order(abs(PAPPA), decreasing = TRUE)])

#     table(p2g > 0)
#     FALSE  TRUE 
#         9   121 
    
    #options(repr.plot.width=5,repr.plot.height=5)
    #hist(PAPPA,breaks = 50) #most >=0.2
    #abline(v=0.2,col='red')

    
    
    #output PAPPA p2g links
#     table(PAPPA > 0.25)

#     FALSE  TRUE 
#        64    66

    outbed <- p2g[p2g > 0.25]

    res.df <- data.frame()
    for(i in seq_along(outbed)){
        v = outbed[i]
        names(v) <- NULL
        id = names(outbed[i])
        id.split <- unlist(strsplit(id,split = ':|-') )
        res.df <- rbind.data.frame(res.df ,data.frame('chr'=id.split[1],'start'=id.split[2],'end'=id.split[3],'id'=paste(id,'-',gene,sep=''),'score'=v,'strand'='+'))

    }

    write.table(res.df,file='peak_gene_link_liger/all.p2g.bed',sep='\t',quote = FALSE,row.names = FALSE,col.names =FALSE,append = TRUE)

}

# ######visualize with liger track
# regnet.filter <- regnet[!duplicated(colnames(regnet)),]
# 299537 x 18190

# makeInteractTrack(regnet.filter, genes.list = 'PAPPA', path_to_coords = 'atac_activity/genebodyandpromoter.GRCh38.all.bed')






#######final save object###
saveRDS('a.placenta.te.rds',object = a.placenta.te)
#a.placenta.te=readRDS('a.placenta.te.rds')

#cluster.all.coord <- readRDS('cluster.all.coord.rds')
cluster.df.add <- readRDS('cluster.df.add.rds')
cluster.df.add.te <- readRDS('cluster.df.add.te.rds')

#all.equal(cluster.all.coord[,c('cluster','UMAP1','UMAP2')], cluster.df.add[,c('cluster','UMAP_1','UMAP_2')]  )


