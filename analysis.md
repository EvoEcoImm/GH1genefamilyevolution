### 1. database buildup

__extract GH1 proteins from dbCAN by Biopython__
extract GH1 hmm from whole hmm database. and run_dbCAN.py with only GH1 prediction.

```bash
##setup run_dbCAN with Anoconda
conda create -n run_dbcan python=3.8 diamond hmmer prodigal -c conda-forge -c bioconda
conda activate run_dbcan
pip install run-dbcan==2.0.11
test -d db || mkdir db
cd db \
    && wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312019.fa.nr && diamond makedb --in CAZyDB.07312019.fa.nr -d CAZy \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt && mv dbCAN-HMMdb-V8.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff

#database and running tools, run_dbcan, diamond, hmmer,hotpep; form the hmmer database, diamond database and hotpep databases.
##tools
wget https://files.pythonhosted.org/packages/ae/87/431b9b651544ead5b0889c39e6fddb9e1dd98e773c2ed5c229e18e063ee6/run_dbcan-2.0.11.tar.gz && tar zxvf run_dbcan-2.0.11.tar.gz
wget http://bcb.unl.edu/dbCAN2/download/Tools/hotpep-python-08-06-2020.tar.gz
wget http://bcb.unl.edu/dbCAN2/download/Tools/hmmscan-parser.gz
#wget http://bcb.unl.edu/dbCAN2/download/Tools/run_dbcan.tar.gz

pip install run-dbcan

#database
wget http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V9.txt
wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa
wget http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07302020.fam.subfam.ec.txt
wget http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07302020.fam-activities.txt
wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm
wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa
wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm
wget http://bcb.unl.edu/dbCAN2/download/Databases/tf.fa
wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm

##filter GH1 for make database
#edit Hotpep family file in Hotpep folder
mv fam_list.txt fam_list_BACKUP.txt && echo "GH1" > fam_list.txt
#get GH1 protein sequences and build up database for diamond
awk 'BEGIN{RS=">";FS="\n"}$1~/\<GH1\>/{print ">"$0}' CAZyDB.07312020.fa |sed '/^$/d' > CAZyDB.GH1.fa && diamond makedb --in CAZyDB.GH1.fa -d CAZy
#get GH1 profile and build up profile for diamond
hmmfetch dbCAN-HMMdb-V9.txt GH1.hmm > GH1.dbCAN.txt && hmmpress GH1.dbCAN.txt 
diamond makedb --in tcdb.fa -d tcdb
hmmpress tf-1.hmm
hmmpress tf-2.hmm
hmmpress stp.hmm
```

### 2. Download genome assemblies and annotations

search genome resources
Refine "Eukaryotes" in Genome of NCBI; search "Arthropoda" ; restriction in: 1)Asseembly level, Chromosome+Scaffold; 2) RefSeq category, reference+representative
in table, select RefSeqFTP not empty; keep all chromosome assembly and one representative genome in rest of each genus:1)level, chromosome > scaffold; 2) modified Date, newest> older 

```r
#!/usr/bin/env Rscript

#download genome resources from NCBI genome database

library(readODS)
library(RCurl)

args = commandArgs(trailingOnly=TRUE)
setwd(".")
print(args[1])
assembly_list<-read_ods(args[1],sheet=2)

dla_rsync<- function(doc){
    doc["GenBank FTP"]<-lapply(doc["GenBank FTP"],function(colsub) sub("ftp","rsync",colsub))
    for (ni in 1:nrow(doc)){
        dirn<-sub(" ","_",doc[ni,"#Organism Name"])
        dir.create(dirn)
        downloadscript_cds<- paste0("rsync --copy-links --times --verbose ",doc[ni,"GenBank FTP"],"/",doc[ni,"Assembly"],"_*cds_from_genomic.fna.gz ./", dirn)
        downloadscript_protein<- paste0("rsync --copy-links --times --verbose ",doc[ni,"GenBank FTP"],"/",doc[ni,"Assembly"],"_*protein.faa.gz ./", dirn)
        downloadscript_gff<- paste0("rsync --copy-links --times --verbose ",doc[ni,"GenBank FTP"],"/",doc[ni,"Assembly"],"_*genomic.gff.gz ./", dirn)
        system(downloadscript_cds)
        system(downloadscript_protein)
        system(downloadscript_gff)
        }
    }

dla_https<- function(doc){
    doc["RefSeq FTP"]<-lapply(doc["RefSeq FTP"],function(colsub) sub("ftp","https",colsub))
    for (ni in 1:nrow(doc)){
        dirn<-sub(" ","_",doc[ni,"#Organism Name"])
        filepref<-gsub(".*/","",doc[ni,"RefSeq FTP"])
        dir.create(dirn)
        downloadscript_cds<- paste0("wget ",doc[ni,"RefSeq FTP"],"/",filepref,"_cds_from_genomic.fna.gz -P ./", dirn)
        downloadscript_protein<- paste0("wget ",doc[ni,"RefSeq FTP"],"/",filepref, "_protein.faa.gz -P ./", dirn)
        downloadscript_gff<- paste0("wget ",doc[ni,"RefSeq FTP"],"/",filepref,"_genomic.gff.gz -P ./", dirn)
        system(downloadscript_cds)
        system(downloadscript_protein)
        system(downloadscript_gff)
        }
    }

setwd("Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq")
dla_https(assembly_list)
```

seperate download Blattella gemanica genome
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/018/175/GCA_003018175.1_Bger_1.1

### 3.predict/extract GH1 from assembly of genomes in Genome/resource/Arthropoda_Refseq/genomes

```bash
#predict GH1 with run_dbCAN.py
ls /media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/ >Genome_specids #get species id
for i in `cat Genome_specids`; do python run_dbcan.py /media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/$i/GC_*_protein.faa protein --out_pre $i --out_dir $i --db_dir ../ ;done 
```

### 4. get details and sequences of predicted proteins.

__R scripts to extract gff from genes__

```r
library(rtracklayer)
library(Biostrings)
library(rentrez)
library(tidyverse)

###This getting cds from genebank probably is not perfect because the protein sequences do not match.
if (FALSE) {
get_gb<-function(gid){ #retrieve genbank format based on gene id from ncbi
    gidsum<-entrez_summary(db="gene",id=gid)
    if(!length(gidsum$genomicinfo)){
        gbcommand<-paste0("efetch -db gene -id ",gid, " -format docsum|xtract -pattern DocumentSummary -block LocationHist -first ChrAccVer ChrStart ChrStop")
        cd_out<-system(gbcommand,intern=T)
        gbcoo<- unlist(strsplit(cd_out,"\t"))
        stopco<- as.numeric(gbcoo[3])+1
        startco<- as.numeric(gbcoo[2])+1
        acc<-gbcoo[1]
        gidgb<-entrez_fetch(db="nuccore", id=acc, rettype="gb", seq_start=startco, seq_stop=stopco)
        } else {
        gidgb<-entrez_fetch(db="nuccore", id=gidsum$genomicinfo$chraccver, rettype="gb", seq_start=gidsum$genomicinfo$chrstart+1, seq_stop=gidsum$genomicinfo$chrstop+1)
        }
    return(gidgb)
    }
}
###This can work with nuc seqs.
get_genes_gff<-function(spec){ #define function to extract and output gene gffs,represent_proteins
    GH1<-read.delim(paste0(spec,"/",spec,"overview.txt"),stringsAsFactors=F)
    GH1s<-GH1[GH1$X.ofTools>=3,]
    if (!nrow(GH1s)){
    genes_gff<-""
    file(paste0(spec,"_GH1.gff3"),"w")
    file(paste0(spec,"_GH1.fna"),"w")
    } else {
    genome_folder_path<-"/media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/"
    gff_fn<-list.files(paste0(genome_folder_path,spec),pattern = "\\.gff$")
    gff_fp<-file.path(genome_folder_path,spec,gff_fn)
    gff<-import(gff_fp,format="gff3")
    if("locus_tag"%in%colnames(mcols(gff))){
        locus_tag_value<-subset(gff,gff$protein_id%in%GH1s$Gene.ID)$locus_tag
        if (any(is.na(locus_tag_value))){
        GH1s_gff<-gff[gff$gene%in%(unique(subset(gff,gff$protein_id%in%GH1s$Gene.ID)$gene))]
        } else {
        GH1s_gff<-gff[gff$locus_tag%in%(unique(locus_tag_value))]
        }
        } else {
        GH1s_gff<-gff[gff$gene%in%(unique(subset(gff,gff$protein_id%in%GH1s$Gene.ID)$gene))]
        }
    con<-file(paste0(spec,"_GH1.gff3"),"w")
    export(GH1s_gff,con,format="gff3")
    close(con)
    genes_gff<- GH1s_gff[GH1s_gff$type=="gene"]
    #if("locus_tag"%in%colnames(mcols(genes_gff))&!any(is.na(genes_gff$locus_tag))){
    #    gene_names<- unlist(lapply(genes_gff$locus_tag, function(x) sub("[A-z]*[A-z]","",x)))
    #    } else {
    #    gene_names<- unlist(lapply(genes_gff$gene, function(x) sub("[A-z]*[A-z]","",x)))
    #    }
    #genegbs<-""
    #for (gname in gene_names){
    #    gidgbs<-get_gb(gname)
    #    genegbs<-paste0(genegbs,gidgbs)
    #    }
    #gbfilename <- paste0(spec,"_GH1s.gb")
    #write(genegbs,file=gbfilename,append=T)
    return(genes_gff)
    }
}
specs<-readLines("Genome_specids") #get species ids
specs_genes<-list()
for (specid in specs){ #except "Culex_quinquefasciatus" and "Pediculus_humanus", manually download, get genbank files and gff3 structures for all the genes
    try(specs_genes[[specid]]<-get_genes_gff(specid))
    }

###############################################################
#Summary genes and gene details
################################################################
#__retrieve summary data__

library(dplyr)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)
specs<-readLines("Genome_specids_selected_tree")
options(ucscChromosomeNames=FALSE)

get_r_summary<-function(gff3){
    sTxDB<-makeTxDbFromGFF(gff3,format="gff3")
    genelist<-transcriptsBy(sTxDB,by="gene")
    #tranlist<-lapply(exonsBy(sTxDB,by="gene"),function(x)length(disjoin(x)))
    #for(i in names(genelist)){genelist[[i]]$nexon<-tranlist[[i]]}
    gened<-as.data.frame(genelist)
    rownames(gened)<-gened$tx_name
    t2nexon<-as.data.frame(sapply(exonsBy(sTxDB,by="tx",use.name=T),function(x)length(x)))
    t2cdsl<-as.data.frame(sapply(cdsBy(sTxDB,by="tx",use.name=T),function(x)sum(width(x))))
    t2cdsn<-as.data.frame(sapply(cdsBy(sTxDB,by="tx",use.name=T),function(x)unique(x$cds_name)))
    g1<-cbind(t2cdsn,t2cdsl[,1][match(rownames(t2cdsn),rownames(t2cdsl))],t2nexon[,1][match(rownames(t2cdsn),rownames(t2nexon))])
    colnames(g1)<- c("proteinname","wcds","nexon")
    genesummary<-merge(gened,g1,by="row.names")[,-1]
    return(genesummary)
    }

for (spec in specs){
    GH1<-read.delim(paste0(spec,"/",spec,"overview.txt"),stringsAsFactors=F)
    GH1s<-GH1[GH1$X.ofTools>=3,]
    nps<-nrow(GH1s)
    if (!nps) {
    geneids<-""
    } else {
    GH1s$ndom<-str_count(GH1s$HMMER,'\\+')+1
    abc<-get_r_summary(paste0(spec,"_GH1.gff3"))
    abcn<-merge(abc,GH1s[,c(1,2,7)],by.x="proteinname",by.y="Gene.ID")
    abcn$Org<-spec
    if (exists("allsummary")){allsummary<-rbind(allsummary,abcn)} else {allsummary<-abcn}
    }
    }
write.table(allsummary,"all_genome_GH1.gene.summary.csv",quote=F,sep="\t",row.names=F)

#plot gene number to tree
library(tidyverse)
library(phytools)
all<-read.csv("all_genome_GH1.gene2protein.summary.csv",stringsAsFactor=F,row.names=1)
summaryall<- all %>% group_by(Org) %>% summarise(proteincount=n(),genecount=n_distinct(group_name)) %>% as.data.frame
spectreetxt<-"((Limulus_polyphemus,((Varroa_jacobsoni,(Dermatophagoides_pteronyssinus,Tetranychus_urticae)),Centruroides_sculpturatus)),(Penaeus_vannamei,(Daphnia_magna,(Folsomia_candida,((Blattella_germanica,(Zootermopsis_nevadensis,Cryptotermes_secundus)),((Frankliniella_occidentalis,((Diaphorina_citri,(Aphis_gossypii,Acyrthosiphon_pisum)),(Nilaparvata_lugens,Cimex_lectularius))),(Pediculus_humanus,((Cephus_cinctus,(Orussus_abietinus,((Nasonia_vitripennis,Copidosoma_floridanum),(Polistes_dominula,((Harpegnathos_saltator,(Atta_colombica,Acromyrmex_echinatior)),(Dufourea_novaeangliae,(Apis_mellifera,Bombus_terrestris))))))),((Agrilus_planipennis,((Onthophagus_taurus,Nicrophorus_vespilloides),(Tribolium_castaneum,((Anoplophora_glabripennis,Diabrotica_virgifera),(Sitophilus_oryzae,Dendroctonus_ponderosae))))),((Plutella_xylostella,((Papilio_machaon,(Danaus_plexippus,Bicyclus_anynana)),(Hyposmocoma_kahamanoa,((Galleria_mellonella,Amyelois_transitella),((Helicoverpa_armigera,Spodoptera_litura),(Bombyx_mandarina,Manduca_sexta)))))),(((Culex_quinquefasciatus,Aedes_aegypti),Anopheles_gambiae),(Drosophila_melanogaster,(Rhagoletis_zephyria,Musca_domestica)))))))))))));"
spectree<-read.tree(text=spectreetxt)
#dotTree(spectree,snumbers[,1:2])
nogene<-data.frame("Org" = c("Limulus_polyphemus","Dermatophagoides_pteronyssinus","Tetranychus_urticae","Centruroides_sculpturatus","Penaeus_vannamei"), "proteincount" = rep(0,5), "genecount" = rep(0,5))
countall<-rbind(summaryall,nogene)
rownames(countall)<-countall[,1]
plotTree.barplot(spectree,countall[,c(3,2)],args.plotTree=list(cex=0.5),args.barplot=list(col=c("blue","red"),legend.text=c("gene","protein"),main="No. of genes and proteins",args.legend=list(bty="n"),beside=TRUE))
#dotTree(spectree,snumbers[,1:2])
dev.print(pdf,"Number_summary_pandg.pdf")
dev.off()

summarynumber<- function(spec){ #return no. of proteins and no. of genes
    GH1<-read.delim(paste0(spec,"/",spec,"overview.txt"),stringsAsFactors=F)
    GH1s<-GH1[GH1$X.ofTools>=3,]
    nps<-nrow(GH1s)
    if (!nps) {
    ngs<-0
    } else {
    GH1s_gff<-import(paste0(spec,"_GH1.gff3"),format="gff3")
    ngs<-length(GH1s_gff[GH1s_gff$type=="gene"])
    }
    stotal<-c(nps,ngs)
    return(stotal)
    }

specs<-readLines("Genome_specids")
specs_genes_su<-list()
for (spec in specs){ specs_genes_su[[spec]]<- summarynumber(spec)} #return the no. of proteins and no. of genes for all species
summarytable<-t(as.data.frame(specs_genes_su))
colnames(summarytable)<-c("protein","gene")
write.csv(summary)
Allgenes<- as(specs_genes,"GRangesList")

#plot gene number to tree
library(tidyverse)
library(phytools)
numbers<-read.csv("summary_p2g.csv",stringsAsFactor=F,row.names=1)
spectreetxt<-"((Limulus_polyphemus,((Varroa_jacobsoni,(Dermatophagoides_pteronyssinus,Tetranychus_urticae)),Centruroides_sculpturatus)),(Penaeus_vannamei,(Daphnia_magna,(Folsomia_candida,((Blattella_germanica,(Zootermopsis_nevadensis,Cryptotermes_secundus)),((Frankliniella_occidentalis,((Diaphorina_citri,(Aphis_gossypii,Acyrthosiphon_pisum)),(Nilaparvata_lugens,Cimex_lectularius))),(Pediculus_humanus,((Cephus_cinctus,(Orussus_abietinus,((Nasonia_vitripennis,(Trichogramma_pretiosum,Copidosoma_floridanum)),(Polistes_dominula,((Harpegnathos_saltator,(Linepithema_humile,(Camponotus_floridanus,(Pogonomyrmex_barbatus,(Solenopsis_invicta,(Atta_colombica,Acromyrmex_echinatior)))))),(Dufourea_novaeangliae,(Megachile_rotundata,(Habropoda_laboriosa,(Apis_mellifera,Bombus_terrestris))))))))),((Agrilus_planipennis,((Onthophagus_taurus,Nicrophorus_vespilloides),(Tribolium_castaneum,((Anoplophora_glabripennis,Diabrotica_virgifera),(Sitophilus_oryzae,Dendroctonus_ponderosae))))),((Plutella_xylostella,((Papilio_machaon,(Danaus_plexippus,Bicyclus_anynana)),(Hyposmocoma_kahamanoa,((Galleria_mellonella,Amyelois_transitella),((Helicoverpa_armigera,Spodoptera_litura),(Bombyx_mandarina,Manduca_sexta)))))),(((Culex_quinquefasciatus,Aedes_aegypti),Anopheles_gambiae),(Drosophila_melanogaster,(Rhagoletis_zephyria,Musca_domestica)))))))))))));"
spectree<-read.tree(text=spectreetxt)
snumbers<-numbers[rownames(numbers)%in%spectree$tip.label,]
#dotTree(spectree,snumbers[,1:2])
plotTree.barplot(spectree,snumbers[,c(2,1)],args.plotTree=list(cex=0.5),args.barplot=list(col=c("blue","red"),legend.text=c("gene","protein"),main="No. of genes/proteins",args.legend=list(bty="n")))
dev.print(pdf,"summary_p2g.pdf")
dev.off()
```

```r
#############################################
##############    R plot genes    ##############
#############################################
library(genoPlotR) #http://genoplotr.r-forge.r-project.org/vignette.php
library(rtracklayer)
library(GenomicFeatures)
library(Gviz)

grange2dnaseg<-function(x){ #define a function to transform grange object to dna_seg object
    gd<-as.data.frame(granges(x))
    gd<-gd[,-which(names(gd)%in%c("width"))]
    colnames(gd)<-c("name","start","end","strand")
    droplevels(gd)
    levels(gd$strand) <- c(levels(gd$strand), 1,-1)
    gd$strand[gd$strand=="-"]<--1 #can use with function to replace
    #linshi$strand<- with(gd,ifelse(strand=="-",-1,1))
    gd$strand[gd$strand=="+"]<-1
    return(dna_seg(gd))
    }

getannotation<-function(x){
    pl<-split(x,x$protein_id)
    pll<-lapply(pl,function(x)as.data.frame(range(x),stringsAsFactors=F))
    ad<-do.call(rbind.data.frame, pll)
    #data.frame(t(sapply(pll,c)))
    #dcast(melt(pll),L1~variable)
    ad$text<-rownames(ad)
    ad$x1<-ad$start
    ad$x2<-ad$end
    ad$rot<- 15
    anno<-as.annotation(ad[,c("text","x1","x2","rot")])
    return(anno)
    }

getxlim<-function(x){
    pl<-split(x,x$protein_id)
    pll<-lapply(pl,function(x)as.data.frame(range(x),stringsAsFactors=F))
    ad<-do.call(rbind.data.frame, pll)
    #data.frame(t(sapply(pll,c)))
    #dcast(melt(pll),L1~variable)
    ad$text<-rownames(ad)
    ad$x1<-ifelse(ad$start>ad$end,ad$start+50,ad$start-50)
    ad$x2<-ifelse(ad$start>ad$end,ad$end-50,ad$end+50)
    ads<-ad[order(ad$x1),]
    lims=unlist(as.vector(t(unique(ads[,c("x1","x2")]))))
    return(lims)
}

sps=readLines("Genome_specids_selected_tree")

gffs<-import(paste0(sps[38],"_GH1.gff3"),format="gff3")
#seqlevels(gff, pruning.mode="coarse") <- seqlevelsInUse(gffs)
gffss<-gffs[gffs$type=="CDS"]
gffssList<-split(gffss,seqnames(gffss))
#gffssList<-split(gffss,gffss$gene)
dna_segs<-lapply(gffssList,grange2dnaseg)
annotations<-lapply(gffssList,getannotation)
xlims<-lapply(gffssList,getxlim)
plot_gene_map(dna_segs,annotations=annotations, xlims=xlims,scale=TRUE, dna_seg_scale=TRUE)
#,dna_seg_labels=filenames,gene_type=) #produce names as protein ids
dev.print(pdf,paste0(sps[38],"_GH1.geno_xlims.gff3.pdf"))
plot_gene_map(dna_segs,annotations=annotations,scale=TRUE, dna_seg_scale=TRUE)
dev.print(pdf,paste0(sps[8],"_GH1.geno.gff3.pdf"))

options(ucscChromosomeNames=FALSE)
sTxDB<-makeTxDbFromGFF(paste0(sps[54],"_GH1.gff3"),format="gff3")
#seqlevels(sTxDB)
#columns(sTxDB)
#keytypes(sTxDB)
tracks<-GeneRegionTrack(sTxDB)
chrosomen<- seqlevels(tracks)
scTrack<-GenomeAxisTrack(scale=500,labelPos="below",exponent=3) #setup a genome scale for each track
grid.newpage()
pushViewport(viewport(layout=grid.layout(length(chrosomen),1)))
for (i in seq_along(chrosomen)){
pushViewport(viewport(layout.pos.col=1, layout.pos.row=i))
plotTracks(list(tracks,scTrack),chromosome=chrosomen[i],groupAnnotation="id",transcriptAnnotation="exon",add=TRUE,main=chrosomen[i],cex.main=0.8)
popViewport(1)
}
dev.print(pdf,paste0(sps[54],"_GH1.gff3.pdf"))


###################################################################################
###############    get representing protein seqs and cds seqs for each gene    ##############
###################################################################################

"""r script"""
library(rtracklayer)
library(tidyverse)
genesummary<-read.csv("all_genome_GH1.gene2protein.summary.csv",header=T,sep=",",stringsAsFactors=F)
sgenesummary<- genesummary %>% group_by(group_name) %>% filter(wcds==max(wcds)) %>% top_n(1) %>% ungroup() %>% as.data.frame
specs<-readLines("Genome_specids")
genome_folder_path<-"/media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/"

get_protein_cds<-function(spec){ #define function to extract and output gene gffs,represent_proteins
    protein_fn<-list.files(paste0(genome_folder_path,spec),pattern = "\\.faa$")
    nuc_fn<-list.files(paste0(genome_folder_path,spec),pattern = "\\.fna$")
    protein_fp<-file.path(genome_folder_path,spec,protein_fn)
    nuc_fp<-file.path(genome_folder_path,spec,nuc_fn)
    specsummary=sgenesummary[sgenesummary$Org==spec,]
    proseqs=readAAStringSet(protein_fp)
    profas=AAStringSet()
    nucseqs=readDNAStringSet(nuc_fp)
    genefas=DNAStringSet()
    for (sid in specsummary$proteinname){
        sidq=nucseqs[grepl(sid,names(nucseqs))]
        names(sidq)=sid
        genefas=c(genefas,sidq)
        pidp=proseqs[grepl(sid,names(proseqs))]
        profas=c(profas,pidp)
        if (len(sidp)%%3){
            print(paste(spec,sid))
            }
    fasfilename <- paste0(spec,"_GH1s.fna")
    profilename <- paste0(spec,"_GH1s.faa")
    writeXStringSet(genefas,fasfilename,format='fasta')
    writeXStringSet(profas,profilename,format='fasta')
    }
}



lapply(specs, get_protein_cds)
```

### 5. run taxonomy text for all proteins by top 10 blastp search from ncbi

```bash
for i in {00..09}; do ~/GH1/taxonomy_contamination_detect_transcriptome.py -i GH1s_genome_protein.$i.faa.alignment.xml -o GH1s_genome_protein.$i.faa.alignment.xml.tax; done;cat *xml.tax > GH1s_genome_protein.faa.alignment.xml.tax; rm GH1s_genome_protein.*.faa.alignment.xml.tax
#manully check the sequences from bacteria.
```

Becareful check the bacterial origin proteins, should check the flanking genes around the protein.

### 6.CAFE analysis (not suitable because of not for genes but for multiple gene families) and ANOVA test predict gene numbers with feeding groups or orders

__preparing tree and analysis__

```r
#1. Download the time phylogeny for part of the species.[from the supplementary file in the paper -Gene content evolution in the arthropods-](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1925-7)

library(phytools)
library(ape)
treefrompaper<- "((((((Ladona_fulva:355.768668,Ephemera_danica:355.768668):71.765861,((Blattella_germanica:148.521245,Zootermopsis_nevadensis:148.521245):262.478755,((((((((Copidosoma_floridanum:107.287992,Trichogramma_pretiosum:107.287992):7.283957,Nasonia_vitripennis:114.571949):66.550479,(((Dufourea_novaeangliae:45.641976,Lasioglossum_albipes:45.641976):16.236509,(((((Apis_mellifera:8.103582,Apis_florea:8.103582):20.296418,((Bombus_impatiens:6.188111,Bombus_terrestris:6.188111):17.229702,Melipona_quadrifasciata:23.417813):4.982187):5.656707,Eufriesea_mexicana:34.056707):10.448774,Habropoda_laboriosa:44.505481):9.118989,Megachile_rotundata:53.624470):8.254016):32.021514,(((((((Atta_colombica:9.937519,Acromyrmex_echinatior:9.937519):23.084663,Solenopsis_invicta:33.022181):.740737,Cardiocondyla_obscurior:35.762918):4.330127,Pogonomyrmex_barbatus:40.093046):9.270144,Camponotus_floridanus:49.363190):3.121819,Linepithema_humile:52.485009):11.488588,Harpegnathos_saltator:63.973597):29.926403):87.222428):8.778331,Orussus_abietinus:189.900759):8.601893,Cephus_cinctus:198.502651):28.694151,Athalia_rosae:227.196802):163.862573,((((((Anoplophora_glabripennis:104.531206,Leptinotarsa_decemlineata:104.531206):32.621188,Dendroctonus_ponderosae:137.152394):27.716131,Tribolium_castaneum:164.868525):51.492868,Onthophagus_taurus:216.361394):26.273864,Agrilus_planipennis:242.635258):133.806773,(((((Bombyx_mandarina:80.077648,Manduca_sexta:80.077648):26.739891,(Heliconius_melpomene:73.030688,Danaus_plexippus:73.030688):33.786851):28.317991,Plutella_xylostella:135.135530):145.951499,Limnephilus_lunatus:281.087029):80.418316,(((Aedes_aegypti:83.961536,Culex_quinquefasciatus:83.961536):80.113960,((Anopheles_gambiae:39.856741,Anopheles_funestus:39.856741):48.584539,Anopheles_albimanus:88.441279):75.634217):146.794275,((((((Lucilia_cuprina:74.112753,Musca_domestica:74.112753):29.696044,Glossina_morsitans:103.808797):37.597288,Ceratitis_capitata:141.406085):16.264027,((Drosophila_pseudoobscura:49.846935,Drosophila_melanogaster:49.846935):25.133061,Drosophila_grimshawi:74.979997):82.690116):120.951839,Mayetiola_destructor:278.621951):19.431902,Lutzomyia_longipalpis:298.053853):12.815918):50.635574):14.936686):14.617344):14.029973,(((((((Halyomorpha_halys:108.538881,Oncopeltus_fasciatus:108.538881):71.274364,Cimex_lectularius:179.813245):47.917855,Gerris_buenoi:227.731100):77.820934,Homalodisca_vitripennis:305.552034):35.712699,(Acyrthosiphon_pisum:306.647247,Pachypsylla_venusta:306.647247):34.617486):30.725536,Frankliniella_occidentalis:371.990270):16.709075,Pediculus_humanus:388.699345):16.390004):5.910651):16.534529):59.131713,Catajapyx_aquilonaris:486.666242):35.725916,((Hyalella_azteca:487.000000,Eurytemora_affinis:487.000000):19.394942,Daphnia_magna:506.394942):15.997216):45.942043,Strigamia_maritima:568.334201):2.019415,((((((Latrodectus_hesperus:86.575583,Parasteatoda_tepidariorum:86.575583):52.883152,Stegodyphus_mimosarum:139.458735):112.347677,Loxosceles_reclusa:251.806412):146.585282,Centruroides_sculpturatus:398.391694):71.123177,(Metaseiulus_occidentalis:391.827067,Ixodes_scapularis:391.827067):77.687804):27.994573,Tetranychus_urticae:497.509444):72.844172);"
spectreetxt<-"((Limulus_polyphemus,((Varroa_jacobsoni,(Dermatophagoides_pteronyssinus,Tetranychus_urticae)),Centruroides_sculpturatus)),(Penaeus_vannamei,(Daphnia_magna,(Folsomia_candida,((Blattella_germanica,(Zootermopsis_nevadensis,Cryptotermes_secundus)),((Frankliniella_occidentalis,((Diaphorina_citri,(Aphis_gossypii,Acyrthosiphon_pisum)),(Nilaparvata_lugens,Cimex_lectularius))),(Pediculus_humanus,((Cephus_cinctus,(Orussus_abietinus,((Nasonia_vitripennis,(Trichogramma_pretiosum,Copidosoma_floridanum)),(Polistes_dominula,((Harpegnathos_saltator,(Linepithema_humile,(Camponotus_floridanus,(Pogonomyrmex_barbatus,(Solenopsis_invicta,(Atta_colombica,Acromyrmex_echinatior)))))),(Dufourea_novaeangliae,(Megachile_rotundata,(Habropoda_laboriosa,(Apis_mellifera,Bombus_terrestris))))))))),((Agrilus_planipennis,((Onthophagus_taurus,Nicrophorus_vespilloides),(Tribolium_castaneum,((Anoplophora_glabripennis,Diabrotica_virgifera),(Sitophilus_oryzae,Dendroctonus_ponderosae))))),((Plutella_xylostella,((Papilio_machaon,(Danaus_plexippus,Bicyclus_anynana)),(Hyposmocoma_kahamanoa,((Galleria_mellonella,Amyelois_transitella),((Helicoverpa_armigera,Spodoptera_litura),(Bombyx_mandarina,Manduca_sexta)))))),(((Culex_quinquefasciatus,Aedes_aegypti),Anopheles_gambiae),(Drosophila_melanogaster,(Rhagoletis_zephyria,Musca_domestica)))))))))))));"

tree<-read.tree(text=treefrompaper)
tree2<-read.tree(text=spectreetxt)
bc<-c("Trichogramma_pretiosum","Linepithema_humile","Camponotus_floridanus","Pogonomyrmex_barbatus","Solenopsis_invicta","Megachile_rotundata","Habropoda_laboriosa")
tree2<-drop.tip(tree2,bc)
tips<-tree$tip.label[!tree$tip.label%in%tree2$tip.label]
bc<-append(tips,c("Centruroides_sculpturatus","Tetranychus_urticae"))
tree<-drop.tip(tree,bc)
write.tree(tree, "Species_dated_tree.newick")

tree<-read.tree(text=treefrompaper)
tree2<-read.tree(text=spectreetxt)
bc<-c("Trichogramma_pretiosum","Linepithema_humile","Camponotus_floridanus","Pogonomyrmex_barbatus","Solenopsis_invicta","Megachile_rotundata","Habropoda_laboriosa")
tree2<-drop.tip(tree2,bc)
tips<-tree$tip.label[!tree$tip.label%in%tree2$tip.label]
bc<-append(tips,c("Centruroides_sculpturatus","Tetranychus_urticae"))
tree<-drop.tip(tree,bc)
write.tree(tree, "Species_dated_tree.newick")
#edit gene numbers
numbers<-read.csv("summary_p2g.csv",stringsAsFactor=F,row.names=1)
snumbers<-numbers[rownames(numbers)%in%tree$tip.label,]
write.(t(snumbers),"gene_numbers.txt")
#edit gene_numbers.txt table to suit cafe analysis.
```

__run cafe analysis, multiple times to check convergency. CAFE is suitable for analysis a large number of families not a single families. And it does not count the gene phylogeny.__

```bash
for i in {1..5}; do cafexp -i gene_numbers.txt -t Species_dated_tree.newick -o result$i; done
```

__compare genes with feeding groups and orders__

```r
#1-separate feeding groups as Figure 1;2-separate feeding groups as their plant cell wall feeding or not.
tr<-read.tree("Species_dated_tree.newick")
tree<-drop.tip(tr,"Daphnia_magna")
genet<-read.csv("gene_n_type_annotated.csv",header=T)
geneto<-genet[match(tree$tip.label,genet$Species_name),]
#test between herbivores and non-herbivores
phylANOVA(tree,geneto$Herbivorous,geneto$GeneN,nsim=1000,posthoc=T,p.adj="holm")

#test separation groups of direct feeding and non-direct feeding on plant celll walls.
phylANOVA(tree,geneto$Plant_cell_wall,geneto$GeneN,nsim=1000,posthoc=T,p.adj="holm")
```

### 6. Dual domain proteins

__#extract protein with multiple domains with "python"__

```python
import pandas as pd
from Bio import SeqIO
import os
with open("Genome_specids_selected_tree",'rt') as file:
    specs=[c.strip('\n') for c in file.readlines()]

def retrieve_dual_domain_protein(spec):
    GH1o=pd.read_table("./"+spec+"/"+spec+"overview.txt")
    GH1os=GH1o[GH1o['#ofTools']>2]
    if any(GH1os['HMMER'].str.contains("\+")):# the "+" means another domains.
        pids=GH1os['Gene ID'].values.tolist()
        genome_folder_path="/media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/"
        p_fn=[file for file in os.listdir(genome_folder_path+spec) if file.endswith("_protein.faa")]
        p_fp=os.path.join(genome_folder_path,spec,p_fn[0])
        pseq=[seq for seq in SeqIO.parse(p_fp,'fasta') if seq.id in pids]
        cds_fn=[file for file in os.listdir(genome_folder_path+spec) if file.endswith(".fna")]
        cds_fp=os.path.join(genome_folder_path,spec,cds_fn[0])
        cds=[seq for seq in SeqIO.parse(cds_fp,'fasta') if seq.id.split("_cds_")[-1].rsplit("_",1)[0] in pids]
        for seq in cds:
            seq.id=seq.id.split("_cds_")[-1].rsplit("_",1)[0]
        return (pseq,cds)

allseqs={spec:retrieve_dual_domain_protein(spec) for spec in specs if retrieve_dual_domain_protein(spec)!=None }
os.mkdir("Dual_domain_protein")
os.chdir("Dual_domain_protein")
for sp,se in allseqs.items():
    SeqIO.write(se[0],sp+".faa",'fasta')
    SeqIO.write(se[1],sp+".fna",'fasta')
```

__extract protein with multiple domains with R__

```r
library(tidyverse)
library(Biostrings)
specs<-readLines("Genome_specids_selected_tree")
for(spec in specs){
GH1os<-read.delim(paste0(spec,"/",spec,"overview.txt"),stringsAsFactors=F) %>% filter(.,X.ofTools>2)%>%mutate(Org=spec)
specre<- if (any(grepl("\\+",GH1os$HMMER))) GH1os else GH1os[0,] # the "+" means another domains.
#genome_folder_path<-"/media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/"
#protein_fn<-list.files(paste0(genome_folder_path,spec),pattern = "\\.faa$")
#protein_fp<-file.path(genome_folder_path,spec,protein_fn)
#specreaa<- readAAStringSet(protein_fp,format="fasta")
#GH1h<-read.delim(paste0(spec,"/",spec,"hmmer.out"),stringsAsFactors=F) %>% filter(., Gene.ID%in%GH1os$Gene.ID)
#specre<-dplyr::filter(GH1,GH1$Gene.ID%in%GH1$Gene.ID[duplicated(GH1$Gene.ID)])%>%mutate(Org=spec)
if (!exists('specsre')) {
    specsre<-specre
    } else {
    specsre<- bind_rows(specsre,specre)
    }
}
write.table(specsre,"./Dual_domain_protein/Dual_GH1_domain_gene.csv",sep="\t",quote=F)
```

### 7. TREE BUILDING for whole GH1s from genomes

```bash
#alignment with AQUA
AQUA.tcl GH1s_genome_protein.faa GH1_genome_AQUA
#build tree with iqtree
iqtree -s GH1s_genome_protein.faa.mafft.rascal -nt 10 -bb 1000 -alrt 1000
```

__format with python to phylip for aligment to use exbasyes__ __*not used*__

```python
import re
import pandas as pd
from Bio import AlignIO
ali=AlignIO.read("GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal",'fasta')
fsp=re.compile(r"\d+\.\d+\_(\w*)\.*\_\w*$")
i=0
ids2fullids={}
for seq in ali:
    i+=1
    sids=fsp.findall(seq.id)[0]
    sid=sids[:7]+f"{i:03d}" if len(sids)>7 else sids+f"{i:03d}"
    ids2fullids[sid]=seq.id
    seq.id=sid
AlignIO.write(ali,"GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.phylip",'phylip')
#df=pd.DataFrame.from_dict(ids2fullids,orient="index")
#df.to_csv("GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.id2phylipid")
pd.Series(ids2fullids).to_csv('GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.id2phylipid',sep="\t",header=False)
```

#run exabayes with 2 runs 4 chains, 0.01 sdsf, LG protein model.

> `mpirun -np 8 exabayes -f GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.phylip -m PROT -r GH1sexabayes -s 258 -c config.nex`

__to retrieve tips *R*__

```r
library(phytools)
library(stringi)
library(rtracklayer)
library(GenomicFeatures)
library(Gviz)
library(genoPlotR) #http://genoplotr.r-forge.r-project.org/vignette.php

tree=read.tree("GH1s_genome_protein.faa.mafft.rascal.contree")
extract_node_id<-fastMRCA(tree, "XP_026292938.1_Frankliniella_o","XP_026273946.1_Frankliniella_o")
extract_clade_tree<- extract.clade(tree, extract_node_id)
extract_seqs_ids<-unlist(lapply(extract_clade_tree$tip.label,function(x) str_extract(x,"^\\w+_\\d+\\.\\d+")))
extract_clade_specs<-unique(unlist(lapply(extract_clade_tree$tip.label,function(x) str_extract(x,"[:alpha:]+_[:alpha:]+$"))))
fp<-"../GH1_from_dbCAN/run_dbCAN"
gff_fn<-list.files(fp,pattern=paste0(extract_clade_specs[0],".*gff3")
gff_fp<-file.path(fp,gff_fn)
options(ucscChromosomeNames=FALSE)
sTxDB<-makeTxDbFromGFF(gff_fp,format="gff3")
cdsset<-cdsBy(sTxDB,by='gene')
scdsset<-cdsset[sapply(cdsset,function(x)any(x$cds_name%in%extract_seqs_ids))]

cdsgrange2dnaseg<-function(x){ #define a function to transform grange object to dna_seg object
    gd<-as.data.frame(granges(x))
    gd<-gd[,-which(names(gd)%in%c("width"))]
    colnames(gd)<-c("name","start","end","strand")
    droplevels(gd)
    levels(gd$strand) <- c(levels(gd$strand), 1,-1)
    gd$strand[gd$strand=="-"]<--1 #can use with function to replace
    #linshi$strand<- with(gd,ifelse(strand=="-",-1,1))
    gd$strand[gd$strand=="+"]<-1
    return(dna_seg(gd))
    }

dna_segs<-lapply(scdsset,cdsgrange2dnaseg)
plot_gene_map(dna_segs,scale=TRUE, dna_seg_scale=TRUE,dna_seg_labels=names(scdsset),gene_type=)
```

__to retrieve tips *Python*__

```python
from ete3 import Tree
from Bio import SeqIO,Seq
GH1t=Tree("GH1s_genome_protein.faa.mafft.rascal.contree")
detach_node=GH1t.get_common_ancestor(["XP_026292938.1_Frankliniella_o","XP_026273946.1_Frankliniella_o"])
#detach_node.detach() #detach
#GH1t_a_d_names=GH1t.get_leaf_names() 
#GH1t.prune([detach_node]) #prune
tids=[i.rsplit('_',2)[0] for i in detach_node.get_leaf_names()]
tseqs=[seq for seq in SeqIO.parse("/home/shulinhe/GH1/GH1_from_dbCAN/run_dbCAN/GH1s_genome_protein.faa",'fasta') if seq.id.rsplit('_',2)[0] in tids]
SeqIO.write(tseqs,"test.fasta",'fasta')
```

### 8.NOTUNG reconsilation of gene tree from species tree

__both tree should be rooted__
__prepare tree tips names__

```r
library(ape)
library(tidyverse)
genetree=read.tree("GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree")
tipnames<-as.data.frame(sapply(genetree$tip.label,function(x) str_replace(x,"_[:alpha:]+$","")))
colnames(tipnames)<-"V1"
tipnames$V1<-str_replace(tipnames$V1,"humanus","Pediculus") %>%str_replace(.,"virgifera","Diabrotica") %>% str_replace(.,"plexippus","Danaus") %>%str_replace(.,"str.","Anopheles")
genetree$tip.label<-tipnames[["V1"]][match(rownames(tipnames),genetree$tip.label)]
write.tree(genetree,"GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified") ##root tree with dendroscope as well
spectreetxt<-"((Limulus,((Varroa,(Dermatophagoides,Tetranychus)),Centruroides)),(Penaeus,(Daphnia,(Folsomia,((Blattella,(Zootermopsis,Cryptotermes)),((Frankliniella,((Diaphorina,(Aphis,Acyrthosiphon)),(Nilaparvata,Cimex))),(Pediculus,((Cephus,(Orussus,((Nasonia,Copidosoma),(Polistes,((Harpegnathos,(Atta,Acromyrmex)),(Dufourea,(Apis,Bombus))))))),((Agrilus,((Onthophagus,Nicrophorus),(Tribolium,((Anoplophora,Diabrotica),(Sitophilus,Dendroctonus))))),((Plutella,((Papilio,(Danaus,Bicyclus)),(Hyposmocoma,((Galleria,Amyelois),((Helicoverpa,Spodoptera),(Bombyx,Manduca)))))),(((Culex,Aedes),Anopheles),(Drosophila,(Rhagoletis,Musca)))))))))))));"
spectree<-read.tree(text=spectreetxt)
write.tree(spectree,"GH1s_genome_species.tree.modified")
```

__run notung__

```bash
mv GH1s_genome_species.tree.modified GH1s_genome_species.tree.modified.nwk
mv GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified.rooted GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified.rooted.nws
java -jar ~/opt/Notung/Notung-2.9.1.5.jar -g GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified.rooted -s GH1s_genome_species.tree.modified --rearrange --speciestag postfix --threshold 85% --bootstraps name --outputdir GH1_reconciletest --log --events --parsable --treestats --progressbar --savepng --saveweakedgespng --homologtablecsv
```

__correlation test with duplication/loss with BUSCO scores__

```r
library(ggpubr)
abc=read.csv("busco_correlation_test.csv",header=T,sep="\t",stringsAsFactors=F)
ggscatter(abc[abc$Color!="",], x="M", y="Loss", add="reg.line",conf.int=T, cor.coef=T,cor.method="pearson",xlab="BUSCO Missing", ylab="Gene Loss")
dev.print(pdf,"Missing2Loss.pdf")
ggscatter(abc[abc$Color!="",], x="D", y="Dup", add="reg.line",conf.int=T, cor.coef=T,cor.method="pearson",xlab="BUSCO Duplication", ylab="Gene Duplicate")
dev.print(pdf,"Dup2Duplication.pdf")
```

### 9. get cds sequences from database based on protein sequences

```python
from Bio import SeqIO
import os
import re

def getsequences(pseq,cdss):
    if len(cdss)%3==1:
        print(cdss.id,"has",1,"more nucleotides in cds.")
        c=cdss[:-1]
    elif len(cdss)%3==2:
        print(cdss.id,"has",2,"more nucleotides in cds. Removed")
        c=cdss[:-2]
    else:
        c=cdss
    if len(pseq)*3>len(c):
        print(pseq.id,"has more aa than cds")


fid=re.compile(r"\_cds\_(\w*\.\d+)\_\d+$")

def retrieve_pcds(specful,ids):
    spec="_".join(specful.split(" ")[:2])
    genome_folder_path="/media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/"
    p_fn=[file for file in os.listdir(genome_folder_path+spec) if file.endswith("_protein.faa")]
    p_fp=os.path.join(genome_folder_path,spec,p_fn[0])
    psseqs=[seq for seq in SeqIO.parse(p_fp,'fasta') if seq.id in ids]
    cds_fn=[file for file in os.listdir(genome_folder_path+spec) if file.endswith("genomic.fna")]
    cds_fp=os.path.join(genome_folder_path,spec,cds_fn[0])
    cdsseqs=[seq for seq in SeqIO.parse(cds_fp,'fasta') if fid.findall(seq.id)[0] in ids]
    for seq in cdsseqs:
        seq.id=fid.findall(seq.id)[0]
    for pid in ids:
        pse=[seq for seq in psseqs if seq.id==pid][0]
        cse=[seq for seq in cdsseqs if seq.id==pid][0]
        getsequences(pse,cse)
    for seq in cdsseqs:
        seq.id=seq.id+"_"+spec
    for seq in psseqs:
        seq.id=seq.id+"_"+spec
    return psseqs,cdsseqs

def groupseqs(spec,seqs):    #retrieve protein ids list for spec
    pids=[seq.description.split(" ")[1] for seq in seqs if spec in seq.description]
    return pids

pseqs=[a for a in SeqIO.parse("GH1s.genome.protein.faa",'fasta')] #get all proteins
specs=set([a.description.split("[")[-1].strip("]") for a in pseqs]) #get species names
spec2pseqi={spec:groupseqs(spec,pseqs) for spec in specs} #species names 2 protein id list
os.mkdir("GH1s_cds")
os.chdir("GH1s_cds")
for specful,idl in spec2pseqi.items():
    speci="_".join(specful.split(" ")[:2])
    print(speci)
    try:
        proteinseq,cdsseq=retrieve_pcds(specful,idl)
        SeqIO.write(proteinseq,speci+"_genome_protein_GH1.faa",'fasta')
        SeqIO.write(cdsseq,speci+"_genome_cds_GH1.fna",'fasta')
    except AttributeError:
        print("Error in processing",speci)

os.chdir("../")     
#manually check to remove the additional nucleotides and aa in files###
############### Modified_cds_note #####################################
```

### 10.Tree building according to orders based on the edited cds/protein files

```python
from Bio import SeqIO
import re
"""#older version
allseq=[seq for seq in SeqIO.parse("GH1s_genome_protein.faa",'fasta')]
order2spec={"Chelicearata":["Limulus_polyphemus","Varroa_jacobsoni","Dermatophagoides_pteronyssinus","Tetranychus_urticae","Centruroides_sculpturatus"],
"Crustaceans":["Penaeus_vannamei","Daphnia_magna"],
"Springtails":["Folsomia_candida"],
"Cockroaches":["germanica","Zootermopsis_nevadensis","Cryptotermes_secundus"],
"Thrips":["Frankliniella_occidentalis"],
"Hemiptera":["Diaphorina_citri","Aphis_gossypii","Acyrthosiphon_pisum","Nilaparvata_lugens","Cimex_lectularius"],
"Psocodea":["humanus_corporis"],
"Hymenoptera":["Cephus_cinctus","Orussus_abietinus","Nasonia_vitripennis","Copidosoma_floridanum","Polistes_dominula","Harpegnathos_saltator","Atta_colombica","Acromyrmex_echinatior","Dufourea_novaeangliae","Apis_mellifera","Bombus_terrestris"],
"Bettles":["Agrilus_planipennis","Onthophagus_taurus","Nicrophorus_vespilloides","Tribolium_castaneum","Anoplophora_glabripennis","virgifera_virgifera","Sitophilus_oryzae","Dendroctonus_ponderosae"],
"Lepitoptera":["Plutella_xylostella","Papilio_machaon","plexippus_plexippus","Bicyclus_anynana","Hyposmocoma_kahamanoa","Galleria_mellonella","Amyelois_transitella","Helicoverpa_armigera","Spodoptera_litura","Bombyx_mandarina","Manduca_sexta"],
"Diptera":["Culex_quinquefasciatus","Aedes_aegypti","str._PEST","Drosophila_melanogaster","Rhagoletis_zephyria","Musca_domestica"]}
for a,b in order2spec.items():
    seqs=[seq for seq in allseq if seq.id.split('_',2)[-1] in b]
    SeqIO.write(seqs,a+"_genome_protein.faa",'fasta')
"""#end of old version

order2spec={"Chelicearata":["Limulus_polyphemus","Varroa_jacobsoni","Dermatophagoides_pteronyssinus","Tetranychus_urticae","Centruroides_sculpturatus"],
"Crustaceans":["Penaeus_vannamei","Daphnia_magna"],
"Springtails":["Folsomia_candida"],
"Cockroaches":["Blattella_germanica","Zootermopsis_nevadensis","Cryptotermes_secundus"],
"Thrips":["Frankliniella_occidentalis"],
"Hemiptera":["Diaphorina_citri","Aphis_gossypii","Acyrthosiphon_pisum","Nilaparvata_lugens","Cimex_lectularius"],
"Psocodea":["Pediculus_humanus"],
"Hymenoptera":["Cephus_cinctus","Orussus_abietinus","Nasonia_vitripennis","Copidosoma_floridanum","Polistes_dominula","Harpegnathos_saltator","Atta_colombica","Acromyrmex_echinatior","Dufourea_novaeangliae","Apis_mellifera","Bombus_terrestris"],
"Bettles":["Agrilus_planipennis","Onthophagus_taurus","Nicrophorus_vespilloides","Tribolium_castaneum","Anoplophora_glabripennis","Diabrotica_virgifera","Sitophilus_oryzae","Dendroctonus_ponderosae"],
"Lepitoptera":["Plutella_xylostella","Papilio_machaon","Danaus_plexippus","Bicyclus_anynana","Hyposmocoma_kahamanoa","Galleria_mellonella","Amyelois_transitella","Helicoverpa_armigera","Spodoptera_litura","Bombyx_mandarina","Manduca_sexta"],
"Diptera":["Culex_quinquefasciatus","Aedes_aegypti","Anopheles_gambiae","Drosophila_melanogaster","Rhagoletis_zephyria","Musca_domestica"]}
for a,b in order2spec.items():
    try:
        sp_p_seqs=[seq for spec in b for seq in SeqIO.parse(spec+"_genome_protein_GH1.faa",'fasta')]
        sp_cds_seqs=[seq for spec in b for seq in SeqIO.parse(spec+"_genome_cds_GH1.fna",'fasta')]
        SeqIO.write(sp_p_seqs,a+"_genome_protein_GH1e.faa",'fasta')
        SeqIO.write(sp_cds_seqs,a+"_genome_cds_GH1e.fna",'fasta')
    except FileNotFoundError as err:
        print(err)
```

```bash
for i in {"Bettles","Cockroaches","Diptera","Hemiptera","Lepitoptera"}; do \
AQUA.tcl GH1s_genome_protein.faa GH1_genome_AQUA; \

#build tree with iqtree
iqtree -s GH1s_genome_protein.faa.mafft.rascal -nt 10 -bb 1000 -alrt 1000
```

```r
library(ape)
library(phytools)
library(tidyverse)
spectreetxt<-"((Limulus,((Varroa,(Dermatophagoides,Tetranychus)),Centruroides)),(Penaeus,(Daphnia,(Folsomia,((Blattella,(Zootermopsis,Cryptotermes)),((Frankliniella,((Diaphorina,(Aphis,Acyrthosiphon)),(Nilaparvata,Cimex))),(Pediculus,((Cephus,(Orussus,((Nasonia,(Trichogramma,Copidosoma)),(Polistes,((Harpegnathos,(Linepithema,(Camponotus,(Pogonomyrmex,(Solenopsis,(Atta,Acromyrmex)))))),(Dufourea,(Megachile,(Habropoda,(Apis,Bombus))))))))),((Agrilus,((Onthophagus,Nicrophorus),(Tribolium,((Anoplophora,Diabrotica),(Sitophilus,Dendroctonus))))),((Plutella,((Papilio,(Danaus,Bicyclus)),(Hyposmocoma,((Galleria,Amyelois),((Helicoverpa,Spodoptera),(Bombyx,Manduca)))))),(((Culex,Aedes),Anopheles),(Drosophila,(Rhagoletis,Musca)))))))))))));"
spectree<-read.tree(text=spectreetxt)

Bettles_tree<- extract.clade(spectree, fastMRCA(spectree, "Agrilus","Dendroctonus"))
write.tree(Bettles_tree,"Bettles_species.tree.modified")
genetree=read.tree("Bettles_genome_protein.faa.muscle.rascal.contree")
tipnames<-as.data.frame(sapply(genetree$tip.label,function(x) str_replace(x,"_[:alpha:]+$","")))
colnames(tipnames)<-"V1"
#tipnames<- rownames_to_column(tipnames,var="rowname") %>%mutate(V2=str_replace_all(V1, "virgifera$", "Diabrotica")) %>% column_to_rownames(.,var="rowname") %>% select(.,"V2")
genetree$tip.label<-tipnames[["V1"]][match(rownames(tipnames),genetree$tip.label)]
write.tree(genetree,"Bettles_genome_protein.faa.muscle.rascal.contree.modified") ##root tree with dendroscope as well

Lepitoptera_tree<-extract.clade(spectree, fastMRCA(spectree, "Plutella","Manduca"))
write.tree(Lepitoptera_tree,"Lepitoptera_species.tree.modified")
genetree=read.tree("Lepitoptera_genome_protein.faa.mafft.contree")
tipnames<-as.data.frame(sapply(genetree$tip.label,function(x) str_replace(x,"_[:alpha:]+$","")))
colnames(tipnames)<-"V1"
#tipnames<- rownames_to_column(tipnames,var="rowname") %>%mutate(V2=str_replace_all(V1, "plexippus$", "Danaus")) %>% column_to_rownames(.,var="rowname") %>% select(.,"V2")
genetree$tip.label<-tipnames[["V1"]][match(rownames(tipnames),genetree$tip.label)]
write.tree(genetree,"Lepitoptera_genome_protein.faa.mafft.contree.modified") ##root tree with

Cockroaches_tree<- extract.clade(spectree, fastMRCA(spectree, "Blattella","Zootermopsis"))
write.tree(Cockroaches_tree,"Cockroaches_species.tree.modified")
genetree=read.tree("Cockroaches_genome_protein.faa.muscle.rascal.contree")
tipnames<-as.data.frame(sapply(genetree$tip.label,function(x) str_replace(x,"_[:alpha:]+$","")))
colnames(tipnames)<-"V1"
genetree$tip.label<-tipnames[["V1"]][match(rownames(tipnames),genetree$tip.label)]
write.tree(genetree,"Cockroaches_genome_protein.faa.mafft.contree.modified") ##root tree with

Hymenoptera_tree<- extract.clade(spectree, fastMRCA(spectree, "Cephus","Bombus"))
Hymenoptera_trees<- drop.tip(Hymenoptera_tree,"Trichogramma")
write.tree(Hymenoptera_trees,"Hymenoptera_species.tree.modified")
genetree=read.tree("Hymenoptera_genome_protein.faa.muscle.contree")
tipnames<-as.data.frame(sapply(genetree$tip.label,function(x) str_replace(x,"_[:alpha:]+$","")))
colnames(tipnames)<-"V1"
genetree$tip.label<-tipnames[["V1"]][match(rownames(tipnames),genetree$tip.label)]
write.tree(genetree,"Hymenoptera_genome_protein.faa.mafft.contree.modified") ##root tree with

Hemiptera_tree<-extract.clade(spectree, fastMRCA(spectree, "Diaphorina","Cimex"))
write.tree(Hemiptera_tree,"Hemiptera_species.tree.modified")
genetree=read.tree("Hemiptera_genome_protein.faa.muscle.contree")
tipnames<-as.data.frame(sapply(genetree$tip.label,function(x) str_replace(x,"_[:alpha:]+$","")))
colnames(tipnames)<-"V1"
genetree$tip.label<-tipnames[["V1"]][match(rownames(tipnames),genetree$tip.label)]
write.tree(genetree,"Hemiptera_genome_protein.faa.muscle.contree.modified") ##root tree with

Diptera_tree<-extract.clade(spectree, fastMRCA(spectree, "Culex","Musca"))
write.tree(Diptera_tree,"Diptera_species.tree.modified")
genetree=read.tree("Diptera_genome_protein.faa.muscle.rascal.contree")
tipnames<-as.data.frame(sapply(genetree$tip.label,function(x) str_replace(x,"_[:alpha:]+$","")))
colnames(tipnames)<-"V1"
tipnames<- rownames_to_column(tipnames,var="rowname") %>%mutate(V2=str_replace_all(V1, "str.$", "Anopheles")) %>% column_to_rownames(.,var="rowname") %>% select(.,"V2")
genetree$tip.label<-tipnames[["V2"]][match(rownames(tipnames),genetree$tip.label)]
write.tree(genetree,"Diptera_genome_protein.faa.muscle.rascal.contree.modified") ##root tree with
```

__gene conversion test__

```bash
~/opt/pal2nal.v14/pal2nal.pl ./GH1inorder/Tree_bestalin/Bettles_genome_protein_GH1e.faa.muscle.edi Bettles_genome_cds_GH1e.fna -output fasta >Bettles_codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl ./GH1inorder/Tree_bestalin/Lepitoptera_genome_protein_GH1e.faa.muscle.edi ./Lepitoptera_genome_cds_GH1e.fna -output fasta >Lepitoptera_codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl Hemiptera_genome_protein_GH1e.faa.mafft Hemiptera_genome_cds_GH1e.fna -output fasta >Hemiptera_codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl Cockroaches_genome_protein_GH1e.faa.muscle.rascal Cockroaches_genome_cds_GH1e.fna -output fasta >Cockroaches_codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl Hymenoptera_genome_protein_GH1e.faa.muscle Hymenoptera_genome_cds_GH1e.fna -output fasta >Hymenoptera__codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl Diptera_genome_protein_GH1e.faa.muscle.rascal.edi Diptera_genome_cds_GH1e.fna -output fasta >Diptera_codon.aln.fa
```

### 11. GH1 positive selection test

`Cockroaches_genome_protein_GH1e.faa Hymenoptera_genome_protein_GH1e.faa Hemiptera_genome_protein_GH1e.faa Bettles_genome_protein_GH1e.faa Diptera_genome_protein_GH1e.faa Lepitoptera_genome_protein_GH1e.faa Folsomia_candida_genome_protein_GH1.faa Pediculus_humanus_genome_protein_GH1.faa Daphnia_magna_genome_protein_GH1.faa > All_genome_protein_GH1e.faa`
`cat Cockroaches_genome_cds_GH1e.fna Hymenoptera_genome_cds_GH1e.fna Hemiptera_genome_cds_GH1e.fna Bettles_genome_cds_GH1e.fna Lepitoptera_genome_cds_GH1e.fna Diptera_genome_cds_GH1e.fna Folsomia_candida_genome_cds_GH1.fna Pediculus_humanus_genome_cds_GH1.fna Daphnia_magna_genome_cds_GH1.fna > All_genome_cds_GH1e.fna`
__get modified gene tree__

```python
from ete3 import Tree
from Bio import SeqIO,Seq
import re
import os
fsp1=re.compile(r"(\d+[\.]*\d+)\_(\w+)\_*\w*$")
fsp2=re.compile(r"(\d+[\.]*\d+)\_(\w+)\_\w*$")
fsp3=re.compile(r"\d+\_\w+\_\w+$")
def get_phylip_name(id):
    if bool(fsp3.search(id)):
        fr=fsp2.findall(id)
        sid=fr[0][1][:3]+fr[0][0][-9:-2]
    else:
        fr=fsp1.findall(id)
        sid=fr[0][1][:3]+fr[0][0][-9:-2]
    return sid


GH1t=Tree("/home/shulinhe/GH1/GH1s_iqtree/GH1s_full/GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified.rooted.swift")

###test positive selection in beetles##
#######################################
detach_node=GH1t.get_common_ancestor(["XP_011556095.1_Plutella","XP_557098.2_Anopheles"])
####test gene clusters with different groups, myronase in two groups beetles and lepidoptera, beta-glucosidase, male/female specific gene.
#detach_node.detach() #detach
#GH1t_a_d_names=GH1t.get_leaf_names() 
#GH1t.prune([detach_node]) #prune
with open("All_codon.aln.fa.rdp5.csv.recombinates",'r') as file:
    recombinates=[a.strip('\n').rsplit('_',2)[0] for a in file.readlines()]

tids=[i.rsplit('_',1)[0] for i in detach_node.get_leaf_names()]
tids_tidy=[a for a in tids if a not in recombinates]
tsps=[seq for seq in SeqIO.parse("All_genome_protein_GH1e.faa",'fasta') if seq.id.rsplit('_',2)[0] in tids_tidy]
#tsps=[seq for seq in SeqIO.parse("All_genome_protein_GH1e.faa",'fasta') if seq.id.rsplit('_',2)[0] in tids]
tscdss=[seq for seq in SeqIO.parse("All_genome_cds_GH1e.fna",'fasta') if seq.id.rsplit('_',2)[0] in tids_tidy]
tspsids=[seq.id.rsplit("_",1)[0] for seq in tsps]
tspsleafs=[ln for ln in detach_node.get_leaf_names() if ln in tspsids]
detach_node.prune(tspsleafs,preserve_branch_length=False)
for seq in tsps:
    seq.id=get_phylip_name(seq.id)

for seq in tscdss:
    seq.id=get_phylip_name(seq.id)

for leaf in detach_node:
    leaf.name=get_phylip_name(leaf.name)

SeqIO.write(tsps,"lepidip_group2_protein.faa",'fasta')
SeqIO.write(tscdss,"lepidip_group2_cds.fna",'fasta')
detach_node.write(format=9, outfile="lepidip_group2_tree.nwk")
```

__align proteins and faas__

```console
mafft beetle_group5_protein.faa > beetle_group5_protein.faa.mafft && muscle -in beetle_group5_protein.faa -out beetle_group5_protein.faa.muscle && ~/opt/rascal1.34/rascal beetle_group5_protein.faa.mafft beetle_group5_protein.faa.mafft.rascal && ~/opt/rascal1.34/rascal beetle_group5_protein.faa.muscle beetle_group5_protein.faa.muscle.rascal
for i in `ls beetle_group5_protein.faa.*` ; do echo $i >> beetle_group5_protein.faa.normd; ~/opt/normd_noexpat/normd $i >> beetle_group5_protein.faa.normd; done

~/opt/pal2nal.v14/pal2nal.pl selected_test_protein.faa.muscle.rascal selected_test_cds.fna -output paml >selected_test_codon.aln.paml
###################codeml to test positive selection in genes####################
##edit "selected_test_tree.nwk" for testing under http://phylotree.hyphy.org/, the output is not rooted, need to add a ";" at the end.
##edit configure files.
~/opt/paml4.9j/bin/codeml test.Sites_model_codeml.ctl
###################hyphy to test positive selection in genes##########
#####################interactive is much better to use################
hyphy meme --alignment selected_test_codon.aln.paml --tree selected_test_tree.nwk --branches All --pvalue 0.1 --output selected_test.meme.json
hyphy fel --alignment selected_test_codon.aln.paml --tree selected_test_tree.nwk --branches All --srv Yes --pvalue 0.1 --output selected_test.fel.json #or test branch label
hyphy fubar --alignment selected_test_codon.aln.paml --tree selected_test_tree.nwk --output selected_test.fel.json
hyphy relax --alignment selected_test_codon.aln.paml --tree selected_test_tree.nwk --mode Classic --models Minimal --test test --output selected_test.fel.json #a significant K>1 would indicate intensified selection on test lineages, and significant K<1 would indicate relaxed selection on test lineages.
```

### plot alignment from full gene tree

```python
#python# to retrieve tips.
from ete3 import Tree
from Bio import SeqIO,Seq
import re
import os
fsp1=re.compile(r"(\d+[\.]*\d+)\_(\w+)\_*\w*$")
fsp2=re.compile(r"(\d+[\.]*\d+)\_(\w+)\_\w*$")
fsp3=re.compile(r"\d+\_\w+\_\w+$")
def get_phylip_name(id):
    if bool(fsp3.search(id)):
        fr=fsp2.findall(id)
        sid=fr[0][1][:3]+fr[0][0][-9:-2]
    else:
        fr=fsp1.findall(id)
        sid=fr[0][1][:3]+fr[0][0][-9:-2]
    return sid


GH1t=Tree("/home/shulinhe/GH1/GH1s_iqtree/GH1s_full/GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified.rooted.swift")

detach_node=GH1t.get_common_ancestor(["XP_011556095.1_Plutella","XP_022825432.1_Spodoptera"])
####test gene clusters with different groups, myronase in two groups beetles and lepidoptera, beta-glucosidase, male/female specific gene.
#detach_node.detach() #detach
#GH1t_a_d_names=GH1t.get_leaf_names() 
#GH1t.prune([detach_node]) #prune
tids=[i.rsplit('_',1)[0] for i in detach_node.get_leaf_names()]
tsps=[seq for seq in SeqIO.parse("GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal",'fasta') if seq.id.rsplit('_',2)[0] in tids]
SeqIO.write(tsps,"test.faa",'fasta')
```
