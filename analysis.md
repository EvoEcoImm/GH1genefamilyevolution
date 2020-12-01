### 1. database buildup
__extract extract GH1 proteins from dbCAN by Biopython__ 
```console
cat << EOF > pyscript.py
#!/usr/bin/python
from Bio import SeqIO

a=[se for se in SeqIO.parse("CAZyDB.07312019.fa","fasta") if se.id.split("|")[-1]=="GH1"]
SeqIO.write(a,"GH1_dbcan.fasta","fasta")

EOF

chmod 755 pyscript.py
./pyscript.py


##################################################################
#get taxonomy for GH1 of dbCAN database, which is NOT THE BEST YET
##################################################################
grep ">" GH1_dbcan.fasta |cut -f1 -d "|"|sed 's/^>//g' > protein.ids
for i in `cat protein.ids`; do \
echo $i >> protein_id2tax; \
esearch -db protein -query $i| \
elink -target taxonomy| \
efetch -format xml| \
xtract -pattern Taxon -tab "," -first TaxId ScientificName \
-group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
-block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
-block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
-block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
-block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
-block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
-block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
-group Taxon -tab "," -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS" \
>> protein_id2tax ; \
done
```
### 2.predict/extract GH1 from assembly of genomes in Genome/resource/Arthropoda_Refseq/genomes
```console
#predict GH1 with run_dbCAN.py
ls /media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/ >Genome_specids #get species id
for i in `cat Genome_specids`; do python run_dbcan.py /media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/$i/GC_*_protein.faa protein --out_pre $i --out_dir $i --db_dir ../ ;done 
```
__R scripts to extract gff from genes__
```R
library(rtracklayer)
library(Biostrings)
library(rentrez)

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


get_genes_gff<-function(spec){ #define function to extract and output gene gffs,represent_proteins
	GH1<-read.delim(paste0(spec,"/",spec,"overview.txt"),stringsAsFactors=F)
	GH1s<-GH1[GH1$X.ofTools>=3,]
	if (!nrow(GH1s)){
	genes_gff<-""
	file(paste0(spec,"_GH1.gff3"),"w")
	file(paste0(spec,"_GH1.gb"),"w")
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
	if("locus_tag"%in%colnames(mcols(genes_gff))&!any(is.na(genes_gff$locus_tag))){
		gene_names<- unlist(lapply(genes_gff$locus_tag, function(x) sub("[A-z]*[A-z]","",x)))
		} else {
		gene_names<- unlist(lapply(genes_gff$gene, function(x) sub("[A-z]*[A-z]","",x)))
		}
	genegbs<-""
	for (gname in gene_names){
		gidgbs<-get_gb(gname)
		genegbs<-paste0(genegbs,gidgbs)
		}
	gbfilename <- paste0(spec,"_GH1s.gb")
	write(genegbs,file=gbfilename,append=T)
	}
	return(genes_gff)
	}

specs<-readLines("Genome_specids") #get species ids
specs_genes<-list()
for (specid in specs){ #except "Culex_quinquefasciatus" and "Pediculus_humanus", manually download, get genbank files and gff3 structures for all the genes
	try(specs_genes[[specid]]<-get_genes_gff(specid))
	}


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
#############################################
##############	R plot genes	##############
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
###############	get representing protein seqs for each gene	##############
###################################################################################
"""r script"""
library(rtracklayer)
specs<-readLines("Genome_specids_selected")
totalgid=data.frame("ID"="","Org"="")
for (spec in specs){ #get genes of each species
	GH1<-read.delim(paste0(spec,"/",spec,"overview.txt"),stringsAsFactors=F)
	GH1s<-GH1[GH1$X.ofTools>=3,]
	nps<-nrow(GH1s)
	if (!nps) {
	geneids<-""
	} else {
	GH1s_gff<-import(paste0(spec,"_GH1.gff3"),format="gff3")
	ID<-GH1s_gff[GH1s_gff$type=="gene"]$ID
	dfg<-as.data.frame(ID)
	dfg$"Org"<-spec
	}
	totalgid<-rbind(totalgid,dfg)
	}
write.csv(totalgid[-1,],"Genome_ids_GH1_list.csv",quote=F,row.names=F)

"""End of R"""

```
