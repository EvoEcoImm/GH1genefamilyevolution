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

#"""End of R"""

```
__python to get represent protein id__
```python
import pandas as pd
from Bio import SeqIO
data=pd.read_csv("Genome_ids_GH1_list.csv")
data['ID']=data['ID'].str.replace("^gene-","")
orgset=set(data['Org'].values.tolist())
spec2ids2={a:data[data.Org==a]['ID'].values.tolist() for a in orgset}

#create species 2 gene ids dictionary
data2list=data.values.tolist()
spec2ids=dict()
for a,b in data2list:
	if b not in spec2ids.keys():
		spec2ids[b]=[a]
	else:
		spec2ids[b].append(a)


def retrieve_p(spec,idlist):
	fa="./"+spec+"_GH1s.gb"
	fq=[seq for seq in SeqIO.parse(fa,'gb')]
	genestr='_'.join(idlist)
	g2p_dict={}
	for seq in fq:
		sfqf=[sf for sf in seq.features if sf.type=="CDS" ]
		slecqf=[sfe for sfe in sfqf if sfe.qualifiers["gene"][0] in genestr]
		g_id=slecqf[0].qualifiers["gene"][0]
		if len(slecqf)==1:
			protein_id=slecqf[0].qualifiers["protein_id"][0]
		else:
			protein_f=slecqf[0]
			for ssf in slecqf[1:]:
				if len(ssf) > len(protein_f):
					protein_f=ssf
			protein_id=protein_f.qualifiers["protein_id"][0]
		g2p_dict[g_id]=protein_id
	return g2p_dict
total={}
for spec,glist in spec2ids.items():
	try:
		total[spec]=retrieve_p(spec,glist)
	except KeyError:
		print(spec,"cds has no gene ids")
	except IndexError:
		print(spec,"has no cds annotation")
for spec in total.keys():
	df=pd.DataFrame.from_dict(total[spec],orient='index')
	df.to_csv(spec+"_GHs_g2p.csv",header=False)

###"""python retrieve protein seqs from database"""###
from Bio import SeqIO,Seq
import os

def retrieve_p_seq(spec):
	genome_folder_path="/media/shulinhe/DATA/Genome_metagenome_transcriptomes/resource/Arthropoda_genome_Chromo_Scaff_ref_repre_RefSeq/"
	p_fn=[file for file in os.listdir(genome_folder_path+spec) if file.endswith("_protein.faa")]
	p_fp=os.path.join(genome_folder_path,spec,p_fn[0])
	id_fn=spec+"_GHs_g2p.csv"
	try:
		with open(id_fn,'r') as file:
			ids=[e.strip('\n').split(',')[-1] for e in file.readlines()]
		pseq=[seq for seq in SeqIO.parse(p_fp,'fasta') if seq.id in ids]
	except FileNotFoundError:
		print(spec,"has no GH1 genes")
		pseq=""
	return pseq

def replaceseqname(seq):
	cid=seq.description.replace('[',"").replace(']',"").split(' ')
	cidn='_'.join([cid[0],cid[-2],cid[-1]])
	seq.id=cidn
	return seq

with open("Genome_specids_selected",'rt') as file:
	specs=[item.strip('\n') for item in file.readlines()]

totalseq=[retrieve_p_seq(spec) for spec in specs]
totalseqs=[replaceseqname(seq) for f in totalseq if f!="" for seq in f]

SeqIO.write(totalseqs,"GH1s_genome_protein.faa",'fasta')

```
__retrieve cds and protein seqs__
```R
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
	#sgenesummary<- genesummary %>% group_by(group_name) %>% filter(wcds==max(wcds)) %>% top_n(1) %>% ungroup() %>% as.data.frame
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
```

### 3. run taxonomy text for all proteins by top 10 blastp search from ncbi
```console
for i in {00..09}; do ~/GH1/taxonomy_contamination_detect_transcriptome.py -i GH1s_genome_protein.$i.faa.alignment.xml -o GH1s_genome_protein.$i.faa.alignment.xml.tax; done;cat *xml.tax > GH1s_genome_protein.faa.alignment.xml.tax; rm GH1s_genome_protein.*.faa.alignment.xml.tax
#manully check the sequences from bacteria.
```
### 4. Dual domain proteins
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
	if any(GH1os['HMMER'].str.contains("\+")):
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

```R
library(tidyverse)
library(Biostrings)
specs<-readLines("Genome_specids_selected_tree")
for(spec in specs){
GH1os<-read.delim(paste0(spec,"/",spec,"overview.txt"),stringsAsFactors=F) %>% filter(.,X.ofTools>2)%>%mutate(Org=spec)
specre<- if (any(grepl("\\+",GH1os$HMMER))) GH1os else GH1os[0,]
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
### 5. TREE BUILDING for whole GH1s from genomes
```console
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
```R
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

### 6.NOTUNG reconsilation of gene tree from species tree
__both tree should be rooted__
__prepare tree tips names__
```R
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
```console
mv GH1s_genome_species.tree.modified GH1s_genome_species.tree.modified.nwk
mv GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified.rooted GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified.rooted.nws
java -jar ~/opt/Notung/Notung-2.9.1.5.jar -g GH1s.genome.protein.faa.mafft.contree.extracted.faa.muscle.rascal.contree.modified.rooted -s GH1s_genome_species.tree.modified --rearrange --speciestag postfix --threshold 85% --bootstraps name --outputdir GH1_reconciletest --log --events --parsable --treestats --progressbar --savepng --saveweakedgespng --homologtablecsv
```

### 7. get cds sequences from database based on protein sequences
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

def groupseqs(spec,seqs):	#retrieve protein ids list for spec
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
### Tree building according to orders based on the edited cds/protein files

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
```console
for i in {"Bettles","Cockroaches","Diptera","Hemiptera","Lepitoptera"}; do \
AQUA.tcl GH1s_genome_protein.faa GH1_genome_AQUA; \

#build tree with iqtree
iqtree -s GH1s_genome_protein.faa.mafft.rascal -nt 10 -bb 1000 -alrt 1000
```
```R
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
```console
~/opt/pal2nal.v14/pal2nal.pl ./GH1inorder/Tree_bestalin/Bettles_genome_protein_GH1e.faa.muscle.edi Bettles_genome_cds_GH1e.fna -output fasta >Bettles_codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl ./GH1inorder/Tree_bestalin/Lepitoptera_genome_protein_GH1e.faa.muscle.edi ./Lepitoptera_genome_cds_GH1e.fna -output fasta >Lepitoptera_codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl Hemiptera_genome_protein_GH1e.faa.mafft Hemiptera_genome_cds_GH1e.fna -output fasta >Hemiptera_codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl Cockroaches_genome_protein_GH1e.faa.muscle.rascal Cockroaches_genome_cds_GH1e.fna -output fasta >Cockroaches_codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl Hymenoptera_genome_protein_GH1e.faa.muscle Hymenoptera_genome_cds_GH1e.fna -output fasta >Hymenoptera__codon.aln.fa
~/opt/pal2nal.v14/pal2nal.pl Diptera_genome_protein_GH1e.faa.muscle.rascal.edi Diptera_genome_cds_GH1e.fna -output fasta >Diptera_codon.aln.fa
```
### 8. GH1 positive selection test
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

