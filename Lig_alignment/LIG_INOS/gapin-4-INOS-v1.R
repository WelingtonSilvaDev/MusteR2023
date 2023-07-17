begin <- Sys.time()
library("tidyverse")
library("bio3d")
library("doMC")

prodpath = "/data/compute/pdb_import"
workpath = getwd()
pdbupper = F

taskline = "######################################################"
emptline = "                                                      "
prodline = "#-#-#-#-#-#-#-#- P R O D U C T I O N  -#-#-#-#-#-#-#-#"
testline = "#-#-#-#-#-#-#-#-       T E S T       #-#-#-#-#-#-#-#-#"
impoline = "                     GENERATING                       "
aligline = "                      ALIGNING                        "
filtline = "                      FILTERING                       "

if (workpath==prodpath){
  bibpath = paste0(prodpath,"/bib/")
  pdbpath = "/data/pdb/"
  print(taskline)
  print(emptline)
  print(prodline)
  print(impoline)
  print(emptline)
  #print(taskline)
}else{
  #bibpath = paste0(workpath,"/bib/")
  bibpath = paste0("../../bib/")
  #pdbpath = "/data/test/pdb/"
  pdbpathin = "pdb-in/"
  pdbpathout = "pdb/"
  pdbpath = pdbpathout
  print(taskline)
  print(emptline)
  print(testline)
  print(filtline)
  print(emptline)
}

print(emptline)
print(paste("workpath:",workpath))
print(paste("bibpath:",bibpath))
print(paste("pdbpath:",pdbpath))
print(emptline)

source(paste0(bibpath,"common.R"))
source(paste0(bibpath,"fun-base-v45.R"))
source(paste0(bibpath,"fun-mass-v45.R"))

a = 0
LOCAL = 1

vdw = set_atom_names()
ttime=c()
b = 1
aanum = c(HOH=Inf,GLY=Inf,ALA=Inf,VAL=2,LEU=3,MET=3,ILE=3,PRO=1,SER=1,THR=2,CYS=1,ASN=3,GLN=4,LYS=4,ARG=6,HIS=5,PHE=6,TYR=7,TRP=9,ASP=3,GLU=4)
meso=c("Q.CD","S.CB","N.CG","E.CD","D.CG","K.CE","T.CB","R.CZ","Y.CZ") 
aacode=c("a","m","p","b","c","i","o","h")
probe = 1.400 #verificar possibilidade de aumentar para 1.5 ou 1.6
nmax.atom = 10000
db.size = c(0)
dba.size = c(0)
dba2.size = c(0)
dbpro.size = c(0)
dbgrp.size = c(0)
dbsel.size = c(0)
dball.size = c(0)
###################################

#### variáveis ####################
ncores = 2
#LIGs = c()
#PDBS = c()
#PDBS = c("bis-isothiourea-iNOSas-Docked-Pose1")
#LIGs = c("UNK,0,HEM,1010")

db.size = c(0)


# prepara para importacao
if(0){
  print(taskline)
  print("READING LIG TABLE...")
  filename = "LIGS/leo-ligs-INOS-v1.csv"
  lig.table = read_delim(file=filename,delim=",")
  decname = "LIGS/heme-decomp-v1.csv"
  dec.table = read_delim(file=decname,delim=",")
  dec.table$atom = str_replace_na(dec.table$atom)
  #scriptname = "gapin-pdb2align.R"
  sep =","
  #for (i in 1:dim(lig.table)[1]){
  #for (i in 1:2){
  PDBS = c()
  LIGs = c()
  TARs = c()
  #for (i in c(1,seq(2,by=3,to=50))){
  for (i in 1:dim(lig.table)[1]){
  #for (i in c(26)){
    #i = 20
    PDBS = c(PDBS,lig.table$pdb_title[i])
    LIGs = c(LIGs,paste0(lig.table[i,2],sep,lig.table[i,3]))
    TARs = c(TARs,paste0(lig.table[i,4],sep,lig.table[i,5]))
  }
}

# Filtra...
if(0){
  print(taskline)
  pre.sufix = ""
  pdb.sufix = ".pdb"
  if (pdbupper){
    PDBS = toupper(PDBS)
  }
  # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=0,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "IMPORTING...")### percentage=101, mensagem de erro.
  ttime = c(proc.time()[3],ttime);print("INITIALIZING DB LIST...")
  if (LOCAL){
    print(paste("WARNING: LOCAL flag enabled, reading local PDB",PDBS[1],"..."))
    copy_from_server(pdbid=PDBS[1],workpath,prodpath,ipaddr,pdbpath)#;readline()
  }else{
    download_pdb(PDBS,pdbpath,sufix=pdb.sufix)
  }
  pdb = list()
  pdb = initialize_db_list(PDBS,pdbpathin,meso,aacode,nmax.atom,workpath,prodpath,ipaddr,
                          vcores=c(ncores=0,ncores=0),pre.sufix=pre.sufix,sufix=pdb.sufix,
                          conform_1_chain=T,ligs=LIGs,tars=TARs,slim="-",pdbpathout=pdbpathout,
                          dec.table=dec.table)
  if (!is.null(pdb[[1]]$error)){
    if (pdb[[1]]$error==2){
      print(paste("WARNING: PDBid",pdb[[1]]$pdbname,"has no interfaces"))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath,"ERROR: No interfaces found")
    }else if (pdb[[1]]$error==1){
      print(paste("WARNING: PDBid",pdb[[1]]$pdbname,"has",pdb[[1]]$tampdb,
                  "atoms and exceeds the limit of:",nmax.atom))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, paste(c("ERROR: maxatom exceeded:",nmax.atom)))
    }else if (pdb[[1]]$error==0){
      print(paste("WARNING: PDB file for",pdb[[1]]$pdbname,"not found"))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "ERROR: nonexistent file")
    }else{
      print(paste("WARNING: some problem with server connection for",pdb[[1]]$pdbname))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "ERROR: server failure")
    }
    print("Halting all!")
    a=0
  }
}

# Alinha...
if (0){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("ALIGNMENT OF TARGET DB LIST...")
  pdb = align_db_list(pdb,pdbpath,workpath,prodpath,vcores=c(ncores=0,ncores=0),pre.sufix=pre.sufix,
                     sufix=pdb.sufix,pdbpathout="pdb-align/")
  pdb.size = object.size(pdb)
  ttime = c(proc.time()[3],ttime);print(paste("Time alignment:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(pdb.size,units="Mb")))
  
}

if(1){
  library("Matrix")
  library("pracma")
  library("cluster")
  library("bio3d")
  library("geometry")
  library("deldir")
  library("caret")
  library("graphkernels")
  library("stringr")
  library("shape")
  library("prodlim")
  library("rARPACK")
  library("stringi")
  library("RColorBrewer")
  library("mongolite")
  library("tnet")
  library("rgl")
  library("plotly")
}

# Importa...
a = 0
if(a){
  PDBS = c()
  for (i in 1:length(pdb)){
    PDBS = c(PDBS,pdb[[i]]$pdbname)
  }
  #PDBS = c(pdb[[1]]$pdbname,pdb[[2]]$pdbname)
  print(taskline)
  pre.sufix = ""
  pdb.sufix = ".pdb"
  if (pdbupper){
    PDBS = toupper(PDBS)
  }
  # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=0,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "IMPORTING...")### percentage=101, mensagem de erro.
  ttime = c(proc.time()[3],ttime);print("INITIALIZING DB LIST...")
  if (LOCAL){
    print(paste("WARNING: LOCAL flag enabled, reading local PDB",PDBS[1],"..."))
    copy_from_server(pdbid=PDBS[1],workpath,prodpath,ipaddr,pdbpath)#;readline()
  }else{
    download_pdb(PDBS,pdbpath,sufix=pdb.sufix)
  }
  db = list()
  db = initialize_db_list(PDBS,pdbpath="pdb-align/",meso,aacode,nmax.atom,workpath,prodpath,ipaddr,
                          vcores=c(ncores=0,ncores=0),pre.sufix=pre.sufix,sufix=pdb.sufix,
                          conform_1_chain=F)
  if (!is.null(db[[1]]$error)){
    if (db[[1]]$error==2){
      print(paste("WARNING: PDBid",db[[1]]$pdbname,"has no interfaces"))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath,"ERROR: No interfaces found")
    }else if (db[[1]]$error==1){
      print(paste("WARNING: PDBid",db[[1]]$pdbname,"has",db[[1]]$tampdb,"atoms and exceeds the limit of:",nmax.atom))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, paste(c("ERROR: maxatom exceeded:",nmax.atom)))
    }else if (db[[1]]$error==0){
      print(paste("WARNING: PDB file for",db[[1]]$pdbname,"not found"))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "ERROR: nonexistent file")
    }else{
      print(paste("WARNING: some problem with server connection for",db[[1]]$pdbname))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "ERROR: server failure")
    }
    print("Halting all!")
    a=0
  }
  
  db.size = object.size(db)
  ttime = c(proc.time()[3],ttime);print(paste("Time dist:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(db.size,units="Mb")))
  #browser()
  
}

# Calcula matriz de areas BSR
# OBS: tem pre-requesito o if anterior
if(a){
  #browser()
  print(taskline)
  # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=10,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "CREATING MATRICES...")
  ttime = c(proc.time()[3],ttime);print("ADDING BSR MATRIX...")
  db = add_bsr_matrix(db,vdw,probe,vcores=c(ncores,ncores))
  db.size = object.size(db)
  ttime = c(proc.time()[3],ttime);print(paste("Time adding bsr:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(db.size,units="Mb")))
  
}

# Expande matrix BSR de ALL para PP e AA
if(a){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("EXPANDING BSR ALL TO AA AND PP...")
  db = expand_bsr_matrix(db,vcores=c(0,0))
  db.size = object.size(db)
  ttime = c(proc.time()[3],ttime);print(paste("Time expanding bsr:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(db.size,units="Mb")))
  
}

# Monta dba a partir de db
if(a){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("BUILDING DBA FROM DB...")
  dba = build_dba_from_db(db,vdw)
  dba.size = object.size(dba)
  ttime = c(proc.time()[3],ttime);print(paste("Time building dba from db:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(dba.size,units="Mb")))
  
}

# Montagem do supercluster: 1a filtragem (descobre elementos conexos)
#TO DO: cercar erros quando AA ou PP for nulo
if(a){
  
  print(taskline)
  # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=30,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "CLUSTERIZATIONS...")
  ttime = c(proc.time()[3],ttime);print("SUPERCLUS 1o FILTER...")
  dba = make_hyper_superclus(dba,icores=3,vcores=c(ncores,ncores,ncores,ncores,0),conex=T,only_a=T,type=1,karpack=4) #karpack=22
  dba.size = object.size(dba)
  ttime = c(proc.time()[3],ttime);print(paste("Time superclus 1o filter:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba:",format(dba.size,units="Mb")))
  
}

# Montagem do supercluster: 2a filtragem (considera elementos conexos acima de cutpar)
if(a){
  
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("SUPERCLUS 2o FILTER...")
  dba2 = make_hyper_filter(dba,icores=0,vcores=c(0),cutpar=10)
  dba2 = make_hyper_superclus(dba2,icores=3,vcores=c(ncores,ncores,ncores,ncores,ncores),conex=F,only_a=F,type=2,karpack=4,lowcut=c(0.5,0.5,21))#;readline()
  dba2 = make_hyper_reorientation(dba2,icores=0,vcores=c(ncores=0))
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time superclus 2o filter:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
  
}

if(0){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("MAKE VISMODEL...")
  # NAO MASS AINDA... 
  #dba2[[1]][[2]] = add_res_information_to_clusters(dba2[[1]][[2]])
  seq.id = seq(2,by=2,length(dba2[[1]]))
  for (i in seq.id ){
    dba2[[1]][[i]] = add_res_information_to_clusters(dba2[[1]][[i]])
  }
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time make vismodel:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}


if(1){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("CONSOLIDATE RES TABLE...")

  seq.id = c(80,82) #Docked-trans-geranylcetone-iNOSas
  #seq.id = c(seq(16,by=2,50)) # top-iNOSas-COX1
  #seq.id = c(seq(16,by=2,24)) # Docked-Caryophyllene-a1-iNOSas-Pose1
  seq.id = c(2,4,12,14) #Crystal-iNOSas-PDB-1NSI-CA-CB-Substrate-ARG/ITU
  m.ligs1 = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))

  seq.id = c(80,82) #Docked-trans-geranylcetone-iNOSas
  #seq.id = c(seq(52,by=2,82)) # top-iNOSas-iNOScs
  #seq.id = c(seq(16,by=2,24)) # Docked-Caryophyllene-a1-iNOSas-Pose1
  seq.id = c(6,8,10) #bis-isothiourea
  m.ligs2 = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  
  seq.id = c(seq(84,by=2,100)) # top-just-iNOSas
  seq.id = c(seq(16,by=2,100))
  m.ligs3 = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))

  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time consolidate res table:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}

if(1){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("PLOT VISMODEL...")
  #m1 = dba2[[1]][[6]]$exp$re[[1]][[1]]$res.table
  p.ligs1 = plot_lig_res_interactions(m.ligs1,scene="scene1",area.lim=100,alpha=0.1,n=1200)
  p.ligs2 = plot_lig_res_interactions(m.ligs2,scene="scene2",area.lim=100,alpha=0.1,n=1200)
  p.ligs3 = plot_lig_res_interactions(m.ligs3,scene="scene3",area.lim=100,alpha=0.1,n=400)
  
  range = c(-8, 8,90, 110,10, 30)
  text = c("Docked trans-geranylcetone iNOSas-Poses",
           "Docked trans-geranylcetone iNOSas-Poses",
           "")
  # text = c("Top iNOSas-COX1 Poses",
  #          "Top iNOSas-iNOScs Poses",
  #          "Top just iNOSas Poses")
  # text = c("Docked Caryophyllene iNOSas Poses",
  #          "Docked Caryophyllene iNOSas Poses",
  #          "")
  text = c("iNOSas ARG ITU crystal (A/B)",
           "Docked iNOSas bis-isothiourea Poses",
           "Docked Top iNOSas Poses")
  p.all = plotly_all_together_ligs(p.ligs1,p.ligs2,p.ligs3,text=text,range=range)
  #p.all = plotly_all_together_ligs(p.ligs1,p.ligs2,text=text,range=range)
  print(p.all)

  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time plot vismodel:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}

if(1){
  library(orca)
  p.sizes = c(1200,600)
  #orca(p.all,"../Images/Caryophillen-iNOSas-v2.png",scale=2,width=p.sizes[1],height=p.sizes[2])
  #orca(p.all,"../Images/iNOSas-bygroup-v2.png",scale=2,width=p.sizes[1],height=p.sizes[2])
  #orca(p.all,"../Images/trans-geranylcetone-iNOSas-v2.png",scale=2,width=p.sizes[1],height=p.sizes[2])
  orca(p.all,"../Images/iNOSas-ARG-ITU-bis-iso-ligs-v2.png",scale=2,width=p.sizes[1],height=p.sizes[2])
}


###################### OLD #################333

if(0){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("CONSOLIDATE RES TABLE...")
  seq.id = c(2,4,12,14)
  m.inos.arg.itu = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  seq.id = c(6,8,10)
  m.inos.bis = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  seq.id = seq(16,by=2,length(dba2[[1]]))
  m.inos.ligs = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  #seq.id = seq(52,by=2,82)
  #m.inos.nos.nos = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  #seq.id = seq(84,by=2,100)
  #m.inos.just.nos = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time consolidate res table:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}

if(0){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("PLOT VISMODEL...")
  #m1 = dba2[[1]][[6]]$exp$re[[1]][[1]]$res.table
  p.inos.arg.itu = plot_lig_res_interactions(m.inos.arg.itu,scene="scene1",area.lim=100,alpha=0.1,n=800)
  p.inos.bis = plot_lig_res_interactions(m.inos.bis,scene="scene2",area.lim=100,alpha=0.1,n=1000)
  p.inos.ligs = plot_lig_res_interactions(m.inos.ligs,scene="scene3",area.lim=100,alpha=0.1,n=200)
  p.all = plotly_all_together_INOSas_ligs(p.inos.arg.itu,p.inos.bis,p.inos.ligs)
                                     
  print(p.all)
  #print(p.inos.arg.itu)
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time plot vismodel:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}


############# COM SEGREGAÇÃO #################3

if(0){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("CONSOLIDATE RES TABLE...")
  seq.id = c(2,4)
  m.inos.arg = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  seq.id = c(2,4,12,14)
  m.inos.arg.itu = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  seq.id = c(6,8,10)
  m.inos.bis = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  seq.id = seq(16,by=2,50)
  m.inos.nos.cox = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  seq.id = seq(52,by=2,82)
  m.inos.nos.nos = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  seq.id = seq(84,by=2,100)
  m.inos.just.nos = consolidate_res_table(dba2,seq.id,noid=c("HCC","HCO","HCR","HFN"))
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time consolidate res table:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}

if(0){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("PLOT VISMODEL...")
  #m1 = dba2[[1]][[6]]$exp$re[[1]][[1]]$res.table
  p.inos.arg = plot_lig_res_interactions(m.inos.arg,scene="scene1",area.lim=100,alpha=0.1,n=1000)
  p.inos.arg.itu = plot_lig_res_interactions(m.inos.arg.itu,scene="scene2",area.lim=100,alpha=0.1,n=1000)
  p.inos.bis = plot_lig_res_interactions(m.inos.bis,scene="scene3",area.lim=100,alpha=0.1,n=1000)
  p.inos.nos.cox = plot_lig_res_interactions(m.inos.nos.cox,scene="scene4",area.lim=100,alpha=0.1,n=300)
  p.inos.nos.nos = plot_lig_res_interactions(m.inos.nos.nos,scene="scene5",area.lim=100,alpha=0.1,n=300)
  p.inos.just.nos = plot_lig_res_interactions(m.inos.just.nos,scene="scene6",area.lim=100,alpha=0.1,n=300)
  #p.inos.nos = plot_lig_res_interactions(m.inos.nos,scene="scene3",area.lim=100,alpha=0.1,n=1000)
  #p.all = subplot(p.inos.nos.cox,p.inos.nos.nos,p.inos.just.nos);print(p.all)
  p.all = plotly_all_together_INOSas(p.inos.arg,p.inos.arg.itu,p.inos.bis,
                                     p.inos.nos.cox,p.inos.nos.nos,p.inos.just.nos)
  print(p.all)
  #print(p.inos.arg.itu)
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time plot vismodel:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}

#################################### ARA/ARO ##############3333

# CONSOLIDATE RES.TABLE ARA/ARO
if(0){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("CONSOLIDATE RES TABLE...")
  #x = pdb[[1]]$atom %>% filter(type =="ATOM", chain == "A") %>% select(resno)
  #y = pdb[[2]]$atom %>% filter(type =="ATOM", chain == "A") %>% select(resno)
  #id.tab = bind_cols(x,y)
  #browser()
  #x = renumerate_res_table(dba2[[1]][[2]]$exp$residues[[1]][[1]]$res.table,adj=-24)
  seq.id = c(2)
  m.cry = consolidate_res_table(dba2,seq.id)
  m.cry = renumerate_res_table(m.cry,adj=-24)
  seq.id = seq(4,by=2,8)
  #browser()
  m.ara = consolidate_res_table(dba2,seq.id)
  seq.id = seq(10,by=2,length(dba2[[1]]))
  m.aro = consolidate_res_table(dba2,seq.id)
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time consolidate res table:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}

# PLOT VISMODEL
# https://plotly.com/r/reference/#mesh3d
# https://cran.r-project.org/web/packages/mvmesh/
# https://rdrr.io/github/obreschkow/cooltools/man/fibonaccisphere.html
# https://plotly.com/r/3d-subplots/#3d-subplot
#https://rpubs.com/bcd/subplot-titles

if(0){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("PLOT VISMODEL...")
  #m1 = dba2[[1]][[6]]$exp$re[[1]][[1]]$res.table
  p.cry = plot_lig_res_interactions(m.cry,scene="scene1",area.lim=300,alpha=0.1)
  p.ara = plot_lig_res_interactions(m.ara,scene="scene2",area.lim=300,alpha=0.1)
  #print(p.ara)
  p.aro = plot_lig_res_interactions(m.aro,scene="scene3",area.lim=300,alpha=0.1)
  #print(p.aro)
  #p.all = subplot(p.ara,p.aro)
  p.all = plotly_all_together(p.cry,p.aro,p.ara)
  print(p.all)
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time plot vismodel:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}

#m1 = dba2[[1]][[2]]$exp$re[[1]][[1]]$res.table
#kk = generate_spheres_mesh(m1 %>% filter(chain == "B") %>% select(x,y,z),n=120)
#plot3d(x=kk$x,y=kk$y,z=kk$z,type="s",size=1)




################ TESTES

if(0){
  
  x = pdb[[1]]$atom %>% filter(type =="ATOM", chain == "A")
  y = pdb[[2]]$atom %>% filter(type =="ATOM", chain == "A")
  
  y = y %>% mutate(eleno=row_number())
  
  y = y %>% mutate(nores = paste0(eleno,resid))
  x = x %>% mutate(nores = paste0(eleno,resid))
  
  x.seq = aa321(pdb[[1]]$atom[pdb[[1]]$calpha,"resid"])
  y.seq = aa321(pdb[[2]]$atom[pdb[[2]]$calpha,"resid"])
  
  seqs <- seqbind(x.seq,y.seq)
  seqaln(seqs)
  
}


if(0){
  
  m = dba2[[1]][[2]]$exp$residues[[1]][[1]]$res.table
  size = 2
  col = "black"
  lab = c("","","")
  add = 0
  lim = c(-50,260)
  
  m1 = m %>% filter(area>500 | resn == "ACD")
  
  temp = round(m1$area,0)
  f = which(temp==max(temp))
  temp[f]=0
  
  fi = m1$resn != "ACD"
  tempi = temp
  tempi[fi] = 0
  
  m1 = m1 %>% mutate(resCi = case_when(
                     resCi == "ACD700" ~ " ",
                     TRUE ~ resCi)
                     )
  
  m1$resCi[f] = "ACD"
  
  m1 = m1 %>% rename(xc = x, yc = y, zc = z) 
  
  #p = plot_ly(x=m$xc, y=m$yc, z=m$zc, type="scatter3d", mode="markers", color=temp)
  #p = p %>% add_trace(text = m$residue, hoverinfo = 'text')
  #print(p)
  t <- list(
    family = "sans serif",
    size = 14,
    color = toRGB("grey50"))
  
  
  x = m1 %>% plot_ly(x= ~xc, y= ~yc, z = ~zc) %>%
    add_markers(
      size = tempi,
      sizes = c(0,2000),
      alpha = 0.1
    ) %>% 
    add_text(
      text= ~resCi,
      textposition = 'middle',
      hoverinfo = "text",
      textfont = list(family = "Roboto Condensed", size = temp/40)#
    ) %>%
    hide_legend()
  
  print(x)
}

if(0){
  
  m = dba2[[1]][[4]]$exp$residues[[1]][[1]]$res.table
  size = 2
  col = "black"
  lab = c("","","")
  add = 0
  lim = c(-50,260)
  resnm = "UNL"
  
  m1 = m %>% filter(area>500 | resn == resnm)
  
  temp = round(m1$area,0)
  f = which(temp==max(temp))
  temp[f]=0
  
  fi = m1$resn != resnm
  tempi = temp
  tempi[fi] = 0
  
  m1 = m1 %>% mutate(resCi = case_when(
    resCi == "UNL1" ~ " ",
    TRUE ~ resCi)
  )
  
  m1$resCi[f] = resnm
  
  m1 = m1 %>% rename(xc = x, yc = y, zc = z) 
  
  #p = plot_ly(x=m$xc, y=m$yc, z=m$zc, type="scatter3d", mode="markers", color=temp)
  #p = p %>% add_trace(text = m$residue, hoverinfo = 'text')
  #print(p)
  t <- list(
    family = "sans serif",
    size = 14,
    color = toRGB("grey50"))
  
  
  x = m1 %>% plot_ly(x= ~xc, y= ~yc, z = ~zc) %>%
    add_markers(
      size = tempi,
      sizes = c(0,2000),
      alpha = 0.1
    ) %>% 
    add_text(
      text= ~resCi,
      textposition = 'middle',
      hoverinfo = "text",
      textfont = list(family = "Roboto Condensed", size = temp/40)#
    ) %>%
    hide_legend()
  
  print(x)
}

  # p = m1 %>% plot_ly(x= ~xc, y= ~yc) %>%
  #     add_markers(
  #       size = tempi,
  #       sizes = c(0,2000),
  #       alpha = 0.1,
  #       hoverinfo = "text"
  #     ) %>%
  #     add_trace(
  #       type="scatter",
  #       mode="text",
  #       text= ~resCi,
  #       textposition = 'middle',
  #       hoverinfo = "text",
  #       textfont = list(family = "Roboto Condensed", size = temp/40)#,
  #       #size = temp,
  #       #sizes = c(1,10)
  #     
  #       #color = temp
  #   ) %>%
  #   hide_legend()

  #print(p)
  

  
  
  ### TESTANDO UM POLYGONS
  # b = m1 %>% filter(resn == "ACD") %>% 
  #     plotly_empty(x = ~xc, y = ~yc)
  # b = b %>% add_polygons(hoverinfo = "none", color = I("black"))
  # print(b)
  
  
  
  #p = plot_ly(x=m$xc, y=m$yc, z=m$zc, text=m$residue,size=2,color=temp)
  #p = m1 %>% plot_ly(x= ~xc, y= ~yc, z= ~zc, text=m1$residue,size=2,color=temp)
  #p = m1 %>% plot_ly(x= ~xc, y= ~yc, mode="text", type="scatter", 
  #                   text= ~resCi, textposition = 'middle right')
  
  
  #p = p %>% add_markers()
  #p = p %>% add_text(textfont = t, textposition = "right")
  #p <- p %>% layout(xaxis = list(range = c(1.6, 3.2)),
  #                      showlegend = FALSE)
  #print(p)
  
  #p = m %>% plot_ly(x ~ xc, y ~ yc, z ~ zc, type="scatter3d",mode="markers",color="black")
  #print(p)
  # %>%
  #add_markers(type="scatter3d", mode="markers")
  
  #plot_ly(x=m[,"xc"], y=m[,"yc"], z=m[,"zc"], type="scatter3d", mode="markers", color=temp)
  
  #ids=plot3d(x=m[,"xc"],y=m[,"yc"],z=m[,"zc"],type="s",size=size,box=F,col=col,
  #           xlab=lab[1],ylab=lab[2],zlab=lab[3],xlim=lim,ylim=lim,zlim=lim,add=add)
  #if (text){
  #  text3d(x=m[,1],y=m[,2],z=m[,3],text=paste(1:dim(m)[1]),adj=adj,cex=textcex)
  #  
  #}
  
  
  #x = dba2[[1]][[2]]
  
  #f = grep("LEU534",x$element_name)
  
  #y = get_res_geom_center_from_atom_list(f,x$element_name,x$pdb$xyz)
  
  #browser()
  
  #xyz1 = x$exp$superclus[[2]]$geomc
  
  #plot_super(exp=x$exp,n=1:4,mfrow=c(2,2))
  #browser()
  
  #multi_plot_model_clus(db=x,krange=4,mfrow=c(1,1),main="",xlab="",ylab="",xlim=c(10,50),ylim=c(10,50))
  #plot_model_clus2(db$kclus[[i]]$clus,clusvis,auxpdb,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
  #                 text.marg=aux.marg, add=add,lcol=lcol,col=col[1:3],dlim=dlim,dres=dres,ri=ri,
  #                 text.cex=text.cex,mtext.cex=mtext.cex,restype=restype,lwd=lwd,lty=lty,font=font)
  
#}


