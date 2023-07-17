# ----------- Backend MusteR 2022 ----------------------
#        __  __           _       _____  
#       |  \/  |         | |     |  __ \ 
#      | \  / |_   _ ___| |_ ___| |__) |
#     | |\/| | | | / __| __/ _ \  _  / 
#    | |  | | |_| \__ \ ||  __/ | \ \ 
#   |_|  |_|\__,_|___/\__\___|_|  \_\
# Created by Carlos Henrique da Silveira and updated by Welington Goncalves Silva
# Last update on march 29, 2023
# My github: https://github.com/WelingtonSilvaDev

backend_MusteR <- function(filename, check_pdb){
library("tidyverse")
library("bio3d")
library("doMC")

a_if = 1;
  
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
  #bibpath = paste0("../bib/")
  bibpath = paste0("bib/")
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

source("show3d/3dmol.R")
source(paste0(bibpath,"common-v2.R"))
source(paste0(bibpath,"fun-base-v46.R"))
source(paste0(bibpath,"fun-mass-v46.R"))
source("main/fun-ligs-v1.R")

a = 0
if(check_pdb) #se check_pdb for TRUE, entao o usuario quer fazer upload LOCAL, entao LOCAL=1
  LOCAL = 1
else
  LOCAL = 0

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

ncores = 4

db.size = c(0)

# prepara para importacao
if(a_if){
  #browser()
  print(taskline)
  print("READING LIG TABLE...")
  #filename = "LIGS/pdb_table_test2.csv" ### ALTEROU
  lig.table = read_delim(file=filename,delim=",")
  scriptname = "gapin-pdb2align.R"
  sep =","
  #for (i in 1:dim(lig.table)[1]){
  #for (i in 1:2){
  LIGs = c()
  PDBS = c()
  #for (i in c(1,seq(2,by=3,to=50))){
  for (i in 1:dim(lig.table)[1]){
    #i = 20
    PDBS = c(PDBS,lig.table$pdb_title[i])
    LIGs = c(LIGs,paste0(lig.table[i,2],sep,lig.table[i,3]))
  }
}


# Filtra...
if(a_if){
  #browser()
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
    download_pdb(PDBS,pdbpathin,sufix=pdb.sufix, destfile=pdbpathin)
  }

  pdb1 = list()

  pdb1 = initialize_db_list(PDBS,pdbpathin,meso,aacode,nmax.atom,workpath,prodpath,ipaddr,
                           vcores=c(ncores=0,ncores=0),pre.sufix=pre.sufix,sufix=pdb.sufix, 
                           conform_1_chain=T,ligs=LIGs,slim="-",pdbpathout=pdbpathout)
  if (!is.null(pdb1[[1]]$error)){
    if (pdb1[[1]]$error==2){
      print(paste("WARNING: PDBid",pdb[[1]]$pdbname,"has no interfaces"))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath,"ERROR: No interfaces found")
    }else if (pdb1[[1]]$error==1){
      print(paste("WARNING: PDBid",pdb1[[1]]$pdbname,"has",pdb1[[1]]$tampdb,
                  "atoms and exceeds the limit of:",nmax.atom))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, paste(c("ERROR: maxatom exceeded:",nmax.atom)))
    }else if (pdb1[[1]]$error==0){
      print(paste("WARNING: PDB file for",pdb1[[1]]$pdbname,"not found"))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "ERROR: nonexistent file")
    }else{
      print(paste("WARNING: some problem with server connection for",pdb1[[1]]$pdbname))
      # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=101,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "ERROR: server failure")
    }
    print("Halting all!")
    a=0
    a_if = 0
  }
}

# Alinha...
if (a_if){
  #x = system2(command = "pymol", args = "-cp Pymol/script-python-v1.py  --  pdb/7NF5.pdb pdb/7NF5.pdb 10", stdout = T);browser()
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("ALIGNMENT OF TARGET DB LIST...")
  pdb1 = align_db_list(pdb1,ids=1:length(pdb1),pdbpath,workpath,prodpath,vcores=c(ncores=0,ncores=0),pre.sufix=pre.sufix,
                      sufix=pdb.sufix,pdbpathout="pdb-align/")
  db.size = object.size(pdb1)
  ttime = c(proc.time()[3],ttime);print(paste("Time alignment:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(db.size,units="Mb")))
  #browser()
}

if(1){
  library(Rpdb)
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
if(a_if){
  #browser()
  PDBS = c()
  for (i in 1:length(pdb1)){
    PDBS = c(PDBS,pdb1[[i]]$pdbname)
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
if(a_if){
  #browser()
  print(taskline)
  # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=10,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "CREATING MATRICES...")
  ttime = c(proc.time()[3],ttime);print("ADDING BSR MATRIX...")
  db = add_bsr_matrix(db,vdw,probe,vcores=c(ncores,ncores))
  db.size = object.size(db)
  ttime = c(proc.time()[3],ttime);print(paste("Time adding bsr:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(db.size,units="Mb")))
  #browser()
}

# Expande matrix BSR de ALL para PP e AA
if(a_if){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("EXPANDING BSR ALL TO AA AND PP...")
  db = expand_bsr_matrix(db,vcores=c(0,0))
  db.size = object.size(db)
  ttime = c(proc.time()[3],ttime);print(paste("Time expanding bsr:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(db.size,units="Mb")))
  
}

# Monta dba a partir de db
if(a_if){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("BUILDING DBA FROM DB...")
  dba = build_dba_from_db(db,vdw)
  dba.size = object.size(dba)
  ttime = c(proc.time()[3],ttime);print(paste("Time building dba from db:",round((ttime[1]-ttime[2])/b,1),"secs"))
  print(paste("Size of db:",format(dba.size,units="Mb")))
  
}

# Montagem do supercluster: 1a filtragem (descobre elementos conexos)
#TO DO: cercar erros quando AA ou PP for nulo
if(a_if){
  icoresvalue = 4 #alterei para 4 quando ncores=4, era 3 antes
  print(taskline)
  # logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=30,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "CLUSTERIZATIONS...")
  ttime = c(proc.time()[3],ttime);print("SUPERCLUS 1o FILTER...")
  dba = make_hyper_superclus(dba,icores=icoresvalue,vcores=c(ncores,ncores,ncores,ncores,0),conex=T,only_a=T,type=1,karpack=22)
  dba.size = object.size(dba)
  ttime = c(proc.time()[3],ttime);print(paste("Time superclus 1o filter:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba:",format(dba.size,units="Mb")))
  
}

# Montagem do supercluster: 2a filtragem (considera elementos conexos acima de cutpar)
if(a_if){
  
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("SUPERCLUS 2o FILTER...")
  dba2 = make_hyper_filter(dba,icores=0,vcores=c(0),cutpar=10)
  dba2 = make_hyper_superclus(dba2,icores=icoresvalue,vcores=c(ncores,ncores,ncores,ncores,ncores),conex=F,only_a=F,type=2,karpack=22,lowcut=c(0.5,0.5,21))#;readline()
  dba2 = make_hyper_reorientation(dba2,icores=0,vcores=c(ncores=0))
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time superclus 2o filter:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
  
}

if(a_if){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("MAKE VISMODEL...")
  # NAO MASS AINDA... 
  #dba2[[1]][[2]] = add_res_information_to_clusters(dba2[[1]][[2]])
  seq.id = seq(2,by=2,length(dba2[[1]]))
  j = 1
  pdb.out = list()
  #browser()
  for (i in seq.id ){ #i = 4;
    dba2[[1]][[i]] = add_res_information_to_clusters(dba2[[1]][[i]])
    #browser()
     # area = 1586 #comentar 
     # area = area - 0.1 
     # browser()
    pdb.out[[j]] = write_molecule_muster_version(pdb1[[j]], dba2[[1]][[i]],
                                                 area.lim = area)
    j = j + 1
  }
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time make vismodel:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}


if(a_if){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("CONSOLIDATE RES TABLE...")
  
  seq.id = c(2)
  m.ligs1 = consolidate_res_table(dba2,seq.id)

  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time consolidate res table:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))
}

if(a_if){
  print(taskline)
  ttime = c(proc.time()[3],ttime);print("PLOT VISMODEL...")
  
  p.ligs1 = plot_lig_res_interactions(m.ligs1,scene="scene1",area.lim= area,alpha=0.1,n=1200)
  # p.ligs2 = plot_lig_res_interactions(m.ligs2,scene="scene2",area.lim=100,alpha=0.1,n=1200)
  # p.ligs3 = plot_lig_res_interactions(m.ligs3,scene="scene3",area.lim=100,alpha=0.1,n=1200)
  print(p.ligs1)
  # mol3d(paste0(lig.table$pdb_title, ".pdb"), lig.table$lig1id, check_style = input$molecule_style)

  
  dba2.size = object.size(dba2)
  ttime = c(proc.time()[3],ttime);print(paste("Time plot vismodel:",round((ttime[1]-ttime[2])/b,1),"segs"))
  print(paste("Size of dba2:",format(dba2.size,units="Mb")))

  
}
#-----------------------------------------END-------------------------------------------------

pdb1 <<- pdb1
dba2 <<- dba2
m.ligs1 <<- m.ligs1
#browser()
p.ligs1 <<- p.ligs1
lig.table <<- lig.table
# fim <<- T
}
