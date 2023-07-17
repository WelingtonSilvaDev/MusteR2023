#rm(list = ls())
begin <- Sys.time()
library("Matrix")
#library("igraph")
library("pracma")
library("cluster")
#library("rgl")
library("bio3d")
library("geometry")
library("deldir")
#library("ppls")
library("caret")
library("graphkernels")
library("stringr")
#library("rPython")
#library("Corbi")
library("shape")
library("prodlim")
library("rARPACK")
library("stringi")
library("RColorBrewer")
library("doMC")
library("mongolite")
library("tnet")
library("tidyverse")

taskline = "######################################################"
emptline = "                                                      "
prodline = "#-#-#-#-#-#-#-#- P R O D U C T I O N  -#-#-#-#-#-#-#-#"
testline = "#-#-#-#-#-#-#-#-       T E S T       #-#-#-#-#-#-#-#-#"
impoline = "                     GENERATING                       "
aligline = "                      ALIGNING                        "

#prodpath = "/data/compute/pdb_import"
prodpath = "/data/compute/pdb_import"
workpath = getwd()
ipaddr = "172.16.86.2"
port = ":27017"
addrpath = paste0(ipaddr,port)
localaddr = "localhost"
localpath = paste0(localaddr,port)
pdbupper = F

# Le da linha de comando

if(0){
	myArgs <- commandArgs(trailingOnly = TRUE)
	PDBS = myArgs[1]
	LIGs = myArgs[2]
	#PDBS = unlist(strsplit(PDBS, split=","))
	#workpath = myArgs[2]
	#LOCAL = as.numeric(myArgs[3])
}

#print(workpath)#;readline()

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
	bibpath = paste0("../bib/")
	#pdbpath = "/data/test/pdb/"
	pdbpath = "pdb/"
	print(taskline)
	print(emptline)
	print(testline)
	print(impoline)
	print(emptline)
	#print(taskline)
}

print(emptline)
#print(taskline)
print(paste("workpath:",workpath))
print(paste("bibpath:",bibpath))
print(paste("pdbpath:",pdbpath))
#print(taskline)
print(emptline)

#bibpath = "/data/compute/pdb_import/bib/"
#bibpath = "/data/compute/test_import/pdb_import/bib/"
#pdbpath = "/data/pdb/"

source(paste0(bibpath,"common.R"))
source(paste0(bibpath,"fun-base-v45.R"))
source(paste0(bibpath,"fun-mass-v45.R"))


#### variáveis ####################
ncores = 2
LIGs = c()
PDBS = c()
PDBS = c("bis-isothiourea-iNOSas-Docked-Pose1")
LIGs = c("UNK,0,HEM,1010")
#LIGs = c("UNL,1,HEM,1")
#PDBS = c("Cedr-8-15-ene-COX-1-dock-Pose1")
#PDBS = c("6LU7","5RGS","7brp")
#PDBS = c("6LU7","7brp")
#PDBS = c("INTZB8C")
#PDBS = c("1PPF","1ACB","1ROR","1HJA")
#PDBS = c("1PPF","1ACB","1ROR")
#PDBS = c("1PPF","1ACB")
#PDBS = c("1PPF","1ACB","1R0R","1TEC")
#PDBS = c("1PPF")
#PDBS = c("1ACB")
#PDBS = c("1ROR")
#PDBS = c("1HJA")
#PDBS = c("1VXA")
#PDBS = c("1PPF","1HJA")
#PDBS = c("1PPF","1TEC","1HJA")
#PDBS = c("1PPF","1HJA","1VXA","1ACB")
#PDBS = c("1PPF","1ACB")
#PDBS = c("1TAW","4SGB","1OYV")
#PDBS = c("1R0R","1TEC")
#PDBS = c("1R0R")
#PDBS = c("1TEC")
#PDBS = c("1VXA","1ACB")
#PDBS = c("1TAW","4SGB","1OYV","1R0R","1TEC")
#PDBS = c("1PPF","1ACB","1TAW","4SGB","1OYV","1R0R","1TEC")
#PDBS = c("1B2V","2O6P")
#PDBS = c("1D8U","1DLW","1KR7","1UT0","1HMD","1A7D","1B71")
#PDBS = c("1PPF","1HJA","1VXA","1ACB","1TAW","4SGB","1OYV","1R0R","1TEC","1B2V","2O6P","1D8U","1DLW","1KR7","1UT0","1HMD","1B71")
#PDBS = c("1HMD")
#PDBS = c("1PPF","1A7D")
#PDBS = c("1A7D")
#PDBS = c("1PPF")
#PDBS = c("1TEC")
#PDBS = c("1PPF","2HHDK")
#PDBS = c("2HHD")
#PDBS = c("1PPF","1ACB")
#PDBS = c("1TAL")
#PDBS = c("2SEC")
#PDBS = c("1PPFNEW")
#PDBS = c("1A22","1BXI","1DAN","1DVF","1FC2","1JRH","1VFB","2PTC","3HFM","1BRSA","1CBWA","1FCCA","1A4YA","1AHWA","1GC1A")

#PDBS=c("1PPF","1ACB","1TEC","1CSE","1MEE","1SBN","3SGB","1R0R")

#PDBS=c("1CHOA")

#PDBS=c("1VXA","1B2V")

#PDBS=c("1PPFZ")

#PDBS=c("1E96","1GRN","1J2J","1KAC","1KTZ","1OC0","1MAH","1YVB","1F34","1GLA","1WQ1","1ZHI","2A9K","2BTF","2HLE","2OOB","3BZD","3DAW","4M76")

#PDBS=c("4R2D")

LOCAL = 1

#"1B8D": exemplo de PDB muito grande
#"1A7D": tem hidrogenios (mas eh cristalografia)

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
#LOCAL=0
#
#prodpath = workpath
#LOCAL = as.numeric(LOCAL)
#print(LOCAL)

a=0

# Le da linha de comando
#if(a){
#	myArgs <- commandArgs(trailingOnly = TRUE)
#	PDBS = myArgs[1]
#	PDBS = unlist(strsplit(PDBS, split=","))
#}

# prepara para importacao
if(0){
  filename = "LIGS/leo-ligs-v1.csv"
  lig.table = read_delim(file=filename,delim=",")
  scriptname = "gapin-pdb2align.R"
  sep =","
  for (i in 1:dim(lig.table)[1]){
  #for (i in 19:19){
    PDBS = c(PDBS,lig.table$pdb_title[i])
    LIGs = c(LIGs,paste0(lig.table[i,2],sep,lig.table[i,3]))
    
    #for (i in 1:1){
    #parn = paste0(lig.table[i,2],sep,lig.table[i,3],sep,lig.table[i,4],sep,lig.table[i,5])
    #command = paste("Rscript",scriptname,lig.table$pdb_title[i],parn)
    #print(paste("RUNNING",command))
    #try(system(command))
    #browser()
  }
}

#browser()

# Importa...
if(0){
	print(taskline)
  pre.sufix = "b"
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
	db = initialize_db_list(PDBS,pdbpath,meso,aacode,nmax.atom,workpath,prodpath,ipaddr,
	                        vcores=c(ncores=2,ncores=2),pre.sufix=pre.sufix,sufix=pdb.sufix,
	                        conform_1_chain=T,ligs=LIGs,slim="-")
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
if(0){
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
if(0){
	print(taskline)
	ttime = c(proc.time()[3],ttime);print("EXPANDING BSR ALL TO AA AND PP...")
	db = expand_bsr_matrix(db,vcores=c(0,0))
	db.size = object.size(db)
	ttime = c(proc.time()[3],ttime);print(paste("Time expanding bsr:",round((ttime[1]-ttime[2])/b,1),"secs"))
	print(paste("Size of db:",format(db.size,units="Mb")))

}

# Monta dba a partir de db
if(0){
	print(taskline)
	ttime = c(proc.time()[3],ttime);print("BUILDING DBA FROM DB...")
	dba = build_dba_from_db(db,vdw)
	dba.size = object.size(dba)
	ttime = c(proc.time()[3],ttime);print(paste("Time building dba from db:",round((ttime[1]-ttime[2])/b,1),"secs"))
	print(paste("Size of db:",format(dba.size,units="Mb")))

}

# Saida CHAIN-CHAIN BSR ALL
if(0){
  print(taskline)
  print("GENERATING ADJACENCY MATRIX FILES... ")
  pdbids = seq(2,by=2,length.out=length(PDBS))
  #browser()
  for (i in pdbids){
    fileall = paste0("graph-out/",dba[[1]][[i]]$matrixname,".csv")
    print(fileall)
    aux = as.data.frame(as.matrix(dba[[1]][[i]]$A$a))
    colnames(aux) = dba[[1]][[i]]$element_name
    #rownames(aux) = colnames(aux)
    write.csv(aux,file=fileall)
  }
}

# Montagem do supercluster: 1a filtragem (descobre elementos conexos)
#TO DO: cercar erros quando AA ou PP for nulo
if(a){

	print(taskline)
	# logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=30,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "CLUSTERIZATIONS...")
	ttime = c(proc.time()[3],ttime);print("SUPERCLUS 1o FILTER...")
	dba = make_hyper_superclus(dba,icores=3,vcores=c(ncores,ncores,ncores,ncores,0),conex=T,only_a=T,type=1,karpack=22)
	dba.size = object.size(dba)
	ttime = c(proc.time()[3],ttime);print(paste("Time superclus 1o filter:",round((ttime[1]-ttime[2])/b,1),"segs"))
	print(paste("Size of dba:",format(dba.size,units="Mb")))

}

# Montagem do supercluster: 2a filtragem (considera elementos conexos acima de cutpar)
if(a){

	print(taskline)
	ttime = c(proc.time()[3],ttime);print("SUPERCLUS 2o FILTER...")
	dba2 = make_hyper_filter(dba,icores=0,vcores=c(0),cutpar=10)
	dba2 = make_hyper_superclus(dba2,icores=3,vcores=c(ncores,ncores,ncores,ncores,ncores),conex=F,only_a=F,type=2,karpack=22,lowcut=c(0.5,0.5,21))#;readline()
	dba2 = make_hyper_reorientation(dba2,icores=0,vcores=c(ncores=0))
	dba2.size = object.size(dba2)
	ttime = c(proc.time()[3],ttime);print(paste("Time superclus 2o filter:",round((ttime[1]-ttime[2])/b,1),"segs"))
	print(paste("Size of dba2:",format(dba2.size,units="Mb")))

}

if(a){

	print(taskline)
	# logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=80,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "SETTING SPOTS...")
	ttime = c(proc.time()[3],ttime);print("PROTEIN OUTPUTS...")
	dbpro = list()
	dbpro = make_hyper_protein_output(dba2,ncores,maxcut=20)
	dbpro.size = object.size(dbpro)
	ttime = c(proc.time()[3],ttime);print(paste("Time protein outputs:",round((ttime[1]-ttime[2])/b,1),"segs"))
	print(paste("Size of dbpro:",format(dbpro.size,units="Mb")))

}

#x=fix_labels("1-1",c(14),c(12,13));print(x)
# define hotspots
if(a){

	print(taskline)
	# logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=95,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "SAVING...")
	ttime = c(proc.time()[3],ttime);print("SETTING HOTSPOTS...")
	dbsel = list()
	dbsel = make_hyper_hotspot(dba2,icores=3,vcores=c(0,ncores,ncores),lowcut=c(0.50,0.50),minclusid=3,maxclusid=21)#lowcut=c(0.30,0.40)
	dbsel.size = object.size(dbsel)
	ttime = c(proc.time()[3],ttime);print(paste("Time hotspot:",round((ttime[1]-ttime[2])/b,1),"segs"))
	print(paste("Size of dbsel:",format(dbsel.size,units="Mb")))

}

if(a){

	print(taskline)
	ttime = c(proc.time()[3],ttime);print("GROUP OUTPUTS...")
	dbgrp = list()
	dbgrp = make_hyper_group_output(dba2,dbsel,maxcut=20,ncores=0)
	dbgrp.size = object.size(dbgrp)
	ttime = c(proc.time()[3],ttime);print(paste("Time group outputs:",round((ttime[1]-ttime[2])/b,1),"segs"))
	print(paste("Size of dbgrp:",format(dbgrp.size,units="Mb")))

}


if(0){
	#clean_database("proteins",workpath,prodpath,addrpath,localpath)
  toMongo(dbpro,"proteins",workpath,prodpath,addrpath,localpath)
	#clean_database("groups",workpath,prodpath,addrpath,localpath)
  toMongo(dbgrp,"groups",workpath,prodpath,addrpath,localpath)
	logger(pdb=PDBS[1], type=1, key=PDBS[1], percentage=100,wp=workpath,pp=prodpath,ap=addrpath,lp=localpath, "FINISHED IMPORT!")
}


###################### VISUALIZADOR #######################

if (0){

	plot_super(dba2[[3]][[2]]$exp,2:7,mfrow=c(2,3),dbgrp[[1]][[2]]$spot_score_bsr_nonpolar_nonpolar)


}

###################### TESTES #######################

if(0){

	pdbid="1PPF"
	x=copy_to_server(pdbid,workpath,prodpath,ipaddr,pdbpath)
	#scp data/pdb/1ACB2.pdb 172.16.86.2:/data/pdb


}



if (0){

	i=3
	#cols = c("black",colorlist$colors)
	cols = NULL
	dxy=c(4,2)
	rotref = diag(4)#;print("0")
	if (0){
		visualize_align_super_cluster(db=dba2[[i]],dball=dball[[i]],pdbid=1,clusid=0,dbhot=dball[[i]]$hot,ids=c(1:4),score=1,cols=cols,mfrow=c(2,2),dxy=dxy,xlim=xlim,ylim=ylim,rotref=rotref,self=F)

	}
	if (0){
		visualize_align_supercluster_prescore(db=dba2[[i]],dball=dball[[i]],dbhot=dball[[i]]$hot,pdbnames=c("1PPF","1TEC"),cols=cols,dxy=dxy, rotref=rotref,self=F)
	}
	if (1){
		pdbid=2
		visualize_align_supercluster_hot(db=dba2[[i]],dball=dbsel[[i]],dbhot=dbsel[[i]]$hot,pdbid=pdbid,cols=cols,dxy=dxy,rotref=rotref)
	}
	#dev.off()
}

if (0){

	#plot_density(dbgrp,"1PPF")
	plot_density(dbgrp,"1BV2")

}


if(0){

	
	x = c(1:10,rep(142,10),rep(147,10),rep(142,10),rep(147,10),rep(142,10))
	z = c(letters[1:10],rep("a",10),rep("b",10),rep("c",10),rep("d",10),rep("e",10))
	xz = data.frame(x,z)
	print(xz)
	xmax = max(xz[,1])
	#k = 0
	d = 100
	#while (k<xmax) k = k + d
	z = unique(xz) 
	y = unique(xz)
	print(y);print(k)
	tam = dim(y)[1]
	yf = factor(y[,1])
	yfl = levels(yf);print(yf)
	
	for (i in yfl){
		yi = which(y[,1]==i);print(yi);readline()
		tamj = length(yi)
		if (tamj>1){
			for (j in 2:tamj){
				dy = y[yi[j],1]+d#;print(dy);readline()
				repeat{
					if (dy %in% y[,1]){
						dy = dy+d
					}else{
						y[yi[j],1] = dy
						break;
					}
				}
				#print(xmax);print(dy);print(y);readline()
			}
		}
		#readline()
	}
	print(z);print(y)
	for (i in 1:tam){
		f = (xz[,1] == z[i,1])&(xz[,2] == z[i,2])
		xz[f,1] = y[i,1]
	}
	print(xz)
	

	;readline()
}

if(0){

	x=c("A.GLY1.O-4-O-A-b","A.GLY1.H1-5-H-A-o","A.GLY1.H1-5-H-A-h")
	m = data.frame(matrix(0,ncol=3,nrow=3))
	m[1,2]=11
	m[2,1]=11
	m[1,3]=13
	m[3,1]=13
	colnames(m)=x

	y=consider_hydrogens(m)
	print(m);print(y)


}


if(0){

	x = read.pdb("/data/pdb/1A7D.pdb")
	y = read.pdb("/data/pdb/1PPF.pdb")
	library("Rpdb")
	z = read.pdb("/data/pdb/1A7D.pdb",ATOM = F,HETATM = T,CRYST1 = F,CONECT = F,TITLE = T,REMARK = F,MODEL = NA)
	detach(package:Rpdb)

	y1 = clean.pdb(y)

}

if(0){

	x=list(a=1,b=2,c=3,d=4,e=5)
	y=c(T,T,F,F,T)
	y=c(F,F,F,F,F)
	y=c(T,T,T,T,T)

	print(x)
	for (i in 5:1){
		if (y[i]){
			x[i]=NULL
		}
	}
	print("---")
	print(x)
}

if(0){

	m = Matrix(c(1.2,6.0,6.0,1.2),ncol=2,byrow=T)
	colnames(m)=c("E.ILE16.N-1","E.ILE16.CA-2")
	x=silveira_romanelli_transformation(m,vdw,probe)


}


if (0){

	x=colnames(db[[1]]$polarity[[1]]$bsr[[1]])
	y=new_codify_polarity3(x,meso,aacode)


}

if (0){

	x=calculate_BSR_matrix2(all,vdw,probe,vcores[3])

}


#for(i in  1:length(PDBS)){
#  distance_matrix = calculate_distance_matrix(PDBS[i],pdbpath,ncores)
  #chamar o metodo do silveira
  
  
  # loop sobre todas as combinações apenas para gravar em arquivo, vai sumir
  #for(i in 1:length(distance_matrix$combinations[1,])){
  #  chain_a_name = distance_matrix$combinations[1,i]
  #  chain_b_name = distance_matrix$combinations[2,i]
  #  
  #  m = Matrix(distance_matrix$dist[[i]])
    
    # aqui chamaria o novo método para cácular tudo passando as matrizes
    # writeMM(m, file=paste0('/Users/biharck/Desktop/', distance_matrix$pdb, '_',chain_a_name, '_' , chain_b_name, '_BSR'))
  #}
  
  
#}

print(taskline)
print(paste("Size all: db:",format(db.size,units="Mb"),
						",dba:",format(dba.size,units="Mb"),
						",dba2:",format(dba2.size,units="Mb"),
						",dbpro:",format(dbpro.size,units="Mb"),
						",dbgrp:",format(dbgrp.size,units="Mb"),
						",dbsel:",format(dbsel.size,units="Mb"),
						",Total:",format(db.size+dba.size+dba2.size+dbpro.size+dbgrp.size+dbsel.size,units="Mb")))
end <- Sys.time() - begin
print(end,1)
print('IMPORT FINISHED')

#vscode

###################### DEPRECIATED #######################

######################################################################################################
# BIHARCK, ESSA PARTE DO (TOMONGODB), ACHO MELHOR VOCÊ FAZER !
# DBHOT E DWEB TEM A SEGUINTE ESTRUTURA:
# dbhot[[polarity]][[pdbids]]][[superclus]]
# ONDE: 
# polarity=1, ALL
# polarity=2, PP
# polarity=3, AA
# SE #PDBS = c("1PPF","1ACB")
# ENTAO:
# pdbids=1, 1PPF
# pdbids=2, 1ACB
# FINALMENTE, SUPERCLUS:
# superclus=1, 3 superclusters
# superclus=2, 4 superclusters
# and so on. 

if (0){
	print(taskline)
  ttime = c(proc.time()[3],ttime)
	ToMongoDB(dbweb=dbweb,collection="alignments")
	ToMongoDB(dbweb=dbhot,collection="hotspots")
	ttime = c(proc.time()[3],ttime);print(paste("Time ToMongoDB:",round((ttime[1]-ttime[2])/b,1),"segs"))
  
}

#Make alignments
if(0){
	print(taskline)
	ttime = c(proc.time()[3],ttime);print("MAKING ALIGNMENTS...")
	dball = list()
	dball = make_hyper_alignment(dba2,dbsel,icores=3,vcores=c(ncores,ncores),lowcut=c(0.50,0.60),minclusid=4,maxclusid=6)
	ttime = c(proc.time()[3],ttime);print(paste("Time alignment:",round((ttime[1]-ttime[2])/b,1),"segs"))
	print(paste("Size of dball:",format(object.size(dball),units="Mb")))
}

#Hotspots outputs
if(0){

	print(taskline)
	ttime = c(proc.time()[3],ttime);print("HOTSPOT OUTPUTS...")
	dbhot = list()
	dbhot = make_hyper_web_output(dba2,dbsel,dbgrp,icores=0,vcores=c(0))
	dbhot.size = object.size(dbhot)
	ttime = c(proc.time()[3],ttime);print(paste("Time hotspot outputs:",round((ttime[1]-ttime[2])/b,1),"segs"))
	print(paste("Size of dbhot:",format(dbhot.size,units="Mb")))

}

#alignment outputs
if(0){

	print(taskline)
	if (exists("dball")){
		ttime = c(proc.time()[3],ttime);print("ALIGNMENT OUTPUTS...")
		dbweb = list()
		dbweb = make_hyper_web_output(dba2,dball,icores=3,vcores=c(0))
		dbweb.size = object.size(dbweb)
		ttime = c(proc.time()[3],ttime);print(paste("Time alignment outputs",round((ttime[1]-ttime[2])/b,1),"segs"))
		print(paste("Size of dbweb:",format(dbweb.size,units="Mb")))
	}

}


