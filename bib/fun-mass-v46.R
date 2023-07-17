#####################################################################
#Biblioteca funções em massa, aplicada sobre todos os complexos
#####################################################################


#	tam = dim(dbareanz[[1]]$Lw$eig$v)[1]
#	filterv=as.factor(dbareanz[[1]]$Lw$eig$v[,tam])
#	levels(filterv) = c(TRUE,TRUE)
#	dbareanzu[[1]] = apply_filter(dbareanz[[1]],as.logical(filterv))


#expid = 4
#k = 4

#filebase = "einzaa-Lw-eig-pam"
#filebase = paste(filebase,k,sep="")
#vertexsufix = make_vertex_sufix(k,c("e","i"))

#dbareanzaa = mass_experiment_02(dbareanzaa,expid,k)
#par(mfrow=c(1,2))
#plot(dbareanzaa[[1]]$exp$pam[[expid]])
#make_mass_pymol_partition(dbareanzaa,1,expid,dirpathin,dirpathout,filebase,edgename,colors,vertexsufix,edpart)


#################################################### HYPER ####################################

###NEW EXP###
make_mass_taball_sasas_areas = function(lote){

	tam = length(lote)
	taball.sasas = c()
	taball.areas = c()
	for (i in 1:tam){
			taball.sasas = rbind(taball.sasas,lote[[i]]$tab.sasas)
			taball.areas = rbind(taball.areas,lote[[i]]$tab.areas)
	}
	res = list()
	res$sasas=taball.sasas
	res$areas=taball.areas
	return(res)
}

###NEW EXP###
make_mass_sasas_areas_exp = function(PDBS,pdbpath,db,dba2,polar_type="all",sufix=".pdb",nmax.atom=6000){

	#print("o")
	pdb = make_mass_pdb_readling(PDBS,pdbpath,sufix=sufix,nmax.atom=nmax.atom)
	#print("o")
	sasas = make_mass_sasas(pdb,db)
	tab.sasas = make_mass_areas_from_sasas(PDBS,sasas,polar_type=polar_type)
	#print("o")
	acons = make_mass_acons(dba2,tam=length(PDBS),polar_type=polar_type)
	tab.areas = make_mass_areas_from_acons(PDBS,acons)
	#print("o")
	res=list()
	res$pdb = pdb
	res$sasas = sasas
	res$tab.sasas = tab.sasas
	res$acons = acons
	res$tab.areas = tab.areas

	return(res)

}


###NEW EXP###
make_mass_acons = function(db,tam,polar_type,pol=1,type=2){

	#tam = length(db)
	aux = list()
	#tam=1

	for (i in 1:tam){
		aux[[i]] = make_acons(db[[pol]][[type*i]],polar_type=polar_type)
		#aux[[i]] = make_areas(dba,i)

	}
	return(aux)
}

###NEW EXP###
make_mass_areas_from_acons = function(pdbids,acons){

	tam = length(db)
	#aux = list()
	auxn = c()
	auxa = c()

	for (i in 1:tam){
		auxn = c(auxn,pdbids[i])
		auxv = make_areas_from_acons(acons[[i]])
		auxa = rbind(auxa,auxv)

	}
	auxr = data.frame(pdbid=auxn,auxa)
	rownames(auxr)=NULL
	return(auxr)

}


###NEW EXP###
make_mass_areas_from_sasas = function(pdbids,sasas,polar_type){

	tam = length(sasas)
	#aux = list()
	auxn = c()
	auxa = c()

	for (i in 1:tam){
		auxn = c(auxn,pdbids[i])
		auxv = make_areas_from_sasas(sasas[[i]],polar_type=polar_type)
		auxa = rbind(auxa,auxv)

	}
	auxr = data.frame(pdbid=auxn,auxa)
	rownames(auxr)=NULL
	return(auxr)

}



###NEW EXP###
make_mass_areas = function(dba,pol=1,type=2){

	tam = length(pdb)
	aux = list()
	#tam=1

	for (i in 1:tam){
		aux[[i]] = make_areas(dba[[pol]][[type*i]],pol)
		#aux[[i]] = make_areas(dba,i)

	}
	return(aux)
}

###NEW EXP###
make_mass_sasas = function(pdb,db,type=2){

	tam = length(pdb)
	aux = list()

	for (i in 1:tam){;print(i)
	#for (i in 3){
		aux[[i]] = make_sasa(pdb[[i]],colnames(db[[i]]$dist[[type]]))

	}
	return(aux)
}


###NEW EXP###
make_mass_pdb_readling = function(pdbids,pdbpath,sufix=".pdb",nmax.atom=6000){

	aux = list()
	tam = length(pdbids)
	#tam=12
	for(i in 1:tam){#;print(i)
		aux[[i]] = make_it_mass_pdb_reading(pdbids[i],pdbpath,sufix=sufix,nmax.atom=nmax.atom)
	}
	return(aux)
}

###NEW EXP###
make_mass_term_par = function(PDBS,sasas){

	auxn = c()
	auxv = c()
	auxa = c()
	tam = length(PDBS)
	for (i in 1:tam){
		auxn = c(auxn,PDBS[i])
		auxv = make_term_par(sasas[[i]])
		auxa = rbind(auxa,auxv)
	}
	auxr = data.frame(pdbid=auxn,auxa)
	rownames(auxr)=NULL
	return(auxr)
	#print(auxr);readline()
}

make_mass_area_par =function(PDBS,areas){

	auxn = c()
	auxv = c()
	auxa = c()
	tam = length(PDBS)
	for (i in 1:tam){
		auxn = c(auxn,PDBS[i])
		auxv = make_areas_par(areas[[i]])
		auxa = rbind(auxa,auxv)
	}
	auxr = data.frame(pdbid=auxn,auxa)
	rownames(auxr)=NULL
	return(auxr)


}

###NEW EXP###
make_mass_table_areas = function(sasas,db,pol=1,type=2){

	tam = length(sasas)
	#tam=1
	aux = list()

	for (i in 1:tam){	
		aux[[i]] = make_table_areas(sasas[[i]],db[[pol]][[type*i]])
	}
	return(aux)
	
}

#build_mass_pdb_colname = function(pdb){

#	aux = list()
#	tam = length(pdb)
#	for(i in 1:tam){
#		aux[[i]] = build_pdb_colname(pdb)
#	}
#	return(aux)

#}

make_it_mass_pdb_reading = function(pdbname,dirpathin,prefix="",sufix=".pdb",file,
                                    maxlines = -1, multi = FALSE, rm.insert = TRUE, 
                                    rm.alt = TRUE, ATOM.only = TRUE, verbose = TRUE,
                                    nmax.atom=5000, slim=NULL){


  slimtitle = ""
	auxpdb=list()
	#pdbname = as.character(unlist(pdbname))
	auxfile = paste(dirpathin,prefix,pdbname,sufix,sep="")#;print(auxfile);readline()
	#browser()
	if (!file.exists(auxfile)){
		auxfile = paste(dirpathin,prefix,toupper(pdbname),sufix,sep="")
	}
	if (!file.exists(auxfile)){
		auxfile = paste(dirpathin,prefix,tolower(pdbname),sufix,sep="")
	}
	if (!is.null(slim)){
	  slimname = unlist(str_split(pdbname,slim))
  	if (length(slimname)>1){
  	  slimname = str_c(str_extract(slimname,"^."),str_extract(slimname,".$"))
  	  slimname = paste0(str_replace_all(slimname, "(.)\\1", "\\1"),collapse = "")
  	  slimname = toupper(slimname)
  	  slimfile = paste(dirpathin,prefix,slimname,sufix,sep="")
  	  if (!file.exists(slimfile)){
  	  #browser()
  	    command = paste0("cp ",auxfile," ",slimfile)
  	    try(system(command))
  	  }
  	  auxfile = slimfile
  	  slimtitle = pdbname
  	  pdbname = slimname
  	  #pdb$pdbname = pdbname
	  }
	}
	#browser()
	if (file.exists(auxfile)){
		#auxpdb = read.pdb(auxfile,ATOM.only=TRUE)
	  library(bio3d)
		auxpdb = bio3d::read.pdb(auxfile,maxlines = maxlines, multi = multi, rm.insert = rm.insert, rm.alt = rm.alt, ATOM.only = ATOM.only, verbose = verbose)
		tampdb = dim(auxpdb$atom)[1]
		if (tampdb>nmax.atom){
			aux = list()
			aux$error = 1
			aux$tampdb = tampdb
			aux$pdbname = pdbname
			return(aux)
		}
		tamchains = length(unique(auxpdb$atom$chain))
		tamhet = length(unique(auxpdb$atom$type))#;print(tamchains);print(tamhet);readline()
		if ((tamchains==1)&(tamhet==1)){
			aux = list()
			aux$error = 2
			aux$pdbname = pdbname
			return(aux)
		}
		library("Rpdb")
		if (tamhet>1){
			extra = read.pdb(auxfile,ATOM = F,HETATM = T,CRYST1 = F,CONECT = F,TITLE = T,REMARK = F,MODEL = NA)
			auxpdb$ligands = extra$atoms
		} else {
			extra = read.pdb(auxfile,ATOM = T,HETATM = T,CRYST1 = F,CONECT = F,TITLE = T,REMARK = F,MODEL = NA)
		}
		extra$title = remove_extra_spaces(extra$title,"TITLE")
		auxpdb$title = extra$title
		#browser()
		if (auxpdb$title==""){
  		if (slimtitle!=""){
  		  auxpdb$title = slimtitle
  		}
		}
		#auxpdb$ligands = extra$atoms
		auxpdb$pdbname = toupper(pdbname)
		auxpdb$chains = unique(auxpdb$atom$chain)
		names(auxpdb$chains) = rep("C",length(auxpdb$chains))
		detach(package:Rpdb)
		auxpdb$atom$elety = gsub("'","",auxpdb$atom$elety) #case 4R2D
    auxpdb$ligands$elename = gsub("'","",auxpdb$ligands$elename) #case 4R2D
		#auxpdb$title=""
		#auxpdb$ligands=""
	}else{
		print(paste("WARNING: file",auxfile,"do not exist"))
		aux = list()
		aux$error = 0
		aux$pdbname = pdbname
		return(aux)
	}
	#browser()
	return(auxpdb)


}

make_hyper_pdb_reading = function(dirname,pdbfile,dirpathin,icores,vcores=0,sep=","){

	pdbids = read.csv(paste(dirname,pdbfile,sep=""),header=F,sep=sep)
	pdbids = as.character(unlist(pdbids))

	tam = length(pdbids)
	aux = list()
	print(paste("Input PDB data..."))

	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_mass_pdb_reading(pdbids[i],dirpathin)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			aux[[i]] = make_it_mass_pdb_reading(pdbids[i],dirpathin)
		}

	}
	return(aux)

}

#make_it_mass_input(dirname,filenames[i],pdbxyz[[i]],meso,aacode,vcores)
make_it_mass_input = function(dirname,filename,pdbxyz,meso,aacode,vcores){

	pdbids = read.csv(paste(dirname,filename,sep=""),header=F,sep=",")
	pdbids=paste(as.matrix(pdbids))
	names(pdbids)=1:length(pdbids)
	#Entrada 
	dba=list()
	#Agora tem opcao de usar matrizes esparsas em matrixtype
	dba = make_mass_input(dirname=dirname,pdbids=pdbids,sufix=".csv",meso=meso,aacode=aacode,ncores=vcores[1],matrixtype="sparse", bin=FALSE,sep=",")#;readline()
		 
	#Elimina linhas/colunas zeros, se houver
	dba = make_mass_only_connected(dba,ncores=vcores[2])#;readline()
	#Guarda em ...$pdb as coordenadas atomicas (PARALELIZADO! Ficou entre 3x a 4x mais rapido...)
	if (length(dba)>1){
		dba = add_mass_pdb_xyz(dba,pdbxyz,ncores=vcores[3],typename=2)
	}else{
		dba = add_mass_pdb_xyz(dba,pdbxyz,ncores=vcores[4],typename=2)
	}
	return(dba)

}

#	dba = make_hyper_mass_input(dirname,filenames,pdbxyz,meso,aacode,icores=3,vcores=c(nc,nc,dc,nc))
make_hyper_mass_input = function(dirname,filenames,pdbxyz,meso,aacode,icores,vcores,sep=","){

	#dba = make_hyper_mass_input(dirname,pdbxyz,meso,aacode,icores=3,vcores=c(nc,nc,dc,nc))
	tam = length(filenames)
	#tam = length(pdbxyz)
	aux = list()
	print(paste("Input hyper data..."))

	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_mass_input(dirname,filenames[i],pdbxyz,meso,aacode,vcores)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			aux[[i]] = make_it_mass_input(dirname,filenames[i],pdbxyz,meso,aacode,vcores)
		}

	}
	return(aux)

}

make_it_superclus = function(dba,vcores,conex,only_a,type=2,karpack,lowcut){

	dba = make_mass_laps(dba,ncores=vcores[1],anorm=TRUE,r=10,only=c(F,T))
	# Faz decomposicao EIG em Lw, agora parcial com pacote rARPACK
	dba = make_mass_decomps(dba,L=FALSE,Lw=TRUE,EIG=TRUE,SVD=FALSE,ncores=vcores[2],karpack=karpack,r=10)
	# Clusterizacao por Spectral Clustering na lista ...$exp$pam
	# type=1: PAM e EIG-Lw
	# conex=T: considera separacao apenas de clusteres perfeitos (elementos nao conexos)
	# Maxrange=15: tenta separar ate no maximo 15 grupos perfeitos
	#if(0){
	if (!only_a){
		dba=make_mass_centrality(dba,ncores=vcores[5])
	}

	if (type==1){
		dba = make_mass_pam_clustering(dba,maxrange=10,type=1,conex=conex,ncores=vcores[3])#maxrange era 15
	# Gera dados supercluster na lista ...$exp$superclus
		dba = make_mass_superclus(dba,ncores=vcores[4],only_a=only_a)
	}else{
		dba = make_mass_pam_superclus(dba,lowcut=lowcut,conex=conex,ncores=vcores[3])
	}
	#}
	return(dba)

}

make_hyper_superclus = function(dba,icores=3,vcores,conex,only_a,type,karpack,lowcut){

	tam = length(dba)
	aux = list()
	print(paste("Hyper superclus..."))
	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_superclus(dba[[i]],vcores,conex,only_a,type,karpack,lowcut)
		}
	}else{
		print(paste("Doing it sequential..."));
		for (i in 1:tam){#;print(i)
		#for (i in 3){
			aux[[i]] = make_it_superclus(dba[[i]],vcores,conex,only_a,type,karpack,lowcut)
		}
	}
	return(aux)
}

make_hyper_filter = function(dba,icores=3,vcores,cutpar){

	tam = length(dba)
	aux = list()
	print(paste("Hyper filter..."))
	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_mass_cluster_filter(dba[[i]],ncores=vcores[1],cutpar=cutpar) 
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			aux[[i]] = make_mass_cluster_filter(dba[[i]],ncores=vcores[1],cutpar=cutpar) 
		}
	}
	return(aux)

}

make_hyper_reorientation = function(dba,icores=3,vcores){

	tam = length(dba)
	aux = list()
	print(paste("Hyper reorientation..."))
	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_mass_superclus_reorientation(dba[[i]],ncores=vcores[1])
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
		#for (i in 3){
			aux[[i]] = make_mass_superclus_reorientation(dba[[i]],ncores=vcores[1])
		}
	}
	return(aux)

}


make_it_mass_hotspot = function(dba,vcores,lowcut,minclusid,maxclusid){

	dbsel = list()
	#dbsel[[i]]=list()
	dbsel = make_mass_superclus_tables(db=dba,dball=dbsel,lowcut=lowcut,minclusid=minclusid,maxclusid=maxclusid,ncores=vcores[1])
	#print(dbsel);readline()
	dbsel = make_mass_self_similarity_assessment(db=dba,dball=dbsel,ncores=vcores[2])#;print("o")
	dbsel = set_mass_hotspot_candidates(db=dba,dball=dbsel,ncores=vcores[3])#;print("o")
	return(dbsel)

}


make_hyper_hotspot = function(dba,icores=3,vcores,lowcut,minclusid,maxclusid){

	tam = length(dba)
	aux = list()
	print(paste("Hyper hotspot..."))
	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_mass_hotspot(dba[[i]],vcores,lowcut,minclusid,maxclusid)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){#;print(paste("i: ",i))
		#for (i in c(3)){
			#print(paste("-------------------------------------->",i))
			aux[[i]] = make_it_mass_hotspot(dba[[i]],vcores,lowcut,minclusid,maxclusid)
		}

	}
	return(aux)
}

make_mass_centrality = function(db,ncores=0){

	tam = length(db)
	aux = list()

	if(ncores){
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			make_centrality(db[[i]])
		}
	}else{
		for (i in 1:tam)		
			aux[[i]] = make_centrality(db[[i]])
	}
	return(aux)
}

make_pair_graph_alignment = function(db,pdbids,type,vcores,lowcut=c(0.50,0.50),cluspar=c(4,2,0.9),inter=c("LIG","CHAIN"),pol=c("ALL","PP","AA")){

	#tam = length(dba)
	aux = list()
	res = list()
	print(paste("Pair aligment..."))
	pdblist = unlist(mapply("[[",db,1)[1,])#;print(pdblist);print(pdbids);print(pdbids);readline()
	pdblist = mapply("[[",strsplit(pdblist,"_"),1)
	str1 = paste0("^",pdbids[1],"$")
	str2 = paste0("^",pdbids[2],"$")
	#pdb1 = grep(pdbids[1],pdblist);print(pdb1)
	#pdb2 = grep(pdbids[2],pdblist);print(pdb2);readline()
	pdb1 = grep(str1,pdblist);print(pdb1)
	pdb2 = grep(str2,pdblist);print(pdb2);readline()

	db1 = db[pdb1]
	db2 = db[pdb2]
	id1 = c(1,which(type[1]==inter),which(type[2]==pol))#;print(ids);readline()
	id2 = c(1,which(type[1]==inter),which(type[2]==pol))
	pretab1 = make_pre_superclus_table(db1,id1,lowcut)
	pretab2 = make_pre_superclus_table(db2,id2,lowcut)
	print(pretab1);print(pretab2);readline()
	tamtabs = dim(pretab1)[1]+dim(pretab2)[1]#;print(tamtabs)

	n = 3
	repeat{
		tab = list()
		tab = try_best_pair_tabs(pdbids,pretab1,pretab2,n)
	
		tab1 = tab$t1
		tab2 = tab$t2
		#print(tab1);print(tab2);print(tab$p1);print(tab$p2)#;readline()

		tabids = expand.grid(1:dim(tab1)[1],1:dim(tab2)[1])
		#print(tab1$group[tabids[,1]]);readline()
		tabids = cbind(tabids,tab1$group[tabids[,1]],tab2$group[tabids[,2]])
		colnames(tabids) = c("id1","id2","grp1","grp2")
		tabids = filter(tabids,grp1==grp2)
		tam = dim(tabids)[1]#;print(tabids);print(tam);readline()
		#print(tab1);print(tab2);print(tab$p1);print(tab$p2);print(tabids);readline()
	
		#tali = cbind(as.character(tab1[tabids[,1],c(1,5:7)]),as.character(tab2[tabids[,2],c(1,5:7)]))
		tali = cbind(tab1[tabids[,1],c(1,6:7)],tab2[tabids[,2],c(1,6:7)])
		#print(tab1);print(tab2);print(tab$p1);print(tab$p2);print(tali);print(tam);readline()
		if (tam){
			break;
		} else {
			n = n + 1
			if (n==tamtabs) break;
		}
	}
	#print(tab1);print(tab2);print(tab$p1);print(tab$p2);print(tali);print(tam);readline()
	if (tam){
		if (vcores[1]){
			print(paste("Doing it parallel with",vcores[1],"cores"))
			registerDoMC(vcores[1])
			aux = foreach(i=1:tam) %dopar% {
				make_graph_alignment(db1,db2,tab1[tabids[i,1],],tab2[tabids[i,2],],vcores)
			}
		}else{
			print(paste("Doing it sequential..."))
			for (i in 1:tam){
			#for (i in c(3)){
				aux[[i]] = make_graph_alignment(db1,db2,tab1[tabids[i,1],],tab2[tabids[i,2],],vcores)
			}

		}
		#print("ok")
		res$tabs = list()
		res$tabs$ngroup = n
		res$tabs$pre1 = tab$p1
		res$tabs$pre2 = tab$p2
		res$tabs$pos1 = tab1
		res$tabs$pos2 = tab2
		res$tabs$aligned = tali
		res$align = aux
	}else{
		print(paste(c("WARNING: it was not possible to align",pdbids,"with types",type),collapse=" "))
	}
	return(res)
}


make_it_mass_alignment = function(dba,dbsel,vcores,lowcut,minclusid,maxclusid){

	dball = list()
	#dball[[i]]=list()
	dball = make_mass_superclus_tables(db=dba,dball=dball,lowcut=lowcut,minclusid=minclusid,maxclusid=maxclusid,ncores=vcores[1])
	#print(dball);readline()
	dball = make_mass_similarity_assessment(db=dba,dball=dball,dbsel=dbsel,id=1,self=F,ncores=vcores[2])
	return(dball)

}


make_hyper_alignment = function(dba,dbsel,icores=3,vcores,lowcut,minclusid,maxclusid){

	tam = length(dba)
	aux = list()
	print(paste("Hyper alignment..."))
	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_mass_alignment(dba[[i]],dbsel[[i]],vcores,lowcut,minclusid,maxclusid)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
		#for (i in c(3)){
			aux[[i]] = make_it_mass_alignment(dba[[i]],dbsel[[i]],vcores,lowcut,minclusid,maxclusid)
		}

	}
	return(aux)
}


make_it_web_output = function(dbsel,dbgrp,vcores){

	tam = length(dba)
	aux=list()
	for (i in 1:tam){
		aux[[i]] = resume_for_web_output(dbsel,dbgrp,id=i)	
	}
	return(aux)
	#resume_for_web_output(dba,dbsel,id)

}

#x = resume_for_web_output(db=dba2[[1]], dball=dbsel[[1]][[1]],id=1);print(x)
make_hyper_web_output = function(dbsel,dbgrp,icores=3,vcores){

	tam = length(dbsel)
	aux = list()
	print(paste("Outputs..."))
	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_web_output(dbsel[[i]],dbgrp[[i]],vcores)
		}
	}else{
		print(paste("Doing it sequential..."))
		#for (i in 1:tam){
		for (i in 1){
			aux[[i]] = make_it_web_output(dbsel[[i]],dbgrp[[i]],vcores)
		}

	}
	return(aux)

}


make_it_protein_output = function(i,dba,maxcut,sep="_",n=2,k=3){

	aux = list()
	#aux$pdb = toupper(substr(dba[[1]][[i]]$pdbname,1,4))
	comb = unlist(strsplit(dba[[1]][[i]]$pdbname,sep))
	aux$pdb = comb[1]
	if (i%%n==1) aux$pdb_type = 2
	else aux$pdb_type = 1
	aux$title = dba[[1]][[i]]$title
	aux$chains = dba[[1]][[i]]$chains
	#comb = unlist(strsplit(dba[[1]][[i]]$pdbname,sep))#;print(comb) 
	aux$combination_chains = comb[2:3]
	aux$rotall = dba[[1]][[i]]$exp$rotref
	aux$rotpp = dba[[2]][[i]]$exp$rotref
	aux$rotaa = dba[[3]][[i]]$exp$rotref
	aux$bsr = data.frame()
	aux$bsr_elements = data.frame()
	aux$bsr_groups = list()
	aux$bsr_groups_polar_polar = list()
	aux$bsr_groups_nonpolar_nonpolar = list()

#round(transform_by_rot_matrix(geomc,rotref),k)
	
	if (!is.null(dba[[1]][[i]]$A)){

		aux$bsr = get_bsr_info_interactions(dba[[1]][[i]]$A$a,dba2[[1]][[i]]$element_name,dba[[1]][[i]]$element_polarity)
		xyz = dba[[1]][[i]]$pdb$xyz
		##xyz = round(transform_by_rot_matrix(xyz,aux$rotall),k)
		aux$bsr_elements = get_bsr_info_elements(dba[[1]][[i]]$element_name,dba[[1]][[i]]$element_polarity,dba[[1]][[i]]$vdw,xyz)
		aux$bsr_groups = get_bsr_info_groups(dba[[1]][[i]]$element_name,dba[[1]][[i]]$exp$pam,maxcut)

	}

	if (!is.null(dba2[[2]][[i]]$A)){
		aux$bsr_groups_polar_polar = get_bsr_info_groups(dba[[2]][[i]]$element_name,dba[[2]][[i]]$exp$pam,maxcut)
	}

	if (!is.null(dba2[[3]][[i]]$A)){
		aux$bsr_groups_nonpolar_nonpolar = get_bsr_info_groups(dba[[3]][[i]]$element_name,dba[[3]][[i]]$exp$pam,maxcut)
	}

	return(aux)

}


make_hyper_protein_output = function(dba,ncores,maxcut){

	tam = length(dba[[1]])

	#tam_pdb = tam/3
	proteins = list()
	aux = list()
	print(paste("Protein outputs..."))
	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_protein_output(i,dba,maxcut)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			aux[[i]] = make_it_protein_output(i,dba,maxcut)
		}

	}

	proteins = re_list(aux)

	return(proteins)

}


make_it_group_output = function(i,dba,dbsel,maxcut,sep="_",n=2){

	aux = list()
	#aux$pdb = toupper(substr(dba[[1]][[i]]$pdbname,1,4))
	comb = unlist(strsplit(dba[[1]][[i]]$pdbname,sep))#;print(comb)
	aux$pdb = comb[1]
	if (i%%n==1) aux$pdb_type = 2
	else aux$pdb_type = 1
	aux$chains = dba[[1]][[i]]$chains
	#comb = unlist(strsplit(dba[[1]][[i]]$pdbname,sep))#;print(comb) 
	aux$combination_chains = comb[2:3]
	aux$rotall = dba[[1]][[i]]$exp$rotref
	aux$rotpp = dba[[2]][[i]]$exp$rotref
	aux$rotaa = dba[[3]][[i]]$exp$rotref
	
	aux$connections_group_bsr = list()
	aux$connections_group_bsr_polar_polar = list()
	aux$connections_group_bsr_nonpolar_nonpolar = list()

	aux$volume_bsr = list()
	aux$volume_bsr_polar_polar = list()
	aux$volume_bsr_nonpolar_nonpolar = list()

	aux$spot_score_bsr = list()
	aux$spot_score_bsr_polar_polar = list()
	aux$spot_score_bsr_nonpolar_nonpolar = list()

	aux$local_score_bsr = list()
	aux$local_score_bsr_polar_polar = list()
	aux$local_score_bsr_nonpolar_nonpolar = list()

	aux$global_score_bsr = list()
	aux$global_score_bsr_polar_polar = list()
	aux$global_score_bsr_nonpolar_nonpolar = list()

	aux$dense_bsr_total = list()
	aux$dense_bsr_total_polar_polar = list()
	aux$dense_bsr_total_nonpolar_nonpolar = list()

	aux$dense_bsr_bipar = list()
	aux$dense_bsr_bipar_polar_polar = list()
	aux$dense_bsr_bipar_nonpolar_nonpolar = list()

	aux$center_bsr = list()
	aux$center_bsr_polar_polar = list()
	aux$center_bsr_nonpolar_nonpolar = list()
	#print(i)
	if (!is.null(dba[[1]][[i]]$A)){
		hot = dbsel[[1]]$hot[[i]]#;print(hot)
		superclus = prune_superclus(dba[[1]][[i]]$exp$superclus,maxcut)
		res =  get_connections_group_bsr(superclus)
		aux$connections_group_bsr = res$conn
		aux$volume_bsr = res$vol
		aux$spot_score_bsr = set_spot_score(length(res$conn),hot)#;print("o")
		aux$local_score_bsr = res$local
		aux$global_score_bsr = res$global
		aux$dense_bsr_total = res$dense_t
		aux$dense_bsr_bipar = res$dense_b
		aux$center_bsr = res$center
	
	}
	if (!is.null(dba[[2]][[i]]$A)){
		hot = dbsel[[2]]$hot[[i]]
		superclus = prune_superclus(dba[[2]][[i]]$exp$superclus,maxcut)
		res =  get_connections_group_bsr(superclus)
		aux$connections_group_bsr_polar_polar = res$conn
		aux$volume_bsr_polar_polar = res$vol
		aux$spot_score_bsr_polar_polar = set_spot_score(length(res$conn),hot)#;print("o")
		aux$local_score_bsr_polar_polar = res$local
		aux$global_score_bsr_polar_polar = res$global
		aux$dense_bsr_total_polar_polar = res$dense_t
		aux$dense_bsr_bipar_polar_polar = res$dense_b
		aux$center_bsr_polar_polar = res$center
	
	}
	if (!is.null(dba[[3]][[i]]$A)){
		hot = dbsel[[3]]$hot[[i]]
		superclus = prune_superclus(dba[[3]][[i]]$exp$superclus,maxcut)
		res =  get_connections_group_bsr(superclus)
		aux$connections_group_bsr_nonpolar_nonpolar = res$conn
		aux$volume_bsr_nonpolar_nonpolar = res$vol#;;print(hot)
		aux$spot_score_bsr_nonpolar_nonpolar = set_spot_score(length(res$conn),hot)#;print("o")
		aux$local_score_bsr_nonpolar_nonpolar = res$local
		aux$global_score_bsr_nonpolar_nonpolar = res$global
		aux$dense_bsr_total_nonpolar_nonpolar = res$dense_t
		aux$dense_bsr_bipar_nonpolar_nonpolar = res$dense_b
		aux$center_bsr_nonpolar_nonpolar = res$center
	}
	
	return(aux)

}

make_hyper_group_output = function(dba,dbsel,maxcut,ncores){

	tam = length(dba[[1]])

	#tam_pdb = tam/3

	group = list()
	aux = list()
	print(paste("Group outputs..."))
	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_group_output(i,dba,dbsel,maxcut)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){#print(i)
		#for (i in 18){
			aux[[i]] = make_it_group_output(i,dba,dbsel,maxcut)
		}

	}

	group = re_list(aux)

	return(group)

}


make_align_output = function(db,sep="_",n=2){

	aux = list()
	if (length(db)){

		tit1 = unlist(strsplit(as.character(db$tabs$pos1[1,1]),sep))
		tit2 = unlist(strsplit(as.character(db$tabs$pos2[1,1]),sep))

		#print(tit1[4]);readline()

		aux$pdb1 = tit1[1]
		aux$pdb2 = tit2[1]
	
		if (tit1[2]=="CHAIN") aux$type = 1
		else aux$type = 2

		aux$combination_chains = tit1[2:3]

		aux$polarity = tit1[4]

		aux$rotself = db$align[[1]]$rotself

		superclus = list()
		#superclus[[1]] = list()
		#superclus[[1]]$a = db$align[[1]]$p1$a
		#superclus[[1]]$geomc = db$align[[1]]$p1$xyz
		#superclus[[1]]$rot = diag(4)
		#superclus[[1]]$score_align = NA
		#superclus[[1]]$score_group = db$align[[1]]$t1[1,6:7]

		#superclus = c(superclus,add_superclus_from_align(db$align))
		superclus = add_superclus_from_align(db$align)#;print("o")

		aux$connections_group_bsr = list()
		aux$connections_group_bsr_polar_polar = list()
		aux$connections_group_bsr_nonpolar_nonpolar = list()

		aux$volume_bsr = list()
		aux$volume_bsr_polar_polar = list()
		aux$volume_bsr_nonpolar_nonpolar = list()

		aux$rot_bsr = list()
		aux$rot_bsr_polar_polar = list()
		aux$rot_bsr_nonpolar_nonpolar = list()

		aux$score_align_bsr = list()
		aux$score_align_bsr_polar_polar = list()
		aux$score_align_bsr_nonpolar_nonpolar = list()

		aux$score_group_bsr = list()
		aux$score_group_bsr_polar_polar = list()
		aux$score_group_bsr_nonpolar_nonpolar = list()

		aux$spot_score_bsr = list()
		aux$spot_score_bsr_polar_polar = list()
		aux$spot_score_bsr_nonpolar_nonpolar = list()

		aux$dense_bsr_total = list()
		aux$dense_bsr_total_polar_polar = list()
		aux$dense_bsr_total_nonpolar_nonpolar = list()

		aux$local_score_bsr = list()
		aux$local_score_bsr_polar_polar = list()
		aux$local_score_bsr_nonpolar_nonpolar = list()

		#print(superclus);readline()

		if (tit1[4]=="ALL"){
			res =  get_connections_group_bsr(super=superclus,original=F)
			aux$connections_group_bsr = res$conn
			aux$volume_bsr = res$vol
			aux$rot_bsr = res$rot
			aux$score_align_bsr = res$score_align
			aux$score_group_bsr = res$score_group
			aux$spot_score_bsr = res$spot
			aux$dense_bsr_total = res$dense
			aux$local_score_bsr = res$qualy

		}

		if (tit1[4]=="PP"){
			res =  get_connections_group_bsr(super=superclus,original=F)
			aux$connections_group_bsr_polar_polar = res$conn
			aux$volume_bsr_polar_polar = res$vol
			aux$rot_bsr_polar_polar = res$rot
			aux$score_align_bsr_polar_polar = res$score_align
			aux$score_group_bsr_polar_polar = res$score_group
			aux$spot_score_bsr_polar_polar = res$spot
			aux$dense_bsr_total_polar_polar = res$dense
			aux$local_score_bsr_polar_polar = res$qualy

		}

		if (tit1[4]=="AA"){
			res =  get_connections_group_bsr(super=superclus,original=F)
			aux$connections_group_bsr_nonpolar_nonpolar = res$conn
			aux$volume_bsr_nonpolar_nonpolar = res$vol
			aux$rot_bsr_nonpolar_nonpolar = res$rot
			aux$score_align_bsr_nonpolar_nonpolar = res$score_align
			aux$score_group_bsr_nonpolar_nonpolar = res$score_group
			aux$spot_score_bsr_nonpolar_nonpolar = res$spot
			aux$dense_bsr_total_nonpolar_nonpolar = res$dense
			aux$local_score_bsr_nonpolar_nonpolar = res$qualy

		}
	}else{
		print(paste("WARNING: dbali is empty"))
	}
	#print(res);readline()
	
	return(aux)


}

#make_align_output = function(db,dbsup,ncores){

#	tam = length(db)

	#tam_pdb = tam/3

#	group = list()
#	aux = list()
#	print(paste("Aligned outputs..."))
#	if (0){
#		print(paste("Doing it parallel with",ncores,"cores"))
#		registerDoMC(ncores)
#		aux = foreach(i=1:tam) %dopar% {
#			make_it_align_output(i,db,dbsup)
#		}
#	}else{
#		print(paste("Doing it sequential..."))
#		for (i in 1:tam){#print(i)
		#for (i in 4){
#			aux[[i]] = make_it_align_output(i,db,dbsup)
#		}

#	}

#	group = re_list(aux)

#	return(group)

#}


#################################################### GENERAL ##################################

# GENERAL: Faz input em massa de dados
#input_general_data = function(filename,matrixname,sufix,typename="Values",row.names=0,col1onull=T,as.matrix=T,bin=F,weight=F,header=T,rnd=2,matrixtype="normal"){

make_mass_input_general_data =function(datapath,filename,typename="Values",row.names=0,col1onull=T,as.matrix=T,bin=F,weight=F,sufix=c(".csv",".mtx"),header=T,rnd=2,matrixtype="normal"){

	#input_general_data = function(filename,matrixname,typename="Values",col1onull=T,as.matrix=T,bin=F,weight=F)
	#make_mass_input = function(dirname,pdbids,prefix,sufix,typename,meso,aacode,verbose=TRUE,bin=FALSE){
	#if (verbose) print("Entrada de dados")
	#print(filename)
	filelist = read.csv(paste(datapath,filename,sufix[1],sep=""),header=F,sep=",")#;print(filelist)
	filelist = paste(as.matrix(filelist))#;print(filelist)
	tam = length(filelist)
	aux = list()
	for (i in 1:tam){
		#print(paste(dirname,prefix,pdbids[1,i],sufix,sep=""))
		#readline()
		#aux[[i]] = input_data(paste(dirname,prefix,pdbids[1,i],sufix,sep=""),as.character(pdbids[1,i]),typename,meso,aacode)
		#aux[[i]] = input_data(paste(dirname,prefix,pdbids[i],sufix,sep=""),pdbids[i],typename,meso,aacode,bin)
		aux[[i]] = input_general_data(filename=paste(datapath,filelist[i],sufix[2],sep=""),matrixname=filelist[i],sufix=sufix[2],typename=typename,row.names=row.names,col1onull=col1onull,as.matrix=as.matrix,bin=bin,weight=weight,header=header,rnd=2,matrixtype=matrixtype)


	}

	return(aux)

}


make_mass_correct_colnames = function(db,s=c("X","K","V")){

	tam = 1:length(db)

	for (i in tam){
		auxc = db[[i]]$colnames
		tamj = 1:length(s)
		for (j in tamj){
			auxc = gsub(s[j],"",auxc)
		}
		db[[i]]$colnames = auxc
	}
	return(db)
}

#GENERAL: filtra db conforme filter
make_mass_general_filter = function(db,rowfilter,colfilter=rowfilter,colall=T){

	tam = length(db)
	aux = list()

	for (i in 1:tam){
		#if (type==1) auxfactor = as.factor(db[[i]]$ap)
		#if (type==2) auxfactor = as.factor(db[[i]]$apg)
		#if (type==3) auxfactor = as.factor(db[[i]]$ape) #definição mais indicada
		#levels(auxfactor) = vfilter
		if (colall){
			aux[[i]] = apply_general_filter(db[[i]],rowfilter[[i]],colfilter=rep(T,length(db[[i]]$colnames)))
		} else {
			aux[[i]] = apply_general_filter(db[[i]],rowfilter[[i]],colfilter[[i]]) 
		}
	}
	return(aux)
}

#GENERAL: constroi vector de filtro com base no IDs faltantes

build_vector_general_filter = function(dbin,dbout,range=1:length(dbin),colind=1){

	aux=list()
	#print(range)
	j=1
	for (i in range){
		#print(range);readline()
		aux1=unique(dbout[[i]]$A$a[,colind])
		aux2=as.numeric(dbin[[i]]$rownames)
		aux[[j]] = !(aux2 %in% aux1)
		j=j+1
	}
	#print(aux);readline()
	return(aux)
}

# GENERAL: Produz em massa as Decomposições EIG e SVD
make_mass_decomps_general = function(db,A,L,Lw,EIG,SVD,centered=1,weight=c()){

	#tam = length(pdbids)
	tam=length(db)
	aux = db

	for (i in 1:tam){

		#ppfaz = laps_and_decomps(ppfaz)
		if (!is.null(aux[[i]]$weight)){
			weight = aux[[i]]$weight
		}
		aux[[i]] = decomps_general(aux[[i]],A,L,Lw,EIG,SVD,centered=centered,weight=weight)

	}

	return(aux)

}

#GENERAl: prepara experimentos PAM
make_mass_init_pam_exp_general = function (db,expid,nvect){

	tam = length(db)
	#print(tam)
	for (i in 1:tam){
		db[[i]]$exp[[expid]]=list()
		db[[i]]$exp[[expid]]$vec=nvect
		db[[i]]$exp[[expid]]$pam=list()
		db[[i]]$exp[[expid]]$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
	}
	return(db)
}


# GENERAL: Faz PAM sobre db para número de clusters k em range
make_mass_pam_general = function(db,range,expid,nvect,type="AS",matrix="us",verbose=F){

	db = make_mass_init_pam_exp_general(db,expid,nvect)
	#print(1)
	#expid = paminit
	for (k in range){
		if (verbose) print(k)
		#db = mass_experiment_02(db,expid,k)
		db = mass_pam_general(db,k,expid,nvect,type,matrix)		
		#pamid = pamid + 1
	}
	return(db)

}

# GENERAL: PAM com k partições
mass_pam_general = function(db,k,expid,nvect,type,matrix){

	tam = length(db)
	#no = n-k+1
	#mlogic = diag(rep(1,n))>0
	#print("ok")
	for (i in 1:tam){
		#expam = list()
		#print(i);print(expid);print(k);readline()
		#db[[i]]$exp[[expid]]=list()
		#db[[i]]$exp[[expid]]$pam=list()
		#db[[i]]$exp[[expid]]$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
		#print("ok")
		if(substr(type,1,1)=="A"){
			if (substr(type,2,2)=="S"){
				#n = dim(db[[i]]$A$svd$v)[1]
				#no = n-k+1
				#print(i)
				if (matrix == "us"){
					aux = pam(db[[i]]$A$svd$us[,nvect],k)
				}
				if (matrix == "u"){
					aux = pam(db[[i]]$A$svd$u[,nvect],k)
				}
				if (matrix == "a"){
					aux = pam(db[[i]]$A$a[,nvect],k)
				}
				#print("ok")
			}
		}

		#exp1 = db[[i]]$exp$vlist[[1]]
		#exp2 = make_vector_indicator_list(aux$clustering)	
		#exp = make_logical_vector_product(exp1,exp2)
		#db[[i]]$exp$vlist[[expid]]=exp
		#print(aux);readline()
		db[[i]]$exp[[expid]]$pam[[k]] = aux
		#print(db[[i]]$exp[[expid]]$pam[[k]]);readline()

		#exp = list()
		#exp = make_vector_indicator_list(db[[i]]$ide)
		#db[[i]]$exp[[expid]]=exp
	}
	return(db)	
}

#GENERAL ERNEST: a partir dados Ernesto
make_mass_stat_res_by_decomp = function(db,k){

	tamdb = 1:length(db)

	for (i in tamdb){
		db[[i]]$exp1=list()
		svd = db[[i]]$A$svd
		db[[i]]$exp1 = stat_res_by_decomp(svd,k)

	}
	return(db)

}

make_mass_ernest_plot = function(db,nc=1:3,changex=T,main="",xlab="",ylab=""){

	taml = length(db)

	tamc = length(nc)
	
	par(mfrow=c(taml,tamc))

	

	for (i in 1:taml){

		for (j in nc){
			auxm = db[[i]]$exp1[[j]]$mr
			if (changex){
				aux1 = db[[i]]$A$a[,1]
				auxm2 = auxm[,2:dim(auxm)[2]]
				auxm = cbind(aux1,auxm2)
			}
			if (j==1){
				main=paste("Orig Norm Spectra sample",i)
			}else{
				main=paste("Rebuild Matrix with autovector: ",1,":",j,sep="")
			}
			plot_ernest(auxm,main=main,xlab=xlab,ylab=ylab)
		}
	}

}

make_mass_ernest_res_plot = function(db,type=1,nc=1:3,ylim=0,changex=T,main="",xlab="",ylab=""){

	taml = length(db)

	tamc = length(nc)
	#print(nc);readline()
	
	if (type==1) par(mfrow=c(taml,tamc))
	if (type==2) par(mfrow=c(taml,tamc))

	for (i in 1:taml){

		for (j in nc){

			if (type==1){
				auxm = db[[i]]$exp1[[j]]$dr
				if (is.null(dim(auxm))){
					auxm = db[[i]]$exp1[[j]]$mr
				}
				if (changex){
					aux1 = db[[i]]$A$a[,1]
					auxm2 = auxm[,2:dim(auxm)[2]]
					auxm = cbind(aux1,auxm2)
				}
				if (j==1){
					main=paste("Orig Norm Spectra sample",i)
				}else{
					main=paste("Residue Plot: ",1,":",j," - ","Original",sep="")
				}
				if (ylim[j]==0){
					plot_ernest(auxm,main=main,xlab=xlab,ylab=ylab)
				}else{
					plot_ernest(auxm,ylim=c(-ylim[j],ylim[j]),main=main,xlab=xlab,ylab=ylab)
				}
			}
			if (type==2){
				#a = db[[1]]$exp1[[10]]$mr
				#b = db[[1]]$exp1[[9]]$mr
				#x= a - b
				#plot_ernest(a)
				#plot_ernest(b)
				#x = cbind(a[,1],x[,2:dim(x)[2]])
				#plot_ernest(x,ylim=c(-0.001,0.001))
				if (j==1){
					auxm = db[[i]]$exp1[[j]]$mr
					main=paste("Orig Norm Spectra sample",i)
					if (changex){
						aux1 = db[[i]]$A$a[,1]
						auxm2 = auxm[,2:dim(auxm)[2]]
						auxm = cbind(aux1,auxm2)
					}
					plot_ernest(auxm,main=main,xlab=xlab,ylab=ylab)
				}else{
					auxm = db[[i]]$exp1[[j]]$cp
					if (changex){
						aux1 = db[[i]]$A$a[,1]
						auxm2 = auxm[,2:dim(auxm)[2]]
						auxm = cbind(aux1,auxm2)
					}
					main=paste("Component: ",j," (",1,":",j,"-",1,":",j-1,")",sep="")
					plot_ernest(auxm,ylim=c(-ylim[j],ylim[j]),main=main,xlab=xlab,ylab=ylab)
				}
			}
		}
	}

}


# GENERAL: View silhouette plot para métrica avg_width em função do número de clusters
# Efeito colateral: retorna matriz com os avg_width
view_mass_pam_silhouette_avg_width_general = function(db,expid,xlim,ylim,vleg,tleg,main="",xlab="",ylab="",type=1,ret=F,mtext="", cex.leg=0.8,cex.pch=0.6,pch=21,mtext.cex=0.8,line=0.3,zlim=-0.05){

	tamdb = 1:length(db)
	#print(tamdb);print(db[[1]]$pdbname);print(db[[2]]$pdbname)
	tampam = length(db[[1]]$exp[[expid]]$pam)
	#print(tampam)
	aux=c()
	auxneg=c()
	auxall = c()
	legtext = c()
	main = paste(main,strcat(paste(db[[1]]$exp[[expid]]$vec)))
	plot(1:tampam,seq(-1,1,length.out=tampam),type="n",ylim=ylim,xlim=xlim,main=main,xlab=xlab,ylab=ylab)#;print("ok");readline()
	#if (mtext=="") mtext("",cex=mtext.cex,line=line)
	#else
	#print(mtext)
	mtext(mtext,cex=mtext.cex,line=line)
	for (i in tamdb){
		#par(mfrow=c(1,1))
		#tampam = length(db[[i]]$exp$pam)
		aux = c()
		auxneg = c()
		for (j in 2:tampam){
			temp = db[[i]]$exp[[expid]]$pam[[j]]$silinfo$widths[,3]
			auxsum = sum(as.numeric(temp<zlim))
			if (auxsum>0) auxneg = c(auxneg,"BLACK")
			else {
				temp1 = as.logical(temp>=zlim)
				temp2 = as.logical(temp<(-zlim))
				temp12 = temp1*temp2
				auxsum = sum(as.numeric(temp12))
				if (auxsum>0) auxneg = c(auxneg,"GRAY")
				else auxneg = c(auxneg,"WHITE")
			}
			if (type==1) aux = c(aux,db[[i]]$exp[[expid]]$pam[[j]]$silinfo$avg.width)
			if (type==2) {
				aux = c(aux,mean(db[[i]]$exp[[expid]]$pam[[j]]$silinfo$clus.avg.widths))
				#print(mean(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths));print(sd(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
			}
			if (type==3) {
				aux = c(aux,sd(db[[i]]$exp[[expid]]$pam[[j]]$silinfo$clus.avg.widths))
				#print(mean(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths));print(sd(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
			}
			
			#plot(i,dbpam[[i]]$silinfo$avg.width,type="b")
		}
		
		lines(2:tampam,aux,lty=i)
		#print(auxneg)
		#readline()
		points(2:tampam,aux,pch=pch,bg=auxneg,cex=cex.pch)
		#legtext = c(legtext,db[[i]]$pdbname)
		legtext = c(legtext,db[[i]]$matrixname)
		auxall = rbind(auxall,aux)
	}
	legend(vleg[1],vleg[2],legtext,lty=tamdb,cex=cex.leg)
	#rownames(auxall)=pdbids[tamdb]
	colnames(auxall)=2:tampam
	if (ret) return(auxall)
}


#	k=5
#	indv=dbv[[1]]$exp[[7]]$pam[[5]]$clustering==k
#	m=dbv[[1]]$A$a[indv,]
#	#view_matrix(m,pcex=1,minv=min(m),maxv=max(m),sym=1,pal=2,main.text="")
#	view_matrix_general(m,paste("grupo ",k,sep=""))


#GENERAL: visualiza matrix
view_mass_matrix_general = function(db,expid,g,krange,mfrow,pcex,x.lab="",y.lab="",prefix="Grupo ",pal=1,rd=2,mt.cex=0.6,mt.line=0.3,mt.prefix="Avg Width: ",saveg=c(),u=c(),sufix=".csv",outname="v1",savetype=2,s=c(),colbyc=T){
	
	par(mfrow=mfrow)

	clus = db$exp[[expid]]$pam[[g]]$clustering
	if (pal==1){
		tcol=unique(clus)
	}
	if (pal==2){
		ncol=rainbow(g)
		tcol = ncol
	}
	m=db$A$a
	vec=db$exp[[expid]]$vec
	matrixname=db$matrixname
	#xdim = dim(m)[2]
	#ydim = dim(m)[1]
	#plot(1,1,type="n",xlim=c(1,xdim),ylim=c(1,ydim),xlab=x.lab,ylab=y.lab)
	#title(main.text,col.main=tcol)
	clusavw=round(db$exp[[expid]]$pam[[g]]$silinfo$clus.avg.widths,rd)
	widths=db$exp[[expid]]$pam[[g]]$silinfo$widths
	rn = as.numeric(rownames(widths))
	o=order(rn)
	widths=as.matrix(widths[o,])
	#print(widths);readline()
	j=1
	for (i in krange){
		indv=clus==i
		#m=db$A$a[indv,]
		aux=m[indv,]
		auxw=widths[indv,]
		auxrownames = as.numeric(db$rownames[indv])
		#print(auxrownames);readline()
		#print(m)f
		if(!is.matrix(aux)){
			#print(1)
			aux = t(as.matrix(aux))
			#auxw = t(as.matrix(auxw))
			#print(m)
		}
		#print(is.matrix(auxw))
		#print(auxw);readline()
		if(is.matrix(auxw)){
			aux=cbind(aux,round(auxw[,3],rd))
		}else{
			#print(aux);print(auxw)
			aux=cbind(aux,round(auxw[3],rd))
			rownames(aux)=which(indv==T)
			#print(which(indv==T));readline()
			#print(indv);print(o)
			#print(o[indv]);readline()
			#print(aux);readline()
		}
		#print(aux);readline()
		#aux=cbind(aux,as.numeric(rownames(auxw)))
		#print(as.numeric(rownames(auxw)));readline()
		xdim = dim(aux)[2]+1
		ydim = dim(aux)[1]
		plot(1,1,type="n",xlim=c(1,xdim),ylim=c(1,ydim),xlab=x.lab,ylab=y.lab)
		title(paste(matrixname," - ",prefix,i,sep=""),col.main=tcol[j])
		if (!is.null(mt.prefix)) {
			mtext(paste(mt.prefix,clusavw[j],"  n: ",sum(indv)," u: ",strcat(paste(vec))),cex=mt.cex,line=mt.line)
		}
		#view_matrix_general(aux,paste(prefix,i,sep=""),tcol[j])
		#print(aux);readline()
		view_matrix_general(aux,auxrownames,pcex=pcex,colbyc=colbyc)
		j=j+1
		if (i %in% saveg){
			vec=db$exp[[expid]]$vec
			if (savetype==1){
				#vec=db$exp[[expid]]$vec
				aux=as.data.frame(aux)
				#print(aux);readline()
				aux=cbind(aux,rep(expid,ydim))
				aux=cbind(aux,rep(strcat(paste(vec)),ydim))
				aux=cbind(aux,rep(g,ydim))
				aux=cbind(aux,rep(i,ydim))
				#TO DO: nao esta appendando...
				#write.csv(aux,file=paste(db$matrixname,"-",outfile,sufix,sep=""),append=TRUE)
				write.table(aux,file=paste(db$matrixname,"-",outname,sufix,sep=""),append=TRUE,col.names=F,sep=",")
				#write.csv(aux,file=paste(db$matrixname,"-e",expid,"-u",u,"-k",g,"-g",i,sufix,sep=""),col.names=FALSE)
			}
			if (savetype==2){
				#print(aux);readline()
				#auxs=data.frame()
				#auxs=as.numeric(rownames(aux))
				#auxs=as.numeric(rownames(aux))
				auxs=auxrownames
				auxs=as.data.frame(auxs)
				#print(auxs);readline()
				auxs=cbind(auxs,rep(s,ydim))
				#print(auxs);readline()
				#aux=cbind(aux,rep(u,ydim))
				#aux=cbind(aux,rep(g,ydim))
				auxs=cbind(auxs,rep(i,ydim))
				auxs=cbind(auxs,aux[,xdim-1])
				auxs=cbind(auxs,rep(strcat(paste(vec)),ydim))
				#print(auxs);readline()
				filename=paste(db$matrixname,"-",outname,sufix,sep="")
				colnames(auxs)=c("id","group","subgroup","avg_width","svd-u")
				#TO DO: nao esta appendando...
				#write.csv(aux,file=paste(db$matrixname,"-",outfile,sufix,sep=""),append=TRUE)
				if(file.exists(filename)){
					print("One file will be rewrite... Are you sure?");readline()
					write.table(auxs,file=filename,append=T,col.names=F,sep=",",row.names=F)
				}else{
					write.table(auxs,file=filename,append=F,col.names=T,sep=",",row.names=F)
				}
				#write.csv(aux,file=paste(db$matrixname,"-e",expid,"-u",u,"-k",g,"-g",i,sufix,sep=""),col.names=FALSE)
			}
		}
	} 

}

eliminate_mass_zero_column = function(db,bkp=T){

	tam = 1:length(db)

	for (i in tam){
		db[[i]] = eliminate_zero_column(db[[i]],bkp=bkp)

	}
	return(db)


}

#GENERAL: anota dados dos grupos
resume_stat_mass_group_general = function(db,id=c(),iexp=c(),ipam=c(),isubg=c()){

	aux = db
	#return(aux)
	#tam = length(db)
	if (is.null(id)){
		id = 1:length(aux)
	}


	for (i in id){
		gid=1
		if (is.null(iexp)){
			exp = 1:length(aux[[i]]$exp)
		}else{
			exp=iexp
		}
		auxa=db[[i]]$A$a
		ids = db[[i]]$A
		rownames = db[[i]]$rownames
		for (j in exp){
			if (is.null(ipam)){
				pam = 2:length(aux[[i]]$exp[[j]]$pam)
			}else{
				pam=ipam
			}
			groups=list()
			groups[[1]]=aux[[i]]$exp[[j]]$pam[[1]]

			for (k in pam){
				#print(i);print(j);print(k);print(is.list(aux[[i]]$exp))
				#return(aux)
				widths=aux[[i]]$exp[[j]]$pam[[k]]$silinfo$widths
				#print(widths);readline()
				rn = as.numeric(rownames(widths))
				o=order(rn)
				widths=as.matrix(widths[o,])
				if (is.null(isubg)){
					subg = 1:k
				}else{
					subg=isubg
				}
				clus = aux[[i]]$exp[[j]]$pam[[k]]$clustering
				#print(clus)
				subs=list()
				for (g in subg){
					#print(g);readline()	
					indv=clus==g
					auxm=auxa[indv,]
					auxw=widths[indv,]
					auxrn=rownames[indv]
					if(!is.matrix(auxm)){
						auxm = t(as.matrix(auxm))
					}
					if(!is.matrix(auxw)){
						auxw = t(as.matrix(auxw))
						rownames(auxw)=rownames[which(indv==T)]
					}else{
						rownames(auxw)=auxrn
					}
					#if(is.matrix(auxw)){
					#	aux=cbind(aux,round(auxw[,3],rd))
					#}else{
					#	aux=cbind(aux,round(auxw[3],rd))
					#	rownames(aux)=which(indv==T)
					#}
					subs[[g]] = list()
					subs[[g]]$finger = finger_print_matrix(auxm)
					subs[[g]]$widths = auxw
					subs[[g]]$gid = gid
					subs[[g]]$mark = F
					gid=gid+1
				}
				groups[[k]] = list()
				groups[[k]]$subs = subs 
					
			}
			aux[[i]]$exp[[j]]$groups = groups
		
		}
		
	}
	return(aux)

}

#GENERAL: seleciona melhores grupos...
select_mass_groups_general = function(db,expid,pamid,srange,one){

	aux=db
	if (is.null(aux[[1]]$savegroups)){
		#print(1)
		aux[[1]]$savegroups=c()
		last=1
	}else{
		auxrn = as.numeric(rownames(aux[[1]]$savegroups))
		tam = length(auxrn)
		last = auxrn(tam)
		print("Saving groups... Are you sure?");readline()
	
	}
	#auxg=aux[[1]]$savegroups
	auxg=c()
	for (i in srange){
		#db$exp[[expid]]$groups[[pamid]]$subs[[i]]$mark=T
		auxgid=aux[[1]]$exp[[expid]]$groups[[pamid]]$subs[[i]]$gid
		auxrow=c(expid,pamid,i,auxgid)#;print(auxgid)
		auxg=rbind(auxg,auxrow)
	}
	#print(auxg)
	#aux[[1]]$savegroups=auxg
	colnames(auxg)=c("expid","pamid","subg","gid")
	if (one){
		rownames(auxg)=rep(last,length(srange))
	}else{#ARRUMAR !!!
		print(last);print(last+length(srange))
		rownames(auxg)=last:(last+length(srange))
	}
	
	aux[[1]]$savegroups=rbind(aux[[1]]$savegroups,auxg)
	return(aux)
}


###################################################### NILMA ##################################
# Agrupa informações gerais sobre os PDBids usados
make_general_input = function(db,pdbids){

	tam = 1:length(db)

	aux = list()
	aux$enzids=c()
	aux$inbids=c()

	for (i in tam){
		aux$enzids = c(aux$enzids,db[[i]]$enzcode)
		aux$inbids = c(aux$inbids,db[[i]]$inbcode)
	}
	aux$enztype = substr(aux$enzids,3,3)
	aux$pdbids = pdbids
	#aux$stat = list()
	#names(aux$enzids) = tam
	#names(aux$inbids) = tam
	#names(aux$pdbids) = tam
	#names(aux$enztype) = tam
	#print(aux);readline()
	return(aux)
}

# Faz entrada em massa da matriz de p-values ks.test()
make_mass_stat_input=function(db,pdbids,verbose=FALSE){

	#dks=list()
	#dks[[1]]=make_stat_input(dba,diagz=1,verbose=FALSE)
	#dks=add_bin_stat_input(dks,pdbids,diagz=1)
	#dks[[3]]=make_stat_input(dba,diagz=0,verbose=FALSE)
	#dks=add_bin_stat_input(dks,pdbids,diagz=0)
	#dks=add_center_stat_input(dks)
	aux=list()
	aux[[1]]=make_stat_input(db,diagz=1,verbose=verbose)
	aux=add_bin_stat_input(aux,pdbids,diagz=1)
	aux[[3]]=make_stat_input(db,diagz=0,verbose=verbose)
	aux=add_bin_stat_input(aux,pdbids,diagz=0)
	aux=add_center_stat_input(aux)
	return(aux)
}

# Insere em massa informações das coordenadas atômicas PDBs
add_mass_pdb_xyz = function(db,pdbxyz,ncores=0,verbose=TRUE,typename=1){

	tam = length(db)

	#aux=db
	aux = list()
	print(paste("Input PDB xyz..."))
	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		
		aux = foreach(i=1:tam) %dopar% {
			#insert_pdb_xyz(db[[i]],dirpathin,typename=typename)
			insert_pdb_xyz(db[[i]],pdbxyz[[i]],typename=typename)
		}

	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			if (verbose) print(paste("Processing pdb xyz for",toupper(db[[i]]$pdbname),"..."))
			#aux[[i]]$pdb = insert_pdb_xyz(db[[i]],dirpathin,typename=typename)
			#aux[[i]] = insert_pdb_xyz(db[[i]],dirpathin,typename=typename)
			aux[[i]] = insert_pdb_xyz(db[[i]],pdbxyz[[i]],typename=typename) 
		}
	}

	return(aux)
}

add_mass_clus_net_info = function(db,net){

	tam = 1:length(db)	

	for (i in tam){
		#db[[i]]$clus = list()
		db[[i]] = add_clus_net_info(db[[i]],net)
	}
	return(db)
}

# Acrescenta informação dos cluster ids por atom conforme k em kclus
add_mass_atom_loop_clustering = function(db,tamclus){

	tam = 1:length(db)
	aux = db
	
	for (i in tam){
		aux[[i]]$kclus = list()
		aux[[i]]$kclus[[1]] = "There is no data for cluster with k = 1"
		for (j in 2:tamclus){
			aux[[i]]$kclus[[j]] = map_atomic_loop_clustering(db[[i]],j)
		}
	}
	return(aux)
	
}


# OBSOLETED: Acrescenta informação dos cluster ids por atom conforme k em kclus
# CUIDADO: se rodar mais de uma vez ele appenda inter mais de uma vez
add_mass_atom_loop_clustering_OLD = function(db,kclus){

	tam = 1:length(db)
	aux = db
	
	for (i in tam){
		aux[[i]]$inter=c(aux[[i]]$inter,map_atomic_loop_clustering(db[[i]],kclus[i]))
	}
	return(aux)
	
}

# Faz associação em massa entre atomos do looping com k cluster definidos em inter
add_mass_loopid_per_clustering = function(db){

	tami = 1:length(db)
	aux = db
	for (i in tami){
		tamj = 2:length(db[[i]]$kclus)
		for (j in tamj){
			aux[[i]]$kclus[[j]]$clus=add_loopid_per_clustering(db[[i]]$inter,db[[i]]$kclus[[j]])
		}
	}
	return(aux)
}


# OBSOLETED: Faz associação em massa entre atomos do looping com k cluster definidos em inter
add_mass_loopid_per_clustering_OLD = function(db){

	tam = 1:length(db)
	aux = db
	for (i in tam){
		aux[[i]]$clus=add_loopid_per_clustering(db[[i]]$inter)
	}
	return(aux)
}

# Acrescenta informação do loop inibitório.
add_mass_inibitory_loop = function(db){

	tam = 1:length(db)
	aux = db
	for (i in tam){
		aux[[i]]$inter=find_inibitory_loop(db[[i]]$name,db[[i]]$ei)
	
	}
	return(aux)
}

# Acrescenta informação do loop inibitório em nível atômico
add_mass_atom_loop_list = function(db){

	tam = 1:length(db)
	aux = list()
	for (i in tam){
		aux[[i]]=add_atom_loop_list(db[[i]])
	
	}
	return(aux)
}

#Insere classificação das cadeias em db conforme códigos em chainids
add_mass_chains_classification = function(db,chainids,pdbids,chain_1="E",chain_2="I"){

	tam = dim(chainids)[1]
	aux = as.matrix(chainids)
	
	for (i in 1:tam){
		#print(i)
		auxline = aux[i,]
		auxpdb = auxline[1]
		auxce = unlist(strsplit(auxline[2],""))
		auxci = unlist(strsplit(auxline[3],""))
		id=which(auxpdb==pdbids)
		if (length(id)==1){
			#print(c(auxpdb,auxce,auxci,id))
			#readline()
			if (is.null(db[[id]]$ei)) db[[id]]$ei=rep("-",length(db[[id]]$ide))
			j=1
			while (auxce[j]!="e") j=j+1
			j=j+1
			while (j<=length(auxce)){
				xe = auxce[j]
				#print(xe)
				#readline()
				#dba[[which(as.vector(chainids[1,1])==pdbids)]]$ide
				fil=db[[id]]$ide==toupper(xe)
				db[[id]]$ei[fil]="E"
				j=j+1
			}
			j=1
			while (auxci[j]!="i") j=j+1
			j=j+1
			while (j<=length(auxci)){
				xi = auxci[j]
				#dba[[which(as.vector(chainids[1,1])==pdbids)]]$ide
				fil=db[[id]]$ide==toupper(xi)
				db[[id]]$ei[fil]="I"
				j=j+1
			}
			#print(db[[id]]$ei)
			#print(db[[id]]$ide)
			#readline()
			#aux = levels(as.factor(db[[i]]$ide))
			#if (length(aux)==2){
		
			#}
		} else {
			print(paste("WARNING: chainids do not have correspondence with pdbids: ",auxpdb))
		}
	}
	return(db)
}

add_mass_geometric_center_per_clustering = function(db){

	tam = 1:length(db)
	aux = db
	
	for (i in tam){
		#print(i)
		aux[[i]] = add_geometric_center_per_clustering(aux[[i]])
	}
	return(aux)

}

add_mass_geometric_center_per_super_clustering = function(db){

	tam = 1:length(db)
	aux = db
	
	for (i in tam){
		#print(i)
		aux[[i]] = add_geometric_center_per_super_clustering(aux[[i]])
	}
	return(aux)
}

tomatrix = function(m){

	if (!is.matrix(m)){
		return(t(as.matrix(m)))
	}else{
		return(m)
	}	


}

make_mass_cluster_filter = function(db,cutpar,ncores=0){

	tam = length(db)
	aux = list()
##	auxr = list()

	print(paste("db filtering..."))
	if (ncores){
		print(paste("Doing it parallel with",ncores,"ncores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			make_all_cluster_filter(db[[i]],cutpar)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			aux[[i]] = make_all_cluster_filter(db[[i]],cutpar)
		}
	}
		#auxf = make_cluster_filter(db[[i]],cutpar)
##		auxr = make_cluster_filter(db[[i]],cutpar)
		#print(auxf);readline()
		#aux[[i]] = apply_filter(db[[i]],auxf)	
##		aux[[i]] = apply_filter(db[[i]],auxr$atomfil)
		#aux[[i]]$mv = auxr$mfil	
##		aux[[i]]$exp$preclus = list()
##		auxlen = length(db[[i]]$exp$superclus)
##		if (!is.null(auxr$mfil)){
##			aux[[i]]$exp$preclus$ids = auxr$mfil
##			a = db[[i]]$exp$superclus[[auxlen]]$a[auxr$mfil,auxr$mfil]
##			aux[[i]]$exp$preclus$a = tomatrix(a)
##			geomc = db[[i]]$exp$superclus[[auxlen]]$geomc[auxr$mfil,]
##			aux[[i]]$exp$preclus$geomc = tomatrix(geomc)
##		}
##	}
	
	return(aux)
}



make_mass_superclus_tables = function(db,dball,lowcut,minclusid,maxclusid,metricpar=c(T,T,T),verbose=T,ncores=0){

	tam = length(db)
	aux = list()


	if (verbose) print(paste("Generating info matrix..."))

	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			make_all_metric_super_clustering(db=db[[i]],id=i,lowcut=lowcut,minclusid=minclusid,maxclusid=maxclusid,metricpar=metricpar)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
		#for (i in 2){
			aux[[i]] = make_all_metric_super_clustering(db=db[[i]],id=i,lowcut=lowcut,minclusid=minclusid,maxclusid=maxclusid,metricpar=metricpar)
		}
	}
	#print(aux);readline()
	aux = do.call(rbind.data.frame, aux)
	colnames(aux) = c("gacc","lacc","verts","edges","idinput")#;print(aux);readline()
	dball$info = aux
	dball = new_make_group_table(db=db,dball=dball,minclusid=minclusid)
	#print(dball);readline()

	return(dball)

}

add_mass_topological_features = function(dball,db,k=2,sep="_"){


	#tam = 1:dim(dball[[1]]$A$a)[1]
	tam = 1:dim(dball[[1]]$A$info)[1]
	#tam = 15
	auxr1 = c()
	auxr2 = c()
	auxr3 = c()
	auxr4 = c()
	auxr5 = c()
	auxr6 = c()
	auxr7 = c()
	auxrn = colnames(dball[[1]]$A$a)
	auxlab = rownames(dball[[1]]$A$info)
	for (j in tam){
		#auxlab = as.character(dball[[1]]$A$info[j,1])
		#print(auxlab[j])
		i = as.numeric(unlist(strsplit(auxlab[j],sep))[1])
		nc = dball[[1]]$A$info[j,3]
		#print(i);print(nc);readline()
		auxr1 = c(auxr1,count_components(dba2[[i]]$exp$superclus[[nc]]$ga))#;print("ok")
		auxr2 = c(auxr2,as.numeric(automorphisms(dba2[[i]]$exp$superclus[[nc]]$ga)$group_size))#;print("ok")
		auxr3 = c(auxr3, centr_betw(dba2[[i]]$exp$superclus[[nc]]$ga,directed=FALSE)$centralization)#;print("ok")
		auxr4 = c(auxr4, centr_clo(dba2[[i]]$exp$superclus[[nc]]$ga,mode="all")$centralization)#;print("ok")
		auxr5 = c(auxr5, diameter(dba2[[i]]$exp$superclus[[nc]]$ga,unconnected=TRUE))#;print("ok")
		auxr6 = c(auxr6, radius(dba2[[i]]$exp$superclus[[nc]]$ga,mode="all"))#;print("ok")
		auxr7 = c(auxr7, transitivity(dba2[[i]]$exp$superclus[[nc]]$ga))#;print("ok")
		
	}#;print("0k")
	auxr3[is.nan(auxr3)] = 0;auxr3 = round(auxr3,k)
	auxr4[is.nan(auxr4)] = 0;auxr4 = round(auxr4,k)
	auxr7[is.nan(auxr7)] = 0;auxr7 = round(auxr7,k)
	auxrn1 = c(auxrn,"nv","ne","comp","auto","betcen","clocen","dia","rad","trans")
	dball[[1]]$A$a = cbind(dball[[1]]$A$a,dball[[1]]$A$info[,3:4],auxr1,auxr2,auxr3,auxr4,auxr5,auxr6,auxr7)
	dball[[1]]$colnames = auxrn1
	auxdiff = length(auxrn1)-length(auxrn)
	dball[[1]]$weight = c(dball[[1]]$weight,rep(1,auxdiff))
	colnames(dball[[1]]$A$a) = auxrn1

	return(dball)
}


get_rot_by_names = function(align,graphname){

	tam = length(align)

	for (i in 1:tam){
		if (align[[i]]$names[2]==graphname){
			return(align[[i]]$rot)
		}
	}
}

plot_mass_new_model_clus = function(db,dball,pdbid,clusid,ids,mfrow,main="",xlab="",ylab="",xlim=NULL,ylim=NULL,sep="_",dxy=c(4,2),cols=NULL,type=c(1,1)){

	graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$score[,1])[1]
	auxlim = get_xylim_from_model(db,dball,graphnames[1],dxy=dxy,sep=sep)
	xlim = auxlim$xlim
	ylim = auxlim$ylim

	par(mfrow=mfrow)
	graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$score[,1])[ids]
	tam = length(graphnames)#;print(tam);readline()
	for (i in 1:tam){
		#auxrot = get_rot_by_names(dball$res[[pdbid]]$alignment[[clusid]]$align,graphnames[i])
		#auxzz = dball$res[[pdbid]]$alignment[[clusid]]$aligndata[[idgrp]]$supxyz
		pars = unlist(strsplit(graphnames[i],sep))
		j = as.numeric(pars[5])
		z = as.numeric(pars[6])
		#print(pars);readline()
		if (type[1]){
			idgrp = which(graphnames[i]==dball$group$fakename)
			rot = dball$res[[pdbid]]$alignment[[clusid]]$aligndata[[idgrp]]$rot
			plot_new_model_clus(db=db[[j]],clus=z,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,rot=rot)
		}
		if (type[2]){
			visualize_align_super_cluster(db=dba2,dball=dball,pdbid=pdbid,clusid=clusid,ids=ids[i],cols=cols,xlim=xlim,ylim=ylim)
			#visualize_align_super_cluster(db=dba2,dball=dball[[i]],ids=c(1:12),cols=cols,mfrow=c(3,4),dxy=dxy,xlim=xlim,ylim=ylim)
#visualize_align_super_cluster(db=dba2,dball=dball[[i]],pdbid=1,clusid=0,ids=c(1:12),cols=cols,mfrow=c(3,4),dxy=dxy,xlim=xlim,ylim=ylim)
		}
	}

}

get_mass_list_graph = function(dball,db,j,ga=T,gd=F,bin=F,sep="_"){

	auxinfo = dball$A$info
	#print(auxinfo)
	auxlab = rownames(auxinfo)
	#print(auxlab)
	tam = 1:dim(auxinfo)[1]
	#auxg = list()
	auxga = list()
	#auxgd = list()
	#print(tam)

	for (k in tam){
		auxnc = auxinfo[k,3]
		i = as.numeric(unlist(strsplit(auxlab[k],sep))[1])
		#print(auxlab[k]);print(i);print(auxnc);print(j);print(k);readline()
		if (ga){
			if (j<auxnc){
				auxga[[k]] = db[[i]]$exp$superclus[[auxnc-j]]$ga
				
			}else{
				auxga[[k]] = db[[i]]$exp$superclus[[1]]$ga
			}
			if (bin){
				E(auxga[[k]])$weight = rep(1,ecount(auxga[[k]]))

			}
			#print(db[[i]]$exp$superclus[[auxnc-j]]$a);readline()
		}
		#if (gd){
		#	auxgd[[i]] = db[[i]]$exp$superclus[[auxnc]]$gd
		#}
		#print(auxg[[i]]);readline()
	}
	#dball$ga = auxga
	#dball$gd = auxgd
	return(auxga)

}

add_mass_kernel_graph = function(dball,db,only_last=0,randome=F,randomg=F,shortest=F,shortestkv=F,k=0,pare=.01,parg=0.5,diag=F,bin=T){

#dball[[1]]$A$a = cbind(dball[[1]]$A$a,diag(round(CalculateExponentialRandomWalkKernel(dball[[1]]$ga,1.2),2)),diag(CalculateShortestPathKernel(dball[[1]]$ga)))

	#dball[[1]] = add_mass_list_graph(dball[[1]],dba2)
	for (j in 0:only_last){
		auxga = get_mass_list_graph(dball=dball[[1]],db=db,j=j,bin=bin)#;print("ok")
		print(length(auxga))
		auxr = c()
		if (randome){
			auxrwe = CalculateExponentialRandomWalkKernel(auxga,pare)
			if (diag){
				auxrk = round(diag(auxrwe),k)
			}else{
				auxrk = round(auxrwe[,1],k)
			}
			auxr = cbind(auxr,RwK=auxrk)
			dball[[1]]$A$a = cbind(dball[[1]]$A$a,auxr)
		}
		if(shortest){
			
			auxsp = CalculateShortestPathKernel(auxga)
			#auxsp = shortest_path_kv_kernel(auxga)
			if (diag){
				auxrk = round(diag(auxsp),k)
			}else{
				auxrk = round(auxsp[,1],k)
			}
			auxr = cbind(auxr,PK=auxrk)
			#dball[[1]]$A$a = cbind(dball[[1]]$A$a,auxr)
		}
		if(shortestkv){
			
			#auxsp = CalculateShortestPathKernel(auxga)
			auxsp = shortest_path_kv_kernel(auxga)
			if (diag){
				auxrk = round(diag(auxsp),k)
			}else{
				auxrk = round(auxsp[,1],k)
			}
			auxr = cbind(auxr,PK=auxrk)
			#dball[[1]]$A$a = cbind(dball[[1]]$A$a,auxr)
		}
		if (randomg){
			auxrwg = CalculateGeometricRandomWalkKernel(auxga,parg)

			if (diag){
				auxrk = round(diag(auxrwg),k)
			}else{
				auxrk = round(auxrwg[,1],k)
			}
			auxr = cbind(auxr,RwK=auxrk)
			#dball[[1]]$A$a = cbind(dball[[1]]$A$a,auxr)
		}

	}
	#if (gd){
	#	auxrk = round(diag(CalculateExponentialRandomWalkKernel(dball$gd,gdpar)),k)
	#	auxr = cbind(auxr,auxrk)
	#	auxrk = round(diag(CalculateShortestPathKernel(dball$gd)),k)
	#	auxr = cbind(auxr,auxrk)
	#	dball$A$a = cbind(dball$A$a,auxr)
	#}
	
	#return(dball)
	#return(dball[[1]]$A$a)
	if (randome) return(round(auxrwe,k))
	if (shortest) return(round(auxsp,k))
	if (shortestkv) return(round(auxsp,k))
	if (randomg) return(round(auxrwg,k))
}

make_mass_super_matrix_old2 = function(db,geopar=NULL,areapar=NULL,metricpar=NULL,betpar=NULL,clopar=NULL,shopar=NULL,
lowcut=0.66,join=F,only_last=T,minpar=4,pre=F,sep="_",list_g=T){

	tam = 1:length(db)
	#tam = 5
	#auxm = db
	auxgeor = c()
	auxarear = c()
	auxmetr = c()
	auxbetr = c()
	auxclor = c()
	auxshor = c()
	auxrn = c()
	#auxncr = c()
	#auxi = make_super_matrix_geocut(db[[i]],geocutpar)
	#auxr = rbind(auxr,auxi)
	#auxgeo = NULL
	#auxarea = NULL
	for (i in tam){
		#print(i)
		if (!is.null(geopar)){
			auxgeo = make_super_matrix(db=db[[i]],type="geo",par=c(geopar,lowcut),id=i,minpar=minpar,pre=pre)
			#print(auxgeo);readline()
			if (only_last){
				auxdim = dim(auxgeo)[1]
				auxgeo = auxgeo[auxdim,]
				#auxgeor = rbind(auxgeor,auxgeo)
			} else if (join){
				#auxdim = dim(auxgeo)[1]
				auxgeo=apply(auxgeo,2,sum)
				#auxgeo=c(auxgeo,nc=auxgeodim)
			}
			auxgeor = rbind(auxgeor,auxgeo)
			
		}
		if (!is.null(areapar)){
			auxarea = make_super_matrix(db=db[[i]],type="area",par=c(areapar,lowcut),id=i,minpar=minpar,pre=pre)
			#print(auxarea);readline()
			if (only_last){
				auxdim = dim(auxarea)[1]
				auxarea = auxarea[auxdim,]
			}else if (join){
				#auxdim = dim(auxarea)[1]
				auxarea=apply(auxarea,2,sum)
				#auxarea=c(auxarea,nc=auxareadim)
			}
			#print(auxarea);readline()
			auxarear = rbind(auxarear,auxarea)
		}
		if (!is.null(metricpar)){
			auxmet = make_super_matrix_metric(db=db[[i]],metricpar=c(metricpar,lowcut),id=i,minpar=minpar,pre=pre)
			if (only_last){
				auxdim = dim(auxmet)[1]
				#print(auxmet);readline()
				auxmet = auxmet[auxdim,]
				if (metricpar[3]){
					auxmet=c(auxmet,nc=auxdim)
				}
			}else if (join){
				auxdim = dim(auxmet)[1]
				#auxmet = c(acc=apply(auxmet,2,sum))
				auxmet = c(acc=apply(auxmet,2,min))
				if (metricpar[3]){
					auxmet=c(auxmet,nc=auxdim)
				}
			}
			auxmetr = rbind(auxmetr,auxmet)
		}
		if (!is.null(betpar)){
			auxbet = make_super_matrix(db=db[[i]],type="bet",par=c(betpar,lowcut),id=i,minpar=minpar,pre=pre)
			if (only_last){
				auxdim = dim(auxbet)[1]
				auxbet = auxbet[auxdim,]
			}else if (join){
				#auxdim = dim(auxbet)[1]
				auxbet=apply(auxbet,2,sum)
				#auxbet=c(auxbet,nc=auxbetdim)
			}
			auxbetr = rbind(auxbetr,auxbet)
		}
		if (!is.null(clopar)){
			auxclo = make_super_matrix(db=db[[i]],type="clo",par=c(clopar,lowcut),id=i,minpar=minpar,pre=pre)
			if (only_last){
				auxdim = dim(auxclo)[1]
				auxclo = auxclo[auxdim,]
			}else if (join){
				#auxdim = dim(auxclo)[1]
				auxclo=apply(auxclo,2,sum)
				#auxclo=c(auxclo,nc=auxclodim)
			}
			auxclor = rbind(auxclor,auxclo)
		}
		if (!is.null(shopar)){
			auxsho = make_super_matrix(db=db[[i]],type="sho",par=c(shopar,lowcut),id=i,minpar=minpar,pre=pre)
			if (only_last){
				auxdim = dim(auxsho)[1]
				auxsho = auxsho[auxdim,]
			}else if (join){
				#auxdim = dim(auxsho)[1]
				auxsho=apply(auxsho,2,sum)
				#auxsho=c(auxsho,nc=auxshodim)
			}
			#print(auxsho);readline()
			auxshor = rbind(auxshor,auxsho)
		}
		#if (!is.null(ncpar)){
			#auxnc = make_super_matrix(db=db[[i]],type="nc",par=c(ncpar,lowcut),id=i,minpar=minpar,pre=pre)
			#print(auxgeo);readline()
			#if (join){
				#auxshodim = dim(auxsho)[1]
			#	auxsho=apply(auxsho,2,sum)
				#auxsho=c(auxsho,nc=auxshodim)
			#}
			
		#	auxncr = rbind(auxncr,auxnc)
		#}
		if (join){
			#if (!is.null(ncpar)){
			#	auxnc = c(nc=auxdim)
			#	auxncr = rbind(auxncr,auxnc)
			#}
			auxrn = c(auxrn,paste(i,db[[i]]$pdbname,i,sep=sep))
			#auxrn = c(auxrn,paste(i,db[[i]]$pdbname,sep=sep))
		}
		#print(auxi);readline()
		#auxr = rbind(auxr,auxgeo)
		
	}
	auxr1 = cbind(auxgeor,auxarear,auxbetr,auxclor,auxshor)
	auxr2 = auxmetr 
	if (join){
		rownames(auxr1) = auxrn
		rownames(auxr2) = auxrn
	}
	auxr = list()
	auxr$a = auxr1
	auxr$info = auxr2 #verificar erro no nc...
	return(auxr)
}

make_mass_super_matrix_old = function(db,geocutpar=NULL,areacutpar=NULL,metricpar=NULL,betpar=NULL,
lowcut=0.66,join=F,minpar=4,pre=F,sep="_"){

	tam = 1:length(db)
	#auxm = db
	auxgeor = c()
	auxarear = c()
	auxmetr = c()
	auxrn = c()
	#auxi = make_super_matrix_geocut(db[[i]],geocutpar)
	#auxr = rbind(auxr,auxi)
	#auxgeo = NULL
	#auxarea = NULL
	for (i in tam){
		#print(i)
		if (!is.null(geocutpar)){
			auxgeo = make_super_matrix_geocut(db=db[[i]],geocutpar=c(geocutpar,lowcut),id=i,minpar=minpar,pre=pre)
			#print(auxgeo);readline()
			if (join){
				auxgeo=apply(auxgeo,2,sum)
			}
			auxgeor = rbind(auxgeor,auxgeo)
		}
		if (!is.null(areacutpar)){
			auxarea = make_super_matrix_areacut(db=db[[i]],areacutpar=c(areacutpar,lowcut),id=i,minpar=minpar,pre=pre)
			if (join){
				auxarea=apply(auxarea,2,sum)
			}
			auxarear = rbind(auxarear,auxarea)
		}
		if (!is.null(metricpar)){
			auxmet = make_super_matrix_metric(db=db[[i]],metricpar=c(metricpar,lowcut),id=i,minpar=minpar,pre=pre)
			if (join){
				auxmet=apply(auxmet,2,sum)
			}
			auxmetr = rbind(auxmetr,auxmet)
		}
		if (!is.null(betcpar)){
			auxmet = make_super_matrix_bet(db=db[[i]],betpar=c(betpar,lowcut),id=i,minpar=minpar,pre=pre)
			if (join){
				auxmet=apply(auxmet,2,sum)
			}
			auxmetr = rbind(auxmetr,auxmet)
		}
		if (join){
			auxrn = c(auxrn,paste(i,db[[i]]$pdbname,i,sep=sep))
		}
		#print(auxi);readline()
		#auxr = rbind(auxr,auxgeo)
		
	}
	#print(auxr);readline()
	#if (!is.null(auxgeor)){
	#	if (!is.null(auxarear)){
	#		auxr = cbind(auxgeor,auxarear)
	#	} else {
	#		auxr = auxgeor
	#	}
	#} else if (!is.null(auxarear)){
	#		auxr = auxarear
	#} else {
	#	print(paste("WARNING: no super matrix calculated. Review cut parameters"))
	#}
	auxr = cbind(auxgeor,auxarear,auxmetr)
	if (join){
		rownames(auxr) = auxrn
	}
	return(auxr)
}


add_mass_res_per_clustering = function(db){

	tam = 1:length(db)
	aux = db
	for (i in tam){
		aux[[i]] = add_res_per_clustering(db[[i]])
	}
	return(aux)
}


add_mass_res_geom_center_per_clustering = function(db){

	tam = 1:length(db)
	aux = db
	for (i in tam){
		aux[[i]] = add_res_geom_center_per_clustering(db[[i]])
	}
	return(aux)
}


add_mass_ellipsoid_per_clustering = function(db,tamclus){

	tam = 1:length(db)
	aux = db
	for (i in tam){
		aux[[i]] = add_ellipsoid_per_clustering(db[[i]],2:tamclus[i])
	}
	return(aux)
}

#Faz varredura em range dos deslocamentos para estimar melhor quantidade de clusters por id
find_mass_better_k_clustering_per_id = function(maw,maxpos,range){

	aux=c()
	for (i in range){
		temp = find_better_k_clustering_per_id(maw,maxpos,i)
		#print(temp)
		#readline()
		aux=cbind(aux,temp)
	}
	colnames(aux)=range
	return(aux)
}

#multi_view_mass_pam_silhouette_avg_width = function(dbaAAu[index],1:length(index),pdbids[index],xlim=c(1,12),ylim=c(0.0,1.05), main=main,xlab=xlab,ylab=ylab,vleg=c(10,1.05),1)

multi_view_mass_pam_silhouette_avg_width = function(db,pdbids,mfrow,tam,type=1,xlim=c(1,12),ylim=c(0.3,1.05), main=main,xlab=xlab,ylab=ylab,vleg=c(10,1.05)){

	if (!is.null(mfrow)) par(mfrow=mfrow)
	#tam = (mfrow[1]*mfrow[2])
	#tamj = 1:mfrow[2]
	tamimpar = tam%%2
	par(mar=c(0,4,4,0),oma=c(2.2,1,1,1))
	
	if (tamimpar) tam=tam-1
	
	i=2
	while (i<=tam){
		index = c(i-1,i)
		view_mass_pam_silhouette_avg_width(db[index],pdbids[index],xlim,ylim,main="",xlab="",ylab,vleg,type,ret=FALSE,mtext="", cex.leg=0.8, cex.pch=0.6,pch=21, mtext.cex=0.8,line=0.3)
		i=i+2
	}
	if (tamimpar) {
		index = c(tam+1)
		view_mass_pam_silhouette_avg_width(db[index],pdbids[index],xlim,ylim,main="",xlab="",ylab,vleg,type,ret=FALSE,mtext="", cex.leg=0.8, cex.pch=0.6,pch=21, mtext.cex=0.8,line=0.3)
	}
}

# View silhouette plot para métrica avg_width em função do número de clusters
# Efeito colateral: retorna matriz com os avg_width
view_mass_pam_silhouette_avg_width = function(db,pdbids,xlim,ylim,main,xlab,ylab,vleg,type,ret=TRUE,mtext="", cex.leg=0.8,cex.pch=0.6,pch=21, mtext.cex=0.8,line=0.3,zlim=-0.05){

	tamdb = 1:length(db)
	#print(tamdb);print(db[[1]]$pdbname);print(db[[2]]$pdbname)
	tampam = length(db[[1]]$exp$pam)
	#print(tampam)
	aux=c()
	auxneg=c()
	auxall = c()
	legtext = c()
	plot(1:tampam,seq(-1,1,length.out=tampam),type="n",ylim=ylim,xlim=xlim,main=main,xlab=xlab,ylab=ylab)#;print("ok");readline()
	#if (mtext=="") mtext("",cex=mtext.cex,line=line)
	#else
	#print(mtext)
	mtext(mtext,cex=mtext.cex,line=line)
	for (i in tamdb){
		#par(mfrow=c(1,1))
		#tampam = length(db[[i]]$exp$pam)
		aux = c()
		auxneg = c()
		for (j in 2:tampam){
			temp = db[[i]]$exp$pam[[j]]$silinfo$widths[,3]
			auxsum = sum(as.numeric(temp<zlim))
			if (auxsum>0) auxneg = c(auxneg,"BLACK")
			else {
				temp1 = as.logical(temp>=zlim)
				temp2 = as.logical(temp<(-zlim))
				temp12 = temp1*temp2
				auxsum = sum(as.numeric(temp12))
				if (auxsum>0) auxneg = c(auxneg,"GRAY")
				else auxneg = c(auxneg,"WHITE")
			}
			if (type==1) aux = c(aux,db[[i]]$exp$pam[[j]]$silinfo$avg.width)
			if (type==2) {
					aux = c(aux,mean(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
					#print(mean(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths));print(sd(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
			}
			if (type==3) {
					aux = c(aux,sd(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
					#print(mean(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths));print(sd(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
			}
			
			#plot(i,dbpam[[i]]$silinfo$avg.width,type="b")
		}
		
		lines(2:tampam,aux,lty=i)
		#print(auxneg)
		#readline()
		points(2:tampam,aux,pch=pch,bg=auxneg,cex=cex.pch)
		legtext = c(legtext,db[[i]]$pdbname)
		auxall = rbind(auxall,aux)
	}
	legend(vleg[1],vleg[2],legtext,lty=tamdb,cex=cex.leg)
	rownames(auxall)=pdbids[tamdb]
	colnames(auxall)=2:tampam
	if (ret) return(auxall)
}

# OBSOLETED: View silhouette plot para métrica avg_width em função do número de clusters
# Efeito colateral: retorna matriz com os avg_width
view_mass_pam_silhouette_avg_width_OLD = function(db,tamdb,pdbids,xlim,ylim,main,xlab,ylab,vleg,type,ret=TRUE,mtext="", cex.leg=0.8,cex.pch=0.6,pch=21, mtext.cex=0.8,line=0.3){

	#tamdb = length(db)
	tampam = length(db[[1]]$exp$pam)
	aux=c()
	auxneg=c()
	auxall = c()
	legtext = c()
	plot(1:tampam,seq(-1,1,length.out=tampam),type="n",ylim=ylim,xlim=xlim,main=main,xlab=xlab,ylab=ylab)
	#if (mtext=="") mtext("",cex=mtext.cex,line=line)
	#else
	#print(mtext)
	mtext(mtext,cex=mtext.cex,line=line)
	for (i in tamdb){
		#par(mfrow=c(1,1))
		#tampam = length(db[[i]]$exp$pam)
		aux = c()
		auxneg = c()
		for (j in 2:tampam){
			auxsum = sum(as.numeric(db[[i]]$exp$pam[[j]]$silinfo$widths[,3]<0))
			if (auxsum>0) auxneg = c(auxneg,"BLACK")
			else auxneg = c(auxneg,"WHITE")
			if (type==1) aux = c(aux,db[[i]]$exp$pam[[j]]$silinfo$avg.width)
			if (type==2) {
					aux = c(aux,mean(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
					#print(mean(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths));print(sd(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
			}
			if (type==3) {
					aux = c(aux,sd(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
					#print(mean(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths));print(sd(db[[i]]$exp$pam[[j]]$silinfo$clus.avg.widths))
			}
			
			#plot(i,dbpam[[i]]$silinfo$avg.width,type="b")
		}
		lines(2:tampam,aux,lty=i)
		#print(auxneg)
		#readline()
		points(2:tampam,aux,pch=pch,bg=auxneg,cex=cex.pch)
		legtext = c(legtext,db[[i]]$pdbname)
		auxall = rbind(auxall,aux)
	}
	legend(vleg[1],vleg[2],legtext,lty=tamdb,cex=cex.leg)
	rownames(auxall)=pdbids[tamdb]
	colnames(auxall)=2:tampam
	if (ret) return(auxall)
}

# Visualiza PAM conforme gráfico padrão do k-medoids no R
view_mass_pam_silhouette = function(db,tam,k,mfrow,type=2,which=2,text=" - Silhouette plot - k="){

	#text = 
	#sep = ": "
	#k = 5
	par(mfrow=mfrow)
	for (i in tam){
		#print(paste(db[[i]]$pdbname,text,k,sep=""))
		view_pam_silhouette(db[[i]],tam=k,which=which,type=type,main=paste(db[[i]]$pdbname,text,sep=""))
		#view_mass_pam_silhouette(dbarea3nzaau[[4]],k,which=2,type=2,main=paste(dbarea3nzaau[[4]]$pdbname,text,k,sep=""))
	}
}


# View silhouette plot conforme padrão no R
view_pam_silhouette = function(db,tam,which,type,main){

	#tam = length(db$exp$pam)

	if (type == 1){
		par(mfrow=c(1,(tam-1)))
		for (i in 2:tam){
		#for (i in tam){
			plot(db$exp$pam[[i]],which.plot=which)
		}
	}
	if (type == 2){
		#par(mfrow=c(1,(tam-1)))
		#for (i in 2:tam){
		for (i in tam){
			plot(db$exp$pam[[i]],which.plot=which,main=paste(main,i))
		}
	}

}

# Faz transformação em massa para calculo da area de contato a partir da distância
make_mass_silveira_romanelli_transformation = function(db,vdw,probe,liminf,p){

	tam = length(db)
	aux = list()
	for (i in 1:tam){
		#aux[[i]] = only_connected(db[[i]],apply(db[[i]]$A$a,1,Norm)!=0) # elimina linhas zero
		#aux[[i]] = apply_filter(db[[i]],apply(db[[i]]$A$a,1,Norm)!=0) # elimina linhas zero
		aux[[i]] = silveira_romanelli_transformation(db[[i]],vdw,probe,liminf,p)
	}
	return(aux)
}

# # OBSOLETED
initialize_mass_pam_clustering_old = function(db,oneg=F){

	tam = length(db)
	#n = dim(db[[i]]$Lw$eig$v)[1]

	for (i in 1:tam){

		if (!is.null(db[[i]]$name)){
			n = length(db[[i]]$name)
		}else{
			print(paste("WARNING: no size dimention found for",db[[i]]$pdbname))
		}
		if (oneg){
			db[[i]]$exp$pam[[1]] = list()
			db[[i]]$exp$pam[[1]]$clustering = rep(1,n) 
		}else{
			db[[i]]$exp$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
		}
	}
}


make_mass_pam_superclus = function(db,lowcut=c(0.5,0.5),conex=F,ncores=0,init=1,oneg=1,type=1,general=F){

	tam = length(db)
	aux=list()

	print(paste("Making PAM SUPER..."))
	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores at the level of db"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			make_pam_superclus(db[[i]],lowcut,ncores=ncores)
		}
	}else{
		for (i in 1:tam){
			#print(i)
			print(paste("Parallelized with",ncores,"cores at the level of pam"))
			aux[[i]] = make_pam_superclus(db[[i]],lowcut,ncores=ncores)
		}
	}
	return(aux)

}

#Faz PAM sobre db para número de clusters k em range
make_mass_pam_clustering = function(db,maxrange,init=1,type=1,oneg=T,conex=F,verbose=T,general=F,ncores=0,partype=1){

	tam = length(db)
	#no = n-k+1
	#mlogic = diag(rep(1,n))>0
	aux=list()

	print(paste("Making PAM..."))
	if (ncores&partype){
		print(paste("Doing it parallel with",ncores,"cores at the level of db"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			if (init==1) db[[i]]$exp$pam = list()
			make_pam_k_scanning(db[[i]],maxrange,type=type,oneg=oneg,conex=conex,general=general,init=init)
		}
	}else{
		for (i in 1:tam){
			if (verbose){
				if (!general){
					print(paste("Making PAM for",db[[i]]$pdbname))
				}else{
					print(paste("Making PAM for",db[[i]]$matrixname))
				}
			}
			if (init==1) db[[i]]$exp$pam = list()
			if (ncores){
				print(paste("Parallelized with",ncores,"cores at the leve of pam"))
				aux[[i]] = make_pam_k_scanning_par(db[[i]],maxrange,type=type,oneg=oneg,conex=conex,general=general,init=init,ncores=ncores)
			}else{
				print(paste("Semi-sequencial"))
				aux[[i]] = make_pam_k_scanning(db[[i]],maxrange,type=type,oneg=oneg,conex=conex,general=general,init=init)
			}

		}
	}
	return(aux)

}

#make_mass_par_pam_clustering = function(db,maxrange,init=1,type=1,oneg=T,conex=F,verbose=T,general=F,ncor=3){

#	tam = length(db)
	#no = n-k+1
	#mlogic = diag(rep(1,n))>0
	#registerDoMC(ncor)
#	cl<-makeCluster(ncor,outfile="debug.txt")
#	registerDoParallel(cl)
#	ptime = system.time({
#	res=foreach(i=1:tam, .export=c("make_pam_k_scanning","proc.time","pam")) %dopar% {
#		print("oi")
#		if (init==1) db[[i]]$exp$pam = list()
		#make_pam_k_scanning(db[[i]])
#		make_pam_k_scanning(db[[i]],maxrange,type=type,oneg=oneg,conex=conex,general=general,init=init)
		#sqrt(i)

#		}
#	})[3]
#	stopCluster(cl)
#	print(ptime)
#	return(res)

#}


testa = function(x){

	return(sum(sqrt(x)))

}

#make_mass_par_test = function(db,ncor=5){

#	registerDoMC(ncor)
#	tam = length(db[[1]]$exp$pam)
#	i = 1:tam
#	ptime <- system.time({
#		res = foreach(i=1:tam, .combine=c) %dopar% {
#			sum(sqrt(db[[1]]$exp$pam[[i]]$clustering)^2)
#			#testa(db[[1]]$exp$pam[[i]]$clustering)
#			#sqrt(i)
#		}
#	})[3]
#	print(ptime)
#	print(res)
#}

# # OBSOLETED:Faz PAM sobre db para número de clusters k em range
make_mass_pam_clustering_old = function(db,range,expidinit,type=1,oneg=F,conex=F){

	#tam = length(db)
	expid = expidinit

	initialize_mass_pam_clustering(db,oneg)

	for (k in range){
		print(k)
		if (type==1){
			db = mass_experiment_02(db,expid,k,oneg,conex)
		}
		if (type==2){
			db = mass_experiment_03(db,expid,k,oneg,conex)
		}
		expid = expid + 1
	}
	return(db)

}

# # OBSOLETED: EXPERIMENTO 02: PAM com k partições, Lw - EIG
mass_experiment_02_old = function(db,expid,k,oneg=T,conex=F,verbose=T){

	tam = length(db)
	#no = n-k+1
	#mlogic = diag(rep(1,n))>0

	for (i in 1:tam){
		#expam = list()
		n = dim(db[[i]]$Lw$eig$v)[1]
		#n = dim(db[[i]]$Lw$svd$us)[1]
		no = n-k+1
		#if (oneg){
		#	db[[i]]$exp$pam[[1]] = list()
		#	db[[i]]$exp$pam[[1]]$clustering = rep(1,n) 
		#}else{
		#	db[[i]]$exp$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
		#}
		#print(paste(db[[i]]$pdbname,n,k))
		if (k>=n){
			if (verbose) print(paste("WARNING:k=",k,"must be less than n=",n,"in",db[[i]]$pdbname))
			#db[[i]]$exp$pam[[expid]] = 0
		}else{
			aux = pam(db[[i]]$Lw$eig$v[,no:n],k)
			#aux = pam(db[[i]]$Lw$svd$us[,no:n],k)
			#exp1 = db[[i]]$exp$vlist[[1]]
			#exp2 = make_vector_indicator_list(aux$clustering)	
			#exp = make_logical_vector_product(exp1,exp2)
			#db[[i]]$exp$vlist[[expid]]=exp
			db[[i]]$exp$pam[[expid]] = aux
			#exp = list()
			#exp = make_vector_indicator_list(db[[i]]$ide)
			#db[[i]]$exp[[expid]]=exp
		}
	}
	return(db)	
}


# # OBSOLETED: EXPERIMENTO 03: PAM com k partições, Lw - SVD
mass_experiment_03_old = function(db,expid,k,oneg=F,verbose=T){

	tam = length(db)
	#no = n-k+1
	#mlogic = diag(rep(1,n))>0

	for (i in 1:tam){
		#expam = list()
		#n = dim(db[[i]]$Lw$eig$v)[1]
		n = dim(db[[i]]$Lw$svd$v)[1]
		no = n-k+1
		#if (oneg){
		#	db[[i]]$exp$pam[[1]] = list()
		#	db[[i]]$exp$pam[[1]]$clustering = rep(1,n) 
		#}else{
		#	db[[i]]$exp$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
		#}
		#db[[i]]$exp$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
		#aux = pam(db[[i]]$Lw$eig$v[,no:n],k)
		#print(no);print(n)
		#print(db[[i]]$Lw$svd$u[,no:n]);readline()
		if (k>=n){
			if (verbose) print(paste("WARNING:k=",k,"must be less than n=",n,"in",db[[i]]$pdbname))
			#db[[i]]$exp$pam[[expid]] = 0
		}else{
			aux = pam(db[[i]]$Lw$svd$v[,no:n],k)#;print("ok")
			#exp1 = db[[i]]$exp$vlist[[1]]
			#exp2 = make_vector_indicator_list(aux$clustering)	
			#exp = make_logical_vector_product(exp1,exp2)
			#db[[i]]$exp$vlist[[expid]]=exp
			db[[i]]$exp$pam[[expid]] = aux

			#exp = list()
			#exp = make_vector_indicator_list(db[[i]]$ide)
			#db[[i]]$exp[[expid]]=exp
		}
	}
	return(db)	
}

set_mass_hotspot_candidates = function(db,dball,ncores=0){


	tam = length(dball$res)#;print(tam)
	aux = list()

	print(paste("Hotspot candidates..."))


	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			set_hotspot_candidates(db,dball$res[[i]])
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
		#for (i in 1){
			aux[[i]]=set_hotspot_candidates(db,dball$res[[i]])
		}
	}
	dball$hot = aux
	return(dball)
	
}


make_mass_self_similarity_assessment = function(db,dball,ncores=0,dbhot=NULL,rotref=diag(4)){

	#pdblist = as.character(unique(dball$group$pdb))### ARRUMAR, mas acho que nao dara problema mais...
	pdbframe = unique(dball$group[c("pdb","empty")])
	pdblist = as.character(pdbframe[,1])
	pdbstate = pdbframe[,2]
	#print(pdbframe)#;readline()
	#print(pdblist)#;readline()
	#print(pdbstate);readline()
	tam = length(pdblist)#;print(tam)

	aux=list()

	print(paste("Computing self graph similarity"))

	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			new_make_similarity_assessment(db=db,dball=dball,id0=i,pdblist=pdblist,self=T)
			
		}
	}else{
		print(paste("Doing it sequential..."))#;print(tam);readline()
		for (i in 1:tam){#;print(i)
		#for (i in 3){
			aux[[i]] = new_make_similarity_assessment(db=db,dball=dball,id0=i,pdblist=pdblist,pdbstate,self=T)
			
		}
	}
	dball$res = aux#;print(length(dball$res));readline()	
	dball = new_consolidate_pdb_self_similarity_assessment(db,dball,pdblist)
	#print(dball);readline()	
	return(dball)

}


make_mass_similarity_assessment = function(db,dball,dbsel,id,self=F,ncores=0){


#	pdbname = toupper(pdblist[idb])
#	auxres$pdbname = pdbname
#	auxgrp = which(dball$group[,"pdb"]==pdbname)

	pdblist = as.character(unique(dball$group$pdb))
	#pdblist = as.character(unique(dball$group$pdb))[1:2]
	pdbname = toupper(pdblist[id])
	groupid = which(dball$group[,"pdb"]==pdbname)
	dball$res = list()
	parg12 = list()
	stat = c()
	print(paste("Computing graph similarity"))

	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		#parg12 = make_queue_for_par(db=db,dball=dball,dbhot=dbsel$hot,groupid=groupid)
		parg12 = make_queue_for_par(db=db,dball=dball,dbhot=NULL,groupid=groupid)
		#return(parg12)
		tam = length(parg12)
		#print(parg12[[1]]$parg1);print(parg12[[1]]$parg2);print(length(parg12));readline()
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
		#aux = foreach(i=11) %do% {
			new_best_align_graph_pdb(parg12[[i]]$parg1,parg12[[i]]$parg2)
		}
		#return(aux)
		dball = resume_similarity_assessment(dball=dball,parg=parg12,score=aux,pdblist=pdblist)
		dball$res[[1]]$parg = parg12
		dball$res[[1]]$groupid = groupid
	
	}else{
		print(paste("Doing it sequential..."))
		#dball$res[[1]] = new_make_similarity_assessment(db=db,dball=dball,id0=id,pdblist=pdblist,dbhot=dbsel$hot,self=self)
		dball$res[[1]] = new_make_similarity_assessment(db=db,dball=dball,id0=id,pdblist=pdblist,dbhot=NULL,self=self)
		dball$res[[1]] = new_consolidate_pdb_similarity_assessment(db,dball$res[[1]],pdblist)
	}
	dball$hot = dbsel$hot
	
	return(dball)

}

#make_mass_similarity_assessment = function(db,dball,id,pre=T,rotref=diag(4),version="pre"){
make_mass_similarity_assessment_old2 = function(db,dball,dbhot=NULL,id=1,rotref=diag(4),tracktime=T,self=F){

	
	if (tracktime) {
		print("TRACKING PROCESS TIME FOR SIMILARITY...")
		t=c()
		b=60
		rtt = 1
		t = round(c(proc.time()[3],t),rtt)
	}

	pdblist = as.character(unique(dball$group$pdb))
	#pdblist = pdblist[1:2]

	if (self) {
		#dball$sel = list()
		tami = 1:length(pdblist)
	}
	else tami = id

	if (is.null(dball$res)){
		dball$res = list()
	}

	if (is.null(id)){
		dball$res = list()
		tami = 1:length(pdblist)
		if (self) print(paste("Computing self graph similarity for all database..."))
		else print(paste("Computing graph similarity for all database..."))
	}else{
		if (self) print(paste("Computing self partial graph similarity"))
		else  print(paste("Computing partial graph similarity"))
	}


	#for (i in 1:tami){
	for (i in tami){
		#print(pdblist[i])
		auxr = list()
		pdbname = toupper(pdblist[i])#;print(pdbname);readline()
		if (self) print(paste("Making SELF similarity alignment for",pdbname))
		else print(paste("Making NEW similarity alignment for",pdbname))
		auxr$pdbname = pdbname
		auxgrp = which(dball$group[,"pdb"]==pdbname)#;print(auxgrp);readline()
		#if (self) auxgrp = auxgrp[1]
		#auxgrp = 1
		auxr$groupid = auxgrp
		auxr$alignment = list()
		rotid =  grep(pdbname,dball$rotref)
		if (length(rotid)!=1){
			print(paste("WARNING: something wrong with rotref for",pdbname,"getting the first..."))
			rotid=rotid[1]
		}
		#rotref = dball$rotref[[rotid]]$rot
		#print(rotid);print(rotref);readline()
		#if (version=="pre"){
			#auxr$alignment = pre_make_similarity_assessment(db=db,dball=dball,id=auxgrp,pdbname=pdbname,rotref=rotref)#;print("ok")
		#}
		#if (version=="new"){
		#if (self){
		auxr$alignment = new_make_similarity_assessment(db=db,dball=dball,dbhot=dbhot,id=auxgrp,pdbname=pdbname,rotref=rotref,w=c(1,1),self=self)#;print("ok")
		#else{
		#	auxr$alignment = new_make_similarity_assessment(db=db,dball=dball,id=auxgrp,pdbname=pdbname,rotref=rotref,w=c(1,1),self=self)
		#}
		#}
		#print(rotid);print(rotref);readline()
		#dball$res[[i]] = make_pdb_similarity_assessment(db=db,dball=dball,pdbname=pdbname,rotref=rotref,pre=pre)
		#dball$res[[i]]$scorend = consolidate_pdb_similarity_assessment(dball$res[[i]]$alignment)
		#auxr = consolidate_pdb_similarity_assessment(db,auxr,pdblist)
		#auxr = new_consolidate_pdb_similarity_assessment(db,dball$res[[i]],pdblist)
 		#auxr = new_consolidate_pdb_similarity_assessment(db,auxr,pdblist)
		#print(auxr);readline()
		if (self){
			dball$res[[i]] = auxr
		}
		else{
			#auxr = new_consolidate_pdb_similarity_assessment(db,auxr,dball$scoreself,pdblist)
			auxr = new_consolidate_pdb_similarity_assessment(db,auxr,pdblist)
			#auxr = new_consolidate_pdb_similarity_assessment(db,dball$res[[i]],dball$scoreself,pdblist)
			dball$res[[i]] = auxr

		}

	}

	if (self){
		dball = new_consolidate_pdb_self_similarity_assessment(db,dball,pdblist)
#		dball$scoreself = new_consolidate_pdb_self_similarity_assessment(db,dball$res,pdblist)
		#dball$scorend = new_consolidate_pdb_self_similarity_assessment(db,dball$res,pdblist)

	}

	if (tracktime){
		t = round(c(proc.time()[3],t),rtt)
		print(paste("PROCESS TIME SIMILARITY FOR",tami,"PDBs:",round((t[1]-t[2])/b,rtt),"minutes"))
	}
	return(dball)

}


make_mass_superclus_reorientation = function(db,ncores=0,refid=5,verbose=T){

	tam = length(db)

	aux = list()

	print(paste("Superclus reorientation"))
	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			#tamclus = length(dba2[[i]]$exp$superclus)
			#if (refid>tamclus) refid = tamclus
			#m = dba2[[i]]$exp$superclus[[refid]]$geomc
			#rotref = get_rot_matrix(source=m,target=diag_x(3),svd=c(T,F),center=c(T,F))
			make_superclus_reorientation(db=db[[i]],rotref=rotref,refid=refid)
		}
	}else{
		print(paste("Doing it sequential"))
		for (i in 1:tam){
			if (verbose) print(paste("Making supercluster reorientation for",db[[i]]$matrixname))
			#tamclus = length(dba2[[i]]$exp$superclus)
			#if (refid>tamclus) refid = tamclus
			#m = dba2[[i]]$exp$superclus[[refid]]$geomc
			#rotref = get_rot_matrix(source=m,target=diag_x(3),svd=c(T,F),center=c(T,F))
			aux[[i]] = make_superclus_reorientation(db=db[[i]],rotref=rotref,refid=refid)
			#aux[[i]]$exp$rotref = rotref				
		}
	}
	return(aux)


}


# Monta agrupamento superclus global a partir de agrupamentos locais em db
make_mass_superclus = function(db,ncores=0,only_last=NULL,only_a=F,ma=NULL,verbose=T,geo=T){

	tam = length(db)

	aux = list()

	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			make_superclus(db[[i]],only_last=only_last,only_a=only_a,ma=ma,geo=geo)
		}
	}else{
		print(paste("Doing it sequential"))
		for (i in 1:tam){

			if (verbose) print(paste("Making supercluster for",db[[i]]$matrixname))

			aux[[i]] = make_superclus(db[[i]],only_last=only_last,only_a=only_a,ma=ma,geo=geo)
			#if (geo){
			#	aux[[i]] = add_geometric_center_per_super_clustering(aux[[i]])
				#aux[[i]] = add_dists_matrix_per_super_clustering(aux[[i]])
			#}
		}
	}
	return(aux)


}

threshold_mass_filter = function(db,k=0.5,diagz=T){

	tam = 1:length(db)

	aux = list()

	for (i in tam){

		aux[[i]] = threshold_filter(db[[i]],k,diagz)

	}
	return(aux)

}


make_mass_cutoff_scanning = function(db,init=0,step=0.01,imax=0.5){

	tam = 1:length(db)

	aux = list()

	for (i in tam){

		aux[[i]] = cutoff_scanning(db[[i]],init=init,step=step,imax=imax)

	}
	return(aux)

}


# Filtra os nós conforme o tipo de resíduo
make_mass_res_type_filter = function(db,vfilter,type,enzaa,inbaa,enz="E",inb="I"){

	tam = length(db)
	aux = list()
	if (type < 4){
		for (i in 1:tam){
			if (type==1) auxfactor = as.factor(db[[i]]$ap)
			if (type==2) auxfactor = as.factor(db[[i]]$apg)
			if (type==3) auxfactor = as.factor(db[[i]]$ape) #definição mais indicada
			levels(auxfactor) = vfilter
			aux[[i]] = apply_filter(db[[i]],as.logical(auxfactor)) 
		}
		return(aux)
	}
	if (type==4){
		for (i in 1:tam){
			x=db[[i]]$ape
			y=db[[i]]$ei
			ve=y==enz
			vi=y==inb
			#print(ve)
			#fi=x[vi]==aa
			#fe=x[ve]!=aa
			fi=x[vi] %in% inbaa
			fe=x[ve] %in% enzaa
			auxfactor=c(fe,fi)
			#print(x[auxfactor])
			#readline()
			aux[[i]] = apply_filter(db[[i]],as.logical(auxfactor))
		}
		return(aux)
	}
}


# OBSOLETED: Filtra os nós conforme o tipo de resíduo
make_mass_res_type_filter_old = function(db,vfilter,type=3){

	tam = length(db)
	aux = list()
	for (i in 1:tam){
		if (type==1) auxfactor = as.factor(db[[i]]$ap)
		if (type==2) auxfactor = as.factor(db[[i]]$apg)
		if (type==3) auxfactor = as.factor(db[[i]]$ape) #definição mais indicada
		levels(auxfactor) = vfilter
		aux[[i]] = apply_filter(db[[i]],as.logical(auxfactor)) 
	}
	return(aux)
}


# Faz transformação em massa eliminando linhas zero
make_mass_only_connected = function(db,ncores=0,verbose=T){

	tam = length(db)
	auxa = list()
	aux = list()
	if (verbose) print(paste("Applying only-connected filter"))

	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			ids = apply(db[[i]]$A$a,1,Norm)!=0
			apply_filter(db[[i]],ids)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			#aux[[i]] = only_connected(db[[i]],apply(db[[i]]$A$a,1,Norm)!=0) # elimina linhas zero
			ids = apply(db[[i]]$A$a,1,Norm)!=0#;print(ids)
			aux[[i]] = apply_filter(db[[i]],ids) # elimina linhas zero
			aux[[i]]$zeros = which(!ids)

		}
	}
	return(aux)
}

# Faz input em massa de dados conforme PDBids em pdbids.
make_mass_input = function(dirname,pdbids,sufix,meso,aacode,ncores=0,matrixtype="normal",prefix="",typename="",verbose=TRUE,bin=FALSE,sep=","){

	if (verbose) print("Entrada de dados")
	tam = length(pdbids)
	aux = list()

	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			input_data(filename=paste(dirname,prefix,pdbids[i],sufix,sep=""),pdbname=pdbids[i],typename=typename,meso=meso,aacode=aacode, sufix=sufix,matrixtype=matrixtype,bin=bin,sep=sep)
		}

	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			#print(paste(dirname,prefix,pdbids[1,i],sufix,sep=""))
			#readline()
			#aux[[i]] = input_data(paste(dirname,prefix,pdbids[1,i],sufix,sep=""),as.character(pdbids[1,i]),typename,meso,aacode)
			#print(paste(dirname,prefix,pdbids[i],sufix,sep=""));print(dirname);print(prefix);print(pdbids[i]);print(sufix)
			aux[[i]] = input_data(filename=paste(dirname,prefix,pdbids[i],sufix,sep=""),pdbname=pdbids[i],typename=typename,meso=meso,aacode=aacode,sufix=sufix,matrixtype=matrixtype,bin=bin,sep=sep)

		}
	}

	return(aux)

}


#Multiplot para curvas de avg_widt em função da quantidade de clusters separados pela maxpos
multiplot_mass_pam_silhouette_avg_width = function(db,pdbids,tam,maxpos,mfrow,main="",xlab="",ylab="",vleg=c(10,1.05), xlim=c(1,12),ylim=c(0.0,1.05),type=1){

	#tam = levels(as.factor(maxpos))
	par(mfrow=mfrow)
	
	for (i in as.numeric(tam)){
	
		filter = maxpos %in% i
		#print(filter)
		#print(sum(filter))
		#readline()
		temp=view_mass_pam_silhouette_avg_width(db[filter],1:sum(filter),pdbids[filter],mtext=paste("maxclus=",i+1),xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,vleg=vleg,type=type)
	}

}

#OBSOLETED: Multiplot para curvas de avg_widt em função da quantidade de clusters separados pela maxpos
multiplot_mass_pam_silhouette_avg_width_old = function(db,maxpos,pdbids,mfrow,main="",xlab="",ylab="",vleg=c(10,1.05), xlim=c(1,12),ylim=c(0.0,1.05),type=1){

	tam = levels(as.factor(maxpos))
	par(mfrow=mfrow)
	
	for (i in tam){
	
		filter = maxpos %in% i
		#print(filter)
		#print(sum(filter))
		#readline()
		temp=view_mass_pam_silhouette_avg_width(db[filter],1:sum(filter),pdbids[filter],mtext=paste("maxclus=",i+1),xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,vleg=vleg,type=type)
	}

}


#add_mass_center_node = function(db){


	#aux = db
#	aux = list()

#	tam = 1:length(db)

#	for (i in tam){
#		aux[[i]] = add_center_node(db[[i]])

#	}
#	return(aux)

#}

#Gera novo db com base no vetor indicador dado em PAM, levando em conta a quantidade de elementos conexos em maxpos 
make_mass_pam_connected_filter = function(db,maxpos,d=1){

	aux = list()
	tam = 1:length(db)
	#print("ok")
	#print(tam)
	for (i in tam){
		#print(maxpos[i])
		if (maxpos[i]>0){
			auxclusv = db[[i]]$exp$pam[[maxpos[i]+d]]$clustering
			#print(auxclusv)
			auxmax = find_max_occurrence(auxclusv)
			#print(auxmax)
			#readline()
			filter = auxclusv==auxmax
			aux[[i]]=apply_filter(db[[i]],filter)
		} else {
			aux[[i]]=db[[i]]
		}
	}
	return(aux)
}


# Produz em massa as Laplacianas L e/ou Lw
make_mass_laps = function(db,ncores=0,anorm=FALSE,r=10,only=c(T,T)){

	#tam = length(pdbids)
	tam=length(db)
	#aux = db
	aux = list()

	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			laps(db[[i]],anorm=anorm,r=r,only=only)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			#ppfaz = laps_and_decomps(ppfaz)
			#aux[[i]] = laps_and_decomps(aux[[i]])
			if (anorm){
				print(paste("Making Laplacians L e La for",db[[i]]$pdbname))
			}else{
				print(paste("Making Laplacians L e Ls for",db[[i]]$pdbname))
			}
			aux[[i]] = laps(db[[i]],anorm=anorm,r=r,only=only)

		}
	}
	return(aux)

}

# Produz em massa as Decomposições EIG e SVD
make_mass_decomps = function(db,L,Lw,EIG,SVD,ncores=0,karpack=0,r=10){

	auxl = c(L,Lw)
	auxln = c("L","Lw")
	auxd = c(EIG,SVD)
	auxdn = c("EIG","SVD")
	#tam = length(pdbids)
	tam=length(db)
	#aux = db
	aux = list()

	if (ncores){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		aux = foreach(i=1:tam) %dopar% {
			decomps(m=db[[i]],L=L,Lw=Lw,EIG=EIG,SVD=SVD,r=r,karpack=karpack)
		}
	}else{
		for (i in 1:tam){
			#ppfaz = laps_and_decomps(ppfaz)
			print(paste("Making",paste(auxdn[auxd],collapse=","),"in",paste(auxln[auxl],collapse=","),"with r=",r,"for",db[[i]]$pdbname))
			aux[[i]] = decomps(m=db[[i]],L=L,Lw=Lw,EIG=EIG,SVD=SVD,r=r,karpack=karpack)#;readline()

		}
	}
	return(aux)

}

# Produz em massa as Laplacianas L e Lw e Decomposições EIG e SVD
make_mass_laps_and_decomps = function(db){

	tam = length(pdbids)
	aux = db

	for (i in 1:tam){

		#ppfaz = laps_and_decomps(ppfaz)
		aux[[i]] = laps_and_decomps(aux[[i]])

	}

	return(aux)

}

# Plot a matrix estatística conforme os vetores singulares definidos em v.
plot_mass_stat_matrix = function(db,tam,v,pdbids,xlab,ylab,main,aps=1,mtext="",SVD=TRUE,cex=0.6,mtext.cex=0.8,line=0.3){
	
	for (i in tam){
		if (SVD){
			#print(v)
			#print(db[[i]]$A$svd$us[,v[2]])
			#readline()
			x = db[[i]]$A$svd$us[,v[1]]
			y = db[[i]]$A$svd$us[,v[2]]
			#plot(db[[i]]$A$svd$us[,v[1]],dks[[i]]$A$svd$us[,v[2]],xlab=xlab,ylab=ylab,main=main,pch="")
			#text(db[[i]]$A$svd$us[,v[1]],dks[[i]]$A$svd$us[,v[2]], pdbids,cex=cex)
			plot(x,y,asp=asp,type="n",xlab=xlab,ylab=ylab,main=main)
			text(x,y,pdbids,cex=cex)
			if (mtext=="") mtext(db[[i]]$type,cex=mtext.cex,line=line)
			else mtext(mtext,cex=mtext.cex,line=line)
		}
	}
}

# Checa quantidade de linhas zeros
make_mass_check_zero_lines = function(dbd,dbd2,pdbids){

	aux=c()
	for (i in 1:length(pdbids)){
		a1 = dim(dbd[[i]]$A$a)[1]
		a2 = dim(dbd2[[i]]$A$a)[1]
		#a3 = dim(dba[[i]]$A$a)[1]
		diff = a1-a2
		#diff2 = a1-a3
		aux = rbind(aux,c(dbd[[i]]$pdbname,a1,a2,diff))#,a3,diff2))
	}
	print(aux)
}

# Faz decomposição EIG ou SVD da matriz estatística em db
make_mass_stat_decomps = function(db,EIG=FALSE,SVD=TRUE){

	#aux = list()
	#aux[[1]]=make_stat_matrix(db=db,center=0,inter=0,p=p)
	tam = length(db)
	for (i in 1:tam){
		#db[[i]]=make_stat_matrix(db,i,p=p)
		db[[i]]$A=make_decomp(db[[i]]$A,EIG=EIG,SVD=SVD)
		#aux[[i]]=make_decomp(aux[[i]]$a,EIG=EIG,SVD=SVD)
	}
	return(db)


}

# OBSOLETED
# transforma uma matrix estatistica (ex: ks.test()) numa matrix centralizada (com diagonal zero, se o caso), com interação para recalcular
# a media geral apos zerar a diagonal (se o caso)
# Adaptada para fazer sem centralização (Normal) e centralizada, sem testar diagonal = 0.
# OBS: função longe do ideal...
make_mass_stat_matrix = function(db,center=0,inter=4,EIG=FALSE,SVD=TRUE,p=3){

	aux = list()
	#aux[[1]]=make_stat_matrix(db=db,center=0,inter=0,p=p)
	for (i in 1:inter){
		#aux[[i]]=make_stat_matrix(db,center,i-2,p=p)
		aux[[i]]=make_stat_matrix(db,i,p=p)
		aux[[i]]$A=make_decomp(aux[[i]]$A,EIG=EIG,SVD=SVD)
		#aux[[i]]=make_decomp(aux[[i]]$a,EIG=EIG,SVD=SVD)
	}
	return(aux)
}


# OBSOLETED: transforma uma matrix estatistica (ex: ks.test()) numa matrix com diagonal zero, com interação para recalcular
# a media geral apos zerar a diagonal
make_mass_stat_matrix_OLD = function(db,center,inter,EIG=FALSE,SVD=TRUE,p=3){

	aux = list()
	for (i in 1:inter){
		aux[[i]]=make_stat_matrix(db,center,i-1,p=p)
		aux[[i]]$A=make_decomp(aux[[i]]$A,EIG=EIG,SVD=SVD)
		#aux[[i]]=make_decomp(aux[[i]]$a,EIG=EIG,SVD=SVD)
	}
	return(aux)
}

#make_mass_matrix_pam_silhouette_profile = function(db){

	

#}


# Gera scripts para visualização dos clusters no Pymol
# COM BUGS !!! 
make_mass_pymol_partition = function(db,tam,expid,dirpathin,dirpathout,filebase,edgename,colorlist,vertexsufix,edpart,sufixfile=".pml"){

	#tam = length(db)
	#print("ok")
	for (i in tam){
		pdbname = db[[i]]$pdbname
		transname = pymol_transname(db[[i]]$name)
		filename = paste(pdbname,filebase,sufixfile,sep="")
		min = min(db[[i]]$weight)
		max = max(db[[i]]$weight)
		if (length(edpart)==1) vpart = c(1,min,max,edpart)
		else vpart = c(0,edpart)
		# BUG: auxlevels poe as cadeias em ordem, e no PDB elas nem sempre estão. Ex: 1VEQ, B->ENZ e A->INB
		# BUG: ainda não trata mais de uma cadeia em ENZ e INB 
		# GAMBIARRA ate achar uma solução definitiva: pega 1a e última posição em ide
		auxlevels = c(db[[i]]$ide[1],db[[i]]$ide[length(db[[i]]$ide)])
		#auxlevels = levels(as.factor(db[[i]]$ide))
		#print(auxlevels)
		lineall = pymol_view_by_vector_indicators(db[[i]]$A$a,db[[i]]$exp$vlist[[expid]],dirpathin,db[[i]]$name,pdbname,filename,transname,edgename,colorlist,vpart,vertexsufix,auxlevels,db[[i]]$ei)
		#print(lineall)
		write(lineall,file=paste(dirpathout,filename,sep=""))
		#print(paste(dirpathout,filename,sep=""))
	}

}


# Prepara vetores indicadores para visualização cadeia - cadeia (enzima - inibidor)
preprocess_for_mass_experiment = function(db,expid,verbose=FALSE){

	tam = length(db)

	#mlogic = diag(rep(1,n))>0

	for (i in 1:tam){
		exp = list()
		#exp = make_vector_indicator_list(db[[i]]$ide,verbose)
		#exp = make_vector_indicator_list(db[[i]]$ei,verbose)
		exp = 1
		#if (verbose) print(paste(db[[i]]$pdbname,levels(as.factor(db[[i]]$ide))))
		db[[i]]$exp$vlist[[expid]]=exp
	}

	return(db)
	
}


# Visualiza matrizes como um mapa de contato em cores: azul (negativos) ao vermelho (positivos)
view_mass_matriz = function(db,tam,type,pcex,text,mfrow,sym=1,pal=2){

	#tam = length(db)
	par(mfrow=mfrow)

	for (i in tam){
		if (type=="A"){
			view_matrix(db[[i]]$A$a,pcex=pcex,minv=min(db[[i]]$A$a),maxv=max(db[[i]]$A$a),sym=sym,pal=pal,main.text=paste(text,db[[i]]$pdbname))
		}
		if (type=="L"){
			view_matrix(db[[i]]$L$a,pcex=pcex,minv=min(db[[i]]$L$a),maxv=max(db[[i]]$L$a),sym=sym,pal=pal,main.text=paste(text,db[[i]]$pdbname))
		}
		if (type=="Lw"){
			view_matrix(db[[i]]$Lw$a,pcex=pcex,minv=min(db[[i]]$Lw$a),maxv=max(db[[i]]$Lw$a),sym=sym,pal=pal,main.text=paste(text,db[[i]]$pdbname))
		}
	}

}

#Mostra estatisticas ks-test e wilcox-rank-sum-test para distribuição de pesos
make_mass_stats_weight_distribution = function(db,tam,type,main,verbose=FALSE,p=2){

	#tam = length(db)
	tamlen = length(tam)
	auxm=matrix(rep(0,tamlen*tamlen),nrow=tamlen,ncol=tamlen)
	#print(length(db))
	auxnames=c()
	#text = "p-value of"
	auxstat=c("ks-test","wilcox-Rank-Sum-test")
	if (verbose) print(paste(main,": p-value for",auxstat[type]))
	#x = 1
	#for (i in tam){
	for (i in 1:tamlen){
		auxnames = c(auxnames,db[[i]]$pdbname)
		#y = 1
		#for (j in tam){
		for (j in 1:tamlen){
			if (type==1){
				#print(i);print(j)
				x1 = db[[i]]$weight
				x2 = db[[j]]$weight
				aux = ks.test(x1,x2)
				auxm[i,j] = round(aux$p.value,p)
				#auxm[x,y] = round(aux$p.value,p)
			}
			if (type==2){
				x1 = db[[i]]$weight
				x2 = db[[j]]$weight
				#print(x1);readline()
				aux = wilcox.test(x1,x2)
				auxm[i,j] = round(aux$p.value,p)
				#auxm[x,y] = round(aux$p.value,p)
			}
			#y = y+1

		}
		#x=x+1
	}
	colnames(auxm) = auxnames
	rownames(auxm) = auxnames
	if (verbose) print(auxm)
	return(auxm)
}


# Plot grafico para a distribuição de pesos dos grafos em db
view_mass_weight_distribution = function(db,ylim,main,xlab,ylab,vleg,cex,col=c(),mtext="",line=0.5){

	tam = length(db)
	legtext = c()

	if (is.null(col)) plot(density(db[[1]]$weight),ylim=ylim,main=main,xlab=xlab,ylab=ylab)
	else plot(density(db[[1]]$weight),ylim=ylim,main=main,xlab=xlab,ylab=ylab,col=col[1])
	mtext(mtext,line=line,cex=cex)
	legtext = c(legtext,db[[1]]$pdbname)
	for (i in 2:tam){
		if (is.null(col)) lines(density(db[[i]]$weight),lty=i)
		else lines(density(db[[i]]$weight),lty=i,col=col[i])
		legtext = c(legtext,db[[i]]$pdbname)
	}
	if (is.null(col)) legend(vleg[1],vleg[2],legtext,lty=1:tam,cex=cex)
	else legend(vleg[1],vleg[2],legtext,lty=1:tam,cex=cex,col=col)

}

# OBSOLETED: Plot grafico para a distribuição de pesos dos grafos em db
view_mass_weight_distribution_OLD = function(db,ylim,main,xlab,ylab,vleg,cex,col=c(),mtext="",line=0.5){

	tam = length(db)
	legtext = c()

	if (is.null(col)) plot(density(db[[1]]$weight),ylim=ylim,main=main,xlab=xlab,ylab=ylab)
	else plot(density(db[[1]]$weight),ylim=ylim,main=main,xlab=xlab,ylab=ylab,col=col[1])
	mtext(mtext,line=line,cex=cex)
	legtext = c(legtext,db[[1]]$pdbname)
	for (i in 2:tam){
		if (is.null(col)) lines(density(db[[i]]$weight),lty=i)
		else lines(density(db[[i]]$weight),lty=i,col=col[i])
		legtext = c(legtext,db[[i]]$pdbname)
	}
	if (is.null(col)) legend(vleg[1],vleg[2],legtext,lty=1:tam,cex=cex)
	else legend(vleg[1],vleg[2],legtext,lty=1:tam,cex=cex,col=col)

}

# Plot grafico para correlação de pesos de arestas entre area e distância
view_mass_weight_dist_correlation = function(dbarea,dbdist,tam,ylim,xlim,main,xlab,ylab,vleg,cex,type,mfrow,codtp){

	#tam = length(dbarea)
	legtext = c()

	par(mfrow=mfrow)
	x=dbarea[[1]]$ide=="i"
	y=dbarea[[1]]$ide=="e"
	acpvector = codify_residue_type(dbarea[[1]]$ap[y],dbarea[[1]]$ap[x],codtp)
	plot(as.vector(dbarea[[1]]$A$a[x,y]),as.vector(dbdist[[1]]$A$a[x,y]),xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,type="p",pch=paste(acpvector),cex=cex)
	if (type == 2) legend(vleg[1],vleg[2],dbarea[[1]]$pdbname,pch=1:tam,cex=cex)

	if (type == 1) legtext = c(legtext,dbarea[[1]]$pdbname)
	i = 2
	while (i <= tam){
		x=dbarea[[i]]$ide=="i"
		y=dbarea[[i]]$ide=="e"
		acpvector = codify_residue_type(dbarea[[i]]$ap[y],dbarea[[i]]$ap[x],codtp)
		if (type == 1) {
			points(as.vector(dbarea[[i]]$A$a[x,y]),as.vector(dbdist[[i]]$A$a[x,y]),pch=paste(acpvector))
			legtext = c(legtext,dbarea[[i]]$pdbname)
		}
		if (type == 2) {
			plot(as.vector(dbarea[[i]]$A$a[x,y]),as.vector(dbdist[[i]]$A$a[x,y]),xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,type="p",pch=paste(acpvector),cex=cex)
			legend(vleg[1],vleg[2],dbarea[[i]]$pdbname,pch=i,cex=cex)
		}
		i = i + 1

	}
	if (type==1) legend(vleg[1],vleg[2],legtext,pch=1:tam,cex=cex)

}

# Checa consistencia das matrizes A, L, Lw em termos de simetria
view_mass_symmetric_check = function(db){

	tam = length(db)

	for (i in 1:tam){
		aux=symmetric_check(db[[i]]$A$a,paste(db[[i]]$pdbname,"- A"),db[[i]]$typename,TRUE)
		aux=symmetric_check(db[[i]]$L$a,paste(db[[i]]$pdbname,"- L"),db[[i]]$typename,TRUE)
		aux=symmetric_check(db[[i]]$Lw$a,paste(db[[i]]$pdbname,"- Lw"),db[[i]]$typename,TRUE)
		print(" ")
	}

}

# Verifica em massa número de componentes conexos usando Lw
view_mass_connected_component_number = function(db,p){

	tam = length(db)
	print("-----------------------------------")

	for (i in 1:tam){
		connected_component_number(db[[i]],p=p)
	}

}

make_mass_similarity_assessment_old = function(db,dball,id,pre=T,rotref=diag(4)){

	pdblist = as.character(unique(dball$group$pdb))

	#tami = length(pdblist)
	tami = id

	if (is.null(id)){
		dball$res = list()
		tami = 1:length(pdblist)
		print(paste("Computing graph similarity for all database..."))
	}else{
		print(paste("Computing partial graph similarity for..."))
	}


	#for (i in 1:tami){
	for (i in tami){
		#print(pdblist[i])
		pdbname = toupper(pdblist[i])#;print(pdbname);readline()
		rotid =  grep(pdbname,dball$rotref)
		if (length(rotid)!=1){
			print(paste("WARNING: something wrong with rotref for",pdbname,"getting the first..."))
			rotid=rotid[1]
		}
		rotref = dball$rotref[[rotid]]$rot
		#print(rotid);print(rotref);readline()
		dball$res[[i]] = make_pdb_similarity_assessment(db=db,dball=dball,pdbname=pdbname,rotref=rotref,pre=pre)
		#dball$res[[i]]$scorend = consolidate_pdb_similarity_assessment(dball$res[[i]]$alignment)

	}
	return(dball)

}

#make_it_mass_input(dirname,filenames[i],pdbxyz[[i]],meso,aacode,vcores)
make_it_mass_input_old = function(i,dirname,pdbxyz,meso,aacode,vcores){

#1PPF_E_I_BSR_ALL
	#pdbids = read.csv(paste(dirname,filename,sep=""),header=F,sep=",")
	#pdbids=paste(as.matrix(pdbids))
	#names(pdbids)=1:length(pdbids)
	if(i==0){
		

	}

	#Entrada 
	dba=list()
	#Agora tem opcao de usar matrizes esparsas em matrixtype
	dba = make_mass_input(dirname=dirname,pdbids=pdbids,sufix=".csv",meso=meso,aacode=aacode,ncores=vcores[1],matrixtype="sparse", bin=FALSE,sep=",")#;readline()
		 
	#Elimina linhas/colunas zeros, se houver
	dba = make_mass_only_connected(dba,ncores=vcores[2])#;readline()
	#Guarda em ...$pdb as coordenadas atomicas (PARALELIZADO! Ficou entre 3x a 4x mais rapido...)
	if (length(dba)>1){
		dba = add_mass_pdb_xyz(dba,dirpathin,ncores=vcores[3],typename=2)
	}else{
		dba = add_mass_pdb_xyz(dba,dirpathin,ncores=vcores[4],typename=2)
	}
	return(dba)

}

# OBSOLETED: Faz input em massa de dados conforme PDBids em pdbids - Old
make_mass_input_old = function(dirname,pdbids,sufixes,typename){

	tam = length(pdbids)
	aux = list()

	for (i in 1:tam){

		aux[[i]] = input_data(paste(dirname,pdbids[1,i],sufixes[1],sep=""),paste(dirname,pdbids[1,i],sufixes[2],sep=""),as.character(pdbids[1,i]),typename)

	}

	return(aux)

}

make_mass_all_super_clustering_old2 = function(db,geopar=NULL,areapar=NULL,metricpar=NULL,betpar=NULL,clopar=NULL,shopar=NULL,
lowcut=c(0.7,0.8),join=F,only_last=0,minpar=4,pre=F,sep="_",continuum=F,minclusid=4,dcut=c(0.2,0.01),verbose=T,maxclusid=6){

	tam = 1:length(db)
	#tam = 8
	#auxm = db
	auxgeor = c()
	auxarear = c()
	auxmetr = c()
	auxbetr = c()
	auxclor = c()
	auxshor = c()
	auxrn = c()
	auxrni = c()
	#auxncr = c()
	#auxi = make_super_matrix_geocut(db[[i]],geocutpar)
	#auxr = rbind(auxr,auxi)
	#auxgeo = NULL
	#auxarea = NULL
	if (verbose) print(paste("Generating info matrix..."))
	for (i in tam){
		#if (verbose) print(paste("Generating info matrix"))
		if (!is.null(geopar)){
			auxgeo = make_all_super_clustering(db=db[[i]],type="geo",par=c(geopar,lowcut),id=i,minpar=minpar, pre=pre,join=join, only_last=only_last,continuum=continuum)
			#auxgeo = make_super_matrix(db=db[[i]],type="geo",par=c(geopar,lowcut),id=i,minpar=minpar,pre=pre)
			#print(auxgeo);readline()
			#if (only_last){
			#	auxdim = dim(auxgeo)[1]
			#	auxgeo = auxgeo[auxdim,]
				#auxgeor = rbind(auxgeor,auxgeo)
			#} else if (join){
				#auxdim = dim(auxgeo)[1]
			#	auxgeo=apply(auxgeo,2,sum)
				#auxgeo=c(auxgeo,nc=auxgeodim)
			#}
			auxgeor = rbind(auxgeor,auxgeo)
			
		}
		if (!is.null(areapar)){
			auxarea = make_all_super_clustering(db=db[[i]],type="area",par=c(areapar,lowcut),id=i,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
			#auxarea = make_super_matrix(db=db[[i]],type="area",par=c(areapar,lowcut),id=i,minpar=minpar,pre=pre)
			#print(auxarea);readline()
			#if (only_last){
			#	auxdim = dim(auxarea)[1]
			#	auxarea = auxarea[auxdim,]
			#}else if (join){
				#auxdim = dim(auxarea)[1]
			#	auxarea=apply(auxarea,2,sum)
				#auxarea=c(auxarea,nc=auxareadim)
			#}
			#print(auxarea);readline()
			auxarear = rbind(auxarear,auxarea)
		}
		if (!is.null(betpar)){
			auxbet = make_all_super_clustering(db=db[[i]],type="bet",par=c(betpar,lowcut),id=i,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
			#auxbet = make_super_matrix(db=db[[i]],type="bet",par=c(betpar,lowcut),id=i,minpar=minpar,pre=pre)
			#if (only_last){
			#	auxdim = dim(auxbet)[1]
			#	auxbet = auxbet[auxdim,]
			#}else if (join){
				#auxdim = dim(auxbet)[1]
			#	auxbet=apply(auxbet,2,sum)
				#auxbet=c(auxbet,nc=auxbetdim)
			#}
			auxbetr = rbind(auxbetr,auxbet)
		}
		if (!is.null(clopar)){
			auxclo = make_all_super_clustering(db=db[[i]],type="clo",par=c(clopar,lowcut),id=i,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
			#auxclo = make_super_matrix(db=db[[i]],type="clo",par=c(clopar,lowcut),id=i,minpar=minpar,pre=pre)
			#if (only_last){
			#	auxdim = dim(auxclo)[1]
			#	auxclo = auxclo[auxdim,]
			#}else if (join){
				#auxdim = dim(auxclo)[1]
			#	auxclo=apply(auxclo,2,sum)
				#auxclo=c(auxclo,nc=auxclodim)
			#}
			auxclor = rbind(auxclor,auxclo)
		}
		if (!is.null(shopar)){
			auxsho = make_all_super_clustering(db=db[[i]],type="sho",par=c(shopar,lowcut),id=i,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
			#auxsho = make_super_matrix(db=db[[i]],type="sho",par=c(shopar,lowcut),id=i,minpar=minpar,pre=pre)
			#if (only_last){
			#	auxdim = dim(auxsho)[1]
			#	auxsho = auxsho[auxdim,]
			#}else if (join){
				#auxdim = dim(auxsho)[1]
			#	auxsho=apply(auxsho,2,sum)
				#auxsho=c(auxsho,nc=auxshodim)
			#}
			#print(auxsho);readline()
			auxshor = rbind(auxshor,auxsho)
		}
		if (!is.null(metricpar)){
			auxmet = make_all_super_clustering(db=db[[i]],type="met",par=c(metricpar,lowcut),id=i,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
			#print(lowcut);readline()
			lastline = dim(auxmet)[1]
			lastgo = F
			if (auxmet[lastline,3]>maxclusid){
				print(paste("WARNING: found the limit of maxclusid=",maxclusid,"for",db[[i]]$pdbname,"Trimming it... "))
				auxmet = auxmet[1:maxclusid,]
			}else{
				#print(auxmet[lastline,3]);minclusid;readline()
				
				if (auxmet[lastline,3]==minclusid){
					print(paste("WARNING: found the limit of minclusid=",minclusid,"for",db[[i]]$pdbname,"Trying overtake it... "))
					lowcuti = lowcut
					lastgo = T
					repeat{
						lowcuti[1]=lowcuti[1]-dcut[1]
						if (lowcuti[1]<0){
							lowcuti[2] = lowcuti[2]-dcut[2]
							lowcuti[1] = lowcut[1]
						}
						if (lowcuti[2]<0){
								print(paste("WARNING: could not ovetake the minclusid for",db[[i]]$pdbname))
								lastgo = F
								break;
						}
						auxmet = make_all_super_clustering(db=db[[i]],type="met",par=c(metricpar,lowcuti),id=i,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
						#print(auxmet);print(lowcuti);readline()
						lastline = dim(auxmet)[1]
						if (auxmet[lastline,3]>minclusid) {
							auxmet = auxmet[1:(minclusid+1),]
							break;
						}
					}
				}
			}
			if (lastgo) print(paste("Done!"))

			#print(auxmet);readline()
			auxmetr = rbind(auxmetr,auxmet)
		}
		#if (!is.null(metricpar)){
		#if (0){
		#	auxmet = make_super_matrix_metric(db=db[[i]],metricpar=c(metricpar,lowcut),id=i,minpar=minpar,pre=pre)
			#auxdim = dim(auxmet)[1]
		#	auxrni = rownames(auxmet)
			#print(auxmet);readline()
		#	if (only_last){
		#		auxdim = dim(auxmet)[1]
				#auxnc = as.numeric(unlist(strsplit(rownames(auxmet)[auxdim],"_"))[3])
				#print(nc);readline()
				#print(auxmet);readline()
				#auxmet = auxmet[auxdim,]
		#		auxrni = rownames(auxmet)
		#		auxdim0 = auxdim-only_last+1
		#		if (auxdim0>0){
		#			auxmet = auxmet[auxdim0:auxdim,]
		#		}else{
		#			auxmet = auxmet[1:auxdim,]
		#		}
		#		if (!is.matrix(auxmet)){
		#			auxmet = t(as.matrix(auxmet))
		#			rownames(auxmet)=auxrni
		#		}
				#auxnc = as.numeric(unlist(strsplit(rownames(auxmet)[auxdim],"_"))[3])
				#print(auxmet);readline()
				#if (metricpar[3]){
				#	auxmet=c(auxmet,nc=auxnc)
				#}
				#auxrn = c(auxrn,paste(i,db[[i]]$pdbname,auxnc,sep=sep))
				#auxrni = rownames(auxmet)
		#	}
		#	if (join){
		#		auxdim = dim(auxmet)[1]
				#print(auxmet);readline()
				#auxmet = c(acc=apply(auxmet,2,sum))
				#auxmet = apply(auxmet,2,min)
		#		auxrni = rownames(auxmet)[auxdim]
				#auxmet = apply(auxmet,2,mean)
				#auxmet = apply(auxmet,2,min)
		#		auxmet = auxmet[auxdim,]
		#		if (!is.matrix(auxmet)){
		#			auxmet = t(as.matrix(auxmet))
		#			rownames(auxmet)=auxrni
		#		}
				#print(auxmet);print(auxrn);readline()
				#if (metricpar[3]){
				#	auxmet=c(auxmet,nc=auxdim)
				#}
				#auxrnj = rownames(auxrn[2])
				#auxrn = c(auxrn,auxrnj)
		#	}
		#	auxrn = c(auxrn,auxrni)
			#print(auxrn);readline()
		#	auxmetr = rbind(auxmetr,auxmet)
		#}
		#if (!is.null(ncpar)){
			#auxnc = make_super_matrix(db=db[[i]],type="nc",par=c(ncpar,lowcut),id=i,minpar=minpar,pre=pre)
			#print(auxgeo);readline()
			#if (join){
				#auxshodim = dim(auxsho)[1]
			#	auxsho=apply(auxsho,2,sum)
				#auxsho=c(auxsho,nc=auxshodim)
			#}
			
		#	auxncr = rbind(auxncr,auxnc)
		#}
		#if (join){
			#if (!is.null(ncpar)){
			#	auxnc = c(nc=auxdim)
			#	auxncr = rbind(auxncr,auxnc)
			#}
		#	auxrn = c(auxrn,paste(i,db[[i]]$pdbname,i,sep=sep))
			#auxrn = c(auxrn,paste(i,db[[i]]$pdbname,sep=sep))
		#}
		#if (only_last){
		#	auxrn = c(auxrn,paste(i,db[[i]]$pdbname,auxnc,sep=sep))

		#}
		#print(auxi);readline()
		#auxr = rbind(auxr,auxgeo)
		
	}
	#print(auxrn);readline()
	auxr1 = cbind(auxgeor,auxarear,auxbetr,auxclor,auxshor)
	#print(auxr1);print(auxmetr);readline()
	auxr2 = auxmetr
	#print(auxrn);readline() 
	#if (join|only_last){
	#	rownames(auxr1) = auxrn
	#	rownames(auxr2) = auxrn
	#}
	#print(auxr1);print(auxr2);readline()
	auxr = list()
	auxr$a = auxr1
	#print(auxr2);readline()
	colnames(auxr2) = c("gacc","lacc","verts","edges","idinput")
	auxr$info = as.data.frame(auxr2)

	return(auxr)
}


make_hyper_group_output_old = function(dba,sep="_"){

	tam = length(dba[[1]])

	#tam_pdb = tam/3

	group = list()

	k=0
	for (i in 1:tam){
		
		if (i%%3==1){
			k = k + 1
			group[[k]] = list()
			j = 1
		}
		group[[k]][[j]] = list()
		group[[k]][[j]]$pdb = toupper(substr(dba[[1]][[i]]$pdbname,1,4))
		if (i%%3==2) group[[k]][[j]]$pdb_type = 1
		else group[[k]][[j]]$pdb_type = 2
		group[[k]][[j]]$chains = dba[[1]][[i]]$chains
		comb = unlist(strsplit(dba[[1]][[i]]$pdbname,sep))#;print(comb) 
		group[[k]][[j]]$combination_chains = comb[2:3]
		
		group[[k]][[j]]$connections_group_bsr = list()
		group[[k]][[j]]$connections_group_bsr_polar_polar = list()
		group[[k]][[j]]$connections_group_bsr_nonpolar_nonpolar = list()

		if (!is.null(dba[[1]][[i]]$A)){
			group[[k]][[j]]$connections_group_bsr =  get_connections_group_bsr(dba[[1]][[i]]$exp$superclus)
		
		}
		if (!is.null(dba[[2]][[i]]$A)){
			group[[k]][[j]]$connections_group_bsr_polar_polar =  get_connections_group_bsr(dba[[2]][[i]]$exp$superclus)
		
		}
		if (!is.null(dba[[3]][[i]]$A)){
			group[[k]][[j]]$connections_group_bsr_nonpolar_nonpolar =  get_connections_group_bsr(dba[[3]][[i]]$exp$superclus)
		
		}
		j = j + 1

	}
	return(group)

}


make_it_web_output_old = function(dba,dbsel,dbgrp,vcores){

	tam = length(dba)
	aux=list()
	for (i in 1:tam){
		aux[[i]] = resume_for_web_output(dba,dbsel,dbgrp,id=i)	
	}
	return(aux)
	#resume_for_web_output(dba,dbsel,id)

}

#x = resume_for_web_output(db=dba2[[1]], dball=dbsel[[1]][[1]],id=1);print(x)
make_hyper_web_output_old = function(dba,dball,dbgrp,icores=3,vcores){

	tam = length(dba)
	aux = list()
	print(paste("Outputs..."))
	if (icores){
		print(paste("Doing it parallel with",icores,"cores"))
		registerDoMC(icores)
		aux = foreach(i=1:tam) %dopar% {
			make_it_web_output(dba[[i]],dbsel[[i]],dbgrp[[i]],vcores)
		}
	}else{
		print(paste("Doing it sequential..."))
		#for (i in 1:tam){
		for (i in 1){
			aux[[i]] = make_it_web_output(dba[[i]],dbsel[[i]],dbgrp[[i]],vcores)
		}

	}
	return(aux)

}

make_align_output_old = function(db,sep="_",n=2){

	aux = list()
	if (length(db)){

		tit1 = unlist(strsplit(as.character(db$tabs$pos1[1,1]),sep))
		tit2 = unlist(strsplit(as.character(db$tabs$pos2[1,1]),sep))

		#print(tit1[4]);readline()

		aux$pdb1 = tit1[1]
		aux$pdb2 = tit2[1]
	
		if (tit1[2]=="CHAIN") aux$type = 1
		else aux$type = 2

		aux$combination_chains = tit1[2:3]

		aux$polarity = tit1[4]

		superclus = list()
		superclus[[1]] = list()
		superclus[[1]]$a = db$align[[1]]$p1$a
		superclus[[1]]$geomc = db$align[[1]]$p1$xyz
		superclus[[1]]$rot = diag(4)
		superclus[[1]]$score_align = NA
		superclus[[1]]$score_group = db$align[[1]]$t1[1,6:7]

		superclus = c(superclus,add_superclus_from_align(db$align))


		aux$connections_group_bsr = list()
		aux$connections_group_bsr_polar_polar = list()
		aux$connections_group_bsr_nonpolar_nonpolar = list()

		aux$volume_bsr = list()
		aux$volume_bsr_polar_polar = list()
		aux$volume_bsr_nonpolar_nonpolar = list()

		aux$rot_bsr = list()
		aux$rot_bsr_polar_polar = list()
		aux$rot_bsr_nonpolar_nonpolar = list()

		aux$score_align_bsr = list()
		aux$score_align_bsr_polar_polar = list()
		aux$score_align_bsr_nonpolar_nonpolar = list()

		aux$score_group_bsr = list()
		aux$score_group_bsr_polar_polar = list()
		aux$score_group_bsr_nonpolar_nonpolar = list()

		#print(superclus);readline()

		if (tit1[4]=="ALL"){
			res =  get_connections_group_bsr(superclus)
			aux$connections_group_bsr = res$conn
			aux$volume_bsr = res$vol
			aux$rot_bsr = res$rot
			aux$score_align_bsr = res$score_align
			aux$score_group_bsr = res$score_group

		}

		if (tit1[4]=="PP"){
			res =  get_connections_group_bsr(superclus)
			aux$connections_group_bsr_polar_polar = res$conn
			aux$volume_bsr_polar_polar = res$vol
			aux$rot_bsr_polar_polar = res$rot
			aux$score_align_bsr_polar_polar = res$score_align
			aux$score_group_bsr_polar_polar = res$score_group

		}

		if (tit1[4]=="AA"){
			res =  get_connections_group_bsr(superclus)
			aux$connections_group_bsr_nonpolar_nonpolar = res$conn
			aux$volume_bsr_nonpolar_nonpolar = res$vol
			aux$rot_bsr_nonpolar_nonpolar = res$rot
			aux$score_align_bsr_nonpolar_nonpolar = res$score_align
			aux$score_group_bsr_nonpolar_nonpolar = res$score_group

		}
	}else{
		print(paste("WARNING: dbali is empty"))
	}
	#print(res);readline()
	
	return(aux)


}

#make_mass_hot_assessment_bad = function(db,dball,rotref=diag(4),tracktime=T){

	
#	if (tracktime) {
#		print("TRACKING PROCESS TIME FOR HOT ASSESSEMENT...")
#		t=c()
#		b=60
#		rtt = 1
#		t = round(c(proc.time()[3],t),rtt)
#	}

#	pdblist = as.character(unique(dball$group$pdb))
	#pdblist = pdblist[1:2]

#	tami = 1:length(pdblist)
#	tami = id

#	if (is.null(dball$res)){
#		dball$res = list()
#	}

#	if (is.null(id)){
#		dball$res = list()
#		tami = 1:length(pdblist)
#		print(paste("Computing graph similarity for all database..."))
#	}else{
#		print(paste("Computing partial graph similarity"))
#	}


	#for (i in 1:tami){
#	for (i in tami){
		#print(pdblist[i])
#		auxr = list()
#		pdbname = toupper(pdblist[i])#;print(pdbname);readline()
#		print(paste("Making hot assessement for",pdbname))
#		auxr$pdbname = pdbname
#		auxgrp = which(dball$group[,"pdb"]==pdbname)#;print(auxgrp);readline()
		#auxgrp = 1
#		auxr$groupid = auxgrp
#		auxr$alignment = list()
#		rotid =  grep(pdbname,dball$rotref)
#		if (length(rotid)!=1){
#			print(paste("WARNING: something wrong with rotref for",pdbname,"getting the first..."))
#			rotid=rotid[1]
#		}
#		rotref = dball$rotref[[rotid]]$rot
		#print(rotid);print(rotref);readline()
		#if (version=="new"){
#		print(auxgrp);readline()
#		auxr$alignment = make_hot_assessment(db=db,dball=dball,id=auxgrp[1],pdbname=pdbname,rotref=rotref,w=c(1,1))#;print("ok")
#		print(auxr$alignment);readline()
		#}
		#print(rotid);print(rotref);readline()
		#dball$res[[i]] = make_pdb_similarity_assessment(db=db,dball=dball,pdbname=pdbname,rotref=rotref,pre=pre)
		#dball$res[[i]]$scorend = consolidate_pdb_similarity_assessment(dball$res[[i]]$alignment)
		#auxr = consolidate_pdb_similarity_assessment(db,auxr,pdblist)
		#auxr = new_consolidate_pdb_similarity_assessment(db,dball$res[[i]],pdblist)
# 		auxr = new_consolidate_pdb_similarity_assessment(db,auxr,pdblist)
		#print(auxr);readline()
#		dball$res[[i]] = auxr

#	}

#	if (tracktime){
#		t = round(c(proc.time()[3],t),rtt)
#		print(paste("PROCESS TIME SIMILARITY FOR",tami,"PDBs:",round((t[1]-t[2])/b,rtt),"minutes"))
#	}
	#return(dball)

#}

