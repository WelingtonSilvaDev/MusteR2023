library("doMC")
library("bio3d")
library("dplyr")   
library("flexclust")


get_chain_combinations = function(pdb,water="WT",join="+"){

	#chains = unique(pdb$atom$chain)
	chains = pdb$chains

#[1] "E"  "I"  "L1" "L2" "L3" "L4" "L5" "L6" "WT"
	tam = length(chains)

	if (tam>1){

		f = chains != water

		if (sum(f)>1){
			chains = chains[f]
			chains = c(chains,water)
		}else{
			print(paste("WARNING: no chain water found in PDB:",pdb$pdbname))
		}
		comb = data.frame(chains[1:2])
		#;print(comb);readline()
		i = 3
		repeat{
			if (i>tam) break;
			join_chain = paste(comb[,i-2],collapse=join)
			comb = cbind(comb,data.frame(c(join_chain,chains[i])))
			i = i + 1
			#;print(comb);readline()
		}
	}else{
		print(paste("WARNING AND STRANGE: this PDB:",pdb$pdbname,"has no interface"))
	}
	colnames(comb)=NULL
	return(comb)

}

#
# Prepara uma matriz somente com os elementos do tipo ATOM
# o retorno é uma lista contendo as cadeias e uma matriz no formato x, y, z
# os elementos estão nomeados como A.GLN19.OE1 nas linhas, e xyz nas colunas
#
prepare_atom_matrix = function(pdb, chain_a_name, chain_b_name, hetatm=F,join="+"){
  
	chain_join=c()
	chain_a_name = (as.vector(unlist(chain_a_name)))
	chain_b_name = (as.vector(unlist(chain_b_name)))
	if (grepl(join,chain_a_name,fixed=T)){#print(chain_a_name)
		chain_join = unlist(strsplit(chain_a_name,join,fixed=T))#;print(chain_a_name);print(chain_join)#;readline()
		tam_a = length(chain_join)
		f = c()
		for (i in 1:tam_a){
			f = c(f,which(pdb$atom$chain %in% chain_join[i]))
		}
		#f = pdb$atom$chain %in% chain_join#;print(sum(f))#;readline()
		chain_a = pdb$atom[f,]#;print(dim(pdb$atom))#;print(dim(chain_a));readline() 
		#chain_a_1 = chain_join[1]
		#chain_a_2 = chain_join[3]#;print(chain_join[3]);readline()
		#chain_a = filter(pdb$atom, chain == chain_a_1 | chain == chain_a_2)
	}else{
		#chain_a = filter(pdb$atom, chain == chain_a_name, type == "ATOM")
		chain_a = filter(pdb$atom, chain == chain_a_name)
	}
	#print(chain_a$resid);readline()
	if (grepl(join,chain_b_name,fixed=T)){
		chain_join = unlist(strsplit(chain_b_name,join,fixed=T))#;print(chain_b_name);print(chain_join)
		tam_b = length(chain_join)
		f = c()
		for (i in 1:tam_b){
			f = c(f,which(pdb$atom$chain %in% chain_join[i]))
		}
		f = pdb$atom$chain %in% chain_join
		chain_b = pdb$atom[f,]  
		#chain_b_1 = chain_join[1]
		#chain_b_2 = chain_join[3]
		#chain_b = filter(pdb$atom, chain == chain_b_1 | chain == chain_b_2)
	}else{
		#chain_b = filter(pdb$atom, chain == chain_b_name, type == "ATOM")
		chain_b = filter(pdb$atom, chain == chain_b_name)#;print(dim(chain_b));readline()
	}
	#print(chain_b$resid);readline() 
  #chain_a = filter(pdb$atom, chain == chain_a_name, type == "ATOM")
  #chain_b = filter(pdb$atom, chain == chain_b_name, type == "ATOM")
	#chain_all = filter(pdb$atom, type == "ATOM")
	chain_all = pdb
  
	#print(chain_a);readline()
  cx = c(chain_a$x, chain_b$x)
  cy = c(chain_a$y, chain_b$y)
  cz = c(chain_a$z, chain_b$z)
  mat = cbind(cx, cy, cz)
	#print(mat[1:3,]);readline()
  
  #Nível de resíduo
  rownames(mat) <- c(paste0(chain_a_name, '.', chain_a$resid, chain_a$resno, '.' ,chain_a$elety, '-', chain_a$eleno), paste0(chain_b_name, '.', chain_b$resid, chain_b$resno, '.' ,chain_b$elety, '-', chain_b$eleno))

	#print(mat[1:3,]);readline()
  
  pdb_list_filter <- list(
		pdbname = pdb$pdbname,
    chain_a = chain_a,
    chain_b = chain_b,
    chain_a_name = chain_a_name,
    chain_b_name = chain_b_name,
    xyz_matrix = mat,
    chain_all = chain_all
  )
  
  return(pdb_list_filter)
}

# calcula a distância euclideana entre duas matrizes - Cadeias A e B por exemplo
# o retorno é uma lista no formato:
#
# 1 - nome do resíduo seguindo a estrutura cadeia A + cadeia B por conta da matriz de adjascencia
# 2 - nome do átomo seguindo a estrutura cadeia A + cadeia B por conta da matriz de adjascencia
# 3 - ordem do resíduo no PDB
# 4 - Nome da cadeia
# 5 - matriz de adjascencia de distâncias entre os átomos, linha por colunas
#
#  Exemplo:
#  atom  aa chain naa     X      Y      Z
#  1  N   ILE     A   1 7.322 85.054 61.124
#
#
get_euclidean_distance = function(pdb_matrix,b_intra=F,join="+"){
  
  chain_a = pdb_matrix$chain_a
  chain_b = pdb_matrix$chain_b

  #tamanho de linhas do vetor da cadeia A, medindo a quantidade de linhas
  lines_a = length(chain_a[,1])
  #tamanho de linhas do vetor da cadeia B, medindo a quantidade de linhas
  lines_b = length(chain_b[,1])
  
  #tamanho de linhas e colunas da matriz
  tam_matriz = lines_a + lines_b
  
  #inicialização da matriz de distâncias
  matriz_adj_distancias = matrix(0,tam_matriz,tam_matriz)
  
  mat_dist = dist2(chain_a[9:11],chain_b[9:11])

  for(i in 1:lines_a){
    for(j in 1:lines_b){        
      matriz_adj_distancias[i,lines_a+j] = mat_dist[i,j]
      matriz_adj_distancias[lines_a+j,i] = mat_dist[i,j]        
    }
  }

	if (b_intra){
		mat_dist = dist2(chain_b[9:11],chain_b[9:11])
		for(i in 1:lines_b){
    	for(j in 1:lines_b){        
      matriz_adj_distancias[lines_a+i,lines_a+j] = mat_dist[i,j]
      matriz_adj_distancias[lines_a+j,lines_a+i] = mat_dist[i,j]        
    	}
  	}

	}
  #;print("o")
  names_chain_a <- c(paste0(chain_a$chain, '.', chain_a$resid, chain_a$resno, '.', chain_a$elety, '-', chain_a$eleno, '-',chain_a$elesy,'-',substr(chain_a$type,1,1)))
  names_chain_b <- c(paste0(chain_b$chain, '.', chain_b$resid, chain_b$resno, '.', chain_b$elety, '-', chain_b$eleno, '-',chain_b$elesy,'-',substr(chain_b$type,1,1)))#;print(names_chain_a);readline()
  names <- c(names_chain_a, names_chain_b)
  
  colnames(matriz_adj_distancias) <- names
  rownames(matriz_adj_distancias) <- names
  
  return(matriz_adj_distancias)
  #return(Matrix(matriz_adj_distancias))
}

calculate_distance_matrix_for = function(j,combinations,pdb,water="WT"){

	chain_a_name = combinations[1,j]
	chain_b_name = combinations[2,j]

	#preparando para o formato do import
	pdb_matrix = prepare_atom_matrix(pdb,chain_a_name,chain_b_name)#;print(chain_b_name)
	
	#print(pdb_matrix)

	# Calculando a distância Euclideana
	#print(chain_a_name);print(chain_b_name)#;readline()
	chain_b_name = as.vector(unlist(chain_b_name))
	#print(chain_b_name!=water);readline()

	if (chain_b_name!=water){
		res = get_euclidean_distance(pdb_matrix,b_intra=F)
		#print("fui")
	}else{
		res = get_euclidean_distance(pdb_matrix,b_intra=T)
		#print("foi")
	}
	
	return(res)
}

unify_dist_matrices = function(dist){

	tam = length(dist)#;print(tam);readline()
	new_dist = dist[[tam]]
	if (tam>1){	
		dist[[tam]]=NULL
		tam = length(dist)#;print(tam);readline()
		while (tam){
			tamd = dim(dist[[tam]])[1]
			new_dist[1:tamd,1:tamd] = dist[[tam]]
			dist[[tam]]=NULL
			tam = length(dist)#;print(tam);readline()
		}
	}
	return(new_dist)
}

calculate_distance_matrix = function(combinations,pdb,ncores,water="WT"){

	dist = list()
	res = list()
	if(ncores){
		registerDoMC(ncores)
		dist = foreach(j=1:length(combinations[1,])) %dopar%{
			calculate_distance_matrix_for(j,combinations,pdb,water)
			}
	}else{
		# Testando mapply no lugar do for. Nao deu diferenca significativa.
		#v = 1:length(combinations[1,])
		#dist = mapply(make_euclidean_distance,v,MoreArgs = list(combinations=combinations,pdb=pdb))
		for (j in 1:length(combinations[1,])){
			dist[[j]] = calculate_distance_matrix_for(j,combinations,pdb,water)
		}
	}
	#;print("0")
	#res[[1]] = Matrix(unify_dist_matrices(dist))
	res[[1]] = unify_dist_matrices(dist)

	return(res)

}

#	mxyz = matrix(pdb$xyz,ncol=3,byrow=T)
#	f = pdbxyz$atom$eleno %in% db$atom_ids
#	auxyz$xyz = mxyz[f,]
#	auxyz$atid = pdbxyz$atom$eleno[f] 

sep_atom_name_id = function(atom_list,sep="[-]",ncol=5){

	
	matrix(unlist(mapply(strsplit,atom_list,MoreArgs = list(split=sep))),ncol=ncol,byrow=T)
	#if(type==1)
	#	matrix(unlist(mapply(strsplit,atom_list,MoreArgs = list(split=sep))),ncol=2,byrow=T)
	#else if(type==2)
	#	matrix(unlist(mapply(strsplit,atom_list,MoreArgs = list(split=sep))),ncol=3,byrow=T)
	#else if(type==3)
	#	matrix(unlist(mapply(strsplit,atom_list,MoreArgs = list(split=sep))),ncol=4,byrow=T)
	#else
	#	matrix(unlist(mapply(strsplit,atom_list,MoreArgs = list(split=sep))),ncol=5,byrow=T)
}

get_pdb_essentials_for = function(dist,pdb,sepn="[-]"){

	res = list()
	#mxyz = matrix(pdb$xyz,ncol=3,byrow=T)
	atom_names = colnames(dist)#;print(head(atom_names))
	atom_frame = sep_atom_name_id(atom_names)#;print(atom_frame);readline()
	#elen = colnames(dist)#;print(elen[1:10]);readline()
	#elen = matrix(unlist(mapply(strsplit,elen,MoreArgs = list(split=sepn))),ncol=2,byrow=T)#;print(elen);readline()
	#res$element_name = elen[,1]
	atom_ids = as.numeric(atom_frame[,2])
	#filter = which(pdb$atom$eleno %in% atom_ids)#;print(filter);readline()
	filter = match(atom_ids,pdb$atom$eleno)
	#res$xyz = mxyz[filter,]
	res$xyz = pdb$atom[filter,9:11]
	res$eleno = pdb$atom$eleno[filter]
	res$chains = pdb$chains
	#res$filter = filter
	#res$atomids = atom_ids
	#res$atom = pdb$atom[filter,]
	#res$atom2 = pdb$atom
	#res$title = pdb$title
	#res$ligands = pdb$ligands
	#res$element_name[filter]
	#res$atom_ids[filter]
	return(res)

}

get_pdb_essentials = function(dist,pdb,ncores){

	tam = length(dist)
	res=list()
	if(ncores){
		registerDoMC(ncores)
		res = foreach(j=1:tam) %dopar%{
			get_pdb_essentials_for(dist[[j]],pdb)
		}
	}else{
		for (j in 1:tam){
			res[[j]] = get_pdb_essentials_for(dist[[j]],pdb)
		}
	}
	return(res)
}


download_pdb = function(PDBS,pdbpath,sufix=".pdb"){
  for(i in 1:length(PDBS)){
    destfile = paste0(pdbpath,PDBS[i],sufix)
    if(!file.exists(destfile)){
      err = try(download.file(url=paste0('https://files.rcsb.org/view/', PDBS[i], sufix),destfile=destfile))
			if (class(err)=="try-error"){
				print(paste("WARNING: PDBid",PDBS[i],"not found"))
				next
			}
    }
  }
}


#expand_combinations_to_ligands = function(combinations,ligands){
#}

add_polarity_for = function(polarity,pdbid){

	m = polarity
	res = list()
	res$ALL = list()
	res$ALL$bsr = m
	#res$ALL$filter = rep(T,dim(m)[1])
	res$PP = list()
	res$AA = list()
	if (dim(m)[1]>1){
		atom_names = colnames(m)#;print(atom_names);readline()
		atom_codes = sep_atom_name_id(atom_names)[,5]#;print(atom_codes);readline()#;print(head(atom_codes));readline()
		#atom_codes = new_codify_polarity3(atom_names,meso,aacode)
		#atom_codes = new_codify_polarity6(atom_names,meso,aacode)
		
		filter = atom_codes=="a"
		if (sum(!filter)>0){
			res$PP$bsr = m[!filter,!filter] ### cercar possibilidade de BSRPP ou BSRAA ser tudo zero
			f = apply(res$PP$bsr,1,Norm)!=0
			res$PP$bsr = res$PP$bsr[f,f]
			#res$PP$filter = !filter
		}else{
			print(paste("WARNING: PDB",pdbid,"has generated empty matrix for PP"))
			#res$PP$bsr = NULL
			res$PP$bsr = Matrix()
		}
		if (sum(filter)>0){
			res$AA$bsr = m[filter,filter]
			f = apply(res$AA$bsr,1,Norm)!=0
			res$AA$bsr = res$AA$bsr[f,f]
			#res$AA$filter = filter
		}else{
			print(paste("WARNING: PDB",pdbid,"has generated empty matrix for AA"))
			#res$AA$bsr = NULL
			res$AA$bsr = Matrix()
		}
	}else{
			print(paste("WARNING: PDB",pdbid,"has generated empty matrix for ALL, ergo for PP and AA"))
			res$PP$bsr = m
			#res$PP$filter = c()
			res$AA$bsr = m
			#res$AA$filter = c()
	}
	return(res)

}

adjust_repeated_ligands = function(pdb,d=100){

	atom = pdb$atom#;print(dim(atom))
	ligands = pdb$ligands
	#print(head(atom))#;readline()
	#print(head(ligands))#;readline()
	f = atom$type == "HETATM"
	#print(dim(atom[f,]))
	#print(head(atom[f,]));readline()
	#print(dim(ligands));print(dim(atom[f,]));readline()
	het.cr = atom[f,c("chain","resno")]
	het.cr1 = unique(het.cr)#;print(het.cr1);readline()#;print(atom[f,c("elety","resid","chain","resno")]);readline() 
	het.cr2 = het.cr1
	#print(head(het.cr1))#;readline()
	tam = dim(het.cr1)[1]
	het.cr.f = factor(het.cr1[,2])
	het.cr.fl = levels(het.cr.f)#;print(het.cr.f)#;readline()
	
	z = het.cr1
	y = het.cr2
	for (i in het.cr.fl){
		yi = which(y[,2]==i)#;print(yi);readline()
		tamj = length(yi)
		if (tamj>1){
			for (j in 2:tamj){
				dy = y[yi[j],2]+d#;print(dy);readline()
				repeat{
					if (dy %in% y[,2]){
						dy = dy+d
					}else{
						y[yi[j],2] = dy
						break;
					}
				}
				#print(xmax);print(dy);print(y);readline()
			}
		}
		#readline()
	}
	#print(head(z));print(tail(z));print(head(y));print(tail(y));readline()
	#print(dim(z));print(dim(y));readline()
	#print(z);print(y);readline()
	for (i in 1:tam){
		fi = (atom[,"resno"] == z[i,2])&(atom[,"chain"] == z[i,1])
		atom[fi,"resno"] = y[i,2]
	}
	#print();readline()
	#;print();readline()

	tamatom = dim(atom[f,])[1]
	tamlig = dim(ligands)[1]

	if (tamatom==tamlig){
		ligands$resid = atom[f,"resno"]#;print(ligands);readline()
	}else{
		print(paste("WARNING: atom and ligand tables do not agree for",pdb$pdbname,";ligand table wont be adjusted"))
	}
	pdb$atom = atom
	pdb$ligands = ligands
	#print(pdb$atom);readline()
	#print(0);readline()
	return(pdb)
}

adapt_pdb_chain_to_ligands =function(pdb,water="HOH",WT="WT",base=9){

	if (!is.null(pdb$ligands)){
		pdb = adjust_repeated_ligands(pdb)#;print("o")

		f = pdb$ligands$resname != water
		lig_vector = c()

		if (sum(f)>0){
	
			ligands = pdb$ligands[f,]#;print(ligands);readline()
			#lig_name_id = ligands[,c("resname","resid")]#;print(lig_name_id)
			lig_name_id = ligands[,c("resname","chainid","resid")]
			lig = unique(lig_name_id)#;print(lig);print(head(pdb$atom));readline()
			tam = dim(lig)[1]
			#if (tam>base){
			#	print(paste("WARNING: number of ligands is greater than 9"))
			#}
			for (i in 1:tam){
				#f = (pdb$atom$resid == lig[i,1])&(pdb$atom$resno == lig[i,2])
				fi = (pdb$atom$resid == lig[i,1])&(pdb$atom$chain == lig[i,2])&(pdb$atom$resno == lig[i,3]) ### mudou f para fi
				#lig_vector = c(lig_vector,paste0("L",i))
				#new_chain_name = as.character(i%%base)
				new_chain_name = paste0("L",i)
				lig_vector = c(lig_vector,new_chain_name)
				#pdb$atom$chain[f] = paste0("L",i)
				pdb$atom$chain[fi] = new_chain_name
				fi = (pdb$ligands$resname == lig[i,1])&(pdb$ligands$chainid == lig[i,2])&(pdb$ligands$resid == lig[i,3])
				#pdb$ligands$chainid[f] = paste0("L",i)
				pdb$ligands$chainid[fi] = new_chain_name
			}
			#pdb$atom = adjust_repeated_ligands(pdb$atom)
			#print(tail(pdb$atom))
			#maxres = max(pdb$atom$resno);print(maxres)#;readline()
			#hetab = filter(pdb$atom,type=="HETATM");print(head(hetab));readline()

		}else{
			print(paste("WARNING: no ligands found in PDB",pdb$pdbname))
		}
		weye = F
		if (sum(!f)>0){
			weye = T
			fw = pdb$atom$resid == water
			pdb$atom$chain[fw] = WT
			fw = pdb$ligands$chainid == water
			pdb$ligands$chainid[fw] = WT
		}else{
			print(paste("WARNING: no waters found in PDB",pdb$pdbname))
		}
		pdb$chains = unique(pdb$atom$chain)
		names(pdb$chains) = rep("C",length(pdb$chains))
		if (length(lig_vector)){
			fg = pdb$chains %in% lig_vector
			names(pdb$chains)[fg] = rep("1",length(lig_vector))
		}
		if (weye){
			fw = pdb$chains %in% WT
			names(pdb$chains)[fw] = "0"
		}
	}else{
		print(paste("WARNING: no hetero atoms found in PDB",pdb$pdbname))
	}
	#print(head(pdb$atom));readline()
	#print(pdb$atom$resno);readline()
	#print(pdb$chains);readline()
	return(pdb)

}

adapt_pdb_chain_to_output = function(pdb,water="HOH",k=3){

	f=which(pdb$atom$type == "HETATM")
	if (length(f)){
		pdb$atom$chain[f] = 1
		f=which(pdb$atom$resid == water)
		if (length(f)){
			pdb$atom$chain[f] = 0
		}else{
			print(paste("WARNING: there is no",water,"to write in adapted PDB",pdb$pdbname))
		}
	}else{
		print(paste("WARNING: there is no HETATM to write in adapted PDB",pdb$pdbname))
	}
	#m = pdb$atom[,c("x","y","z")]
	#rot = get_rot_matrix(source=m,target=diag_x(3),svd=c(T,F),center=c(T,F))
	#m = round(transform_by_rot_matrix(m,rot),k)
	#pdb$atom[,c("x","y","z")] = m

	return(pdb)

}


write_pdb = function(pdb,pdbpath,force=F,sufix=".pdb"){
	
	res = c()
	filename = paste0(pdbpath,pdb$pdbname,2,sufix)
	if (!file.exists(filename)){	
		write.pdb(pdb,file=filename)
	} else{
		if (force){
			print(paste("WARNING:", filename,"for PDB",pdb$pdbname,"already exists, but overwrinting is set"))
			write.pdb(pdb,file=filename)
		}else{
			print(paste("WARNING:", filename,"for PDB",pdb$pdbname,"already exists"))
		}
	}
	return(res)
}

dist_version_for_ligands = function(dist,chains){

	chain_filter = chains[names(chains)=="C"]#;print(chain_filter)
	atom_names = colnames(dist)
	
	tam = length(chain_filter)
	f=c()
	if (tam>1){
		for (i in 1:tam){
			f = c(f,which(grepl(paste0("^",chain_filter[i],"\\."),atom_names)))#;print(f);readline()
		}
		#print(dim(dist))
		dist[f,f] = 0
		#print(dim(dist));readline()
	}else{
		print(paste("WARNING: this PDB seems to have only one chain"))
	
	}
	return(dist)
	
}

rectify_colnames = function(dist,chains,lig="^L[0-9]+\\.",water="^WT\\."){

	atom_names = colnames(dist)
	atom_names = gsub(lig,"1.",atom_names)
	atom_names = gsub(water,"0.",atom_names)
	
	colnames(dist) = atom_names

	return(dist)
	
}

add_polarity_in_colnames = function(dist,meso,aacode,sep="-"){


	atom_names = colnames(dist)
	#atom_codes = new_codify_polarity6(atom_names,meso,aacode)
	atom_codes = new_codify_polarity7(atom_names,meso,aacode)
	atom_names = paste0(atom_names,sep,atom_codes)#;print(atom_names);readline()
	colnames(dist) = atom_names
	return(dist)

}

consider_hydrogens = function(dist,def=100){

	atom_names = colnames(dist)#;print(atom_names)
	id = which(grepl("h$",atom_names))#;print(id)

	if (length(id)){
		print(paste("WARNING: hydrogens were detected. Removing them"))
		ids = which(dist>0,arr.ind=T)#;print(ids)#;print(head(ids));readline()
		f = ids[,1] %in% id#; print(f)#print(length(f));print(dim(dist));readline()
		ids = tomatrix(ids[f,])#; print(ids);readline()
		#return(ids)
		dist[tomatrix(ids[,c("row","col")])]=def#;print(ids[,c("row","col")])
		dist[tomatrix(ids[,c("col","row")])]=def#;print(ids[,c("col","row")]);readline()
	}
	return(dist)

}


#mapply(contact_area,v,MoreArgs = list(r1=r1,r2=r2,probe=probe))
initialize_db_list_for = function(pdb_name,pdbpath,meso,aacode,nmax.atom,workpath,prodpath,ipaddr,ncores,rm.h=T,sufix=".pdb"){
  
	#combinations = data.frame()
  #file = paste0('/data/pdb/',pdb_name, '.pdb')
  
  #Lendo os PDBs
  #pdb <- read.pdb(file,maxlines = -1, multi = FALSE, rm.insert = FALSE, rm.alt = TRUE, ATOM.only = FALSE, verbose = TRUE)
	#pdb <- read.pdb(file,maxlines = -1, multi = FALSE, rm.insert = FALSE, rm.alt = TRUE, ATOM.only = TRUE, verbose = TRUE)
	#print(pdb_name)
	pdb = make_it_mass_pdb_reading(pdb_name,pdbpath,nmax.atom=nmax.atom,sufix=sufix)
	#return(pdb)
	if (!is.null(pdb$error)){
		return(pdb)
	}
	#;print("0")
	pdb2 = adapt_pdb_chain_to_ligands(pdb)

	#print(pdb2$atom);readline()

	write_pdb(adapt_pdb_chain_to_output(pdb2),pdbpath,force=T)

	#readline()
	#write.pdb(pdb2,file=paste0(pdbpath,pdb_name,2,sufix))
	#library("Rpdb")
	#write.pdb(pdb2,file=paste0(pdbpath,pdb_name,"2a",sufix))
	#detach(package:Rpdb)

  #calculado as combinações possíveis de cadeias
	combinall = list()
  combinall[[1]] = get_chain_combinations(pdb2)#;print(pdb2$atom$chain);print(combinall[[1]]);readline()

	#combinations = expand_combinations_to_ligands(combinations,pdb)
	
	dist = list()
	dist = calculate_distance_matrix(combinall[[1]],pdb2,ncores,"WT")#;print("0")
	dist[[1]] = add_polarity_in_colnames(dist[[1]],meso,aacode)#;print(head(colnames(dist[[1]])));readline()
	#print(head(dist[[1]]));readline()

	if (rm.h){
		dist[[1]] = consider_hydrogens(dist[[1]])#;print(head(dist[[1]]));readline()
	}

	distall =list()


	info = list()
	infoall = list()
	info = get_pdb_essentials(dist,pdb2,ncores)#;print(pdb_name)
	infoall[[1]] = info[[1]]

	#dist[[1]] = adjust_dist_waters(dist[[1]],info[[1]])

	dist[[1]] = rectify_colnames(dist[[1]],info[[1]]$chains)#;print(colnames(dist[[1]]));readline()

	dist = dist_version_for_ligands(dist[[1]],infoall[[1]]$chains)

	distall[[1]] = Matrix(dist)
	combinations = data.frame(c("ANY","LIG"))	

	#combinations = data.frame(c("ANY","ANY"))

	#if(1){
	
	chain_number = length(pdb$chains)
	combinations = cbind(combinations,c("CHAIN","CHAIN"))
	if (chain_number>1){ 
		combinall[[2]] = get_chain_combinations(pdb)#;print(combinall)#;readline()
		dist = calculate_distance_matrix(combinall[[2]],pdb,ncores,"WT")#;print("0")
		dist[[1]] = add_polarity_in_colnames(dist[[1]],meso,aacode)
		if (rm.h){
			dist[[1]] = consider_hydrogens(dist[[1]])#;print(head(dist[[1]]));readline()
		}
		info = get_pdb_essentials(dist,pdb,ncores)
		#dist[[1]] = adjust_dist_waters(dist[[1]],info[[1]])
		distall[[2]] = Matrix(dist[[1]])
		infoall[[2]] = info[[1]]
	}else{
		print(paste("WARNING: no chain-chain found to",pdb$pdbname))
		combinall[[2]] = data.frame()
		distall[[2]] = Matrix()
		infoall[[2]] = data.frame()
	}
	#}
	#combinations = cbind(combinations,c("ANY","LIG"))
	
	#dist = dist_version_for_ligands(distall[[1]],infoall[[1]]$chains)
	#combinall[[3]] = combinall[[1]]
	#distall[[1]] = Matrix(distall[[1]])
	#distall[[3]] = Matrix(dist)
	#infoall[[3]] = infoall[[1]]

	colnames(combinations)=NULL
  dist_pdb <- list(
    pdbid = pdb_name,
		title = pdb$title,
		combinall = combinall,
    combinations = combinations,
		#combinot = c(),
    dist = distall,
		info = infoall,
		ligands = pdb$ligands,
		polarity = list(
			#ALL = list(
			#	bsr = list()
			#)			
		)
  )
  
	#workpath=prodpath ### FOR DEBUGGING
	if (workpath==prodpath){
		res = copy_to_server(pdb_name,workpath,prodpath,ipaddr,pdbpath)
		if (length(res)){
			pdb=list()
			pdb$error = -1
			pdb$pdbname = pdb_name
			return(pdb)
		}
		pdb_name = paste0(pdb_name,2)
		res = copy_to_server(pdb_name,workpath,prodpath,ipaddr,pdbpath)
		if (length(res)){
			pdb=list()
			pdb$error = -1
			pdb$pdbname = pdb_name
			return(pdb)
		}
	}

  return(dist_pdb)
}

initialize_db_list = function(PDBS,pdbpath,meso,aacode,nmax.atom,workpath,prodpath,ipaddr,vcores=c(0,0),sufix=".pdb"){

	tam = length(PDBS)
	db = list()
	aux = list()
	if (vcores[1]){
		print(paste("Doing it parallel with",vcores[1],"cores"))
		registerDoMC(vcores[1])
		db = foreach(i=1:tam) %dopar% {
			initialize_db_list_for(PDBS[i],pdbpath,meso,aacode,nmax.atom,workpath,prodpath,ipaddr,vcores[2],sufix=sufix)
		}
	}else{
		print(paste("Doing it sequential..."))
		for(i in 1:tam){
			db[[i]] = initialize_db_list_for(PDBS[i],pdbpath,meso,aacode,nmax.atom,workpath,prodpath,ipaddr,vcores[2],sufix=sufix)	
		}
	}

	#print(db);readline()

	db = validate_db(db)

	return(db)

}

validate_db = function(db){

	tam = length(db)

	for (i in 1:tam){
		#print(db[[i]])
		if (!is.null(db[[i]]$error)){
			return(db[i])
		}
	}
	return(db)

}

add_bsr_matrix_for = function(db,vdw,probe,ncores){

	tam = length(db$dist)

	#db$polarity[[1]]$bsr = list()

	if ((ncores)&(tam>1)){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		db$polarity = foreach(i=1:tam) %dopar% {
			silveira_romanelli_transformation(db$dist[[i]],vdw,probe)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			db$polarity[[i]] = silveira_romanelli_transformation(db$dist[[i]],vdw,probe)
		}
	}

	#db = adjust_bsr_waters(db)
	#elimina combinacoes sem interface
	#db = remove_empty_bsr_matrix(db)
	
	#db$polarity[[2]] = only_water_on_interface(db$polarity[[2]])

	return(db)

}

add_bsr_matrix = function(db,vdw,probe,vcores=c(0,0)){

	tam = length(db)
	if (vcores[1]){
		print(paste("Doing it parallel with",vcores[1],"cores"))
		registerDoMC(vcores[1])
		db = foreach(i=1:tam) %dopar% {
			add_bsr_matrix_for(db[[i]],vdw,probe,vcores[2])
		}
	}else{
		print(paste("Doing it sequential..."))
		for(i in 1:tam){
		#for(i in 4){
			db[[i]] = add_bsr_matrix_for(db[[i]],vdw,probe,vcores[2])	
		}
	}
	return(db)

}

expand_bsr_matrix_for = function(db,ncores){

	tam = length(db$polarity)

	#db$polarity[[1]]$bsr = list()

	if ((ncores)&(tam>1)){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		db$polarity = foreach(i=1:tam) %dopar% {
			add_polarity_for(db$polarity[[i]],db$pdbid)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){;print(i)
			db$polarity[[i]] = add_polarity_for(db$polarity[[i]],db$pdbid)
		}
	}
	return(db)

}

expand_bsr_matrix = function(db,vcores=c(0,0)){

	tam = length(db)#;tam=1
	if (vcores[1]){
		print(paste("Doing it parallel with",vcores[1],"cores"))
		registerDoMC(vcores[1])
		db = foreach(i=1:tam) %dopar% {
			expand_bsr_matrix_for(db[[i]],vcores[2])
		}
	}else{
		print(paste("Doing it sequential..."))
		for(i in 1:tam){#;print(i)
		#for(i in 4){
			db[[i]] = expand_bsr_matrix_for(db[[i]],vcores[2])	
		}
	}
	return(db)

}


build_dba_from_db = function(db,vdw,sep="_",vcores=c(0,0),type="BSR",pol=c("ALL","PP","AA"),matrixtype="normal"){

	tami = length(db)
	dba = list()
	n=1
	print(paste("Doing it sequential..."))
	for (i in 1:tami){
		#print(paste("i:",i))
		tamj = length(db[[i]]$polarity)
		combina = db[[i]]$combinations
		for (j in 1:tamj){
			#print(paste("j:",j))
			tamk = length(db[[i]]$polarity[[j]])
			dba[[n]]=list()
			data = list()
			for (k in 1:tamk){
				#print(paste("k:",k))
				data[[k]] = list()
				combina1 = as.vector(unlist(db[[i]]$combinations[1,j]))
				combina2 = as.vector(unlist(db[[i]]$combinations[2,j]))
				data[[k]]$matrixname = paste0(db[[i]]$pdbid,sep,combina1,sep,combina2,sep,type,sep,pol[k])
				data[[k]]$matrixtype = matrixtype#
				#data[[k]]$pdbname = db[[i]]$pdbid
				data[[k]]$pdbname = paste0(db[[i]]$pdbid,sep,combina1,sep,combina2)
				data[[k]]$title = db[[i]]$title
				if (dim(db[[i]]$polarity[[j]][[k]]$bsr)[1]>1){
					data[[k]]$A = list()
					#filter = db[[i]]$polarity[[j]][[k]]$filter
					m = db[[i]]$polarity[[j]][[k]]$bsr
					atom_frame = sep_atom_name_id(colnames(m))#;print(head(atom_frame));readline()
					data[[k]]$element_name = atom_frame[,1]
					data[[k]]$element_polarity = atom_frame[,5]
					#data[[k]]$atom_ids = db[[i]]$info[[j]]$atom_ids[filter]
					data[[k]]$A$a = m
					colnames(data[[k]]$A$a) = NULL
					rownames(data[[k]]$A$a) = NULL
					data[[k]]$exp = list()
					filter = db[[i]]$info[[j]]$eleno %in% atom_frame[,2]
					data[[k]]$pdb = list()
					data[[k]]$pdb$xyz = db[[i]]$info[[j]]$xyz[filter,]
					chain_types = names(db[[i]]$info[[j]]$chains)#;print(chain_types)#;readline()
					f = which(chain_types =="C")
					chain_list = db[[i]]$info[[j]]$chains[f]
					if ("1" %in% chain_types) chain_list = c(chain_list,"1")#;print(chain_list)
					if ("0" %in% chain_types) chain_list = c(chain_list,"0")#;print(chain_list);readline()
					data[[k]]$chains = chain_list
					names(data[[k]]$chains) = NULL
					data[[k]]$vdw = get_vdw(atom_name=atom_frame[,1],atom_list=atom_frame[,3],vdw=vdw)
					#print(data[[k]]$vdw);readline()

#get_vdw = function(atom_name,atom_list,vdw,df=vdw["MEAN"],nf="X"){

				}else{
					print(paste("WARNING: no BSR AA or PP matrices were found for",db[[i]]$pdbid))
				}
				#;print(data);readline()
			}
			dba[[n]] = data
			n=n+1
		}
	}
	#return(dba)
	all = lapply(dba,"[[",1)
	pp = lapply(dba,"[[",2)
	aa = lapply(dba,"[[",3)

	res = list()
	res[[1]] = list()
	res[[1]] = all
	res[[2]] = list()
	res[[2]] = pp
	res[[3]] = list()
	res[[3]] = aa

	#return(dba)
	return(res)	

}

build_dbsup_from_dbgrp = function(db,trans=F,sep="_"){

	dbsup = list()

	tami = length(db)

	for (i in 1:tami){
	#for (i in 16){	
		dbsup[[i]] = list()
		tamj = length(db[[i]])
		for (j in 1:tamj){
			pdbname = paste0(c(db[[i]][[j]]$pdb,db[[i]][[j]]$combination_chains),collapse=sep)#;print(pdbname);readline()
			print(paste("Converting",pdbname))#;print(i);print(j)
			dbsup[[i]][[j]] = list()
			dbsup[[i]][[j]]$pdb = pdbname
			dbsup[[i]][[j]]$clus = list()
			dbsup[[i]][[j]]$clus[[1]] = list()
			#if(0){
			if (length(db[[i]][[j]]$connections_group_bsr)){
				conn = db[[i]][[j]]$connections_group_bsr
				local = db[[i]][[j]]$local_score_bsr
				global = db[[i]][[j]]$global_score_bsr
				spot = db[[i]][[j]]$spot_score_bsr
				dense = db[[i]][[j]]$dense_bsr_total
				center = db[[i]][[j]]$center_bsr
				#if (trans) rot = db[[i]][[j]]$rotall
				#else rot = NULL
				rot = db[[i]][[j]]$rotall
				dbsup[[i]][[j]]$clus[[1]] = list()
				dbsup[[i]][[j]]$clus[[1]]$rotself = rot
				res = get_info_connections_group(conn,local,global,spot,dense,center)
				dbsup[[i]][[j]]$clus[[1]]$nclus = res
			}#;print("o")
			dbsup[[i]][[j]]$clus[[2]] = list()
			if (length(db[[i]][[j]]$connections_group_bsr_polar_polar)){
				conn = db[[i]][[j]]$connections_group_bsr_polar_polar
				local = db[[i]][[j]]$local_score_bsr_polar_polar
				global = db[[i]][[j]]$global_score_bsr_polar_polar
				spot = db[[i]][[j]]$spot_score_bsr_polar_polar
				dense = db[[i]][[j]]$dense_bsr_total_polar_polar
				center = db[[i]][[j]]$center_bsr_polar_polar
				#if (trans) rot = db[[i]][[j]]$rotpp
				#else rot = NULL
				rot = db[[i]][[j]]$rotpp
				dbsup[[i]][[j]]$clus[[2]] = list()
				dbsup[[i]][[j]]$clus[[2]]$rotself = rot
				res = get_info_connections_group(conn,local,global,spot,dense,center)
				dbsup[[i]][[j]]$clus[[2]]$nclus = res
			}#;print("o")
			#}
			dbsup[[i]][[j]]$clus[[3]] = list()
			if (length(db[[i]][[j]]$connections_group_bsr_nonpolar_nonpolar)){
				conn = db[[i]][[j]]$connections_group_bsr_nonpolar_nonpolar#;print(head(conn))
				local = db[[i]][[j]]$local_score_bsr_nonpolar_nonpolar#;print(local)
				global = db[[i]][[j]]$global_score_bsr_nonpolar_nonpolar#;print(global)
				spot = db[[i]][[j]]$spot_score_bsr_nonpolar_nonpolar#;print(spot);readline()
				dense = db[[i]][[j]]$dense_bsr_total_nonpolar_nonpolar
				center = db[[i]][[j]]$center_bsr_nonpolar_nonpolar
				#if (trans) rot = db[[i]][[j]]$rotaa
				#else rot = NULL
				rot = db[[i]][[j]]$rotaa
				dbsup[[i]][[j]]$clus[[3]] = list()
				dbsup[[i]][[j]]$clus[[3]]$rotself = rot
				res = get_info_connections_group(conn,local,global,spot,dense,center)
				dbsup[[i]][[j]]$clus[[3]]$nclus = res
			}#;print("o")
		}
	}
	return(dbsup)
}

clean_database = function(collection,workpath,prodpath,addrpath,localpath,prefix="mongodb://"){

	print("W A R N I N G: CLEANING DATA BASE PROCESS !!!")
	if (workpath==prodpath){
		print(paste("PATH:",workpath,"ADDRESS:",addrpath))
	}else{
		print(paste("PATH:",workpath,"ADDRESS:",localpath))
	}
	print("Are you sure it is to clean mongo? (y/n)")
	res = readline()
	if (res=="y"){
		print("Are you REALLY sure??? (y/n)")
		res = readline()
		if (res=="y"){
			### USAR getwd ###
			#library(mongolite)
			if (workpath==prodpath){
				urlpath = paste0(prefix,addrpath)
				print(paste("Cleaning mongodb at",urlpath))
				#db_collection = mongo(url= urlpath, collection = collection, db = "ppi") 
				#db_collection$remove('{}')
			}else{
				urlpath = paste0(prefix,localpath)
				print(paste("Cleaning mongodb at",urlpath))
				db_collection = mongo(url= urlpath, collection = collection, db = "ppi")
				db_collection$remove('{}') 
			}
		}
	}
}

toMongo = function(db,collection,workpath,prodpath,addrpath,localpath,prefix="mongodb://"){
  
  library(mongolite)

	if (workpath==prodpath){
			urlpath = paste0(prefix,addrpath)  	
	}else{
			urlpath = paste0(prefix,localpath)
	}

	db_collection = mongo(url= urlpath, collection = collection, db = "ppi") # create connection, database and collection
	#db_collection = mongo(url= "mongodb://localhost:27017", collection = collection, db = "ppi") # create connection, database and collection  
	#db_collection = mongo(url= "mongodb://172.16.86.2:27017", collection = collection, db = "ppi") # create connection, database and collection

  if(collection == "proteins"){
		print(paste("Saving collection",collection,"to mongo at:",urlpath))
    # para cada pdb
    for(i in 1: length(db)){
      
      # verifica se já existe o PDB. Isso tb é verificado na UI antes de chegar no R
      query = paste0('{"pdb.0" : "', toupper(db[[i]][[1]]$pdb), '" } ')
      proteins_collection <- db_collection$find(query)
      
      #pode gravar pois não existe nada no banco e não vai gravar elementos duplicados
      if(length(proteins_collection)==0){
        # chain chain
        if(length(db[[i]][[2]]$bsr_groups_nonpolar_nonpolar) == 0){
          db_collection$insert(db[[i]][[1]])
          db_collection$insert(db[[i]][[2]])
        }else{
          db_collection$insert(db[[i]][[2]])
          db_collection$insert(db[[i]][[1]])
        }
      }
    }  
  }else if(collection == "groups"){
		print(paste("Saving collection",collection,"to mongo at:",urlpath))
    for(i in 1: length(db)){
      for(j in 1: length(db[[i]])){
        # verifica se já existe o PDB. Isso tb é verificado na UI antes de chegar no R
        query = paste0('{"pdb.0" : "', toupper(db[[i]][[j]]$pdb), '", "combination_chains.0" : "',dbgrp[[i]][[j]]$combination_chains[1],'", "combination_chains.1" : "',dbgrp[[i]][[j]]$combination_chains[2],'" } ')
        groups_collection <- db_collection$find(query)
        
        #pode gravar pois não existe nada no banco e não vai gravar elementos duplicados
        if(length(groups_collection)==0){
          db_collection$insert(db[[i]][[j]])
        }
      }
      
    }
  }else if(collection == "alignments"){
		print(paste("Saving collection",collection,"to mongo at:",urlpath))
    # verifica se já existe o PDB. Isso tb é verificado na UI antes de chegar no R
    query = paste0('{"pdb1.0" : "', db$pdb1, '", "pdb2.0" : "',db$pdb2,'", "polarity" : "',db$polarity,'" } ')
    alignments_collection <- db_collection$find(query)
    
    #pode gravar pois não existe nada no banco e não vai gravar elementos duplicados
    if(length(alignments_collection)==0){
      db_collection$insert(db)
    }
  }
}

fromMongo = function(collection,pdbs,workpath,prodpath,addrpath,localpath,prefix="mongodb://"){
  
  db_grp <- list()
  
	if (workpath==prodpath){
			urlpath = paste0(prefix,addrpath)  	
	}else{
			urlpath = paste0(prefix,localpath)
	}

  #library(mongolite)
	print(paste("Open mongo connection for",collection,"at:",urlpath))
  db_collection = mongo(url= urlpath, collection = collection, db = "ppi") # create connection, database and collection
  for(k in 1 : length(pdbs)){
		print(paste("Loading PDB:",pdbs[k]))
    query = paste0('{"pdb.0":"',pdbs[k],'"}')
    group <- db_collection$find(query)
    grp <- list()
    db_group <- list()
		tamgroup = length(group$pdb)#;print(tamgroup)
		if (tamgroup){
		  for(i in 1 : tamgroup){
		  #for(i in 1 : length(group$pdb)){
		    grp$pdb = group$pdb[[i]]
		    grp$pdb_type = group$pdb_type[[i]]
		    grp$chains = group$chains[[i]]
		    grp$combination_chains = group$combination_chains[[i]]
		    grp$rotall = group$rotall[[i]]
		    grp$rotpp = group$rotpp[[i]]
		    grp$rotaa = group$rotaa[[i]]

		    grp$connections_group_bsr = group$connections_group_bsr[[i]]
		    grp$connections_group_bsr_polar_polar = group$connections_group_bsr_polar_polar[[i]]
		    grp$connections_group_bsr_nonpolar_nonpolar = group$connections_group_bsr_nonpolar_nonpolar[[i]]

		    grp$volume_bsr = group$volume_bsr[[i]]
		    grp$volume_bsr_polar_polar = group$volume_bsr_polar_polar[[i]]
		    grp$volume_bsr_nonpolar_nonpolar = group$volume_bsr_nonpolar_nonpolar[[i]]

		    grp$spot_score_bsr = group$spot_score_bsr[[i]]
		    grp$spot_score_bsr_polar_polar = group$spot_score_bsr_polar_polar[[i]]
		    grp$spot_score_bsr_nonpolar_nonpolar = group$spot_score_bsr_nonpolar_nonpolar[[i]]

		    grp$local_score_bsr = group$local_score_bsr[[i]]
		    grp$local_score_bsr_polar_polar = group$local_score_bsr_polar_polar[[i]]
		    grp$local_score_bsr_nonpolar_nonpolar = group$local_score_bsr_nonpolar_nonpolar[[i]]

		    grp$global_score_bsr = group$global_score_bsr[[i]]
		    grp$global_score_bsr_polar_polar = group$global_score_bsr_polar_polar[[i]]
		    grp$global_score_bsr_nonpolar_nonpolar = group$global_score_bsr_nonpolar_nonpolar[[i]]

		    grp$dense_bsr_total = group$dense_bsr_total[[i]]
		    grp$dense_bsr_total_polar_polar = group$dense_bsr_total_polar_polar[[i]]
		    grp$dense_bsr_total_nonpolar_nonpolar = group$dense_bsr_total_nonpolar_nonpolar[[i]]

		    grp$dense_bsr_bipar = group$dense_bsr_bipar[[i]]
		    grp$dense_bsr_bipar_polar_polar = group$dense_bsr_bipar_polar_polar[[i]]
		    grp$dense_bsr_bipar_nonpolar_nonpolar = group$dense_bsr_bipar_nonpolar_nonpolar[[i]]

				grp$center_bsr = group$center_bsr[[i]]
		    grp$center_bsr_polar_polar = group$center_bsr_polar_polar[[i]]
		    grp$center_bsr_nonpolar_nonpolar = group$center_bsr_nonpolar_nonpolar[[i]]
		   
		    db_group[[i]] = grp
		  }
		}else{
			print(paste("WARNING: it was not possible to reload",pdbs[k],"from Mongo!"))
		}
    db_grp[[k]] = db_group
  }
  return(db_grp)
  
}


logger = function(pdb, type,key, percentage, wp,pp,ap,lp, ...){

	prefix="mongodb://"
	if (wp==pp){
			urlpath = paste0(prefix,ap)  	
	}else{
			urlpath = paste0(prefix,lp)
	}
  
	#print(urlpath);readline()

	description = c("IMPORT","ALIGNMENT")

  log <- list()
  log$pdb = pdb
  log$key = key
  log$type = description[type]
  log$action = paste(...)
  log$percentage = percentage
  log$time = Sys.time()
  
  library(mongolite)
  db_collection = mongo(url= urlpath, collection = "logs", db = "ppi")
  db_collection$insert(log)
  
}

verify_file_in_server = function(pdbid,ipaddr,pdbpath,command="scp",sufix=".pdb"){
	
	pdbpath2 = pdbpath
	pdbpath2 = '/data/pdb/'
	pdbpath2 = paste0(pdbpath2,pdbid,sufix)
	args = c('. &> /dev/null')
	args = c(paste0(ipaddr,":",pdbpath2),args)
	command_path=paste(command,args[1],args[2])
	print(command_path)#;readline()
	res = system2(command,args=args,stdout=T,stderr=T)
	print(res);readline()
	#scp 172.16.86.2:/data/pdb/1PPFa.dpb ./ >& /dev/null


}

#rsync --chmod=u+rwx,g+rwx,o+rwx /path/to/file server:/path/to/file

copy_from_server = function(pdbid,workpath,prodpath,ipaddr,pdbpath,command="rsync",sufix=".pdb",user="pingauser"){

	#verify_file_in_server(pdbid,ipaddr,pdbpath)

	#pdbid="1PPFa" ### FOR DEBUGGING ###
	res=c()
	if (command=="rsync"){
		#prodpath=workpath ### FOR DEBUGGING ###
		if (workpath==prodpath){
				pdbpath2 = pdbpath
				#pdbpath2 = "/data/pdb/" ### FOR DEBUGGING ###
				args = '--chmod=644'
				args = c(args,paste0(user,"@",ipaddr,":",pdbpath2,pdbid,sufix),paste0(pdbpath))
				command_path = paste(c(command,args),collapse=" ")
				print(paste("Copying file:",command_path))
				res = system2(command,args=args,stdout=T,stderr=T)
				if (length(res)){
					print(paste("WARNING: something wrong in copying file:",command_path))
					print(res)#;readline()
				}
		}

	}
	return(res)
	#scp data/pdb/1ACB2.pdb 172.16.86.2:/data/pdb
	#scp remote-host:~/myfile ./ >& /dev/null
	#scp 172.16.86.2:/data/pdb/1PPF.dpb ./ >& /dev/null

}


copy_to_server = function(pdbid,workpath,prodpath,ipaddr,pdbpath,command="rsync",sufix=".pdb",user="pingauser"){

	#verify_file_in_server(pdbid,ipaddr,pdbpath)

	#pdbid="1PPFa" ### FOR DEBUGGING ###
	res=c()
	if (command=="scp"){
		#prodpath=workpath ### FOR DEBUGGING ###
		if (workpath==prodpath){
				pdbpath2 = pdbpath 
				#pdbpath2 = "/data/pdb/" ### FOR DEBUGGING ###
				args = c(paste0(pdbpath,pdbid,sufix),paste0(user,"@",ipaddr,":",pdbpath2))
				command_path = paste(command,args[1],args[2])
				print(paste("Copying file:",command_path))
				res = system2(command,args=args,stdout=T,stderr=T)
				if (length(res)){
					print(paste("WARNING: something wrong in copying file:",command_path))
					print(res)
				}
		}
		return(length(res))
	}
	if (command=="rsync"){
		#prodpath=workpath ### FOR DEBUGGING ###
		if (workpath==prodpath){
				pdbpath2 = pdbpath
				#pdbpath2 = "/data/pdb/" ### FOR DEBUGGING ###
				args = '--chmod=644'
				args = c(args,paste0(pdbpath,pdbid,sufix),paste0(user,"@",ipaddr,":",pdbpath2))
				command_path = paste(c(command,args),collapse=" ")
				print(paste("Copying file:",command_path))
				res = system2(command,args=args,stdout=T,stderr=T)
				if (length(res)){
					print(paste("WARNING: something wrong in copying file:",command_path))
					print(res)#;readline()
				}
		}

	}
	return(res)
	#scp data/pdb/1ACB2.pdb 172.16.86.2:/data/pdb
	#scp remote-host:~/myfile ./ >& /dev/null
	#scp 172.16.86.2:/data/pdb/1PPF.dpb ./ >& /dev/null

}

###################### OLDS #######################


logger_old2 = function(pdb, key, type, percentage, ...){
  description = c("IMPORT","ALIGNMENT")

  action = paste(...)
  
  if(key != ''){
    key = pdb
  }
  
  log <- list()
  log$pdb = pdb
  log$key = key
  log$type = description[type]
  log$action = action
  log$percentage = percentage
  log$time = Sys.time()
  
  library(mongolite)
  db_collection = mongo(url= "mongodb://localhost:27017", collection = "logs", db = "ppi")
  db_collection$insert(log)
  
}


logger_old = function(pdb, key=pdb, type, percentage, ...){
  description = c("IMPORT","ALIGNMENT")
  action = paste(...)
  
  log <- list()
  log$pdb = pdb
  log$key = key
  log$type = description[type]
  log$action = action
  log$percentage = percentage
  log$time = Sys.time()
  
  library(mongolite)
  db_collection = mongo(url= "mongodb://localhost:27017", collection = "logs", db = "ppi")
  db_collection$insert(log)
  
}

prepare_atom_matrix_olds = function(pdb, chain_a_name, chain_b_name, hetatm=F,join="\\+"){
  
	chain_a_name = (as.vector(unlist(chain_a_name)))
	chain_b_name = (as.vector(unlist(chain_b_name)))
	if (grepl(join,chain_a_name)){#print(chain_a_name)
		chain_join = unlist(strsplit(chain_a_name,join))
		chain_a_1 = chain_join[1]
		chain_a_2 = chain_join[3]#;print(chain_a_1)
		chain_a = filter(pdb$atom, chain == chain_a_1 | chain == chain_a_2)
	}else{
		#chain_a = filter(pdb$atom, chain == chain_a_name, type == "ATOM")
		chain_a = filter(pdb$atom, chain == chain_a_name)
	}
	if (grepl(join,chain_b_name)){
		chain_join = unlist(strsplit(chain_b_name,join))
		chain_b_1 = chain_join[1]
		chain_b_2 = chain_join[3]
		chain_b = filter(pdb$atom, chain == chain_b_1 | chain == chain_b_2)
	}else{
		#chain_b = filter(pdb$atom, chain == chain_b_name, type == "ATOM")
		chain_b = filter(pdb$atom, chain == chain_b_name)#;print("o")
	} 
  #chain_a = filter(pdb$atom, chain == chain_a_name, type == "ATOM")
  #chain_b = filter(pdb$atom, chain == chain_b_name, type == "ATOM")
	#chain_all = filter(pdb$atom, type == "ATOM")
	chain_all = pdb
  
	if(hetatm){
		chain_all = filter(pdb$atom, type == "ATOM")
		chain_a_hetatm = filter(pdb$atom, chain == chain_a_name, type == "HETATM")
		chain_b_hetatm = filter(pdb$atom, chain == chain_b_name, type == "HETATM")
		
		chain_a = rbind(chain_a,chain_a_hetatm)
		chain_b = rbind(chain_b,chain_b_hetatm)
		
		#chain_all = filter(pdb$atom, type == "ATOM")
		chain_all_hetatm = filter(pdb$atom, type == "HETATM")
		
		chain_all = rbind(chain_all, chain_all_hetatm)
  }
	#print(chain_a);readline()
  cx = c(chain_a$x, chain_b$x)
  cy = c(chain_a$y, chain_b$y)
  cz = c(chain_a$z, chain_b$z)
  mat = cbind(cx, cy, cz)
	#print(mat[1:3,]);readline()
  
  #Nível de resíduo
  rownames(mat) <- c(paste0(chain_a_name, '.', chain_a$resid, chain_a$resno, '.' ,chain_a$elety, '-', chain_a$eleno ), paste0(chain_b_name, '.', chain_b$resid, chain_b$resno, '.' ,chain_b$elety, '-', chain_b$eleno ))

	#print(mat[1:3,]);readline()
  
  pdb_list_filter <- list(
    chain_a = chain_a,
    chain_b = chain_b,
    chain_a_name = chain_a_name,
    chain_b_name = chain_b_name,
    xyz_matrix = mat,
    chain_all = chain_all
  )
  
  return(pdb_list_filter)
}


#
# dado um vetor, essa função retorna todas as possibilidades combinadas 2 a 2
#
get_chain_combinations_old = function(pdb,water="WT"){
  
  chains = unique(pdb$atom$chain);print(chains)#;readline()
  
  if(length(chains) < 2)
    return(chains)
  
	comb = combn(chains, 2)

	wt.ids = which(comb[2,] == water)
	new = paste(comb[1,wt.ids],collapse="+")
	new = data.frame(c(new,water))

	comb = as.data.frame(comb)
	#print(comb);print(wt.ids);readline()
	comb[,wt.ids] = NULL
	
	comb = cbind(comb,new)
	colnames(comb)=NULL
	rownames(comb)=NULL
	#print(comb);readline()
	
	return(comb)	

  #return(combn(chains, 2))
}



#new_codify_polarity(data$element_name,meso,aacode)

###################### TESTES #######################


import_to_test_dbgrp = function(PDBS,pdbpath,meso,aacode,vdw,ncores){

	download_pdb(PDBS)
	db = list()
	db = initialize_db_list(PDBS,pdbpath,meso,aacode,vcores=c(ncores,ncores))
	db = add_bsr_matrix(db,vdw,probe,vcores=c(ncores,ncores))
	db = expand_bsr_matrix(db,vcores=c(0,0))
	dba = build_dba_from_db(db,vdw)
	dba = make_hyper_superclus(dba,icores=3,vcores=c(ncores,ncores,ncores,ncores),conex=T,only_a=T)
	dba2 = make_hyper_filter(dba,icores=0,vcores=c(0),cutpar=10)
	dba2 = make_hyper_superclus(dba2,icores=3,vcores=c(0,0,ncores,ncores),conex=F,only_a=F)
	dba2 = make_hyper_reorientation(dba2,icores=3,vcores=c(ncores))
	dbsel = list()
	dbsel = make_hyper_hotspot(dba2,icores=3,vcores=c(0,ncores,ncores),lowcut=c(0.50,0.50),minclusid=3,maxclusid=10)
	dbgrp = list()
	dbgrp = make_hyper_group_output(dba2,dbsel,ncores)
	return(dbgrp)

}


remove_empty_bsr_matrix = function(db){

	tam = length(db$dist)#;print(tam)
	if (tam>1){
		f = rep(F,tam)
		#for (i in 1:tam){
		for (i in tam:1){
			if (dim(db$polarity[[i]])[1]==0){
				print(paste(c("WARNING: PDB:",db$pdbid,"has a combination chain:",db$combinations[,i],"that has no BSR interface"),collapse=" "))
				db$polarity[[i]]=NULL		
				db$info[[i]]=NULL
				f[i]=T
			}
		}
		#print(f);readline()
		if (sum(f)>0){
			#db$polarity[[f]]=NULL
			#db$info[[f]]=NULL
			#print(data.frame(db$combinations[,f]));readline()
			if (sum(f)==1){
				#db$combinot = t(tomatrix(db$combinations[,f]))
				db$combinot = data.frame(db$combinations[,f])
				colnames(db$combinot)=NULL
			}else{
				#db$combinot = tomatrix(db$combinations[,f])
				db$combinot = db$combinations[,f]
			}
			if (sum(!f)==1){
				#db$combinations = t(tomatrix(db$combinations[,!f]))
				db$combinations = data.frame(db$combinations[,!f])
				colnames(db$combinations)=NULL
			}else if (sum(!f)>1){
				#db$combinations = tomatrix(db$combinations[,!f])
				db$combinations = db$combinations[,!f]
			}else{
				print(paste("WARNING: PDB:",db$pdbid,"seems to have no BSR interface. It will be put away from chosen list"))
				return(NULL)
			}
		}
	}else if (dim(db$polarity[[1]])[1]==0){
			print(paste("WARNING: PDB:",db$pdbid,"seems to have no BSR interface. This PDB will be put away from chosen list"))
			return(NULL)
	}
	return(db)
}



adjust_bsr_waters = function(db,target="chain",zero=0.1){

	#print(db$combinations);readline()
	id = which(apply(apply(db$combinations,2,"==",target),2,sum)==2);print(id)

	if (length(id)==1){
		bsr = db$polarity[[id]]
		chains = db$info[[id]]$chains
		tam = length(chains)
		atom_names = colnames(bsr)
		waters = list()
		if (tam>1){
			for (i in 1:tam){
				ids = which(grepl("HOH",atom_names)&(grepl(paste0("^",chains[i],"\\."),atom_names)))
				if (length(ids)){
					waters[[i]] = ids
				}else{
					print(paste("WARNING: it was not found waters to adjust BSR for chain",chains[i],"in",db$pdbid))
					waters[[i]] = c()
				}
			}
			for (i in 1:(tam-1)){
				wa = waters[[i]]
				wb = c()
				for (j in (i+1):tam){
					wb = c(wb,waters[[j]])
				}
				subsr1 = round(bsr[wa,wa],1);print(colnames(subsr1))
				subsr2 = round(bsr[wa,c(wa,wb)],1);print(colnames(subsr2))#;readline()
				print(subsr1);print(subsr2)#;readline()
				subsr1.lines = apply(apply(subsr1,2,">",zero),1,sum);print(subsr1.lines)
				subsr2.lines = apply(apply(subsr2,2,">",zero),1,sum);print(subsr2.lines)
				subsr21.lines = subsr2.lines - subsr1.lines;print(subsr21.lines)
				id.zero = subsr21.lines==0;print(id.zero)
				subsr1[id.zero,]=0
				#subsr1[,id.zero]=0
				;print(subsr1)#;readline()
				#print(colnames(subsr));print(subsr);readline()
			}
		}else{
			print(paste("WARNING: PDB",db$pdbid,"has only one chain",chains[i],". No water BSR adjusts do do"))
		}
		
		#print(dim(m));readline()
	}else if (length(id)==0){
		print(paste("WARNING: PDB:",db$pdbid,"has no chains to adjust waters"))
	}else if (length(id)>1){
		print(paste("WARNING: it is not possible to adjust BSR waters in",db$pdbid))
	}
	return(db)

}



adjust_dist_waters = function(dist,info){

		chains = info$chains
		tam = length(chains)
		#print(chains);print(colnames(dist));readline()
		atom_names = colnames(dist)
		
		if (tam>1){
			for (i in 1:(tam-1)){
				w_a = which(grepl("HOH",atom_names)&(grepl(paste0("^",chains[i],"\\."),atom_names)))#;print(w_a);readline()
				tam_w = length(w_a)
				if (length(tam_w)){
					mat_dist = dist2(info$xyz[w_a,],info$xyz[w_a,])
					tam_w = length(w_a)
					for(j in 1:tam_w){
						for(k in 1:tam_w){
							dist[w_a[j],w_a[k]] = mat_dist[j,k]
							dist[w_a[k],w_a[j]] = mat_dist[j,k]
						}
					}
				}else{
					print(paste("WARNING: adjustments of dist waters couldn't be done..."))
				}
			}
		}
		return(dist)
}



enqueue_dist_matrix = function(all){

	queue = list()

	tami = length(all)
	k = 1
	for (i in 1:tami){
		tamj = length(all[[i]]$dist)
		for (j in 1:tamj){
			queue[[k]] = all[[i]]$dist[[j]]
			k=k+1
		}
	}
	return(queue)

}

dequeue_dist_matrix = function(all,queue){

	tami = length(all)
	k = 1
	for (i in 1:tami){
		tamj = length(all[[i]]$dist)
		all[[i]]$bsr = list()
		for (j in 1:tamj){
			all[[i]]$bsr[[j]] = queue[[k]]
			k=k+1
		}
	}
	return(all)

}

calculate_BSR_matrix2 = function(all,vdw,probe,ncores){

	dist = list()
	dist = enqueue_dist_matrix(all)
	tam = length(dist)
	bsr = list()

	if ((ncores)&(tam>1)){
		registerDoMC(ncores)
		bsr = foreach(i=1:tam) %dopar% {
			silveira_romanelli_transformation(dist[[i]],vdw,probe)
		}
	}else{
		for (i in 1:tam){
			bsr[[i]] = silveira_romanelli_transformation(dist[[i]],vdw,probe)
		}
	}
	all = dequeue_dist_matrix(all,bsr)

}


add_bsr_matrix_for_old = function(db,vdw,probe,ncores){

	tam = length(db$dist)

	#db$polarity[[1]]$bsr = list()

	if ((ncores)&(tam>1)){
		print(paste("Doing it parallel with",ncores,"cores"))
		registerDoMC(ncores)
		db$polarity = foreach(i=1:tam) %dopar% {
			silveira_romanelli_transformation(db$dist[[i]],vdw,probe)
		}
	}else{
		print(paste("Doing it sequential..."))
		for (i in 1:tam){
			db$polarity[[i]] = silveira_romanelli_transformation(db$dist[[i]],vdw,probe)
		}
	}
	#elimina combinacoes sem interface
	#db = remove_empty_bsr_matrix(db)
	
	#db$polarity[[2]] = only_water_on_interface(db$polarity[[2]])

	return(db)

}

toMongo_old = function(db, collection){
  
  #library(mongolite)
	print(paste("Open connection to MongoDB"))
  db_collection = mongo(url= "mongodb://localhost:27017", collection = collection, db = "ppi") # create connection, database and collection
  
  if(collection == "proteins"){
    # para cada pdb
    for(i in 1: length(db)){
      # chain chain
      if(length(db[[i]][[2]]$bsr_groups_nonpolar_nonpolar) == 0){
        db_collection$insert(db[[i]][[1]])
        db_collection$insert(db[[i]][[2]])
      }else{
        db_collection$insert(db[[i]][[2]])
        db_collection$insert(db[[i]][[1]])
      }
      
    }  
  }else{
		
    for(i in 1: length(db)){
			print(paste("Writing",db[[i]][[1]]$pdb))
      for(j in 1: length(db[[i]])){
        db_collection$insert(db[[i]][[j]])
      }
      
    }
  }
}


toMongo_old2 = function(db, collection){
  
  library(mongolite)
  db_collection = mongo(url= "mongodb://localhost:27017", collection = collection, db = "ppi") # create connection, database and collection
  
  if(collection == "proteins"){
      
    # para cada pdb
    for(i in 1: length(db)){
      
      # verifica se já existe o PDB. Isso tb é verificado na UI antes de chegar no R
      query = paste0('{"pdb.0" : "', toupper(db[[i]][[1]]$pdb), '" } ')
      proteins_collection <- db_collection$find(query)
      
      #pode gravar pois não existe nada no banco e não vai gravar elementos duplicados
      if(length(proteins_collection)==0){
        # chain chain
        if(length(db[[i]][[2]]$bsr_groups_nonpolar_nonpolar) == 0){
          db_collection$insert(db[[i]][[1]])
          db_collection$insert(db[[i]][[2]])
        }else{
          db_collection$insert(db[[i]][[2]])
          db_collection$insert(db[[i]][[1]])
        }
      }
    }  
  }else if(collection == "groups"){
    for(i in 1: length(db)){
      for(j in 1: length(db[[i]])){
        # verifica se já existe o PDB. Isso tb é verificado na UI antes de chegar no R
        query = paste0('{"pdb.0" : "', toupper(db[[i]][[j]]$pdb), '", "combination_chains.0" : "',dbgrp[[i]][[j]]$combination_chains[1],'", "combination_chains.1" : "',dbgrp[[i]][[j]]$combination_chains[2],'" } ')
        groups_collection <- db_collection$find(query)
        
        #pode gravar pois não existe nada no banco e não vai gravar elementos duplicados
        if(length(groups_collection)==0){
          db_collection$insert(db[[i]][[j]])
        }
      }
      
    }
  }else if(collection == "alignments"){
    db_collection$insert(db)
  }
}

adapt_pdb_chain_to_ligands_old =function(pdb,water="HOH",WT="WT",base=9){


	f = pdb$ligands$resname != water
	lig_vector = c()

	if (sum(f)>0){
	
		ligands = pdb$ligands[f,];print(ligands);readline()
		lig_name_id = ligands[,c("resname","resid")]#;print(lig_name_id)
		lig = unique(lig_name_id)#;print(lig);readline()
		tam = dim(lig)[1]
		#if (tam>base){
		#	print(paste("WARNING: number of ligands is greater than 9"))
		#}
		for (i in 1:tam){
			f = (pdb$atom$resid == lig[i,1])&(pdb$atom$resno == lig[i,2])
			#lig_vector = c(lig_vector,paste0("L",i))
			#new_chain_name = as.character(i%%base)
			new_chain_name = paste0("L",i)
			lig_vector = c(lig_vector,new_chain_name)
			#pdb$atom$chain[f] = paste0("L",i)
			pdb$atom$chain[f] = new_chain_name
			f = (pdb$ligands$resname == lig[i,1])&(pdb$ligands$resid == lig[i,2])
			#pdb$ligands$chainid[f] = paste0("L",i)
			pdb$ligands$chainid[f] = new_chain_name
		}

	}else{
		print(paste("WARNING: no ligands found in PDB",pdb$pdbname))
	}
	weye = F
	if (sum(!f)>0){
		weye = T
		f = pdb$atom$resid == water
		pdb$atom$chain[f] = WT
		f = pdb$ligands$chainid == water
		pdb$ligands$chainid[f] = WT
	}else{
		print(paste("WARNING: no waters found in PDB",pdb$pdbname))
	}
	pdb$chains = unique(pdb$atom$chain)
	names(pdb$chains) = rep("C",length(pdb$chains))
	if (length(lig_vector)){
		f = pdb$chains %in% lig_vector
		names(pdb$chains)[f] = rep("1",length(lig_vector))
	}
	if (weye){
		f = pdb$chains %in% WT
		names(pdb$chains)[f] = "0"
	}
	#print(pdb$chains);readline()
	return(pdb)

}

#if (0){
#	print(chain_a$chain);print(pdb_matrix$chain_a_name)
#	if (grepl(join,pdb_matrix$chain_a_name,fixed=T)){
#		chain_join = unlist(strsplit(chain_a_name,join,fixed=T))
#		w_a = which((chain_a$resid=="HOH")&(chain_a$chain %in% chain_join))
#		if (sum(diff(w_a)>1)){
#			print(paste("WARNING: no contiguous waters found in",pdb_matrix$pdbname,"-",pdb_matrix$chain_a_name))
#		}
#	}else{
#		w_a =  which((chain_a$resid=="HOH")&(chain_a$chain==pdb_matrix$chain_a_name))
#	}
#	lines_w_a = length(w_a);print(w_a);readline()
#	mat_dist = dist2(chain_a[w_a,9:11],chain_a[w_a,9:11])
#	for(i in 1:lines_w_a){
#  	for(j in 1:lines_w_a){        
#    matriz_adj_distancias[(w_a[1]-1)+i,(w_a[1]-1)+j] = mat_dist[i,j]
#    matriz_adj_distancias[(w_a[1]-1)+j,(w_a[1]-1)+i] = mat_dist[i,j]        
#  	}
#	}#;print("o");readline()
#	if (grepl(join,pdb_matrix$chain_b_name,fixed=T)){
#		chain_join = unlist(strsplit(chain_b_name,join,fixed=T))
#		w_b = which((chain_b$resid=="HOH")&(chain_b$chain %in% chain_join))
#		if (sum(diff(w_b)>1)){
#			print(paste("WARNING: no contiguous waters found in",pdb_matrix$pdbname,"-",pdb_matrix$chain_b_name))
#		}
#	}else{
#		w_b =  which((chain_b$resid=="HOH")&(chain_b$chain==pdb_matrix$chain_b_name))#;print(lines_a+w_b[1]-1);readline()
#	}
#	lines_w_b = length(w_b);print(w_b);readline()
#	mat_dist = dist2(chain_b[w_b,9:11],chain_b[w_b,9:11])
#	for(i in 1:lines_w_b){
#  	for(j in 1:lines_w_b){        
#    matriz_adj_distancias[(lines_a+w_b[1]-1)+i,(lines_a+w_b[1]-1)+j] = mat_dist[i,j]
#   matriz_adj_distancias[(lines_a+w_b[1]-1)+j,(lines_a+w_b[1]-1)+i] = mat_dist[i,j]        
#  	}
#	}
#}


