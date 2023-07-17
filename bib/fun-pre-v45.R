
transrotation = function(rotall, rot, rot0){

	auxrotall = rotall
	auxrot = matrix(rot,ncol=4,byrow=T)
	auxrotall[1:3,1:3] = auxrot[1:3,1:3] %*% auxrotall[1:3,1:3] 
	#auxrotall[,4] = auxrot[,4]
	#auxrotall[4,] = auxrot[4,]
	auxrotall[,4] = rot0[,4]
	auxrotall[4,] = rot0[4,]
	#print(auxrot)
	return(auxrotall)

}

recover_transformation = function(out,id,r=5){

	auxrot = out[id[2]:id[3]]
	auxrot = str_c(auxrot,collapse="")
	#print(nchar(auxrot))
	auxrot = str_sub(auxrot,2,(nchar(auxrot)-1))
	auxrot = round(as.numeric(str_split(auxrot,",")[[1]]),r)
	return(auxrot)
}


align_graph_pdb = function(command,flags,script,args,idline=c(32,34,49),dline=19,tam=10,v1=".pml",v2="a.pml",r=1,dr=0.001){
#align_graph_pdb = function(command,flags,script,args,idline=c(36,38,53),dline=19,tam=10,v1=".pml",v2="a.pml",r=1,dr=0.001){

	auxrmsd = c()
	auxwin = c()
	auxrot = c()
	auxlist = list()
	
	#print(command)#;readline()
	#print(args);readline()
	args1 = c(flags,paste(script,v1,sep=""),args)
	#print(args)#;readline()
	auxout = system2(command,args=args1,stdout=T)
	auxline = unlist(strsplit(auxout[idline[1]]," "))
	#print(auxout)#;readline()
	#print(auxline);readline()
	auxlist[[1]] = list()
	if (auxline[1]=="RMSD"){
		auxnum = as.numeric(auxline[2])
		auxrot = recover_transformation(auxout,idline)
		#print(auxrot)#;readline()
		idline2 = idline+dline
		auxline = unlist(strsplit(auxout[idline2[1]]," "))
		if (auxline[1]=="RMSD"){
			auxnum2 = as.numeric(auxline[2])
			auxrot2 = recover_transformation(auxout,idline2)
			auxrot2a = matrix(auxrot2,ncol=4,byrow=T)[1:3,1:3]
			#print(auxrot2a)#;readline()
			#if (auxrot2a==diag(3)) {
			auxtest = sum(abs(auxrot2a-diag(3)))
			#print(auxtest);readline()
			if (auxtest<dr){
				#print("one");
				auxlist[[1]]$rmsd = auxnum
				auxlist[[1]]$rot = auxrot
				auxlist[[1]]$sym = 0
			}else{
				#print("two");
				i=1
				#print(auxrot);
				args2 = c(flags,paste(script,v2,sep=""),args)
				auxout = system2(command,args=args2,stdout=T)
				#print(auxout);readline()
				#idline = idline+8
				auxrotall = diag(4)
				auxline = unlist(strsplit(auxout[idline[1]]," "))
				auxnum = as.numeric(auxline[2])
				auxrot = recover_transformation(auxout,idline)
				auxrotall = matrix(auxrot,ncol=4,byrow=T)
				#auxrotall = transrotation(auxrotall,auxrot)
				#transform_by_rot_matrix(y,matrix(r,ncol=4,byrow=T))
				#print(auxrotall);readline()
				auxlist[[i]]=list()
				auxlist[[i]]$rmsd = auxnum
				auxlist[[i]]$rot = as.vector(t(auxrotall))
				auxlist[[i]]$sym = i
				#auxrotall2 = diag(4)
				idline = idline+dline
				auxrot0 = auxrotall
				for (i in 2:tam){
					#print(idline);readline()
					auxline = unlist(strsplit(auxout[idline[1]]," "))
					#print(auxline)
					auxnum = as.numeric(auxline[2])
					auxrot = recover_transformation(auxout,idline)
					auxrotall = transrotation(auxrotall,auxrot,auxrot0)
					#print(auxrotall);readline()
					#Revelou-se caro demais este teste:
					#auxrotall2 = round(transrotation(auxrotall2,auxrot),r)
					#print(auxrot);print(auxrotall2)#;readline()
					#if (sum(diag(auxrotall2))>=dr){
					#	print("convergiu!");readline()
					#}
					#print(matrix(auxrot,ncol=4,byrow=T));readline()
					#auxrot = matrix(auxrot,ncol=4,byrow=T)
					#auxrotall[1:3,1:3] = auxrot[1:3,1:3] %*% auxrotall[1:3,1:3] 
					#auxrotall[,4] = auxrot[,4]
					#auxrotall[4,] = auxrot[4,]
					#print(auxrot)
					#print(auxrotall);readline()
					auxlist[[i]]=list()
					auxlist[[i]]$rmsd = auxnum
					auxlist[[i]]$rot = as.vector(t(auxrotall))
					auxlist[[i]]$sym = i
					idline = idline+dline
					
				}
			}
			#print(auxrot2a);readline()
			#auxrot = auxout[idline[2]:idline[3]]
			#auxrot = str_c(auxrot,collapse="")
			#print(nchar(auxrot))
			#auxrot = str_sub(auxrot,2,(nchar(auxrot)-1))
			#auxrot = as.numeric(str_split(auxrot,",")[[1]])
			#print(auxrot);readline()
			#auxrotlist = c(auxrotlist,auxrot) 
			#auxrotlist[[k]] = auxrot;k=k+1
			#auxlist$rmsd = auxnum
			#auxlist$rot = auxrot
		}else{
			print(paste("WARNING: something wrong in cealign"))
		}
	}else{
		#print(paste("WARNING: RMSD line not found in cealign output for",args))
		auxlist[[1]]$rmsd = 0
		auxlist[[1]]$rot = as.vector(diag(4))
		auxlist[[1]]$sym = 0	
	}
	return(auxlist)

}


pre_align_graph_pdb = function(name1,name2,id1,id2,command="pymol",script="script-align-v2",workdir="../Pymol/",minwin=3,sufix=".pdb",r=5,guide=0){

	script = paste(workdir,script,sep="")
	#command = paste(workdir,command,sep="")
	name1 = paste(workdir,name1,sufix,sep="")
	name2 = paste(workdir,name2,sufix,sep="")
	auxrmsd = c()
	auxwin = c()
	auxrot = c()
	auxsym = c()
	#auxrotlist = list()
	auxlist = list()

	auxa=c()

	i=minwin		
	#i = 13
	flags = c("-cq")
	#args = c("-cp",script,name1,name2,i,guide) #ultimo se guide=0 (todos os atomos), guide=1 so CA
	args = c(name1,name2,i,guide)
	#auxlist = align_graph_pdb(command,args)
	#print(command);print(flags);print(script);print(args);readline()


	while(1){

		auxlist = align_graph_pdb(command,flags,script,args)
		#if (auxlist[[1]]$rmsd == 0) break;
		#print(auxlist);print(unlist(auxlist));readline()
		auxlen = length(auxlist)
		if (auxlen>1){
			for (j in 1:auxlen){
				auxr = c(id1=id1,id2=id2,win=i,unlist(auxlist[j]))
				#print(auxr);readline()
				auxa = rbind(auxa,auxr)
			}

		}else{
			#if (auxlist[[1]]$rmsd == 0) break;
			auxr = c(id1=id1,id2=id2,win=i,unlist(auxlist))
			auxa = rbind(auxa,auxr)
			if (auxlist[[1]]$rmsd == 0) break;
		}
		#auxr=unlist(auxlist)
		#print(length(auxlist));readline()
		#auxa = rbind(auxa,auxr)
		#print(auxa);readline()
		args[3] = i
		i=i+1

	}
	#print(auxa);readline()
	return(auxa)
}

par_alignments = function(group,pdb,filename){

	#pdbname = as.character(pdb)
	tamj = which(group$pdb==pdb)
	tamk = dim(group)[1]
	auxa = c()

	for (j in tamj){
		auxnamei = as.character(group[j,"fakename"])			
		for (k in 1:tamk){
			#j=5
			auxnamek = as.character(group[k,"fakename"])
			#print(auxnamei);print(auxnamek)
			#print(sprintf("%s - %s: %.3f - %.3f",auxnamei,auxnamek,j,k));readline()
			auxr=pre_align_graph_pdb(auxnamei,auxnamek,j,k)
			#print(auxr[dim(auxr)[1],])
			#print(dim(auxr)[1]);readline()
			auxa = rbind(auxa,auxr)
			#print(sprintf("%s - %s: %i - %i - %i",auxnamei,auxnamek,j,k,dim(auxa)[1]))#;readline()
		}
	}
	#filename = paste(outdir,toupper(pdbname),sep,sufix,sep="")
	#print(filename);readline()
	write.table(auxa,filename,row.names=F,col.names=T,sep=",")

}

put_to_pre_process = function(queue,n,pdblist,outdir,command,par,prescript,filenamein,sufixout){

	i = queue$i
	while(n){
		#print(i)	
		if (i>queue$tami) {
			print(paste("No more PDB to queue. Finishing..."))
			break
		}
		pdbi = pdblist[i]
		pdbname = as.character(pdbi)
		filenameout = paste(outdir,toupper(pdbname),sufixout,sep="")
		if (!file.exists(filenameout)){
			#print(is.numeric(n))
			if (n){
				print(paste("Initializing alignment precomputation for",pdbname))
				t = proc.time()[3]
				#dt = round((t2-t1)/b,k)
				preline = data.frame(pdb=filenameout,t0=t,t1=0,dt=0)
				queue$pretab = rbind(queue$pretab,preline)
				#print(filequeue)
				args1 = paste(par,prescript,pdbname,filenamein,filenameout)#;print(args1);readline()
		    		x=system2(command,args=args1,wait=F)
				n = n - 1
				#i = i + 1
			} 	
		}
		i = i + 1
	}
	queue$i = i
	return(queue)
}

#pre_computation_of_alignments = function(group,generalpath,version){

#	auxall=list()
	#tamk = dim(dball$group)[1]#;tam=3
	#tam = c(1)
	#pdblist = levels(dball$group$pdb)
#	pdblist = unique(group$pdb)
#	tami = length(pdblist)
	#print(pdblist);readline()
	#auxa = c()

	#for (i in 1:tami){
#	for (i in 1){
#		pdbi = pdblist[i]
		#pdbname = as.character(pdbi)
		#tamj = which(dball$group$pdb==pdbi)
		#auxa = c()
#		print(paste("Pre-computation alignment for",pdbname))
#		par_alignments(group,pdbi)
		#for (j in tamj){
		#	auxnamei = as.character(dball$group[j,"fakename"])			
		#	for (k in 1:tamk){
				#j=5
		#		auxnamek = as.character(dball$group[k,"fakename"])
				#print(auxnamei);print(auxnamek)
				#print(sprintf("%s - %s: %.3f - %.3f",auxnamei,auxnamek,j,k));readline()
		#		auxr=pre_align_graph_pdb(auxnamei,auxnamek,j,k)
				#print(auxr[dim(auxr)[1],])
				#print(dim(auxr)[1]);readline()
		#		auxa = rbind(auxa,auxr)
		#		print(sprintf("%s - %s: %i - %i - %i",auxnamei,auxnamek,j,k,dim(auxa)[1]))#;readline()
		#	}
		#}
		#filename = paste(outdir,toupper(pdbname),sep,sufix,sep="")
		#print(filename);readline()
		#write.table(auxa,filename,row.names=F,col.names=T,sep=",")
		#print(filename);readline()
	
#	}

#}


