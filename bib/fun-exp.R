
get_info_connections_group_exp_old = function(conn,local,global,spot,dense,center){

	tam = length(conn)
	res=list()

	for (i in 1:tam){
		#i=5
		res[[i]] = list()#;print(conn[[i]]);readline()
		res[[i]]$a = conn_to_adj_matrix(conn[[i]])
		res[[i]]$xyz = conn_to_xyz_matrix(conn[[i]])#;print(res);readline()
		res[[i]]$sensitivity = local[[i]]
		res[[i]]$accuracy = global[[i]]
		#res[[i]]$spot = list()
		#res[[i]]$spot$score = spot$score[[i]]#;print(spot$score[[i]]);readline()
		#res[[i]]$spot$colors = spot$colors[[i]]
		res[[i]]$spot = spot[[i]]$score
		res[[i]]$dense = dense[[i]]
		res[[i]]$center = center[[i]]
		#res[[i]]$ga = graph.adjacency(res[[i]]$a,weighted=TRUE,mode="undirected",diag=FALSE)
		#vertex_attr(res[[i]]$ga) = list(label=rep(1,i))
		#print(res);readline()
	}

	return(res)

}



build_dbsup_from_dbgrp_exp_old = function(db,trans=F,sep="_"){

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
				dense = db[[i]][[j]]$dense_bsr_total#;print(spot);print(dense);readline()
				center = db[[i]][[j]]$center_bsr
				#if (trans) rot = db[[i]][[j]]$rotall
				#else rot = NULL
				rot = db[[i]][[j]]$rotall
				dbsup[[i]][[j]]$clus[[1]] = list()
				dbsup[[i]][[j]]$clus[[1]]$rotself = rot
				res = get_info_connections_group_exp(conn,local,global,spot,dense,center)#;print(res[[5]]);readline()
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
				res = get_info_connections_group_exp(conn,local,global,spot,dense,center)
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
				res = get_info_connections_group_exp(conn,local,global,spot,dense,center)
				dbsup[[i]][[j]]$clus[[3]]$nclus = res
			}#;print("o")
		}
	}
	return(dbsup)
}


complete_dbsup_exp = function(dba2,dbsup){

	tami = length(dbsup)

	for (i in 1:tami){
		tamj = length(dbsup[[i]][[2]]$clus)
		for (j in 1:tamj){
			tamk = length(dbsup[[i]][[2]]$clus[[j]]$nclus)
			ele = dba2[[j]][[2*i]]$element_name
			pol = dba2[[j]][[2*i]]$element_polarity
			cen = dba2[[j]][[2*i]]$centrality[,"n.closeness"]
			a = dba2[[j]][[2*i]]$A$a
			dbsup[[i]][[2]]$clus[[j]]$element_name = ele
			dbsup[[i]][[2]]$clus[[j]]$element_polarity = pol
			dbsup[[i]][[2]]$clus[[j]]$a = a
			dbsup[[i]][[2]]$clus[[j]]$center = cen
			for (k in 1:tamk){
				clusids = dba2[[j]][[2*i]]$exp$pam[[k]]$clustering
				dbsup[[i]][[2]]$clus[[j]]$nclus[[k]]$clusids = clusids 
			}
		}	
	}
	return(dbsup)

}


make_data_exp = function(dbsup,filename,inter=2,type=1,aanum=NULL,bb=c("N","C","O","CA","CB"),sep="_"){


	tami = length(dbsup)
	tab = read.csv(file=filename)
	tab$PDB = toupper(tab$PDB)
	tab$Residue = aa123(as.character(tab$Residue))
	tab14 = apply(mapply(paste,tab[,1:4]),1,paste,collapse="_")
	tab = data.frame(tab,Ele14=tab14)
	o = order(tab$Ele14)
	tab = tab[o,]
	tab.lev = levels(tab$Ele14)

	dbexp=list()
	
	for (i in 1:tami){
		dbexp[[i]] = list()
		pdbid = unlist(strsplit(dbsup[[i]][[inter]]$pdb,sep))[1]#;print(pdbid);readline()
		
		tamj = length(dbsup[[i]][[inter]]$clus)
		#auxr = c()
		for(j in 1:tamj){ # j em polarity
			dbexp[[i]][[j]] = list()
			ele = dbsup[[i]][[inter]]$clus[[j]]$element_name
			#tamk = length(ele)
			eletab = data.frame(sep_atom_name_id(ele,sep="[.]",ncol=3))
			eletab = data.frame(eletab,X4=eletab$X3)#;print(head(eletab));readline()
			#colnames(eletab) = c("chain","resn","atom")
			resn = substr(eletab[,2],1,3)
			resi = as.numeric(str_extract(eletab[,2],"[[:digit:]]+"))
			eletab$X2 = resn
			eletab$X3 = resi
			colnames(eletab)=c("Chain","Residue","Sequence","Atom")
	
			eletab = data.frame(PDB=pdbid,eletab)
			ele14 = apply(mapply(paste,eletab[,1:4]),1,paste,collapse="_")
			#eletab = data.frame(eletab,PCRS=ele14,Dg=ele14)
			eletab = data.frame(eletab,PCRS=ele14)
			ele15 = apply(eletab[,c("PCRS","Atom")],1,paste,collapse="-")
			eletab = data.frame(eletab,PCRSA=ele15)
			eletab = data.frame(eletab,Dg=ele14)
			eletab.lev = levels(eletab$Dg)#;print(eletab.lev);readline()
			#tab14 = apply(mapply(paste,tab[,1:4]),1,paste,collapse="_")
			#tab = data.frame(tab,Ele14=tab14)
			v = rep(0,length(eletab.lev))
			f21 = tab.lev %in% eletab.lev#;print(tab.lev);print(eletab.lev);print(f21);readline()
			g = tab$Gobs[f21]#;print(g);readline()
			f12 = eletab.lev %in% tab.lev
			v[f12] = g
			levels(eletab$Dg)=v
			eletab$Dg = as.numeric(as.vector(eletab$Dg))
			f = eletab$Atom %in% bb
			eletab$Dg[f]=0

			#print(head(eletab[,"PCRS"]));readline()

			ele.pcrs = data.frame(levels(eletab[,"PCRS"]))
			eletab$Residue = factor(eletab$Residue)
			ele.res =  data.frame(levels(eletab[,"Residue"]))#;print(ele.res)

			res.join = unlist(apply(ele.res,1,filter_sum,eletab,"Residue","Residue","aa",aanum))#; print(res.join);readline()
			res.join = res.join+rnorm(length(res.join),sd=0.001)
			eletab = data.frame(eletab,Len=eletab$Residue)
			levels(eletab$Len) = res.join
			eletab$Len = round(as.numeric(as.vector(eletab$Len)),0)
			#print(eletab);readline()

			#res.join = apply(ele.pcrs,1,filter_sum,eletab,"PCRS","PCRS","count")#;print(res.join);readline()
			#res.join = res.join+rnorm(length(res.join),sd=0.001)
			#eletab = data.frame(eletab,Len=eletab$PCRS)
			#levels(eletab$Len) = res.join
			#eletab$Len = round(as.numeric(as.vector(eletab$Len)),0)
			#print(head(eletab,30));readline()

			#Dg = as.numeric(as.vector(eletab$Dg))

			eletab = data.frame(eletab, DgM = round(eletab$Dg/eletab$Len,3))
			#print(eletab);readline()

			a = dbsup[[i]][[inter]]$clus[[j]]$a
			eletab = data.frame(eletab,Area=round(apply(a,1,sum),2))
	
			res.join = apply(ele.pcrs,1,filter_sum,eletab,"PCRS","Area","sum")
			res.join = res.join+rnorm(length(res.join),sd=0.001)
			eletab = data.frame(eletab,AreaRS=eletab$PCRS)
			levels(eletab$AreaRS) = round(res.join,2)
			eletab$AreaRS = as.numeric(as.vector(eletab$AreaRS))

			res.join = apply(ele.pcrs,1,filter_sum,eletab,"PCRS","Area","mean")
			res.join = res.join+rnorm(length(res.join),sd=0.001)
			eletab = data.frame(eletab,AreaRM=eletab$PCRS)
			levels(eletab$AreaRM) = round(res.join,2)
			eletab$AreaRM = as.numeric(as.vector(eletab$AreaRM))

			cen = round(100*dbsup[[i]][[inter]]$clus[[j]]$center,2)

			eletab = data.frame(eletab,Center=cen)
			res.join = apply(ele.pcrs,1,filter_sum,eletab,"PCRS","Center","sum")
			res.join = res.join+rnorm(length(res.join),sd=0.001)
			eletab = data.frame(eletab,CenteRS=eletab$PCRS)
			levels(eletab$CenteRS) = round(res.join,2)
			eletab$CenteRS = as.numeric(as.vector(eletab$CenteRS))

			res.join = apply(ele.pcrs,1,filter_sum,eletab,"PCRS","Center","mean")
			res.join = res.join+rnorm(length(res.join),sd=0.001)
			eletab = data.frame(eletab,CenteRM=eletab$PCRS)
			levels(eletab$CenteRM) = round(res.join,2)
			eletab$CenteRM = as.numeric(as.vector(eletab$CenteRM))
			
			tamk = length(dbsup[[i]][[inter]]$clus[[j]]$nclus)

			if (type[1]){
				for (k in 1:tamk){
					clusid = dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$clusids
				
					eletab = data.frame(eletab,Qualy=clusid)#;print(head(eletab));readline()
					tamele = dim(eletab)[2]
					elename = paste0("Qualy",k)
					colnames(eletab)[tamele]=elename
					qualy = 100*dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$sensitivity
					qualy = qualy+rnorm(length(qualy),sd=0.001)
					eletab[,elename]=factor(eletab[,elename])
					levels(eletab[,elename])=qualy
					eletab[,elename] = as.double(as.character(eletab[,elename]))

				}#;print("o")
			}
			if (type[2]){
				for (k in 1:tamk){
					clusid = dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$clusids
				
					eletab = data.frame(eletab,Spot=clusid)#;print(head(eletab));readline()
					tamele = dim(eletab)[2]
					elename = paste0("Spot",k)
					colnames(eletab)[tamele]=elename
					spot = 100*dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$spot
					spot = spot+rnorm(length(spot),sd=0.001)
					eletab[,elename]=factor(eletab[,elename])
					levels(eletab[,elename])=spot
					eletab[,elename] = as.double(as.character(eletab[,elename]))

				}#;print("o")
			}
			if (type[3]){
				for (k in 1:tamk){
					clusid = dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$clusids
				
					eletab = data.frame(eletab,Dense=clusid)#;print(head(eletab));readline()
					tamele = dim(eletab)[2]
					elename = paste0("Dense",k)
					colnames(eletab)[tamele]=elename
					dense = dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$dense
					dense = dense+rnorm(length(dense),sd=0.001)
					eletab[,elename]=factor(eletab[,elename])
					levels(eletab[,elename])=dense
					eletab[,elename] = as.double(as.character(eletab[,elename]))

				}#;print("o")
			}
			if (type[4]){
				for (k in 1:tamk){
					clusid = dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$clusids
				
					eletab = data.frame(eletab,Center=clusid)#;print(head(eletab));readline()
					tamele = dim(eletab)[2]
					elename = paste0("Center",k)
					colnames(eletab)[tamele]=elename
					center = dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$center
					center = center+rnorm(length(center),sd=0.001)
					eletab[,elename]=factor(eletab[,elename])
					levels(eletab[,elename])=center
					eletab[,elename] = as.double(as.character(eletab[,elename]))

				}#;print("o")
			}
			if (type[5]){
				for (k in 1:tamk){
					clusid = dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$clusids
					eletab = data.frame(eletab,Clusid=clusid)#;print(head(eletab));readline()
					tamele = dim(eletab)[2]
					elename = paste0("Clusid",k)
					colnames(eletab)[tamele]=elename
				}
			}
			if(type[6]){
				for (k in 1:tamk){
					#print(head(eletab))
					elename1 = paste0("Clusid",k)
					elename2 = paste0("AreaC",k)
					#clusid = dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$clusids							
					ele.a = data.frame(levels(factor(eletab[,elename1])))#;print(levels(factor(eletab[,elename1])));readline()
					colnames(ele.a)=elename1
					#print(ele.a)
					ele.join = apply(ele.a,1,filter_sum,eletab,elename1,"Area")
					ele.join = ele.join+rnorm(length(ele.join),sd=0.001)
					eletab = data.frame(eletab,AreaC=eletab[,elename1])#;print(ele.join);readline()
					tamele = dim(eletab)[2]
					colnames(eletab)[tamele]=elename2
					eletab[,elename2] = factor(eletab[,elename2])
					#print(levels(eletab[,elename2]))
					levels(eletab[,elename2])=ele.join
					eletab[,elename2] = round(as.numeric(as.vector(eletab[,elename2])),2)
					#;print(eletab[,elename2]);readline()
				}
			}
			if(type[7]){
				for (k in 1:tamk){
					elename1 = paste0("Clusid",k)
					elename2 = paste0("ClusA",k)
					ele.a = data.frame(levels(factor(eletab[,elename1])))#;print(ele.a);readline()
					ele.join = diag(dbsup[[i]][[inter]]$clus[[j]]$nclus[[k]]$a)
					ele.join = ele.join+rnorm(length(ele.join),sd=0.001)
					eletab = data.frame(eletab,AreaC=eletab[,elename1])
					tamele = dim(eletab)[2]
					colnames(eletab)[tamele]=elename2
					eletab[,elename2] = factor(eletab[,elename2])
					levels(eletab[,elename2])=ele.join
					eletab[,elename2] = round(as.numeric(as.vector(eletab[,elename2])),2)
					#ele.join
				}
			}
			#f=!(colnames(eletab) %in% paste("Clusid",1:20,sep=""))
			#eletab = eletab[,f]
			#print(eletab);readline()
			#dbexp[[1]][[1]]$Dg=as.numeric(as.vector(dbexp[[1]][[1]]$Dg))
			eletab$Dg = as.numeric(as.vector(eletab$Dg))
			dbexp[[i]][[j]]$tab = eletab
			#return(eletab)
			#f = which((eletab$Chain %in% tab$Chain)&(eletab$Residue %in% tab$Residue)&(eletab$Sequence %in% tab$Sequence))
			#print(f);readline()
			#return(eletab)
			#x=list()
			#x$ele = eletab
			#x$tab = tab
			#return(x)
			#print(head(eletab));print(tail(eletab));readline()
			
		}

	}
	return(dbexp)
	
}

sum_nan = function(v){

	return(as.logical(sum(is.nan(v))))


}


view_exp_corr = function(tab,n=300,pol=c("ALL","PP","AA"),sep="-",cex=0.8){


	tam = length(tab)
	par(mfrow=c(3,3))
	for (i in 1:tam){

		tab1 = filter(tab[[i]],Dg!=0)

		#head(tab1);readline()

		cols = colorRampPalette(c("brown","red", "yellow", "white"))(n = n)
		idcol = round(mapply(renorm_interval,tab1[,7],MoreArgs=list(a=c(-1,5),b=c(1,n),inv=T)),0)
		col = cols[idcol]

		pdb = tab1[1,1]
		chains = unique(tab1$Chain)
		chains = paste(chains,collapse=":")

		x = tab1$Area
		y = tab1$Center
		xlab = "Area"
		ylab = "Center"
		main = paste(pdb,sep,chains,sep,pol[i])
		plot(x,y,type="p",bg=col,pch=21,xlab=ylab,ylab=ylab,main=main)
		mtext("per atom",cex=cex)

		x = tab1$AreaRS
		y = tab1$CenteRS
		xlab = "AreaRS "
		ylab = "CenterRS"
		main = paste(pdb,sep,chains,sep,pol[i])
		plot(x,y,type="p",bg=col,pch=21,xlab=ylab,ylab=ylab,main=main)
		mtext("sum per residue",cex=cex)

		x = tab1$AreaRM
		y = tab1$CenteRM
		xlab = "AreaRM "
		ylab = "CenterRM"
		main = paste(pdb,sep,chains,sep,pol[i])
		plot(x,y,type="p",bg=col,pch=21,xlab=ylab,ylab=ylab,main=main)
		mtext("mean per residue",cex=cex)
	}
	if (0){
		tab2 = filter(tab[[3]],Dg!=0)
		x = tab2$AreaRM
		y = tab2$CenteRM
		z = tab2$Dg
		#plot3d(x,y,z,size=1,col=col,aspect="iso",box=FALSE,type="s")
		plot3d(x,y,z,size=1,col=col,aspect=T,box=FALSE,type="s")
	}
	#return(tab1)

}

view_dendro_plot = function(db,id,par,type=1,n=300){

	tab = db[[id[1]]][[id[2]]]$tab

	tab1 = filter(tab,Dg>par)
	#print(head(tab1))
	rownames(tab1) = apply(tab1[,c("PCRS","Atom")],1,paste,collapse="-")	
	#pdbn = levels(tab1$PCRS)
	#tam = length(pdbn)
	#cols = colorRampPalette(c("brown","red", "yellow", "white"))(n = n)
	#idcol = round(mapply(renorm_interval,1:tam,MoreArgs=list(a=c(0,tam),b=c(1,n),inv=T)),0)
	#col = pdbn
	#levels(col)=cols[idcol]

	#pdbn.f = factor(pdbn[,1],ordered=F)
	#tam = length(levels(pdbn.f))
	#cols = colorRampPalette(c("brown","red", "yellow", "white"))(n = n)
	#idcol = round(mapply(renorm_interval,1:tam,MoreArgs=list(a=c(0,tam),b=c(1,n),inv=T)),0)
	#col = pdbn.f
	
	#print(head(x));readline()
	f = colnames(tab1) %in% paste("Clusid",1:15,sep="")
	tab1 = tab1[,f]
	if (type==1){
		#tab2 = melt(t(tab1))
		#px = ggplot(tab2, aes(Var1, Var2)) + geom_tile(aes(fill = value)) + geom_text(aes(label = round(value, 1))) + scale_fill_gradient(low = "white", high = "red")
		#print(px)
		tab2 = melt(t(tab1))
		xf = factor(tab2$value)
		tamxf = length(levels(xf))
		levels(xf) = letters[1:tamxf]
		tab2 = data.frame(tab2,label=xf)
		px = ggplot(tab2, aes(Var1, Var2)) + geom_tile(aes(fill = value)) + geom_text(aes(label = label),size=3) + scale_fill_gradient(low = "white", 		high = "red")
		print(px)
		return(db)
	}
	if(type==2){

		tab1 = data.frame(apply(tab1,2,tofactor,letters))#;print(data.frame(x));readline()
		#dd = tab1 %>% daisy(metric="gower") %>% hclust(method="mcquitty") %>% as.dendrogram
		dd = tab1 %>% daisy(metric="gower") %>% hclust(method="ward.D2") %>% as.dendrogram
		ddata = dendro_data(dd)
		p2 = ggplot(segment(ddata)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
		#print(p2)
		#dd %>% rect.dendrogram(k=3,border=1)
		#readline()
		labs = label(ddata) 
		pdbn = matrix(unlist(strsplit(as.character(labs$label),"-")),ncol=2,byrow=T)

		#levels(col) = cols[idcol]
		labs$group = pdbn[,1]
		p2 <- p2 + geom_text(data=label(ddata), aes(label=label, x=x, y=0, colour=labs$group),hjust = "right",nudge_y=-0.05,size=2.5) + scale_y_continuous(expand = c(.4, .1))
		p2 <- p2 + coord_flip() + theme(legend.position="none")
		print(p2)
		hclu = list()
		hclu$dd = dd
		hclu$ddata = ddata
		hclu$labs = labs
		db[[id[1]]][[id[2]]]$hclu = hclu
		return(db)
	}

}

dentro_analysis_exp = function(db,id,h=NULL,k=NULL,minclus=10){

	tab = db[[id[1]]][[id[2]]]$tab 
	dd = db[[id[1]]][[id[2]]]$hclu$dd
	if (!is.null(h)){
		cutt = cutree(dd,h=h)
	}
	if (!is.null(k)){
		cutt = cutree(dd,k=k)
	}
	cuttab = sort(table(cutt),decreasing=T) 
	tam  = length(cuttab)
	cutnam = names(cuttab)
	i = 1
	#res = list()
	res = data.frame()
	repeat{
		#print(cuttab[i])
		if (i>tam){
			break;
		}
		if (cuttab[i]<minclus){
			break;
		}

		#res[[i]] = list()
		idi = cutnam[i]
		f = cutt == idi
		cutlab = names(cutt[f])
		#print(cutlab);readline()
		#res[[i]]$labels = cutlab
		#aux = data.frame(aux,nodes=cutlab)
		f = tab$PCRSA %in% cutlab
		aux = data.frame(tab[f,c("PDB","Chain","Residue","Atom","PCRS","PCRSA","Dg","DgM","Area","Center")])
		aux = data.frame(aux,Clusid=i)
		res = rbind(res,aux)#;print(res);readline()
		i = i + 1
	}
	db[[id[1]]][[id[2]]]$res = res
	return(db)
}

make_multi_dentro_tab = function(db,pol,k,d=1.9,r=2){

	tami = length(db)
	Dgt = c()
	aux = data.frame()
	pdbs = c()
	for (i in 1:tami){
		res = db[[i]][[pol]]$res
		pdbs = c(pdbs,as.character(res$PDB[1]))
		f = res$Dg > d
		#sumdg = sum(abs(res$DgM[f]))#;print(sumdg);readline()
		sumdg = sum(f)
		tamj = k
		Dgm = c()
		for (j in 1:tamj){
			f = res$Clusid == j
			dbj = sum(res$Dg[f]>d)
			#Dgm = c(Dgm,round(dbj/sumdg,r))
			Dgm = c(Dgm,dbj)
		}
		Dgt =  rbind(Dgt,Dgm)
		#names(Dgm) = paste("Clusid",1:k,sep="")
	}
	aux = as.data.frame(Dgt)
	rownames(aux) = pdbs
	colnames(aux) = paste("Clusid",1:k,sep="")
	return(aux)
	#db$multi = aux
	#return(db)
	#print(aux);readline()
}

find_min_sensi = function(nclus,k){

	tam = length(nclus)

	for (i in 2:tam){
		mins = min(nclus[[i]]$sensitivity)
		if (mins<k){
			return(i-1)
		}
	}
	return(tam)

}


prepair_logistic_tab3_notie = function(db,sup,pol,k,j,lim=2){


	tabno = prepair_logistic_tab2_notie(db,sup,pol=pol,k=k,j=j,type=1,lim=-100)
	#print(head(tab));readline()
	#tami = unique(tab$Id)
	tami = dim(tabno)[1]

	dg4 = c()
	dg2 = c()
	dg0 = c()

	for (i in 1:tami){
	#for (i in 2){
		id = tabno[i,"Id"]#;print(id)
		tab = db[[id]][[pol]]$tab
		#coltab = colnames(tab)
		n=find_min_sensi(sup[[id]][[2]]$clus[[pol]]$nclus,k)#;print(n)
		PDB = tabno[i,"PDB"]#;print(PDB)
		Clusidi = tabno[i,"Clusid1"]#;print(Clusidi)
		ClusAAi = tabno[i,"ClusAA1"]#;print(ClusAAi)
		ClusAn = paste0("ClusA",n)#;print(ClusAn)
		Clusidn = paste0("Clusid",n)#;print(Clusidn)
		tabi = tab %>% filter_at(vars(PDB),all_vars(. == PDB)) 
		tabi = tabi %>% filter_at(vars(ClusAn),all_vars(. == ClusAAi))
		tabi = tabi %>%	filter_at(vars(Clusidn),all_vars(. == Clusidi))
		#tabi = tabi %>%	filter_at(vars("Clusid3"),all_vars(. == 2))
		#print(head(tabi));readline()
		#filter_at(vars(ClusAn),all_vars(. == ClusAA1)) %>%
		#filter_at(vars(Clusidn),all_vars(. == Clusid1)) #%>%
		#return(tabi)
		#print(head(tabi));readline()
		tabi1 = tabi %>% filter_at(vars("Dg"),all_vars(. >= 2*lim))
		dg4 = c(dg4,ifelse(dim(tabi1)[1],1,0))#;print(dg4)
		tabi2 = tabi %>% filter_at(vars("Dg"),all_vars(. > lim))
		dg2 = c(dg2,ifelse(dim(tabi2)[1],1,0))#;print(dg2)
		tabi3 = tabi %>% filter_at(vars("Dg"),all_vars(. <= lim)) %>% filter_at(vars("Dg"),all_vars(. != 0)) 
		dg0 = c(dg0,ifelse(dim(tabi3)[1],1,0))#;print(dg0);readline()

		#filter_at(vars("Db"),all_vars(. == Clusid1)) %>%
		#print(tabi);readline()
		#print(head(tabi));readline()
	}
	tabno = data.frame(tabno,Dg4=dg4,Dg2=dg2,Dg0=dg0)
	return(tabno)
}

prepair_logistic_tab2_notie = function(db,sup,pol,k,j,type=1,lim=1.9){

	tami = length(db)
	res = c()
	#print(lim)
	for (i in 1:tami){
	#for (i in 3){
		
		tab = db[[i]][[pol]]$tab
		#print(tab[1,1])
		tabi = tab[,"PDB"]
		tabi = data.frame(Id=i,tabi)
		if (type==1){
			tabi = data.frame(tabi,Dgb=factor(ifelse(tab$Dg>lim,1,0)))
		}
		if (type==2){
			tabi = data.frame(tabi,Dgb=factor(ifelse((tab$Dg<=lim)&(tab$Dg!=0),1,0)))
		}
		#;print(head(tabi));readline()
		colnames(tabi)[2]="PDB"
		n=find_min_sensi(sup[[i]][[2]]$clus[[pol]]$nclus,k)#;print(n)
		colnamen = paste0("ClusA",n)
		clua = tab[,colnamen]

		tabi = data.frame(tabi,cluaa=clua)
		lentabi = length(colnames(tabi))
		colnamej = paste0("ClusAA",j)
		colnames(tabi)[lentabi] = colnamej

		colnamen = paste0("Clusid",n)
		cluid = tab[,colnamen]

		tabi = data.frame(tabi,cluid=cluid)
		lentabi = length(colnames(tabi))
		colnamej2 = paste0("Clusid",j)
		colnames(tabi)[lentabi] = colnamej2
		#print(head(tabi));readline()
		#print(head(tabi[,c("Dgb",colnamej,colnamej2)]));readline()
		tabi = unique(tabi[,c("Id","PDB","Dgb",colnamej,colnamej2)])
		#return(tabi)
		#print(head(tabi));readline()
		clua = tabi[,colnamej]
				
		clua = scale(clua)

		tabi = data.frame(tabi,clua=clua)#;print(head(tabi));readline()
		lentabi = length(colnames(tabi))
		colnamej = paste0("ClusA",j)
		colnames(tabi)[lentabi] = colnamej

		tabi = data.frame(tabi,cluan=n)
		lentabi = length(colnames(tabi))
		colnamej = paste0("ClusAn",j)
		colnames(tabi)[lentabi] = colnamej
		#tabi = data.frame(tabi,n=n)
		#print(head(tabi));readline()
		res = rbind(res,tabi)
	}
	return(res)
}

prepair_logistic_tab2 = function(db,sup,pol,k,lim=1.9,type=1){

	tami = length(db)
	res = c()
	for (i in 1:tami){
		
		tab = db[[i]][[pol]]$tab
		#print(tab[1,1])
		tabi = tab[,"PCRSA"]
		if (type==1){
			tabi = data.frame(tabi,Dgb=factor(ifelse(tab$Dg>lim,1,0)))
		}
		if (type==2){
			tabi = data.frame(tabi,Dgb=factor(ifelse((tab$Dg<=lim)&(tab$Dg!=0),1,0)))
		}
		colnames(tabi)[1]="PCRSA"
		tamk = length(k)
		for (j in 1:tamk){
			n=find_min_sensi(sup[[i]][[2]]$clus[[pol]]$nclus,k[j])#;print(n);readline()
			colnamen = paste0("ClusA",n)
			clua = tab[,colnamen]

			tabi = data.frame(tabi,cluaa=clua)
			lentabi = length(colnames(tabi))
			colnamej = paste0("ClusAA",j)
			colnames(tabi)[lentabi] = colnamej
			
			#if (noties){
			#	colnamen = paste0("Clusid",n)
			#	cluid = tab[,colnamen]

			#	tabi = data.frame(tabi,cluid=cluid)
			#	lentabi = length(colnames(tabi))
			#	colnamej2 = paste0("Clusid",j)
			#	colnames(tabi)[lentabi] = colnamej2
				#print(head(tabi));readline()
				#print(head(tabi[,c("Dgb",colnamej,colnamej2)]));readline()
			#	tabi = unique(tabi[,c("Dgb",colnamej,colnamej2)])
				#print(head(tabi));readline()
			#	clua = tabi[,colnamej]
				#tabi = unique(c())
	
			#}

			#clua = clua/sum(unique(clua))#;print(clua);readline()
			#clua = clua/max(clua)
			
			clua = scale(clua)

			tabi = data.frame(tabi,clua=clua)#;print(head(tabi));readline()
			lentabi = length(colnames(tabi))
			colnamej = paste0("ClusA",j)
			colnames(tabi)[lentabi] = colnamej

			tabi = data.frame(tabi,cluan=n)
			lentabi = length(colnames(tabi))
			colnamej = paste0("ClusAn",j)
			colnames(tabi)[lentabi] = colnamej
			#tabi = data.frame(tabi,n=n)
			print(head(tabi));readline()
		}
		res = rbind(res,tabi)
	}
	return(res)

}

prepair_logistic_tab = function(db,pol,k,lim=1.9){

	tami = length(db)
	res = c()
	for (i in 1:tami){
		tab = db[[i]][[pol]]$tab
		tabi = tab[,"PCRSA"]
		tabi = data.frame(tabi,Dgb=factor(ifelse(tab$Dg>lim,1,0)))
		colnames(tabi)[1]="PCRSA"
		for (j in k){
			colnamek = paste0("ClusA",j)
			clua = tab[,colnamek]
			clua = clua/sum(unique(clua))#;print(clua);readline()
			tabi = data.frame(tabi,clua=clua)#;print(head(tabi));readline()
			lentabi = length(colnames(tabi))
			colnames(tabi)[lentabi] = colnamek
			#print(head(tabi));readline()
		}
		res = rbind(res,tabi)
	}
	return(res)

}

make_multi_dentro_tab2 = function(db,pol,k,d=2,r=2,type="energy"){

	tami = length(db)
	Dgt1 = c()
	Dgt2 = c()
	pdbs = c()
	for (i in 1:tami){
		tab = db[[i]][[pol]]$tab
		pdbi = as.character(tab$PDB[1])
		pdbi = toupper(substr(pdbi,1,4))
		chaini = paste(as.character(unique(tab$Chain)),collapse=":")
		pdbi = paste(pdbi,chaini)
		pdbs = c(pdbs,pdbi)
		#f = tab$Dg > d
		#sumdg = sum(abs(res$DgM[f]))#;print(sumdg);readline()
		#sumdg = sum(f)
		tamj = k
		Dgm1 = c()
		Dgm2 = c()
		clun = paste("Clusid",k,sep="")
		#colid = which(colnames(res)==clun)
		for (j in 1:tamj){
			if (type=="area"){
				#d=-100
				f = tab[,clun] == j
				dbj1 = sum(tab$Area[f])
				#dbj2 = sum((tab$Dg[f]<d) & (tab$Dg[f]!=0))
				dbj2 = dbj1
			}
			if (type=="total"){
				d=-100
				f = tab[,clun] == j
				dbj1 = sum(tab$Dg[f]>=d)
				#dbj2 = sum((tab$Dg[f]<d) & (tab$Dg[f]!=0))
				dbj2 = dbj1
			}
			if (type=="freq"){
				f = tab[,clun] == j
				dbj1 = sum(tab$Dg[f]>=d)
				dbj2 = sum((tab$Dg[f]<d) & (tab$Dg[f]!=0))
			}
			if (type=="energy"){
				tabj = filter(tab,Dg>=d)
				f = tabj[,clun] == j
				dbj1 = sum(tabj$DgM[f])
				tabj = filter(tab,Dg<d & Dg!=0)
				f = tabj[,clun] == j
				dbj2 = sum(tabj$DgM[f])			
			}
			if (type=="aarea"){
				tabj = filter(tab,Dg>=d)
				f = tabj[,clun] == j
				dbj1 = sum(tabj$Area[f])
				tabj = filter(tab,Dg<d & Dg!=0)
				f = tabj[,clun] == j
				dbj2 = sum(tabj$Area[f])			
			}
			#Dgm = c(Dgm,round(dbj/sumdg,r))
			Dgm1 = c(Dgm1,dbj1)
			Dgm2 = c(Dgm2,dbj2)
		}
		Dgt1 =  rbind(Dgt1,Dgm1)
		Dgt2 =  rbind(Dgt2,Dgm2)
		#names(Dgm) = paste("Clusid",1:k,sep="")
	}
	aux1 = as.data.frame(Dgt1)
	aux2 = as.data.frame(Dgt2)
	rownames(aux1) = pdbs
	#colnames(aux1) = paste("id",1:k,sep="")
	colnames(aux1) = paste("id",1:k,sep="")
	rownames(aux2) = pdbs
	colnames(aux2) = paste("id",1:k,sep="")
	aux=list()
	aux$hot = aux1
	aux$cold = aux2
	return(aux)
	#db$multi = aux
	#return(db)
	#print(aux);readline()
}

view_spot_clus = function(tab){

	
	


}

make_dentro_clu_stat = function(db,ids){

	res = db[[ids[1]]][[ids[2]]]$res

	res$Clusid = as.factor(res$Clusid)
	
	#p1=ggplot(res,aes(x=Clusid,y=Area)) + geom_boxplot()
	#p2=ggplot(res,aes(x=Clusid,y=Center)) + geom_boxplot()

	p1 = qplot(Area,Center,data=res)
	p2 = p1 + facet_wrap( ~ Clusid, nrow = 3)

	grid.arrange(p1,p2,nrow=1)

	x11()

	resf = filter(res,Dg!=0)
	
	p1f = qplot(Area,Center,data=resf)
	p2f = p1f + facet_wrap( ~ Clusid, nrow = 3)

	grid.arrange(p1f,p2f,nrow=1)

	#print(p)

	#par(mfrow=c(1,2))
	#boxplot(Area ~ Clusid, res)
	#boxplot(Center ~ Clusid, res)

	#resf = filter(res,Dg!=0)
	#x11()
	#p=ggplot(x[c("Dg","Clusid")],aes(x=Clusid,y=Dg))+geom_point()
	#print(p)

	#boxplot(Dg ~Clusid, resf)
	

}

view_svd_exp = function(tab,n=300){


	tab1 = filter(tab,Dg>1.5)

	#names = tab1[,c("PCRS","Atom")]
	names = apply(tab1[,c("PCRS","Atom")],1,paste,collapse="-")#;print(names);readline()

	cols = colorRampPalette(c("brown","red", "yellow", "white"))(n = n)
	idcol = round(mapply(renorm_interval,tab1[,7],MoreArgs=list(a=c(-1,5),b=c(1,n),inv=T)),0)
	col = cols[idcol]


	if(1){
		notcol = c("PDB","Chain","Residue","Sequence","Atom","PCRS","Dg","Area","AreaRS","Center","CenteRS","AreaRM","CenteRM")
		#notcol = c(notcol,paste("Qualy",c(1:13,15),sep=""))
		#notcol = c(notcol,paste("Spot",c(1:13,15),sep=""))
		#notcol = c(notcol,paste("Dense",c(1:13,15),sep=""))
		#notcol = c(notcol,paste("Center",c(1:13,15),sep=""))

		f = colnames(tab1) %in% notcol

		tab1 = tab1[,!f]
		print(colnames(tab1))#;print(head(tab1));readline()
		tab1 = scale(tab1,center=T,scale=T)
		v = apply(tab1,2,sum_nan)
		tab1 = tab1[,!v]
		#print(tab1);readline()

		#print(head(tab1));print(tail(tab1));readline()

		#par(mfrow=c(1,2))
		#plot(tab$Area,tab$Center,bg=col,pch=21)#;readline()
		#plot(tab$AreaR,tab$CenteR,bg=col,pch=21);readline()

		#tam = dim(tab)[2]#;print(dim(tab))
		#tab1 = matrix( as.numeric( as.matrix (tab[,8:tam])   )  ,ncol=tam-8+1)#;print(dim(tab1));readline()
		#print(head(tab1));readline()
		#tam = dim(tab1)[2]
		#tab22 = scale(tab1[,2:2],center=T,scale=T)#;print(tab22);readline()
		#tab1a = scale(tab1[,1:tam],center=T,scale=T)

		#v = apply(tab1a,2,sum_nan)#;print(v);readline()
		#tab1a = tab1a[,!v]
		#tam = dim(tab1a)[2]
		#tab1 = cbind(tab1[,1],tab1a)### acrescenta Dg
		#tab1 = cbind(tab22,tab1a)### scale only area
		#tab1 = tab1a

		#print(head(tab1));readline()

		tab1.svd = svd(tab1)
		us = tomatrix(tab1.svd$u) %*% diag(tab1.svd$d)
		tab1.svd$us = us
	
		#print(cbind(idcol,tab[,"Dg"]));readline()
		x=us[,1]
		y=us[,2]
		z=us[,3]
		w=us[,4]
		par(mfrow=c(2,2))
		plot(x,y,type="p",bg=col,pch=21)
		plot(x,z,type="p",bg=col,pch=21)
		plot(y,z,type="p",bg=col,pch=21)
		tabf = filter(tab,Dg!=0)
		#plot(tabf$Dg,tabf$AreaR)
		#plot(x,w,type="p",bg=col,pch=21)
		plot(tab1.svd$d)
		#plot3d(x,y,z,size=1,col=col,aspect="iso",box=FALSE,type="s")
		plot3d(x,y,z,size=1,col=col,aspect=T,box=FALSE,type="p")
		text3d(x,y,z,names,cex=0.6,col=col)
		x = list()
		x$tab = tab
		x$tabn = tab1
		x$svd = tab1.svd
	
		return(x)
	}
	
	#print(head(tab1));readline()

	#tab1c = scale(tab1,center=T,scale=F)
	#tab1c = matrix(as.numeric(as.matrix(dbexp[[1]][[1]][1:5,7:9])),ncol=3)
	#tab1c.svd = svd(tab1c)
	#plot(tab1c.svd.

}




filter_sum = function(v,m,n,cn,ret="sum",va=NULL){
	

	#mf = filter(m,n==v)
	if (ret=="aa"){
		f = names(va) == v
		if (sum(f)){
			return(va[f])
		}else{
			return(Inf)	
		}
	}else{
		mf = filter_at(m, vars(n),any_vars(. == v))
		if (ret=="sum"){
			return(sum(mf[,cn]))
		}
		if (ret=="mean"){
			return(mean(mf[,cn]))
		}
		if (ret=="count"){
			return(length(mf[,cn]))
		}
	}

}

plot_corr_exp = function(x1,x2,x3){

	x1a = filter(x1$tab,Dg!=0)
	x2a = filter(x2$tab,Dg!=0)
	x3a = filter(x3$tab,Dg!=0)
	par(mfrow=c(1,3))
	plot(x1a$Dg,x1a$AreaR);print(cor(x1a$Dg,x1a$AreaR))
	plot(x2a$Dg,x2a$AreaR);print(cor(x2a$Dg,x2a$AreaR))
	plot(x3a$Dg,x3a$AreaR);print(cor(x3a$Dg,x3a$AreaR))

}


tofactor = function(v,a){

	vf = factor(v)
	#tam = length(levels(vf))
	vf.l = as.numeric(levels(vf))#;print(as.numeric(vf.l))
	levels(vf) = a[vf.l]#;print(a[as.numeric(vf.l)]);readline()
	return(vf)

}

############### TESTES #################


make_multi_dentro_tab_old = function(db,pol,k,r=2){

	tami = length(db)
	Dgt = c()
	aux = data.frame()
	pdbs = c()
	for (i in 1:tami){
		res = db[[i]][[pol]]$res
		pdbs = c(pdbs,as.character(res$PDB[1]))
		f = res$DgM!=0
		sumdg = sum(abs(res$DgM[f]))#;print(sumdg);readline()
		tamj = k
		Dgm = c()
		for (j in 1:tamj){
			f = res$Clusid == j
			Dgm = c(Dgm,round(sum(abs(res$DgM[f]))/sumdg,r))
		}
		Dgt =  rbind(Dgt,Dgm)
		#names(Dgm) = paste("Clusid",1:k,sep="")
	}
	aux = as.data.frame(Dgt)
	rownames(aux) = pdbs
	colnames(aux) = paste("Clusid",1:k,sep="")
	db$multi = aux
	return(db)
	#print(aux);readline()
}

view_svd_exp_old = function(tab,n=300){


	cols = colorRampPalette(c("brown","red", "yellow", "white"))(n = n)
	idcol = round(mapply(renorm_interval,tab[,7],MoreArgs=list(a=c(-1,5),b=c(1,n),inv=T)),0)
	col = cols[idcol]

	#par(mfrow=c(1,2))
	#plot(tab$Area,tab$Center,bg=col,pch=21)#;readline()
	#plot(tab$AreaR,tab$CenteR,bg=col,pch=21);readline()

	tam = dim(tab)[2]#;print(dim(tab))
	tab1 = matrix( as.numeric( as.matrix (tab[,8:tam])   )  ,ncol=tam-8+1)#;print(dim(tab1));readline()
	#print(head(tab1));readline()
	tam = dim(tab1)[2]
	#tab22 = scale(tab1[,2:2],center=T,scale=T)#;print(tab22);readline()
	tab1a = scale(tab1[,1:tam],center=T,scale=T)

	v = apply(tab1a,2,sum_nan)#;print(v);readline()
	tab1a = tab1a[,!v]
	#tam = dim(tab1a)[2]
	#tab1 = cbind(tab1[,1],tab1a)### acrescenta Dg
	#tab1 = cbind(tab22,tab1a)### scale only area
	tab1 = tab1a

	print(head(tab1));readline()

	tab1.svd = svd(tab1)
	us = tomatrix(tab1.svd$u) %*% diag(tab1.svd$d)
	tab1.svd$us = us
	
	#print(cbind(idcol,tab[,"Dg"]));readline()
	x=us[,1]
	y=us[,2]
	z=us[,3]
	w=us[,4]
	par(mfrow=c(2,2))
	plot(x,y,type="p",bg=col,pch=21)
	plot(x,z,type="p",bg=col,pch=21)
	tabf = filter(tab,Dg!=0)
	plot(tabf$Dg,tabf$AreaR)
	#plot(x,w,type="p",bg=col,pch=21)
	plot(tab1.svd$d)
	plot3d(x,y,z,size=1,col=col,aspect="iso",box=FALSE,type="s")
	x = list()
	x$tab = tab
	x$tabn = tab1
	x$svd = tab1.svd
	
	return(x)
	
	#print(head(tab1));readline()

	#tab1c = scale(tab1,center=T,scale=F)
	#tab1c = matrix(as.numeric(as.matrix(dbexp[[1]][[1]][1:5,7:9])),ncol=3)
	#tab1c.svd = svd(tab1c)
	#plot(tab1c.svd.

}















