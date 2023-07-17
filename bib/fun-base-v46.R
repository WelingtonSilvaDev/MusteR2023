##########################################################
# Biblioteca de funções básicas, geralmente para estudo de
# um complexo em separado
##########################################################

make_areas_from_acons = function(acons){

	areas.polar = round(sum(acons$polar[,2])/2,0)
	areas.apolar = round(sum(acons$apolar[,2])/2,0)
	areas.total = areas.polar + areas.apolar

	res = c(areas.apolar,areas.polar,areas.total)
	names(res) = c("areas.apolar","areas.polar","areas.total")

	return(res)

}

make_areas_from_sasas = function(sasas,polar_type,zero=-0.001,pos=T){

	
	sasas.diff = filter(sasas,diff<zero)

	sasas.apolar = filter(sasas.diff,pol=="a")
	if (polar_type=="all"){	
		sasas.polar = filter(sasas.diff,pol!="a")
	}else{
		f = sasas.diff$pol %in% polar_type
		sasas.polar = sasas.diff[f,]#;print(sasas.polar);readline()
	}

	if (pos){
		dasa.apolar = (-1)*round(sum(as.numeric(sasas.apolar$diff)),0)
		dasa.polar = (-1)*round(sum(as.numeric(sasas.polar$diff)),0)
	}else{
		dasa.apolar = round(sum(as.numeric(sasas.apolar$diff)),0)
		dasa.polar = round(sum(as.numeric(sasas.polar$diff)),0)
	}
	dasa.total = dasa.apolar + dasa.polar
	

	res = c(dasa.apolar,dasa.polar,dasa.total)
	names(res) = c("dasa.apolar","dasa.polar","dasa.total")

	return(res)

}

###NEW EXP###
make_areas_par = function(areas,classic=T,k=1){

	if (classic){
		f = areas$classic
		areas.apolar = areas[f,]
		areas.polar = areas[!f,]
	}else{
		areas.apolar = filter(areas,pol=="a")
		areas.polar = filter(areas,pol!="a")
	}
	areas.apolar = round(sum(as.numeric(areas.apolar$area))/2,0)
	areas.polar = round(sum(as.numeric(areas.polar$area))/2,0)
	areas.total = areas.apolar + areas.polar

	res = c(areas.apolar,areas.polar,areas.total)
	names(res) = c("area.apolar","area.polar","area.total")
	return(res)

}

###NEW EXP###
make_term_par = function(sasas,a=1.88,b=-1.09,c=log(298.15/385.15),zero=-0.001,classic=T,k=1,sep="_"){

	#print(head(sasas));readline()
	sasas.diff = filter(sasas,diff<zero)
	#print(head(sasas.diff));readline()

	if (classic){
		f = sasas.diff$classic
		#if(1){
		#atomlist = data.frame(t(data.frame(strsplit(sasas.diff$ids,"_"))),sasas.diff$pol)#;print(head(atomlist,200));readline()
		#fs = grepl("S",atomlist[,4])
		#apocar = c("a","m","c","i")
		#fa = atomlist[,6] %in% apocar 
		#f = fa&(!fs)
		sasas.apolar = sasas.diff[f,]#;print(sasas.apolar);readline()
		sasas.polar = sasas.diff[!f,]#;print(sasas.polar);readline()
		#}
		#sasas.apolar = filter(sasas.diff,pol=="a"|pol=="m"|pol=="c"|pol=="i")		
		#sasas.polar = filter(sasas.diff,pol=="p"|pol=="b"|pol=="o"|pol=="h")
	}else{
		sasas.apolar = filter(sasas.diff,pol=="a")
		sasas.polar = filter(sasas.diff,pol!="a")
	}
	#print(head(sasas.apolar));readline()

	dasa.apolar = round(sum(as.numeric(sasas.apolar$diff)),0)
	dasa.polar = round(sum(as.numeric(sasas.polar$diff)),0)#;print(dasa.apolar);print(dasa.polar);readline()

	dcp = round(a*dasa.apolar + b*dasa.polar,k)
	
	ds = round(dcp*c,0)

	res = c(ds,dcp,dasa.apolar,dasa.polar)
	names(res) = c("ds","dcp","dasa.apolar","dasa.polar")
	return(res) 

}

###NEW EXP###
make_sasa = function(pdb,atomlist,hetatom=F,sep="-",rm.h=T){

	#print(atomlist);readline()
	atomtab=data.frame(t(data.frame(strsplit(atomlist,sep)))[,4:5])#;print(atomtab);readline()
	if (rm.h){
		f = atomtab[,2] != "h"
		atomtab = atomtab[f,]
	}
	#print(head(pdb$atom));readline()
	rownames(atomtab)=NULL
	colnames(atomtab)=c("type","pol")
	if (rm.h){
		f = pdb$atom$elety != "H"
		pdb$atom = pdb$atom[f,] 
	}
	#;print(head(atomtab));readline()
	#print(head(pdb$atom));readline()
	if (!hetatom){
		atoms = filter(pdb$atom,type=="ATOM")
		atomtab = filter(atomtab,type=="A")
	}else{
		atoms = pdb$atom
	}
	#print(atoms[,1:5]);readline()
	chains = unique(atoms$chain)
	tam = length(chains)
	sasa.sep = data.frame()
	sasa.join = FreeSASA.diff(atoms)[,1:2]
	sasa.all = data.frame(sasa.join)
	
	for (i in 1:tam){
		
		atomc = filter(atoms,chain==chains[i])
		sasac = FreeSASA.diff(atomc)[,1:2]
		sasa.sep = rbind(sasa.sep,sasac)

	}
	#print(head(sasa.sep));readline()
	sasa.diff = sasa.all[,2]-sasa.sep[,2]
	#classic = 
	sasa.all = data.frame(sasa.all,sasa.sep[,2],sasa.diff,atomtab$pol)
	colnames(sasa.all)=c("ids","join","sep","diff","pol")

	newatomlist = data.frame(t(data.frame(strsplit(sasa.all$ids,"_"))),sasa.all$pol)#;print(head(atomlist,50));readline()
	f = new_to_classic_polarities(newatomlist[,4],newatomlist[,6])
	#fs = grepl("S",newatomlist[,4])
	#apocar = c("a","m","c","i")
	#fa = newatomlist[,6] %in% apocar 
	#f = fa&(!fs)
	sasa.all = data.frame(sasa.all,classic=f)
	#print(head(sasa.all));readline()
	return(sasa.all)	

	#print(chains);readline()

}

new_to_classic_polarities = function(atomlist,newpol,apocar = c("a","m","c","i")){

	fs = grepl("S",atomlist)
	fa = newpol %in% apocar 
	f = fa&(!fs)

	return(f)

}

make_acons = function(db,het="UNK",sep="_",polar_type="all"){
	
### TO DO: 
	atab = t(mapply(splitatom,db$element_name,MoreArgs = list(typename=2)))
	#;print(tail(atab));readline()
	atab[,2] = aa123(atab[,2])
	atab = atab[,c(2,3,1,4)]#;print(tail(atab));readline()
	f = atab[,1]!=het
	atab = atab[f,]
	pol = db$element_polarity[f]
	anames = apply(atab,1,paste,collapse=sep)
	names(anames)=NULL
	#print(head(atab));readline()
	#print(head(anames));readline()

	ar = db$A$a[f,f]

	fp = pol=="a"	
	as = apply(ar[fp,fp],1,sum)#;print(as);readline()
	acon.a = data.frame(anames[fp],round(as,0),pol[fp])
	colnames(acon.a)=c("ids","area","pol")

	if (polar_type=="all"){
		fp = !fp
	}else{
		fp = pol %in% polar_type
	}
	as = apply(ar[fp,fp],1,sum)#;print(as);readline()
	acon.p = data.frame(anames[fp],round(as,0),pol[fp])
	colnames(acon.p)=c("ids","area","pol")
	#acon = cbind(anames,round(as,k))#;print(acon);readline()
	#f = new_to_classic_polarities(atab[,4],pol)
	#acon = data.frame(anames,round(as,0),pol,f)
	#colnames(acon)=c("ids","area","pol","classic")
	#print(acon);readline()
	res = list()
	res$apolar = acon.a
	res$polar = acon.p

	return(res)

}

make_areas = function(db,pol,het="UNK",sep="_"){


	atab0 = t(mapply(splitatom,db$element_name,MoreArgs = list(typename=2)))
	#;print(tail(atab));readline()
	atab0[,2] = aa123(atab0[,2])
	atab0 = atab0[,c(2,3,1,4)]#;print(tail(atab));readline()
	f = atab0[,1]!=het

	### NEW ###
	atab1 = atab0[f,]
	pol = db$element_polarity[f]

	anames = apply(atab1,1,paste,collapse=sep)
	names(anames)=NULL
	#print(head(atab));readline()
	#print(head(anames));readline()
	
	as = apply(db$A$a[f,f],1,sum)#;print(as);readline()
	#print(dim(atab));print(length(pol));print(length(as));print(length(anames));readline()
	#acon = cbind(anames,round(as,k))#;print(acon);readline()
	fp = new_to_classic_polarities(atab1[,4],pol)
	acon1 = data.frame(anames,round(as,0),pol,fp)
	colnames(acon1)=c("ids","area","pol","classic")
	rownames(acon1)=NULL
	#print(head(acon));readline()
	### NEW ###

	### CLASSIC
	#fc = new_to_classic_polarities(atab[,4],pol)
	if (pol==2) fp = !fp
	fc = f&fp
	atab2 = atab0[fc,]
	anames = apply(atab2,1,paste,collapse=sep)
	pol = db$element_polarity[fc]
	as = apply(db$A$a[fc,fc],1,sum)
	acon2 = data.frame(anames,round(as,0),pol)
	colnames(acon2)=c("ids","area","pol")
	rownames(acon2)=NULL
	res = list()
	res$new = acon1
	res$classic = acon2
	### CLASSIC

	#return(acon)
	return(res)

}

###NEW EXP###
make_table_areas = function(sasa,db,sep="_",k=0,water="UNK"){

	atab = t(mapply(splitatom,db$element_name,MoreArgs = list(typename=2)))
	#;print(tail(atab));readline()
	atab[,2] = aa123(atab[,2])
	atab = atab[,c(2,3,1,4)]#;print(tail(atab));readline()
	f = atab[,1]!=water
	atab = atab[f,]
	anames = apply(atab,1,paste,collapse=sep)
	names(anames)=NULL
	#print(anames);readline()

	as = apply(db$A$a[f,f],1,sum)#;print(as);readline()
	#acon = cbind(anames,round(as,k))#;print(acon);readline()
	acon = data.frame(anames,round(as,k))#;print(acon);readline()

	stab = t(as.data.frame(strsplit(sasa[,1],sep)))[,1:4]
	snames = apply(stab,1,paste,collapse=sep)
	names(snames)=NULL
	rownames(stab)=NULL
	

	f = snames %in% anames
	sasa = sasa[f,]

	#print(sasa);print(acon);readline()

	astab = cbind(sasa,acon)
			
	return(astab)

	#return(atab)		
	#print(head(anames));print(head(snames))
	#print(head(sasa));print(head(acon))
	#print(tail(sasa));print(tail(acon))
	#print(dim(sasa));print(dim(acon))
	#;readline()

}

###NEW EXP###
#build_pdb_colname = function(pdb){

	



#}

is.splitable = function(x,s,n=2){

	y = unlist(strsplit(x,s))#;print(typeof(y));readline()
	if (length(y)==n){
		return(y)
	}else{
		return(NULL)
	}	

}

char2double = function(x,s,r=c(">","<","~"," ")){

	y = is.splitable(x,s)#;print(y)
	if (!is.null(y)){
		tam = length(r)
		ya = y[1]
		yb = y[2]
		#noise = F
		for (i in 1:tam){
			yi = is.splitable(ya[1],r[i])#;print(typeof(yi));readline()
			if (!is.null(yi)){
				res = as.double(paste(yi[2],"e",yb,sep=""))
				return(res)
				#noise = T
				#break;
			}
		}
		#if (!noise)	
		res = as.double(paste(ya,"e",yb,sep=""))
		return(res)
		#print(res);readline()
		#ya = is.splitable(y[1],r[1])#;print(ya)
		#if (is.null(ya)){
		#	yb = is.splitable(y[1],r[2])
		#	if (is.null(yb)){
		#		res = as.double(paste(y[1],"e",y[2],sep=""))
		#		return(res)
		#	}else{
		#		res = as.double(paste(yb[2],"e",y[2],sep=""))
		#		return(res)
		#	}
		#}else{
		#	res = as.double(paste(ya[2],"e",y[2],sep=""))
		#	return(res)
		#}
	}else{
		return(NA)
	}
}

ToDouble = function(m,n,s=c("E","e"),verbose=F){

	for (i in n){
		x = as.vector(unlist(m[,i]))
		#x = ">1.5E-04"
		tam = length(x)
		z = c()
		for (j in 1:tam){#;print(x[j])
			if (typeof(x[j])=="character"){
				z[j] = char2double(x[j],s[1])
				if (is.na(z[j])){
					z[j] = char2double(x[j],s[2])
					if (is.na(z[j])){
						if (verbose) print(paste("i:",i,"j: ",j,"=",x[j],"not double"))
					}
				}
				#m[j,i] = z[j]
			}
			else if (typeof(x[j])=="double"){
				z[j] = x[j]

			} else {
				z[j] = m[j,i]
				if (verbose) print(paste("i:",i,"j: ",j,"=",x[j],"indefined type"))
			}
			#;print(z[j]);readline()	
		}
		m[,i] = z
	}
	return(m)
}

filter_skempi = function(db){

	db$m=ToDouble(db$m,n=c(7:10,15:26))
	db$m[,14] = as.numeric(str_extract(db$m[,14],"[[:digit:]]+")) #extract only numbers
	#db$m[,14] = unlist(db$m[,14])
	#print(typeof(db$m[1,14]))
	#rownames(db$m) = NULL
	return(db)

}


make_skempi_energy_table = function(db,r=8.314,k = 0.000239006,d=2){

	#tam = dim(db$m)[1]
	m = db$m[,c(1:3,11:13)]
	#for (i in 1:tam){
	dgmut1 = round((r*as.numeric(db$m[,14])*log(as.numeric(db$m[,8])))*k,d)
	dgwil1 = round((r*as.numeric(db$m[,14])*log(as.numeric(db$m[,10])))*k,d)
	ddg1 = round(as.numeric(dgmut1) - as.numeric(dgwil1),d)

	dgmut2 = round(as.numeric(db$m[,23])-as.numeric(db$m[,14])*as.numeric(db$m[,25])/1000,d)
	dgwil2 = round(as.numeric(db$m[,24])-as.numeric(db$m[,14])*as.numeric(db$m[,26])/1000,d)
	ddg2 = round(as.numeric(dgmut2) - as.numeric(dgwil2),d)

	m = cbind(m,dgmut1,dgwil1,ddg1,dgmut2,dgwil2,ddg2)

	colnames(m) = c("PDB","MUT_PDB","MUT_SKE","REF","PROT1","PROT2","DGMUT1","DGWIL1","DDG1","DGMUT2","DGWIL2","DDG2")

	db$tab = m

	return(db)

	#}
	

}


filter_single_mutations = function(db,mut="A",ncol=2,pattern="[[:digit:]]+"){

	tab = db$tab

	tam = dim(db$tab)[1]
	ids = c()

	for (i in 1:tam){#print(db$tab[i,ncol]);readline()
		x = unlist(strsplit(as.character(db$tab[i,ncol]),pattern));
		if (length(x)==2){
			if (x[2]==mut){
				ids = c(ids,i)
			}
		}
	}
	return(db$tab[ids,])
}


select_first_neighboors = function(db,ids,pamid,clusid){


	res = c()
	for (i in ids){

		x = ego(db$g,2,i)
		x=as.numeric(x[[1]])
		f=db$exp$pam[[pamid]]$clustering[x]==clusid
		x=x[f]
		y=db$rownames[x]#;print(y[1:10])
		f=!grepl("[^0-9]", enc2utf8(y))
		#x=x[f]
		y=y[f]
		res = c(res,y)
	}
	return(unique(res))
}


supercluster_graph_analysis = function(db,id,clusid,textcloud="",cloudlim=300,tamnodes=100,coord=NULL,graphplot=F,cloudplot=F,colortype="Reds"){

	aux = list()
	id=which(db$exp$pam[[id]]$clustering==clusid)
	names = db$colnames[id]
	aux$a = db$A$a[id,id]
	aux$n = db$matrixname
	aux$g = graph.adjacency(aux$a,weighted=T,mode="undirected",diag=F)
	dg = degree(aux$g)
	if(0){
		id=dg!=0
		names=names[id]
		aux$a = aux$a[id,id]
		#aux$n = db[[1]]$matrixname
		aux$g = graph.adjacency(aux$a,weighted=T,mode="undirected",diag=F)
	}
	#auxd = diag(aux$a)
	if (graphplot){
		auxd=tamnodes
		#auxyz = db[[1]]$exp$superclus[[i]]$geomc
		#coord = NULL
		#vlabel = 1:dim(aux$a)[1]
		vlabel = names
		elabel = E(aux$g)$weight
		dg = degree(aux$g)#;print(dg)
		colorid = as.factor(round(mapply(renorm_interval,dg,MoreArgs=list(b=c(1,9),a=c(0,max(dg)))),0))#;print(colorid)#;readline()
		#vbre = brewer.pal(n = 9, name = 'YlOrRd')
		vbre = brewer.pal(n = 9, name = colortype)
		levels(colorid) = vbre[as.numeric(levels(colorid))]#;print(colorid)#;readline()
		vcolor = as.character(colorid)
		sizef = 1
		xlim=NULL
		ylim=NULL
		sub =" "
		graph_plot(aux,labcex=1.0,size=sizef*sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,coordi=coord, ewidth=0.2,ecolor="gray75",lcolor="yellow",xlim=xlim,ylim=ylim,sub=sub)
	}

#	at = articles.title
#	a = db[[1]]$A$a
#	rn = db[[1]]$rownames
#	cn = db[[1]]$colnames
#	pamc = db[[1]]$exp$pam
	if (cloudplot){
		textcloud = textcloud[1:cloudlim,]
		#textcloud[,1] = exp(2*log(sqrt(textcloud[,1])))
		print("tagcloud...")
		colorid = as.factor(round(mapply(renorm_interval,textcloud[,1],MoreArgs=list(b=c(3,9),a=c(min(textcloud[,1]),max(textcloud[,1])))),0))
		#vbre = brewer.pal(n = 9, name = 'YlOrBr')
		#vbre = brewer.pal(n = 9, name = 'Reds')
		vbre = brewer.pal(n = 9, name = colortype)
		levels(colorid) = vbre[as.numeric(levels(colorid))]
		tagcolor = as.character(colorid)
		#tagcloud(textcloud[,2],textcloud[,1],col=tagcolor)
		wordcloud(textcloud[,2], textcloud[,1], max.words=300,col=tagcolor,ordered.colors=T,random.order=F)
	}

}

which.pam = function(pam,i,j){

	which(pam[[i]]$clustering == j)

}

title2text = function(db,id,clusid,titles,minchar=3){

	rn = db$rownames
	at = titles

	x = rn[which.pam(db$exp$pam,id,clusid)]
	#print(x)
	f=!grepl("[^0-9]", enc2utf8(x))#;readline()
	y = as.numeric(x[f])
	#print(at[y])
	text = c()
	tam = length(y)

	print("concat...")
	for (j in 1:tam){
	#for (j in 3){
		text = paste0(text,enc2utf8(at[y[j]]))
	}
	text = tolower(text)
	text = unlist(strsplit(text, "\\W+"))
	textab = data.frame(table(text))#;print(textab)
	temp = textab[,1]
	textab[,1] = textab[,2]
	textab[,2] = temp
	print("nchar...")
	textid = as.numeric(mapply(nchar,as.character(textab[,2])))#;print(textid)
	textab = textab[textid>minchar,]#;print(textab)
	print("order...")
	remwords = c("that","with","from","this","into","while","were","other","some","when","here","both","such","been","have","they")
	remwords = c(remwords,"than","many","much","1021","1007","1016","most","their","which","over","these","those","also","about")
	remwords = c(remwords,"more","within","after","under","well","1186","1109","1371","however","further","each","without","major")
	remwords = c(remwords,"each","early","where","large","high","between","among","several","will","going","using","used","1093","1092")
	remwords = c(remwords,"various","only","next","small","big","them","highly","1002")
	textab = remove_words(textab,remwords,ord=F)#;readline()
	o = order(textab[,1],decreasing=T)#;print(textab)
	textab = textab[o,]
	return(textab)
	#textcloud = textab[1:300,]
	#textcloud[,1] = textcloud[,1]


}

plot_cluster_quality = function(db,n,v=NULL){

	x=c();for(i in 1:n) x=c(x,db$exp$superclus[[i]]$accall[1])
	plot(x,type="l")
	y=c();for(i in 1:n) y=c(y,mean(db$exp$superclus[[i]]$accbyclass[,1]))
	lines(y,lty=2)
	z=c(1);for(i in 2:n) z=c(z,db$exp$pam[[i]]$silinfo$avg.width)
	lines(z,lty=3)
	w=c();for(i in 1:n) w=c(w,median(db$exp$superclus[[i]]$accbyclass[,1]))
	lines(w,lty=4)
	if (!is.null(v))abline(v=v,lty=5)

}

make_tab_cluster_decoy = function(db,id,isca){

	#id = c(12,29,125,180)
	#isca = c("3972","200","1016","1049","2595","370","1142")

	out = list()
	k = 1
	for (j in id){
		rowid = c()
		clusid = c()
		for (i in isca){#print(i)
			s = paste("^",i,"$", sep="")
			g.id = grep(s,enc2utf8(db$rownames))#;print(g.id);readline()
			if (length(g.id)==1){
				rowid = c(rowid,g.id)
				clusid = c(clusid,db$exp$pam[[j]]$clustering[g.id])
			}else{
				print("error")
			}
		}
		#print(clusid);readline()
		out[[k]] = clusid
		k=k+1
	}
	#print(out);readline()
	out = data.frame(out)
	colnames(out)=id
	res = cbind(isca,rowid,out)
	return(res)

}

klogistic = function(x,xo=0,k=1,L=1,inv=1){

	aux = L/( 1+exp( -k*(x-xo) ) )

	if (inv){
		return(1-aux)
	}else{
		return(aux)
	}

}


make_stat_score = function(dball,w){

	tam = length(dball$aligndata)
	auxd = c()
	auxdn = c()
	auxgv = c()
	auxnv = c()
	auxwv = c()
	auxsv = c()
	auxsc = c()

	for (i in 1:tam){
		auxd=c(auxd,dball$aligndata[[i]]$score$dist)
		auxdn=c(auxdn,dball$aligndata[[i]]$score$distn*w[1])
		auxgv=c(auxgv,dball$aligndata[[i]]$score$degv*w[2])
		auxnv=c(auxnv,dball$aligndata[[i]]$score$nodv*w[3])
		auxwv=c(auxwv,dball$aligndata[[i]]$score$weiv*w[4])
		auxsv=c(auxsv,dball$aligndata[[i]]$score$scov)
		auxsc=c(auxsc,dball$aligndata[[i]]$score$scor)
	}
	
	dball$stat = list()
	dball$stat$dist = auxd
 	dball$stat$distn = auxdn
	dball$stat$degv = auxgv
	dball$stat$nodv = auxnv
	dball$stat$weiv = auxwv
	dball$stat$scov = auxsv
	dball$stat$scor = auxsc
	return(dball)

}

adjust_ca_in_xyz = function(xyz,div=9){

	auxyz = xyz

	tam = dim(auxyz)[1]

	auxi = seq(1,tam,by=div)
	#print(auxi)

	for (i in auxi){
		temp = auxyz[i,]
		auxyz[i,] = auxyz[i+1,]
		auxyz[i+1,] = temp
	}
	return(auxyz)

}

matrix2pdb = function(xyz,remark="",type=2,div=9){

#REMARK 500 STANDARD TABLE:
#ATOM      1  CA  XXX E  16      25.182  49.029  12.446  1.00  3.00           C

	atom="ATOM"
	ca="CA"
	chain="E"
	res="XXX"
	c="C"
	r=3
	rem="REMARK 500 "
	auxs=c()
	auxr=c()

	tam = dim(xyz)[[1]]
	auxr = paste(rem,remark,sep="")

	if (type==2){
		#print(xyz)
		xyz = adjust_ca_in_xyz(xyz)
		#print(xyz);readline()
		ca = c("N","CA","C","O","CB","CD","CG","OE1","NE2")
		ca = rep(ca,tam/div)
		c = c("N","C","C","O","C","C","C","O","N")
		c = rep(c,tam/div)
		kr = 1
	}
	#print(ca);print(c);readline()
	for (k in 1:tam){
		if (type==1){
			#auxs = sprintf("%s%6.0f%4s%5s%2s%4.0f%12.3f%8.3f%8.3f%6.2f%6.2f%12s",atom,k,ca,res,chain,k,xyz[k,1],xyz[k,2],xyz[k,3],1,r,c)
			auxs = sprintf("%s%6.0f%4s%5s%2s%4.0f%12.3f%8.3f%8.3f%6.2f%6.2f%12s",atom,k,ca,res,chain,k,xyz[k,1],xyz[k,2],xyz[k,3],1,r,c)
		}
		if (type==2){
			auxs = sprintf("%s%7.0f  %-4s%-4s%-2s%3.0f%12.3f%8.3f%8.3f%6.2f%6.2f%12s",atom,k,ca[k],res,chain,kr,xyz[k,1],xyz[k,2],xyz[k,3],1,r,c[k])
			if (k%%div==0){
				kr=kr+1
			}
		}
		auxr = rbind(auxr,auxs)
		#print(auxr);readline()
	}
#	auxs = sprintf("%s%6.0f%4s%5s%2s%4.0f%12.3f%8.3f%8.3f%6.2f%6.2f%12s",atom,1,ca,res,chain,1,25.182,49.029,12.45,1,3,c)
#	auxs = rbind(auxs,sprintf("%s%6.0f%4s%5s%2s%4.0f%12.3f%8.3f%8.3f%6.2f%6.2f%12s",atom,100,ca,res,chain,100,5.182,0.029,-20.445,1,3,c))
	#print(auxr);readline()
	rownames(auxr)=NULL
	return(auxr)
}

increase_matrix_to_min_random = function(xyz,minxyz,min=-1.0,max=1.0){

	tam = dim(xyz)[1]
	auxr = xyz
	tamr = tam
	j=1
	while(tamr<minxyz){
			for (i in j:tam){
				auxu = runif(3,min=min,max=max)
				auxr = rbind(auxr,auxr[i,]+auxu)
			}
			#print(auxr);readline()
			tamr = tamr + tam			
	}
	return(auxr)


}

incubic = function(xyz,l=1,n=3){


	auxm = ntobin(n)
	#print(xyz);print(auxm)
	for (i in 1:n){
		auxi = auxm[,i]==0
		#print(auxi)
		auxm[auxi,i] = xyz[i]-l
		auxm[!auxi,i] = xyz[i]+l
	}
	auxm = rbind(xyz,auxm)
	rownames(auxm)=NULL	
	#print(auxm);readline()
	
	return(auxm)
	
}

increase_matrix_to_cubic = function(xyz,minxyz,l=1){

	tam = dim(xyz)[1]
	auxr = c()
	tamr = tam
	j=1
	#while(tamr<minxyz){
	#print(xyz)
	for (i in j:tam){
		auxi = incubic(xyz[i,])
		auxr = rbind(auxr,auxi)
		#auxu = runif(3,min=min,max=max)
		#auxr = rbind(auxr,auxr[i,]+auxu)
	}
		#print(auxr);readline()
			#tamr = tamr + 8

	if (dim(auxr)[1]<minxyz){
		print(paste("WARNING: matrix length less than minxyz"))
	}
	#print(auxr);readline()			
	return(auxr)
	

}

new_make_group_table = function(db,dball,type=1,minclusid=3,minpdb=16,verbose=T,seps="_",workdir="Pymol/",sufix=".pdb",filename=NULL){

	if (is.null(dball$info$fakename)){
		auxinfo = dball$info
		#print(auxinfo);readline()
		tami = dim(auxinfo)[1]
		#print(tami)
		auxfilepdb = c()

		if (verbose) print(paste("Generating new pdbs names "))

		for (i in 1:tami){
			#tamj = length(db[[i]]$exp$superclus)
			j = as.numeric(auxinfo[i,3])#;j=1
			k = as.numeric(auxinfo[i,5])#;print(j);print(i);print(k);print(minclusid);readline()
			#print(is.data.frame(as.numeric(j)))
			#print(j<minclusid);readline()
			if ((j<minclusid)&(j>0)){
				#print("0k1")
			#	if (verbose) print(paste("WARNING: PDB",db[[k]]$pdbname,"has superclus length",j,"less than minclusid",minclusid))
				auxfilepdb = c(auxfilepdb,"none")#;print("0k")
			}else{
			#for (j in tamj){
				#print("0k2")
				auxname = paste(db[[k]]$matrixname,seps,k,seps,j,sep="")
			#	auxfilefake = paste(workdir,auxname,sufix,sep="")
				#print(auxfilefake);readline()
			#	if (file.exists(auxfilefake)){
			#		if (verbose) print(paste("WARNING: fake PDB",auxfilefake,"detected. It will not be overwritten"))
			#	}else{ 			
			#		auxyz = db[[k]]$exp$superclus[[j]]$geomc
					#print(auxyz)
			#		if (type==1){
			#			auxyz = increase_matrix_to_min_random(auxyz,minxyz=minpdb)
			#		}
			#		if (type==2){
			#			auxyz = increase_matrix_to_cubic(auxyz,minxyz=minpdb)
			#		}
					#print(auxyz);readline()		
					#auxname = paste(db[[k]]$matrixname,seps,k,seps,j,sep="")
			#		auxpdb = matrix2pdb(auxyz,remark=toupper(auxname))
					#auxpdbname = paste(auxname,sufix,sep="")
					#auxpdbname = auxname
			#		write.table(file=auxfilefake,auxpdb,quote=F,row.names=F,col.names=F)
			#	}
				auxfilepdb = c(auxfilepdb,auxname)
				#print(auxpdb);readline()
			}
		}
		#print(auxfilepdb);readline()
		auxcol = colnames(auxinfo)
		auxinfo = data.frame(auxinfo,auxfilepdb)
		colnames(auxinfo)=c(auxcol,"fakename")
		dball$info = auxinfo#;print(auxinfo);readline()
		}
	else{
		print("column fakename already exists; so doing nothing with fakepdbs")
	}
	#print(dball$info);print(tami)
	auxv = dball$info[,6]!="none"
	if (sum(auxv)==0){
		#print(typeof(as.data.frame(dball$info)))
		dball$info$fakename = as.character(dball$info$fakename)#;print(dball$info);readline()
		j = as.numeric(auxinfo[tami,3])
		k = as.numeric(auxinfo[tami,5])
		auxname = paste(db[[k]]$matrixname,seps,k,seps,j,sep="")
		dball$info[tami,6]=auxname#;print(dball$info);readline()
		auxv = dball$info[,6]!="none"
	}
	#print(dball$info);print(auxv)
	auxinf = dball$info[auxv,]#;print(auxinf[,6])
	auxr = pdbsplit(auxinf[,6])#;print(auxinf[,6]);readline()
	auxr = data.frame(auxinf[,c(1:2)],auxr,auxinf[,c(6)])
	colnames(auxr) = c("gacc","lacc","pdb","chain1","chain2","method","polarity","id","nclus","fakename")
	auxnull = rep(0,dim(auxr)[1])
	f = auxr$nclus == 0
	auxnull[f] = 1
	auxr = data.frame(auxr,empty=auxnull)
	#print(auxr);readline()
	dball$group = auxr
	if (!is.null(filename)){
		filename = paste(outdir,toupper(filename),sep,sufix,sep="")
		write.csv(auxr,filename)
	}
	#print(dball);readline()
	
	return(dball)

}






whichmin = function(v){

	id = which(v==min(v))[1]
	res = rep(0,length(v))
	res[id]=v[id]

	return(res)

}

retvi = function(v,i,k){

	return(v[i]==k)

}

is_valid_change_cel = function(i,m1,m2,zero=0.001){

	msg = "WARNING: something strange in make matrix symmetrics:"

	id = which(m1[i,]>zero)#;print(paste("valid:",i,id))
	id.len = length(id)
	if (id.len==1){
		if (m1[i,id]!=m2[id,i]){
			return(id)
		}else{
			return(0)
		}
	}else if (id.len==0){
		return(1)
	}else{
		print(paste(msg,1,"a"))
	}

}


make_matrix_symmetrics = function(mx,my,an,zero=0.001){

	#tam = dim(an)[1]
	msg = "WARNING: something strange in make matrix symmetrics:"
	yes = 0
	no = 0
	#print("make")
	while(length(an)){
		i=an[[1]][1]
		j=an[[1]][2]
		k=an[[1]][3]#;print(mx);print(my);readline()
		if (k==1){
			id = is_valid_change_cel(i=j,m1=my,m2=mx)
			if(id){
				my[j,i] = mx[i,j]
				my[j,id] = 0
				an[[1]] = NULL
				#yes = yes + 1
			}else{
				mx[i,j] = 0
				an[[1]] = NULL
				#no = no + 1
			}
		}else
		if (k==2){
			id = is_valid_change_cel(i=j,m1=mx,m2=my)			
			if(id){
				mx[j,i] = my[i,j]
				mx[j,id] = 0
				an[[1]] = NULL
				#yes = yes + 1				
			}else{
				my[i,j] = 0
				an[[1]] = NULL
				#no = no + 1
			}			
		}else{
			print(paste(msg,3))#;readline()
		}
	
	}
	res = list()
	res$mx = mx
	res$my = my
	#if(yes==no){
	#	res$go = T
	#}else{
	#	res$go = F
	#}
	return(res)

}


#ATENÇÃO: mx,my,mz devem estar arrendodadas!

asymetric_index = function(v,mx,my,zero=0.001){

	if (mx[v[1],v[2]]!=my[v[2],v[1]]){
		if ( (mx[v[1],v[2]]>zero)&(my[v[2],v[1]]>zero) ){
			v = c(v,0)
		}else if (mx[v[1],v[2]]>zero){
			v = c(v,1)
		}else{
			v = c(c(v[2],v[1]),2)
		}
		return(v)
	}
}

find_asymetric_nodes = function(mx,my,zero=0.001){

	#print("o")

	#tam = dim(mx)[1]

	#ids = which(diag(tam)<2,arr.ind=T);print(ids)#;readline()
	ids = which(!(mx<0),arr.ind=T)#;print(ids)#;readline()
	res = apply(ids,1,asymetric_index,mx,my)#;print(res)
	res = res[which(!sapply(res,is.null))]#;print(res);readline()
	
	return(res)
	
}

pasterow = function(v,inv=T,sep="-"){
	
	if (inv){
		return(paste(v[2],v[1],sep=sep))
	}else{
		return(paste(v[1],v[2],sep=sep))
	}

}

find_nearest_nodes = function(xyz1,xyz2,k=3,zero=0.001,verbose=F,force=F){

	msg = "WARNING: something strange in make matrix symmetrics:"

	ml = list()
	#print(xyz1);print(xyz2)
	d12 = round(dist.xyz(xyz1,xyz2),k)#;print(d12)
	d21 = round(dist.xyz(xyz2,xyz1),k)#;print(d21)
	d12[d12<=zero]=2*zero#;print(mx)  ### ATENCAO PARA O <= 
	d21[d21<=zero]=2*zero#;print(my)  ### ATENCAO PARA O <=

	#print(d12);print(d21);readline()
	d12.diag = diag(d12)
	if (all(d12.diag<=zero)){  ### ATENCAO PARA O <=
		dx = d12.diag
		n = which(d12<zero,arr.ind=T)#;print(cbind(n,dx));readline()
		n = apply(n,1,pasterow)
		names(dx) = n
		#print(dx);readline()
		return(dx)
	}

	mx=t(tomatrix(apply(d12,1,whichmin)))#;print(mx)#;readline()
	my=t(tomatrix(apply(d21,1,whichmin)))#;print(my);readline()


	#v = my!=0
	#mz = mx
	#mz[v] = my[v]#;print(mz)#;readline()

	if (0){
		gx=graph.adjacency(mx,weighted=TRUE,mode="directed",diag=T)
		gy=graph.adjacency(my,weighted=TRUE,mode="directed",diag=T)
		x11()
		plot(gx);plot(gy,add=T,vertex.color="lightblue")#;readline()
	}
	nogo = 0
	an = find_asymetric_nodes(mx,my)#;print(an)#;readline()
	if (!is.null(an)){
		ml = make_matrix_symmetrics(mx,my,an)#;print(ml)#;readline()		
		mx = ml$mx
		my = ml$my
		#go = ml$go
		an = find_asymetric_nodes(mx,my)#;print(an);readline()
		if (!is.null(an)){
			print(paste(msg,4))
		}else{
			#print(paste("Done!"))
			if (!force){
				nogo = sum(mapply("<",apply(mx,1,sum),zero))#;print(nogo);readline()
			}
			if (0){
				gx=graph.adjacency(mx,weighted=TRUE,mode="directed",diag=T)
				gy=graph.adjacency(my,weighted=TRUE,mode="directed",diag=T)
				x11()
				plot(gx);plot(gy,add=T,vertex.color="lightblue")#;readline()
			}
		}				
	}

	#print(mx);print(my);readline()
	if (!nogo){
		mx=t(mx)#;print(mx)
		dx = mx[mx>zero]#;print(dx);readline()
		n = which(mx>zero,arr.ind=T)#;print(cbind(n,dx));readline()
		n = apply(n,1,pasterow)
		names(dx) = n
		#print(dx);readline()	
		return(dx)
	}else{
		if (verbose) print(paste("WARNING: it was not possible to find a symmetric nearest distance vector between the subgraphs"))
		return(NULL)
	}
}

expand_to_nearest_global_nodes = function(dist,xyz1,xyz2,tid,sid,k=3,sep="-",retract=T,zero=0.001){

	msg = "WARNING: something strange in make matrix symmetrics:"

	ml = list()

	#print(xyz1);print(xyz2)

	tam1 = dim(xyz1)[1]
	tam2 = dim(xyz2)[1]

	tida = 1:tam1
	sida = 1:tam2

	tid1 = which(!(tida %in% tid))#;print(tid1)
	sid1 = which(!(sida %in% sid))#;print(sid1)#;readline()

	#id0 = sort(c(tid1,sid1));print(id0)

	if (length(tid1)&length(sid1)){
		#d12 = round(dist.xyz(xyz1[tid1,],xyz2[sid1,]),k);print(d12)
		#d21 = round(dist.xyz(xyz2[sid1,],xyz1[tid1,]),k);print(d21)#;readline()
		#mx=t(tomatrix(apply(d12,1,whichmin)));print(mx)#;readline()
		#my=t(tomatrix(apply(d21,1,whichmin)));print(my)#;readline()
		
		#print(xyz1[tid1,])
		#print(xyz2[sid1,])

		dista=find_nearest_nodes(xyz1[tid1,],xyz2[sid1,],force=T)#;print(dista)#;print(paste("dista:",dista));print(dista)#;readline()
		labels = fix_labels(labels=names(dista),id1=tid1,id2=sid1)#;print(labels)#;readline()
		labels_old = fix_labels_old(labels=names(dista),id1=tid1,id2=sid1)
		#print(labels);print(labels_old);readline()	


		names(dista)=labels#;print(dista);readline()

		d.r = c(dist,dista)

		d.m=matrix(as.numeric(unlist(strsplit(names(d.r),sep))),ncol=2,byrow=T)

		d.m = cbind(d.m,d.r)
	
		#print(d.m);readline()

		d.id = which(!(tida %in% d.m[,1]))#;print(d.id)
		d.id.len = length(d.id)
		if (d.id.len){
			d.id = t(tomatrix(d.id))
			d.id = cbind(d.id,tomatrix(mapply(rep,c(0,1/0),d.id.len)))#;print(d.id.len);print(d.id);print(tomatrix(mapply(rep,c(0,1/0),d.id.len)))
			d.m = rbind(d.m,d.id)#;print(d.m);readline()
		}
		#if (length(d.id)){
		#	d.m = rbind(d.m,c(d.id,0,1/0));print(d.m);readline()
		#}
		o = order(d.m[,1])
		d.m = d.m[o,]
	
	
		#;print(d.m);readline()

		d.d = d.m[,3]#;print(d.d);
		labels = apply(d.m[,1:2],1,pasterow,inv=F)
		names(d.d) = labels

		return(d.d);
	}else{
		return(NULL)
	}
	
	#print(d.d)

	#;readline()



	#d12 = round(dist.xyz(xyz1,xyz2),k);print(d12)
	#d21 = round(dist.xyz(xyz2,xyz1),k);print(d21);readline()

	#d12.0 = round(dist.xyz(xyz1[id0,],xyz2[id0,]),k)
	#d21.0 = round(dist.xyz(xyz2[id0,],xyz1[id0,]),k)

	#print(d12);print(d21);print(d12.0);print(d21.0);print(tid);print(sid);print(dist);readline()


}

verify_edge_superimpostion = function(dist,e1,e2,xyz1,xyz2,sep="-"){

	#print(dist);print(e1);print(e2);readline()
	auxc=as.numeric(unlist(strsplit(dist,sep)))
	#print(auxc)
	x = matrix(auxc,ncol=2,byrow=T)
	#print(x)#;readline()
	#tam1 = dim(e1)[1]
	tam1 = length(dist)
	rr = rep(0,tam1)
	rw = rr
	rc = rr
	#auxr=list()
	
	for (i in 1:tam1){
		xi = x[i,1]#;print(paste("i:",i,"xi:",xi))
		#v = e1[,1]==xi;print(paste("i:",i,"xi:",xi,"v:",v))
		#if (all(!v)) break;
		#me1 = e1[v,]
		me1 = e1[e1[,1]==xi,] #conexoes de xi
		#print(me1)#;readline()#;print(is.matrix(me1))
		if (!is.matrix(me1)){
			me1 = t(as.matrix(me1))
		}
		#print(me1)#;readline()
		tame1 = dim(me1)[1]
		if (tame1>0){
			for (j in 1:tame1){
				#print(paste("--- j: ",j))
				ok = F
				c1 = x[me1[j,1],2]#;print(c1)#;readline()#;print(c1)
				c2 = x[me1[j,2],2]#;print(paste(c1,"-",c2))
				if (c2){
					#c2 = x[me1[j,2],2];print(me1[j,2])#;print(c2);readline()

					c12 = c(c1,c2)#;print(c12);readline()
					me2 = e2[,1:2]
					if (!is.matrix(me2)){
						me2 = t(as.matrix(me2))
					}
					e2id1 = row.match(c12,me2)#;print(e2id1);readline()
					if (!is.na(e2id1)){
						ok=T
						e2id = e2id1
					}
					c21 = c(c2,c1)
					e2id2 = row.match(c21,me2)
					if (!is.na(e2id2)){
						ok=T
						e2id = e2id2
					}#;print(ok);readline()
					#if (!( is.na(row.match(c12,e2[,1:2]))&is.na(row.match(c21,e2[,1:2])) )){ #de Morgan...
					if (ok){
						#print(c12);readline()
						#n = n + 1
						#print(e2)
						rww = min(me1[j,3],e2[e2id,3])/max(me1[j,3],e2[e2id,3])
						#print(me1[j,]);print(e2[e2id,]);print(rww);readline()
						v = xyz1[me1[j,1],]-xyz1[me1[j,2],]
						w = xyz2[e2[e2id,1],]-xyz2[e2[e2id,2],]
						xyz.cos = cos_vec(v,w)
						#print(xyz1);print(xyz2);print(v);print(w);print(xyz.cos);readline()
						#rww = min(me1[j,3],e2[e2id,3])/max(me1[j,3],e2[e2id,3])
						rr[me1[j,1]]=rr[me1[j,1]]+1
						rr[me1[j,2]]=rr[me1[j,2]]+1
						rw[me1[j,1]]=rw[me1[j,1]]+rww
						rw[me1[j,2]]=rw[me1[j,2]]+rww
						rc[me1[j,1]]=rc[me1[j,1]]+xyz.cos
						rc[me1[j,2]]=rc[me1[j,2]]+xyz.cos
						#print(rr);print(rw);readline()
						#print(me1[j,1]);print(me1[j,2]);print(rr);readline()
					}
					#print(c12);print(c21);print(rr);readline()
				}
			}
		}

	}
	auxr=list()
	auxr$rr = rr
	auxr$rw = rw
	auxr$rc = rc 
	#print(rr);readline()
	#return(rr)
	#print(auxr);readline()
	return(auxr)
}


maxmin = function(v){

	#print(v)#;readline()
	vmax = max(v)#;print(vmax)
	vmin = min(v)#;print(vmin)
	res = 1-(vmax-vmin)#;print(res);readline()
	return(res)

}

vector_matching_pair = function(v1,v2,type=1){

	if (type==1){
		rr = v1/v2
		rr[is.nan(rr)]=1
		rr[rr>1]=1/rr[rr>1]

		return(rr)
	}
	if (type==2){

		v = rbind(v1,v2)
		x = apply(v,2,maxmin)


	}

}

make_score_matching_pair = function(dist,v1,v2,sort=F,sep="-",type=1){


	v1a = v1
	v2a = v2
	#print(dist);print(v1);print(v2)#;readline()
	auxc=as.numeric(unlist(strsplit(dist,sep)))
	x = matrix(auxc,ncol=2,byrow=T)#;print(x);print(v1);print(v2);readline()
	v = x[,2]!=0
	#x = x[v,]
	#v1 = v1[v]
	#v2 = v2[v]
	#da1 = diag(a1)
	#da2 = diag(a2)	

	tam1 = length(dist)

	#names(da1) = x[
	ea1 = rep(0,tam1)
	ea2 = rep(1,tam1)
	ea1[v] = v1[x[v,1]]
	ea2[v] = v2[x[v,2]]

	#print(x);print(ea1);print(ea2);readline()

	names(ea1) = x[,1]
	names(ea2) = x[,2]
	
	
	sa1 = ea1
	sa2 = ea2
	#print(sa1);print(sa2);readline()
	if (sort){
		sa1 = sort(ea1)
		sa2 = sort(ea2)
	}
	#print(auxc);print(x);print(v1);print(v2);print(ea1);print(ea2);print(sa1);print(sa2);readline()
	
	#print(sa1);print(sa2)
	#rr = sa1/sa2

	ro = order(names(sa1),decreasing=F)
	sa1 = sa1[ro]
	sa2 = sa2[ro]
	
	rr = vector_matching_pair(sa1,sa2,type)

	#rr = sa1/sa2#;print(rr)
	#rr[is.nan(rr)]=1 ### CUIDADO! Verificar consequencias no score...
	#print(names(sa1));
	#print(sa1);print(sa2);readline()
	names(rr)=paste(names(sa1),names(sa2),sep=sep)
	#print(rr);readline()
	#rr[rr>1]=1/rr[rr>1]
	#print(rr);readline()
	return(rr)
}

pull_edgeweight_to_nodes = function(a,e,w){

	tama = length(a)
	rr = rep(0,tama)
	
	if (!is.matrix(e)){
		e=t(as.matrix(e))
	}
	tame = dim(e)[1]
	#print(a);print(e);print(w);readline()

	if (!is.null(w)){
		for (i in 1:tame){
			e1 = e[i,1]
			e2 = e[i,2]
			s12 = a[e1]+a[e2]
			r1 = (a[e1]/s12)*w[i]
			r2 = (a[e2]/s12)*w[i]
			#print(a[e1]);print(a[e2]);print(s12);print(w[i]);print(r1);print(r2)#;readline()
			rr[e1]=rr[e1]+r1
			rr[e2]=rr[e2]+r2
			#print(rr);readline()
		}
	}
	return(rr)

}

#my_gsub = function(labels,)


fix_labels = function(labels,id1,id2,sep="-"){


	tam = length(labels)
	res = c()

	for (i in 1:tam){
		pos = as.numeric(unlist(strsplit(labels[i],sep)))#;print(pos);readline()
		res =c(res,paste0(id1[pos[1]],sep,id2[pos[2]]))#;print(res);readline()
	}
	return(res)

}


fix_labels_old = function(labels,id1,id2){
#id1 = menor; id2=maior

	#labels0=labels
	#print(id1);print(id2)
	#print(labels)
	labels = matrix(labels)
	tam1 = length(id1)
	tam2 = length(id2)
	pos1 = matrix(1:tam1)
	pos2 = matrix(1:tam2)
	id1 = matrix(id1)
	id2 = matrix(id2)
	#print(pos1);print(pos2);print(id1);print(id2);readline()

	#id1 = apply(id1,1,paste0,"+-")
	id1 = paste0(id1[,1],"+-")
	#id2 = stri_reverse(apply(id2,1,paste0,"+"))
	id2 = paste0("+",id2[,1])

	#pos1 = apply(pos1,1,paste0,"-")
	pos1 = paste0(pos1[,1],"-")
	#pos2 = stri_reverse(apply(pos2,1,paste0,"-"))
	pos2 = paste0("-",pos2[,1])
	
	#print(labels);print(pos1);print(pos2);print(id1);print(id2);readline()

	for (i in 1:tam1){
		pattern = pos1[i]#;print(pattern)
		replacement = id1[i]#;print(replacement)
		aux = apply(labels,1,gsub,pattern=pattern,replacement=replacement)
		labels = matrix(aux)
		#print(labels);readline()
	}
	for (i in 1:tam2){
		pattern = pos2[i]#;print(pattern)
		replacement = id2[i]#;print(replacement)
		aux = apply(labels,1,gsub,pattern=pattern,replacement=replacement)
		labels = matrix(aux)
		#print(labels);readline()
	}
	labels = as.vector(labels)
	labels = str_replace(labels,"\\+\\+","\\-")
	#if (id2[1]=="+2"){
	#	print(labels0);print(pos1);print(pos2);print(id1);print(id2);print(labels);
	#	print("---");readline()
	#}
	#print(labels)
	return(labels)
}

#make_score_graph_superimposition(xyz1,xyz2,g1i,g2i,a1i,a2i,tid,sid,rot,w=w,type="local")

#make_score_graph_superimposition = function(dist,g1,g2,a1,a2,id1,id2,type="global",w=c(1,1,1,1,1),consol=3,sep="-",k=2,xo=4,comp=T,edgw=T){
#make_score_graph_superimposition = function(xyz1,xyz2,g1,g2,a1,a2,id1,id2,rot,type="global",w=c(1,1,1,1,1),consol=4,sep="-",k=2,xo=5,comp=F,edgw=F){

#score0.local = make_score_graph_superimposition(dist0i,g1i,g2i,a1i,a2i,xyz1i,xyz2t0i,w=w)
make_score_graph_superimposition = function(dist,g1,g2,a1,a2,xyz1,xyz2,h1=NULL,h2=NULL,w=c(1,1,1,1,1),consol=5,sep="-",k=2,xo=5,comp=F,edgw=F){


	#print(dist);readline()
	auxsc = list()
	auxc = names(dist)#;print("o")

	#print(dist);print(auxc);readline()

	#grau dos nos
	d1 = degree(g1)#;print(d1);readline() 
	d2 = degree(g2)

	# soma dos pesos das arestas de no divido pelo grau
	# TO DO: testar possibilidade de ser "soma"
	#f1 = strength(g1)/d1
	#f2 = strength(g2)/d2
	#f1[is.nan(f1)]=0
	#f2[is.nan(f2)]=0

	#print(strength(g1));print(strength(g2))
	#print(f1);print(f2)
	#print(d1);print(d2);readline()
	#print(f1/d1);print(f2/d2);readline()

	tam1 = length(d1)
	tam2 = length(d2)

	auxrd = c()
	auxra = c()
	auxrf = c()
	auxre = c()
	auxrb = c()
	auxrw = c()
	auxs = c()
	auxa = c()
	auxe = c()
	auxd = c()
	auxw = c()
	auxb = c()
	auxcos = c()
	auxh = c()
	#auxsa = c()

	#print(paste("tam1:",tam1));print(paste("tam2:",tam2));readline()

	#if (tam1<=tam2){ # sempre g1 < g2
	#if (1){ # sempre g1 < g2
		#auxd1 = d1 # graus do grafo g1
		#auxd2 = d2 # graus do grafo g2
		auxa1 = a1 # labels dos nos g1 
		auxa2 = a2 # labels dos nos g2
		#auxf1 = f1 # labels das arestas normalizados por grau de g1
		#auxf2 = f2 # labels das arestas normalizados por grau de g2
		auxg1 = g1 # grafo g1
		auxg2 = g2 # grafo g2
		tam = tam1 # quantidade de nos
		auxyz1 = xyz1
		auxyz2 = xyz2
	#} #else {
#		auxd1 = d2
#		auxd2 = d1
#		auxa1 = a2
#		auxa2 = a1
#		auxf1 = f2
#		auxf2 = f1
#		auxg1 = g2
#		auxg2 = g1
#		tam = tam2
#		auxyz1 = xyz2
#		auxyz2 = xyz1
#	}
	#auxt1 = ends(auxg1,E(auxg1))
	#auxt2 = ends(auxg2,E(auxg2))
	auxt1 = ends(auxg1,E(auxg1))
	auxt2 = ends(auxg2,E(auxg2))
	auxg1w = E(auxg1)$weight
	auxg2w = E(auxg2)$weight
	if (is.null(auxg1w)){
		auxg1w = as.matrix(rep(0,dim(auxt1)[1]))
	}
	if (is.null(auxg2w)){
		auxg2w = as.matrix(rep(0,dim(auxt2)[1]))
	}
	#auxg1w = as.matrix(auxg1w)
	#auxg2w = as.matrix(auxg2w)
	auxt1 = cbind(auxt1,auxg1w)
	auxt2 = cbind(auxt2,auxg2w)#;print(auxc)
	#print(auxc);print(auxt1);print(auxt2);print(is.matrix(auxt2))
	auxel = verify_edge_superimpostion(auxc,auxt1,auxt2,auxyz1,auxyz2)#;print(auxel)#;readline()
	#auxe = auxel$rr
	#auxw = auxel$rw
	#auxcos = auxel$rc/tam
	#auxcos = auxel$rc/d1
	#auxcos[is.nan(auxcos)]=0
	auxcos =  auxel$rc

	#print(auxel);readline()
	
	#names(auxe)=auxc
	#names(auxw)=auxc
	names(auxcos)=auxc
	#print(auxe);print(auxw)#;readline()
	#auxw = auxw/auxe
	#auxw[is.nan(auxw)]=0
	#print(auxw);readline()
	#auxmaxd = klogistic(dist,k=2,xo=3.5) # decaimento logistico...
	#auxmaxd = klogistic(dist,k=1,xo=3) # decaimento logistico...
	#auxmaxd = klogistic(dist,k=2,xo=4) # decaimento logistico...
	auxmaxd = klogistic(dist,k=k,xo=xo) # decaimento logistico...
	#auxmaxd = klogistic(dist,k=2,xo=3) # decaimento logistico...
	#auxmaxd = klogistic(dist,k=2,xo=2.5) # decaimento logistico...
	#auxmaxd = klogistic(dist,k=2,xo=2.8) # decaimento logistico...
	#auxes = sum(auxe)#;print(auxd1)
	#if (auxes) {
		#auxre = auxe/auxes
	#	auxre = auxe/auxd1
	#	auxre[is.nan(auxre)]=0
	#}else{
	#	auxre = auxe
	#}
	#names(auxre)=auxc

	#auxs1 = sum(auxa1)
	#auxs2 = sum(auxa2)
	#print(auxa1);print(auxa2);
	#plot(auxg1)
	#print(E(auxg1)$weight);print(E(auxg2)$weight);readline()
	auxw1 = round(pull_edgeweight_to_nodes(diag(auxa1),auxt1,E(auxg1)$weight),k)
	auxw2 = round(pull_edgeweight_to_nodes(diag(auxa2),auxt2,E(auxg2)$weight),k)
	#print(auxw1);print(auxw2);readline()
	#auxax = make_score_matching_pair(auxc,diag(auxa1),diag(auxa2),sort=T)
	#auxdx = make_score_matching_pair(auxc,auxd1,auxd2,sort=F)
	#print(auxw1);print(auxw2);print(diag(auxa1));print(diag(auxa2));print(diag(auxa1)+auxw1);print(diag(auxa2)+auxw2)
	#auxa = make_score_matching_pair(auxc,diag(auxa1),diag(auxa2),sort=T)
	#print(auxa1);print(auxa2);readline()
	#auxa = make_score_matching_pair(auxc,diag(auxa1)+auxw1,diag(auxa2)+auxw2,sort=T)
#;print(auxc);readline()
	auxa = make_score_matching_pair(auxc,diag(auxa1)+auxw1,diag(auxa2)+auxw2,sort=F)#;print(auxa);readline()#testando intracadeia
	#auxa = make_score_matching_pair(auxc,diag(auxa1),diag(auxa2),sort=F)
	#if (comp) auxb = make_score_matching_pair(auxc,diag(auxa1)+auxw1,diag(auxa2)+auxw2,sort=T)#compensacao
	#print(auxa);readline()
	#print(auxd1);print(auxd2);readline()
	#auxd = make_score_matching_pair(auxc,auxd1,auxd2,sort=F)
	#print(auxd1);print(auxd2);print(auxd);readline()
	#auxw = make_score_matching_pair(auxc,auxw1,auxw2,sort=T)
	#print(auxw);readline()
	##for (i in 1:tam){

	if (!is.null(h1)&!is.null(h2)){
		auxh = make_score_matching_pair(auxc,h1,h2,sort=F,type=2)#;print(auxh);readline()
		auxrh = auxh*auxmaxd
	}else{
		auxh = 0
		auxrh = 0
	}

	#auxvt = compare_vectors(auto1.

	#auxre = auxre*auxmaxd # Edge superimposition
	auxra = auxa*auxmaxd # area
	#auxra = auxa 
	#auxrd = auxd*auxmaxd # degree
	auxrcos = auxcos*auxmaxd
	
	#auxrcos = auxcos
	#if (edgw) auxrw = auxw*auxmaxd
	#if (comp) auxrb = auxb*auxmaxd
	#auxrw = auxw*auxmaxd
	#auxs = w[1]*auxmaxd+w[2]*auxre+w[3]*auxrd+w[4]*auxra+w[5]*auxrw
	#auxs = w[1]*auxmaxd+w[2]*auxre+w[3]*auxrd+w[4]*auxra
	#auxmaxp = klogistic(dist,k=k,xo=xo-1)
	#auxs = w[1]*auxmaxd+w[2]*auxre+w[4]*auxra
	#print(w);readline()
	##auxs = 2*auxra + 1*auxrcos
	w = c(2,1,0)
	auxs = w[1]*auxra + w[2]*auxrcos + w[3]*auxrh#;print(auxs)
	#auxs = auxmaxd + auxrcos
	#print(auxs);readline()
	#print(auxc);readline()
	#auxsc$auto1 = auto1.svd$d
	#auxsc$auto2 = auto2.svd$d
	#v = auto1.svd$d
	#w = auto2.svd$d
	#cos_t=vector_matching_pair(v,w)
	#print(v);print(w)
	#cos_t = round(dot(w,v)/(Norm(w)*Norm(v)),10)#;print(cos_t);readline()
	#auxsc$cos = cos_t
	auxsc$dist = dist # distancia euclideana nos
	auxsc$distn = c(auxmaxd) # distancia logistica 	
	#auxsc$distp = c(auxmaxp) # distancia logistica mais restrita
	#auxsc$edsi = c(auxe) # edge superimposition
	#auxsc$edsv = c(auxre) # edge superimposition
	auxsc$nodi = c(auxa) # diff nos (areas)
	auxsc$nodv = c(auxra) # diff nos (areas)
	auxsc$cosi = c(auxcos)
	auxsc$cosv = c(auxrcos)
	auxsc$hoti = c(auxh)
	auxsc$hotv = c(auxrh)
	auxsc$scov = c(auxs) # score final
	#if (type=="local"){
	#	labels = fix_labels(names(auxsc$dist),id1,id2)
	#	names(auxsc$dist)=labels
	#	names(auxsc$distn)=labels
	#	names(auxsc$edsi)=labels
	#	names(auxsc$edsv)=labels
	#	names(auxsc$nodi)=labels
	#	names(auxsc$nodv)=labels
	#	names(auxsc$scov)=labels
	#}
	#print(auxsc);readline()
	#auxsc$edsi = auxe # edge superimposition
	#auxsc$edsv = auxre # edge superimposition
	#if (edgw){
	#	auxsc$edwi = auxw # diff nos (areas) compensacao
	#	auxsc$edwv = auxrw # diff nos (areas) compensacao
	#}
	#auxsc$degi = auxd # diff graus
	#auxsc$degv = auxrd # diff graus
	#auxsc$weii = auxw # diff peso arestas
	#auxsc$weiv = auxrw # diff peso arestas
	#auxsc$nodi = auxa # diff nos (areas)
	#auxsc$nodv = auxra # diff nos (areas)
	#if (comp){
	#	auxsc$noci = auxb # diff nos (areas) compensacao
	#	auxsc$nocv = auxrb # diff nos (areas) compensacao
	#}
	#auxsc$weiv = auxrf # diff pesos normalizados graus
	#auxsc$scov = auxs # score final
	
	if (consol==1)	auxsc$scor = median(auxs)#+auxe
	if (consol==2)	auxsc$scor = mean(auxs)#+auxe
	if (consol==3)	auxsc$scor = sum(auxs)#+auxe ### dependente No de nós
	if (consol==4)	{
		auxsc$scor = c(sum(auxs)+sum(auxsc$cos),sum(auxsc$distn),sum(auxsc$edsv),sum(auxsc$nodv))
		names(auxsc$scor) = c("total","distn","edsv","nodv")
	}
	if (consol==5)	{
		auxsi = w[1]*auxsc$nodi+w[2]*auxsc$cosi+w[3]*auxsc$hoti
		auxsc$scoi = c(sum(auxsi),w[1]*sum(auxsc$nodi),w[2]*sum(auxsc$cosi),w[3]*sum(auxsc$hoti))
		#auxsc$scoi = c(sum(auxsi),sum(auxsc$distn),sum(auxsc$cosi))
		names(auxsc$scoi) = c("total","nodi","cosi","hoti")
		#auxsc$scor = c(sum(auxs),sum(auxsc$nodv),sum(auxsc$cosv))
		auxsc$scor = c(sum(auxs),w[1]*sum(auxsc$nodv),w[2]*sum(auxsc$cosv),w[3]*sum(auxsc$hotv))#;print(sum(auxsc$hotv));readline()
		names(auxsc$scor) = c("total","nodv","cosv","hotv")
	}

	#auxsc$scor = sum(auxsc$cos)

	#auxsc$degs = sum(auxr)
	#print(auxsc);readline()
	return(auxsc) 
	

}

saveall = function(){

	save.image();savehistory()

}

get_colors = function(col){

	return(colors()[grep(col,colors())])

}


#reduz sempre ao 1o quadrante
cos_vec = function(w,v,k=10){

	Nw = Norm(w)
	Nv = Norm(v)

	if (!Nw | !Nv){
		print(paste("WARNING: cos_vec has one vector with Norm 0"))
		return(0)
	}else{
		return(abs(round(dot(w,v)/(Nw*Nv),k)))#;print(cos_t)
	}
}

roundto = function(n,bin=0.5){

	tam = length(n)
	auxr = c()

	for (i in 1:tam){

		n1 = floor(n[i])
		n2 = ceiling(n[i])
	
		if (n1!=n2){
			x = seq(n1,n2,by=bin)
			diff = abs(n[i]-x)
			id=which(min(diff)==diff)
			#print(x);print(diff);print(x[id[length(id)]])
			#return(x[id[length(id)]])
			auxr = c(auxr,x[id[length(id)]])
		}else{
			#return(n)
			auxr=c(auxr,n[i])
		}
	}
	return(auxr)
}


#make_rot_from_matrices = function(vsid,vtid,xyzs,xyzt,svd.v,rotagain=T){ #2a versao
make_rot_from_matrices = function(vsid,vtid,xyzs,xyzt,pdb1="",pdb2="",rotagain=T){


	if (pdb1!=pdb2){
	
		xyzs = xyzs[vsid,] #source 
		#xyz2 = xyz2[v2,] #target
		xyzt = xyzt[vtid,]

		#print(xyz1);print(xyz2)
		#xyzt.svd = svd(scale(xyzt,s=F))
	
		rotlist = get_rot_matrix(target=xyzt,source=xyzs,base=T)
		rot1 = rotlist$rot
		svd.v = rotlist$svd.t$v

		if (rotagain){
			r180.1 = rot180(svd.v[,1])
			r180.2 = rot180(svd.v[,2])
			rot2 = rot1
			rot3 = rot1
			rot4 = rot1
			rot2[1:3,1:3] = rot1[1:3,1:3] %*% r180.1[1:3,1:3] # rotacao em x (1o autovetor)
			rot3[1:3,1:3] = rot2[1:3,1:3] %*% r180.2[1:3,1:3] # equivale rotacao em z (3o autovetor)
			rot4[1:3,1:3] = rot3[1:3,1:3] %*% r180.1[1:3,1:3] # equivale rotacao em y (2o autovetor)

			#x = rbind(as.vector(rot1),as.vector(rot2))
			#print(x);readline()
			#print(c(as.vector(rot1),as.vector(rot2),as.vector(rot3),as.vector(vtid),as.vector(vsid)));readline()
			rot1 = c(as.vector(rot1),as.vector(vtid),as.vector(vsid))
			rot2 = c(as.vector(rot2),as.vector(vtid),as.vector(vsid))
			rot3 = c(as.vector(rot3),as.vector(vtid),as.vector(vsid))
			rot4 = c(as.vector(rot4),as.vector(vtid),as.vector(vsid))
			#return(c(as.vector(rot1),as.vector(rot2),as.vector(rot3),as.vector(vtid),as.vector(vsid)))
			return(c(rot1,rot2,rot3,rot4))
		}else{
			return(as.vector(rot1))
		}
	}else{
		rot1 = diag(4)
		rot1 = c(as.vector(rot1),as.vector(vtid),as.vector(vsid))	
		return(rot1)
	}
	#print(rot);readline()
	#return(rot)

}

make_rot_comb = function(targetid,sourceids,xyzt,xyzs,pdb1,pdb2){

	#reslist = list()
	
	#xyzt = xyzt[targetid,] # 1a versao
	### xyzt.svd = svd(scale(xyzt,s=F)) ### EM DUVIDA !!!

	#reslist$align = apply(sourceids,2,make_rot_from_matrices,xyzs,xyzt,xyzt.svd$v[,1])
	#res = apply(sourceids,2,make_rot_from_matrices,xyzs,xyzt,xyzt.svd$v[,1:2]) # 1a versao
	#res = apply(sourceids,2,make_rot_from_matrices,targetid,xyzs,xyzt,xyzt.svd$v[,1:2]) # 2a versao
	res = apply(sourceids,2,make_rot_from_matrices,targetid,xyzs,xyzt,pdb1,pdb2)
	
	#print(res);readline()
	return(res)

}

generate_all_new_alignments = function(target,source,pdb1,pdb2,minclus=4,maxclus=5,frame=16){

	#print(target);print(source)
	tamt = dim(target)[1]
	tams = dim(source)[1]
	tamin = min(tamt,tams)

	if (tamin<minclus) minclus = tamin
	if (tamin<maxclus) maxclus = tamin

	res = list()
	j = 1

	for (i in minclus:maxclus){
#	for (i in 4){

		targetids = combn(tamt,i)
		sourceids = combn(tams,i)

		#1OYV - 1PPF
		#targetids = matrix(c(1,3,4,5,6))
		#sourceids = matrix(c(1,2,4,5,6))

		#1PPF - 1OYV
		#targetids = matrix(c(1,2,4,5,6))
		#sourceids = matrix(c(1,3,4,5,6))

		#1PPF - 1OYV
		#targetids = matrix(c(1,2,3,4))
		#sourceids = matrix(c(1,2,3,4))

		#targetids = combn(tamt,minclus)
		#sourceids = combn(tams,minclus)

		#print(combt);print(combs);readline()

		#tamct = dim(tamt)[2]
		#tamcs = dim(tams)[2]

		res0=apply(targetids,2,make_rot_comb,sourceids,target,source,pdb1,pdb2)

		#print(res0[1:10,]);readline()
		#print(res0[1:10,]);readline()

		framei = frame+2*i

		res0 = as.vector(res0)
		res[[j]] = list()
		res[[j]]$rotmatrix = matrix(res0,ncol=framei,byrow=T)
		res[[j]]$minclus = i

		#print(res);readline()

		#res[[j]] = as.vector(res0)
		j = j + 1
		#print(dim(res));readline()

		#return(as.vector(res))
	}
	return(res)
}


#make_rot_comb = function(targetid,sourceids,xyzt,xyzs,pdb1,pdb2){
generate_all_new_alignments_par = function(target,source,pdb1,pdb2,ncores=0,minclus=4,maxclus=4,frame=16,maxspace=350000){ #max=550000

	print(paste("Generating alignment space..."))
	#print(target);print(source)
	tamt = dim(target)[1]#;print(tamt)
	tams = dim(source)[1]#;print(tams)

	tamin = min(tamt,tams)
	if (tamin<2){
		print(paste("WARNING: one of the graphs",pdb1," - ",pdb2,"is too small to align"))
		return(NULL)
	}

	if (tamin<minclus) minclus = tamin
	#if (tamin<maxclus) maxclus = tamin

	if (minclus>2){
		rotspace = choose(tamt,minclus)*choose(tams,minclus)*4#;print(rotspace);print(minclus)
		repeat{
			if (rotspace>maxspace){
				minclus = minclus - 1
				if (minclus==2) break;
			}else{
				break;
			}
			rotspace = choose(tamt,minclus)*choose(tams,minclus)*4
		}
	}
	#;print(rotspace);print(minclus);readline()

	#print(minclus);print(maxclus);readline()

	res = list()
	aux = list()
	k = 1
	#ncores = 10
	#for (i in minclus:maxclus){
	for (i in minclus){
#	for (i in 4){

		targetids = combn(tamt,i)#;print(targetids)#;readline()
		sourceids = combn(tams,i)#;print(sourceids);readline()

		#targetids = matrix(c(1,2,3,4,1,2,3,5,1,2,4,5),ncol=3)
		#sourceids = matrix(c(1,2,3,4,1,2,3,5,1,2,4,5),ncol=3)

		#print(targetids)
		#;print(sourceids);readline()

		tam = dim(targetids)[2]

		if ((ncores)&(tam>1)){
			print(paste("Doing it parallel with",ncores,"cores"))
			registerDoMC(ncores)
			aux = foreach(j=1:tam) %dopar% {
				make_rot_comb(targetids[,j],sourceids,target,source,pdb1,pdb2)
			}
		}else{
			print(paste("Doing it sequential..."))
			for (j in 1:tam){
				aux[[j]] = make_rot_comb(targetids[,j],sourceids,target,source,pdb1,pdb2)
			}
		}
		#res0=apply(targetids,2,make_rot_comb,sourceids,target,source,pdb1,pdb2)
		framei = frame+2*i
		aux =	unlist(aux)
		
		#x1 = matrix(aux,ncol=framei,byrow=T)
		#x2 = matrix(res0,ncol=framei,byrow=T)

		#res0 = as.vector(res0)
		res[[k]] = list()
		res[[k]]$rotmatrix = matrix(aux,ncol=framei,byrow=T)#;print(dim(res[[k]]$rotmatrix));readline()
		res[[k]]$minclus = i
		k = k + 1
		#print(x1);print(x2);
		#print(sum(apply(x1-x2,1,Norm)));
		#readline()
	}
	return(res)

}

score_kernel = function(g1,g2,type="WL"){

	gl = list()
	gl[[1]] = g1
	gl[[2]] = g2
	if (type=="WL"){
		if (all(degree(g1)==0)|all(degree(g2)==0)){
			ker.local = 0
		}else{
			ker = CalculateWLKernel(gl,5)#;print(ker)
			ker.diag = sum(diag(ker))
			if (ker.diag){
				ker.local = ((2*ker[1,2])/ker.diag)
			}else{
				print(paste("WARNING: zero kernel found"))#;readline()
				ker.local = 0
			}
		}

	}
	return(ker.local)

}

score_svd = function(disti,dist,xyz1,xyz2,sep="-",k=2,xo=5,r=3,zero=0.001){

	#print(disti);print(dist)

	if (length(dist)>=length(disti)){
		ddi = matrix(as.numeric(unlist(strsplit(names(disti),sep))),ncol=2,byrow=T)
		dd = matrix(as.numeric(unlist(strsplit(names(dist),sep))),ncol=2,byrow=T)#;print(ddi);print(dd)
		v1 = dd[,1] %in% ddi[,1]#;print(v1)
		dk = round(klogistic(dist,k=k,xo=xo),r)
		v2 = dk>zero#;print(v2)
		v = v1 | v2 #;print(v)
		xyz1.svd = svd(scale(xyz1[dd[v,1],],s=F))
		xyz2.svd = svd(scale(xyz2[dd[v,2],],s=F))
		svd.global = vector_matching_pair(xyz1.svd$d,xyz2.svd$d)
		svd.global.sum = sum(svd.global)
	}else{
		print(paste("WARNING: something wrong in score_svd"))
	}
	return(svd.global.sum)
}

get_hot_score = function(hot,tam){

	#tama1 = dim(a1i)
	hot.n = sapply(hot,dim)
	hot.id = which(tam==hot.n[1,])
	hot.score = hot[[hot.id]]$score
	return(hot.score)

}


plot_super = function(exp,n,mfrow,col=NULL,trans=1,title="1o"){

  #browser()
  
	clus=c()
	par(mfrow = mfrow)
	tam = length(n)
	tamclus = length(exp$superclus)
	if (trans) xyz1 = exp$superclus[[n[tam]]]$geomc
	else xyz1 = exp$superclus[[n[1]]]$geomc0
	limxyz=get_lim_xyz(xyz1,xyz1)
	
	
	for (i in n){

		if (i>tamclus) break

		super = exp$superclus[[i]]
		aux = list()
		aux$n = title

		aux$a = super$a
		auxd = diag(aux$a)
		aux$g = super$ga #g1
		#print(aux);readline()
		if (trans) coord = super$geomc
		else coord = super$geomc0
		vlabel = 1:length(auxd)
		#vcolor = "white"
		if (is.null(col)){
			vcolor = rep("white",length(auxd))
		}else{
			vcolor = col[[i]]$colors
		}
		#vcolor[tid0]="skyblue4"
		elabel = E(aux$g)$weight
		elabel = NULL
		sub=""
		#if (!is.null(dist0i)) sub=paste(round(base0,2),round(score0.local$scor[1],2),round(score0.global$scor[1],2))
		#xlim=c(-35,-10)
		#ylim=c(50,90)
		#xlim=c(min(super$geomc[,1]),max()
		xlim = limxyz$xlim
		ylim = limxyz$ylim
		sizef=10
		
    	graph_plot(aux,labcex=0.7,size=sizef*sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,coordi=coord, 
	ewidth=0.2,ecolor="gray75",lcolor="yellow",xlim=xlim,ylim=ylim,sub=sub)#;readline()

		clus = cbind(clus,exp$pam[[n[i]]]$clustering)

	}
	#print(clus)
	#print(table(apply(clus,1,paste,collapse="")))	
	
}


plot_g = function(parg,limxyz,title="1o",col=NULL,sub="",dx=1.1, vlabel = 1:dim(parg$a)[1], 
                  ecolor="gray75"){

  #browser()
	aux = list()
	aux$n = title

	aux$a = parg$a
	auxd = diag(aux$a)
	aux$g = parg$g #g1
	coord = parg$xyz
	#vlabel = 1:length(auxd)
	#vcolor = "white"
	if (is.null(col)){
		vcolor = rep("white",length(auxd))
	}else{
		vcolor = col
	}
	#vcolor[tid0]="skyblue4"
	elabel = E(aux$g)$weight
	elabel = NULL
	#sub=""
	xlim = limxyz$xlim
	ylim = limxyz$ylim
	#xmin = dx*min(coord[,1])
	#xmax = dx*max(coord[,1])
	#ymin = dx*min(coord[,2])
	#ymax = dx*max(coord[,2])
	#xlim = c(xmin,xmax)
	#ylim = c(ymin,ymax)
	#if (!is.null(dist0i)) sub=paste(round(base0,2),round(score0.local$scor[1],2),round(score0.global$scor[1],2))
	#xlim=c(0,60)
	#ylim=c(40,120)
	#xlim=c(-35,-10)
	#ylim=c(50,90)
	#xlim=c(-3,50)
	#ylim=c(40,100)
	sizef=10

	graph_plot(aux,labcex=0.7,size=sizef*sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,coordi=coord, 
ewidth=0.2,ecolor=ecolor,lcolor="yellow",xlim=xlim,ylim=ylim,sub=sub)#;readline()


}

#paste_t12 = function(v,t1,t2){

#	cbind(t1

#}

plot_from_db2 = function(db,tam=NULL,mfrow,trans=T,r=2){

	par(mfrow=mfrow)
	par(cex.main=0.9)
	tam  = length(db$score_list)

	base = db$base_vector
	base = base[base>0]
	o = order(base,decreasing=T)
	
	p1 = db$p1
	p2 = db$p2
	p2$xyz = db$score_list[[o[1]]]$xyz2t
	rotself = db$rotself
	if (trans){
		p1$xyz = round(transform_by_rot_matrix(p1$xyz,rotself),r+1)
		p2$xyz = round(transform_by_rot_matrix(p2$xyz,rotself),r+1)
	}
	limxyz=get_lim_xyz(p1$xyz,p2$xyz)
	tit1 = as.character(db$t1[1,1])
	tit2 = as.character(db$t2[1,1])
	sub1 = paste(c(round(db$t1[1,7],r)),collapse=" ")
	col1 = db$p1$spot$colors
	plot_g(p1,limxyz,tit1,col1)
	mtext(sub1,line=0,cex=0.6)
	for (i in o){
				
		p2 = db$p2
		p2$xyz = db$score_list[[i]]$xyz2t
		if (trans){
			p2$xyz = round(transform_by_rot_matrix(p2$xyz,rotself),r+1)
		}
		sub2 = paste(c(round(db$score_list[[i]]$score$scor[1],r),round(db$t2[1,7],r)),collapse=" ")
		col2 = db$p2$spot$colors
		plot_g(p2,limxyz,tit2,col2)
		mtext(sub2,line=0,cex=0.6)

	}

}


plot_from_db = function(db,tam=NULL,mfrow,r=2){

	par(mfrow=mfrow)
	par(cex.main=0.9)
	if (is.null(tam)) tam = 1:length(db)#;print(tam)

	for (i in tam){
		p1 = db[[i]]$p1
		p2 = db[[i]]$p2
		limxyz=get_lim_xyz(p1$xyz,p2$xyz)
		tit1 = as.character(db[[i]]$t1[1,1])
		tit2 = as.character(db[[i]]$t2[1,1])
		sub1 = paste(c(round(db[[i]]$score$scor[1],r),round(db[[i]]$t1[1,7],r)),collapse=" ")#;print(sub1)
		sub2 = paste(c(round(db[[i]]$score$scor[1],r),round(db[[i]]$t2[1,7],r)),collapse=" ")
#		plot_g(p1,limxyz,db[[i]]$t1,db[[i]]$score$scor[1:3]);plot_g(p2,limxyz,db[[i]]$t2,db[[i]]$score$scor[1:3])
		plot_g(p1,limxyz,tit1)
		mtext(sub1,line=0,cex=0.6)
		plot_g(p2,limxyz,tit2)
		mtext(sub2,line=0,cex=0.6)
	}
}

get_lim_xyz = function(xyz1,xyz2,dx=0.1){

	xyz = rbind(xyz1,xyz2)
	xmin = min(xyz[,1])-dx*min(xyz[,1])
	xmax = max(xyz[,1])+dx*min(xyz[,1])
	ymin = min(xyz[,2])-dx*min(xyz[,2])
	ymax = max(xyz[,2])+dx*min(xyz[,2])
	xlim = c(xmin,xmax)
	ylim = c(ymin,ymax)
	res = list()
	res$xlim = xlim
	res$ylim = ylim
	return(res)

}

plot_density = function(dbgrp,pdbid,sep="_"){

	j=0
	for (i in 1:length(dbgrp)){
		if (dbgrp[[i]][[1]][[1]]==pdbid){
			j=i
			break;
		}
	}
	#print(j);readline()
	if (j>0){
		par(mfrow=c(3,2))
		#par(cex.main=1)
		x = mapply(mean,dbgrp[[j]][[1]]$dense_bsr_bipar)
		if (length(x)){
			plot(x,type="l",main=paste0(c(dbgrp[[j]][[1]]$pdb,dbgrp[[j]][[1]]$combination_chains,"ALL"),collapse=sep))
			x = mapply(mean,dbgrp[[j]][[1]]$dense_bsr_total)
			lines(x,lty=2)
		}else{
			plot(0,type='n',main=paste0(c(dbgrp[[j]][[1]]$pdb,dbgrp[[j]][[1]]$combination_chains,"ALL"),collapse=sep))
		}
		x = mapply(mean,dbgrp[[j]][[2]]$dense_bsr_bipar)
		if (length(x)){
			plot(x,type="l",main=paste0(c(dbgrp[[j]][[2]]$pdb,dbgrp[[j]][[2]]$combination_chains,"ALL"),collapse=sep))
			x = mapply(mean,dbgrp[[j]][[2]]$dense_bsr_total)
			lines(x,lty=2)
		}else{
			plot(0,type='n',main=paste0(c(dbgrp[[j]][[2]]$pdb,dbgrp[[j]][[2]]$combination_chains,"ALL"),collapse=sep))
		}

		x = mapply(mean,dbgrp[[j]][[1]]$dense_bsr_bipar_polar_polar)
		if (length(x)){
			plot(x,type="l",main=paste0(c(dbgrp[[j]][[1]]$pdb,dbgrp[[j]][[1]]$combination_chains,"PP"),collapse=sep))
			x = mapply(mean,dbgrp[[j]][[1]]$dense_bsr_total_polar_polar)
			lines(x,lty=2)
		}else{
			plot(0,type='n',main=paste0(c(dbgrp[[j]][[1]]$pdb,dbgrp[[j]][[1]]$combination_chains,"PP"),collapse=sep))
		}
		x = mapply(mean,dbgrp[[j]][[2]]$dense_bsr_bipar_polar_polar)
		if (length(x)){
			plot(x,type="l",main=paste0(c(dbgrp[[j]][[2]]$pdb,dbgrp[[j]][[2]]$combination_chains,"PP"),collapse=sep))
			x = mapply(mean,dbgrp[[j]][[2]]$dense_bsr_total_polar_polar)
			lines(x,lty=2)
		}else{
			plot(0,type='n',main=paste0(c(dbgrp[[j]][[2]]$pdb,dbgrp[[j]][[2]]$combination_chains,"PP"),collapse=sep))
		}
		x = mapply(mean,dbgrp[[j]][[1]]$dense_bsr_bipar_nonpolar_nonpolar)
		if (length(x)){
			plot(x,type="l",main=paste0(c(dbgrp[[j]][[1]]$pdb,dbgrp[[j]][[1]]$combination_chains,"AA"),collapse=sep))
			x = mapply(mean,dbgrp[[j]][[1]]$dense_bsr_total_nonpolar_nonpolar)
			lines(x,lty=2)
		}else{
			plot(0,type='n',main=paste0(c(dbgrp[[j]][[1]]$pdb,dbgrp[[j]][[1]]$combination_chains,"AA"),collapse=sep))
		}
		x = mapply(mean,dbgrp[[j]][[2]]$dense_bsr_bipar_nonpolar_nonpolar)
		if (length(x)){
			plot(x,type="l",main=paste0(c(dbgrp[[j]][[2]]$pdb,dbgrp[[j]][[2]]$combination_chains,"AA"),collapse=sep))
			x = mapply(mean,dbgrp[[j]][[2]]$dense_bsr_total_nonpolar_nonpolar)
			lines(x,lty=2)
		}else{
			plot(0,type='n',main=paste0(c(dbgrp[[j]][[2]]$pdb,dbgrp[[j]][[2]]$combination_chains,"AA"),collapse=sep))
		}
	}

}

built_rmsd_list = function(rmsd.list,db,i,j,k,sep="-",plotg=T){

	p1 = dbali[[i]]$align[[j]]$p1
	p2 = dbali[[i]]$align[[j]]$p2
	p21 = p2
	p21$xyz = dbali[[i]]$align[[j]]$score_list[[k]]$xyz2t
	dist = dbali[[i]]$align[[j]]$score_list[[k]]$score$dist
	print(dist)
	d.m = matrix(as.numeric(unlist(strsplit(names(dist),sep))),ncol=2,byrow=T)
	f = d.m[,2]!=0
	d.m = d.m[f,]
	#print(d.m);readline()
	v1 = as.vector(scale(p1$xyz[d.m[,1],],scale=F))
	v2 = as.vector(scale(p2$xyz[d.m[,2],],scale=F))
	v21 = as.vector(scale(p21$xyz[d.m[,2],],scale=F))
	#print(d.m[,1]);print(d.m[,2]);readline()
	rmsd.list$vn = c(rmsd.list$vn,rmsd(v1,v2))
	rmsd.list$va = c(rmsd.list$va,rmsd(v1,v21))
	if (plotg){
		limxyz = get_lim_xyz(p1$xyz,p21$xyz)
		plot_g(p1,x,"","black")
		x11()
		plot_g(p21,x,"","black")
	}
	return(rmsd.list)
}


new_best_self_graph_pdb = function(parg1,parg2,rotref=diag(4),w=c(1,1,1,1,1),minclus=4,maxclus=4,r=5,comp=F,frame=16,sep="-",zero=0.001){

	res = list()
	xyz1=parg1$xyz#;print(xyz1)
	xyz2=parg2$xyz#;print(xyz2)
	g1=parg1$g
	g2=parg2$g
	a1=parg1$a
	a2=parg2$a
	tam1 = dim(xyz1)[1]#;print(tam1)
	tam2 = dim(xyz2)[1]#;print(tam2)
	#par(mfrow=c(1,2));plot_g(parg1);plot_g(parg2)
	pdbname1=parg1$pdbname
	pdbname2=parg2$pdbname
	#print(parg1);print(parg2)
	h1 = NULL
	h2 = NULL
	ker.global = score_kernel(g1,g2)
	#hot1=parg1$hot#;print(hot1);
	#hot2=parg2$hot#;print(hot2);readline()
	dist = find_nearest_nodes(xyz1,xyz2,force=T)#;print(dist)#;readline()
	dist.n = names(dist)
	dist.df = data.frame(strsplit(dist.n,"[-]"))#;print(dist.df)
	tid0 = as.numeric(as.character(unlist(dist.df[1,])));names(tid0)=NULL#;print(tid0)
	sid0 = as.numeric(as.character(unlist(dist.df[2,])));names(sid0)=NULL#;print(sid0)#;readline()
	if (length(dist)!=tam1){
		#par(mfrow=c(1,2));plot_g(parg1);plot_g(parg2)
		print(paste(c("WARNING: dist does not complete in",pdbname1," - ",pdbname2,":",dist),collapse=" "))
		print("Trying to correct it")#;readline()
	#}
	#if (sum(diff(tid0)!=1)){
		##a1i = a1[tid0,tid0]
		##a2i = a2[sid0,sid0]
		##g1i = graph.adjacency(a1i,weighted=TRUE,mode="undirected",diag=FALSE)
		##g2i = graph.adjacency(a2i,weighted=TRUE,mode="undirected",diag=FALSE)
		##xyz1i = xyz1[tid0,]
		##xyz2i = xyz2[sid0,]
		#p1=list();p1$a=a1i;p1$g=g1i;p1$xyz=xyz1i
		#p2=list();p2$a=a2i;p2$g=g2i;p2$xyz=xyz2i
		#par(mfrow=c(1,2));plot_g(p1);plot_g(p2)
		##dist2 = find_nearest_nodes(xyz1i,xyz2i)#;print(dist2)#;readline()
		##labels = fix_labels(labels=names(dist),id1=tid0,id2=sid0)
		##names(dist2) = labels;#;print(dist2);readline()
		#par(mfrow=c(1,2));limxyz=get_lim_xyz(parg1$xyz,parg2$xyz);plot_g(parg1,limxyz);plot_g(parg2,limxyz)
		#print(dist);print(dim(xyz1));print(dim(xyz2));print(tid0);print(sid0);
		dist3 = expand_to_nearest_global_nodes(dist,xyz1,xyz2,tid0,sid0)#;print(dist3);readline()
		if (length(dist3)!=tam1){
			print(paste(c("WARNING: unfortunately, expanded dist3 does not correct it:",dist3),collapse=" "))
			#print(xyz1);print(xyz2)
			#print(dist);print(dist2);print(dist3);readline()
		}else{
			print(paste(c("YEAP! it was done with expanded dist3:",dist3),collapse=" "))
			#print(dist);print(dist2);print(dist3)#;readline()
		}#;print(dist3)
		score.local = make_score_graph_superimposition(dist3,g1,g2,a1,a2,xyz1,xyz2,h1,h2,w=w)#;print(score.local);readline()
	}else{
		score.local = make_score_graph_superimposition(dist,g1,g2,a1,a2,xyz1,xyz2,h1,h2,w=w)#;print(score.local);readline()
	}
	svd.local = 0 
	svd.local.sum = sum(svd.local)#;print(svd.local.sum)
	#local.total = score0.local$scor[1]
	score.local$scor = c(score.local$scor,svd=svd.local.sum)#,ker=ker.local)
	score.local$scoi = c(score.local$scoi,svd=svd.local.sum)#,ker=ker.local)
	score.local$scor[1] = (score.local$scor[1]+svd.local.sum)#+ker.local)
	score.local$scoi[1] = (score.local$scoi[1]+svd.local.sum)
	
	score.local$tid = tid0
	score.local$sid = sid0
	score.global = score.local
	base = score.local$scor[1]
	base = 2*base
	res$imax = 1
	res$score = list()
	res$score$local = score.local
	res$score$global = score.global
	res$score$scot = base*(1+ker.global)
	res$score$ker = ker.global
	
	rot = diag(4)	 
	res$rot = round(rot,r)
	#xyz2t = transform_by_rot_matrix(xyz2,rot)#;print(rotref)
	#if (!is.null(rotref)){
	#	res$supxyz = transform_by_rot_matrix(xyz2t,rotref)
	#}else{
	res$supxyz = xyz2
	#}
	res$names = c(parg1$auxname,parg2$auxname)
	return(res)

}

save_best_scores = function(best,res0){

	base0 = res0$score$scor[1]
	id = best$baseid
	best$score_list[[id]] = res0
	best$base_vector[id] = base0
	id = which.min(best$base_vector)
	best$baseid = id
	best$base = best$base_vector[id]

	return(best)

}

initialize_best_list = function(n){

	best = list()
	best$base = 0
	best$baseid = 1
	best$base_vector = rep(0,n)
	best$score_list = vector("list",n)
	return(best)

}

scan_rotations = function(j,vbin,rotmatrix,framer,framet,frames,p1,p2,n=8){

	rangej = which(vbin==j)

	xyz1 = p1$xyz
	a1 = p1$a
	g1 = p1$g

	xyz2 = p2$xyz
	a2 = p2$a
	g2 = p2$g

	
	best = list()
	best = initialize_best_list(n)#;print(best);readline()
	base = best$base

	#rangej=55649
	for (i in rangej){#;print(i)
	#for (i in 5){
		rot0 = rotmatrix[i,framer]#;print(rot0)
		tid0 = rotmatrix[i,framet]#;print(tid0)
		sid0 = rotmatrix[i,frames]#;print(sid0)
		xyz1i = xyz1[tid0,]
		xyz2i = xyz2[sid0,]
		#xyz1i = scale(xyz1i,s=F)
		#xyz2i = scale(xyz2i,s=F)
		xyz2t0i = transform_by_rot_matrix(xyz2i,rot0)#;print("o")
		xyz2t0 = transform_by_rot_matrix(xyz2,rot0)
		#p2$xyz = xyz2t0;par(mfrow=c(1,2));limxyz=get_lim_xyz(p1$xyz,p2$xyz);plot_g(p1,limxyz);plot_g(p2,limxyz)
		dist0i = find_nearest_nodes(xyz1i,xyz2t0i)#;print(xyz1i);print(xyz2i);print(xyz2t0i);readline()
		#print(dist0i);readline()
		if (!is.null(dist0i)){
			labels = fix_labels(labels=names(dist0i),id1=tid0,id2=sid0)#;print(labels)
			dist1i = dist0i
			names(dist1i) = labels;
			#print(dist1i)#;readline()
			dist = expand_to_nearest_global_nodes(dist1i,xyz1,xyz2t0,tid0,sid0)
			#print(dist)#;readline()
			if(!is.null(dist)){
				local = F	
				score0 = make_score_graph_superimposition(dist,g1,g2,a1,a2,xyz1,xyz2t0)
				#print(score0);readline()
			}else{
				local = T
				print(paste("WARNING: it was impossible expand dist. Making only local dist..."))
				#xyz1i = xyz1[tid0,]
				#xyz2t0i = xyz2t0[sid0,]
				a1i = a1[tid0,tid0]
				a2i = a2[sid0,sid0]
				g1i = graph.adjacency(a1i,weighted=TRUE,mode="undirected",diag=FALSE)
				g2i = graph.adjacency(a2i,weighted=TRUE,mode="undirected",diag=FALSE)
				score0 = make_score_graph_superimposition(dist0i,g1i,g2i,a1i,a2i,xyz1i,xyz2t0i)
				#print(score0);readline()
			}
			base0 = score0$scor[1]#;print(base0)
			if (base0 > base){
				res0 = list()
				res0$jmax = i
				res0$score = score0
				if (local) res0$score = correct_score_labels(res0$score,tid0,sid0)
				res0$rot = rot0
				res0$xyz2t = xyz2t0
				best = save_best_scores(best,res0)
				base = best$base
				#print(best$base_vector);readline()
			}
		}else{
			#print(paste("WARNING: it was impossible to do a symmetric align with dist"))
		}
	}
	return(best)
}

correct_score_labels = function(score,tid,sid){

	labels = fix_labels(labels=names(score$dist),id1=tid,id2=sid)
	names(score$dist) = labels
	names(score$distn) = labels
	names(score$nodi) = labels
	names(score$nodv) = labels
	names(score$cosi) = labels
	names(score$cosv) = labels
	#names(score$hoti) = labels
	#names(score$hotv) = labels
	names(score$scov) = labels
	return(score)

}

make_graph_alignment = function(db1,db2,t1,t2,vcores,minclus=4,maxclus=4,r=5,comp=F,frame=16,sep="-",zero=0.001,n=8){

	#print(t1);print(t2)#;readline()

	#print(paste())
	pdbname1 = as.character(t1[1,"pdbname"])#;print(pdbname1)
	xyz1 = db1[[t1[1,"pdbid"]]][[t1[1,"inter"]]]$clus[[t1[1,"pol"]]]$nclus[[t1[1,"nclus"]]]$xyz
	a1 = db1[[t1[1,"pdbid"]]][[t1[1,"inter"]]]$clus[[t1[1,"pol"]]]$nclus[[t1[1,"nclus"]]]$a
	g1 = graph.adjacency(a1,weighted=TRUE,mode="undirected",diag=FALSE)
	spot1 = db1[[t1[1,"pdbid"]]][[t1[1,"inter"]]]$clus[[t1[1,"pol"]]]$nclus[[t1[1,"nclus"]]]$spot
	dense1 = db1[[t1[1,"pdbid"]]][[t1[1,"inter"]]]$clus[[t1[1,"pol"]]]$nclus[[t1[1,"nclus"]]]$dense
	rotself = db1[[t1[1,"pdbid"]]][[t1[1,"inter"]]]$clus[[t1[1,"pol"]]]$rotself
	qualy1 = db1[[t1[1,"pdbid"]]][[t1[1,"inter"]]]$clus[[t1[1,"pol"]]]$nclus[[t1[1,"nclus"]]]$sensitivity

	pdbname2 = as.character(t2[1,"pdbname"])#;print(pdbname2)
	xyz2 = db2[[t2[1,"pdbid"]]][[t2[1,"inter"]]]$clus[[t2[1,"pol"]]]$nclus[[t2[1,"nclus"]]]$xyz
	a2 = db2[[t2[1,"pdbid"]]][[t2[1,"inter"]]]$clus[[t2[1,"pol"]]]$nclus[[t2[1,"nclus"]]]$a
	g2 = graph.adjacency(a2,weighted=TRUE,mode="undirected",diag=FALSE)
	spot2 = db2[[t2[1,"pdbid"]]][[t2[1,"inter"]]]$clus[[t2[1,"pol"]]]$nclus[[t2[1,"nclus"]]]$spot
	dense2 = db2[[t2[1,"pdbid"]]][[t2[1,"inter"]]]$clus[[t2[1,"pol"]]]$nclus[[t2[1,"nclus"]]]$dense
	qualy2 = db2[[t2[1,"pdbid"]]][[t2[1,"inter"]]]$clus[[t2[1,"pol"]]]$nclus[[t2[1,"nclus"]]]$sensitivity
	#rotself2 = db2[[t2[1,"pdbid"]]][[t2[1,"inter"]]]$clus[[t2[1,"pol"]]]$rotself
	#print(paste("Initiating align between",pdbname1," - ",pdbname2))

	p1 = list()
	p1$xyz = xyz1
	p1$a = a1
	p1$g = g1
	p1$spot = spot1
	p1$dense = dense1
	p1$qualy = qualy1

	p2 = list()
	p2$xyz = xyz2
	p2$a = a2
	p2$g = g2
	p2$spot = spot2
	p2$dense = dense2
	p2$qualy = qualy2

	best = list()
	best = initialize_best_list(n)#;print(best);readline()
	base = best$base

	#print(xyz1);print(xyz2);readline()

	rotlist = generate_all_new_alignments_par(target=xyz1,source=xyz2,minclus=minclus,maxclus=maxclus,frame=frame,pdb1=pdbname1,pdb2=pdbname2, ncores=vcores[2])

	if (is.null(rotlist)){
		return(NULL)
	}

	tamlist = length(rotlist)#;print(tamlist);readline()
	#print(rotlist);readline()
	##### ker.global = score_kernel(g1,g2)#;print("o") ####
	score.local = NULL
	base=0
	bin = 10
	print(paste("Initiating align between",pdbname1," - ",pdbname2))

	for (i in 1:tamlist){### tamlist sempre 1 ###
	#for (i in 2){
		#print("o")
		minclus = rotlist[[i]]$minclus
		framer = 1:frame
		framet = (frame+1):(frame+minclus)
		frames = (frame+minclus+1):(frame+2*minclus)
		#frame = frame+2*minclus
		rotmatrix = rotlist[[i]]$rotmatrix#;print(rotmatrix)
		tamatrix = dim(rotmatrix)[1]#;print(tamatrix);readline()
		print(paste("Alignment space:",tamatrix,"Minclus:",minclus))
		if (bin>tamatrix) bin = tamatrix
		aux = list()
		v = 1:tamatrix#;print(tamatrix)
		if (bin>1){
			vbin = cut(v,bin)
			levels(vbin) = 1:bin#;print(table(vbin))
		}else{
			vbin = bin
		}
		#bin=7
		if ((vcores[3])&(bin>1)){
			print(paste("Doing it parallel with",vcores[3],"cores"))
			registerDoMC(vcores[3])
			aux = foreach(j=1:bin) %dopar% {
				#ids = which(vbin==j)
				#;print(j)
				scan_rotations(j,vbin,rotmatrix,framer,framet,frames,p1,p2)
			}
		}else{
			print(paste("Doing it sequential..."))
			#v = 1:tamatrix#;print(tamatrix)
			#vbin = cut(v,bin)
			#levels(vbin) = 1:bin#;print(vbin)
			for (j in 1:bin){
			#for (j in 7){
				#ids = which(vbin==j)#;print(ids);readline()
	  		aux[[j]] = scan_rotations(j,vbin,rotmatrix,framer,framet,frames,p1,p2)
				#print(aux[[j]]);readline()
			}
		}#;print("0")
		tambin = length(aux)
		base = 0
		for (k in 1:tambin){
			tamlist = length(aux[[k]]$score_list)#;print(tamlist)
			for (m in 1:tamlist){
				base0 = aux[[k]]$score_list[[m]]$score$scor[1]#;print(aux[[k]]$score_list)
				if (!is.null(base0)){
					if (base0 > base){
						res0 = list()
						res0 = aux[[k]]$score_list[[m]]
						#res0$jmax = i
						#res0$score = score0
						#if (local) res0$score = correct_score_labels(res0$score,tid0,sid0)
						#res0$rot = rot0
						#res0$xyz2t = xyz2t0
						best = save_best_scores(best,res0)
						base = best$base
						#print(best$base_vector);readline()
					}
				}
			}
			#base0 = aux[[k]]$score$scor[1]
			#if (base0 > base){
			#	ikjmax = c(i,k,aux[[k]]$jmax)
			#	rot = aux[[k]]$rot
			#	score = aux[[k]]$score
			#	base = base0
			#	xyz2t = aux[[k]]$xyz2
			#}
		}
	}
	#print(score$dist);
	#print(score$scor);
	##res = list()
	##res$ikjmax = ikjmax
	##res$score = score
	##res$rot = rot
	#res$xyz2 = xyz2t
	##res$t1 = t1
	##res$t2 = t2
	##res$p1 = p1
	##p2$xyz = xyz2t
	##res$p2 = p2
	#par(mfrow=c(1,2));limxyz=get_lim_xyz(p1$xyz,p2$xyz);plot_g(p1,limxyz);plot_g(p2,limxyz)
	#readline()
	best$t1 = t1
	best$t2 = t2
	best$p1 = p1
	best$p2 = p2
	best$minclus = minclus
	best$rotself = rotself
	#best$rotself2 = rotself2
	print(paste("Finish align between",pdbname1," - ",pdbname2))
	return(best)
}


# NOVO: metodo de similiridade usuando alinhamento proprio!!! ;-) in memory
#new_best_align_graph_pdb = function(xyz1,xyz2,g1,g2,a1,a2,pdbname1,pdbname2,hot1,hot2,rotref,w=c(1,1,1,1,1),minclus=4,maxclus=4,r=5,comp=F,frame=16,sep="-",zero=0.001){

new_best_align_graph_pdb = function(parg1,parg2,rotref=diag(4),w=c(1,1,1,1,1),minclus=4,maxclus=4,r=5,comp=F,frame=16,sep="-",zero=0.001){

	#print(g1);print(g2)
	#print(framer);print(framet);print(frames);readline()
	#print(parg1);print(parg2);readline()

	xyz1=parg1$xyz
	xyz2=parg2$xyz
	g1=parg1$g
	g2=parg2$g
	a1=parg1$a
	a2=parg2$a
	pdbname1=parg1$pdbname
	pdbname2=parg2$pdbname
	hot1=parg1$hot#;print(hot1);
	hot2=parg2$hot#;print(hot2);readline()

	res = list()

	tam1 = dim(xyz1)[1]#;print(tam1)
	tam2 = dim(xyz2)[1]#;print(tam2)

	#xyz1 = scale(xyz1,s=F)
	#xyz2 = scale(xyz2,s=F)

	rotlist = generate_all_new_alignments(target = xyz1,source = xyz2,minclus=minclus,maxclus=maxclus,frame=frame,pdb1=pdbname1,pdb2=pdbname2)#;print(rotlist);readline()
	#print("o")
	base=0
	tamlist = length(rotlist)#;print(tamlist);readline()
	#xyz1.svd = svd(scale(xyz1,s=F))
	#print("o");print(g1)
	ker.global = score_kernel(g1,g2)#;print("o")
	score.local = NULL
	#ker.global = 0

	for (i in 1:tamlist){
	#for (i in 2){
		#print("o")
		minclus = rotlist[[i]]$minclus
		framer = 1:frame
		framet = (frame+1):(frame+minclus)
		frames = (frame+minclus+1):(frame+2*minclus)
		#frame = frame+2*minclus
		rotmatrix = rotlist[[i]]$rotmatrix#;print(rotmatrix)
		tamatrix = dim(rotmatrix)[1]#;print(tamatrix);readline()
		for (j in 1:tamatrix){
		#for (j in 11){
		#for (j in 1:4){
			#print(paste(i,"-",j));
			rot0 = rotmatrix[j,framer]#;print(rot)
			tid0 = rotmatrix[j,framet]
			sid0 = rotmatrix[j,frames]
			if(0){
			print(tid0);print(sid0)
			view3d_with_base_axis(m=xyz1[tid0,],centralize=T,new=T,text=T)
			xyz2t0 = transform_by_rot_matrix(xyz2,rot0)
			view3d_with_base_axis(m=xyz2t0[sid0,],col="blue",k=0.5,add=T,centralize=T,text=T)#;readline()
			}
			xyz1i = xyz1[tid0,]
			xyz2i = xyz2[sid0,]
			#xyz1i = scale(xyz1i,s=F)
			#xyz2i = scale(xyz2i,s=F)
			xyz2t0i = transform_by_rot_matrix(xyz2i,rot0)#;print("o")
			xyz2t0 = transform_by_rot_matrix(xyz2,rot0)
			dist0i = find_nearest_nodes(xyz1i,xyz2t0i)#;print(dist0i)#;readline()
			#labels = fix_labels(labels=names(dist0i),id1=tid0,id2=sid0);print(labels)
			if (!is.null(dist0i)){
				#if(1){
				a1i = a1[tid0,tid0]
				a2i = a2[sid0,sid0]
				g1i = graph.adjacency(a1i,weighted=TRUE,mode="undirected",diag=FALSE)
				g2i = graph.adjacency(a2i,weighted=TRUE,mode="undirected",diag=FALSE)#;print(dist0i);readline()
				#par(mfrow=c(1,2));plot(g1i);plot(g2i)
				#ker.local = length(tid0)*score_kernel(g1i,g2i)
				#ker.local = 0
				if ( (!is.null(hot1)) & (!is.null(hot2)) ){
					h1 = get_hot_score(hot1,dim(a1i))#;print(h1)
					h2 = get_hot_score(hot2,dim(a2i))#;print(h2);readline()
				}else{
					h1 = NULL
					h2 = NULL
				}
				score0.local = make_score_graph_superimposition(dist0i,g1i,g2i,a1i,a2i,xyz1i,xyz2t0i,h1,h2,w=w)#;print(score0.local);readline()
				xyz1i.svd = svd(scale(xyz1i,s=F))
				xyz2i.svd = svd(scale(xyz2t0i,s=F))
				svd.local = vector_matching_pair(xyz1i.svd$d,xyz2i.svd$d) ### OTIMIZAR NO FUTURO, CALCULANDO JUNTO COM ROT
				svd.local.sum = sum(svd.local)#;print(svd.local.sum)
				#local.total = score0.local$scor[1]
				score0.local$scor = c(score0.local$scor,svd=svd.local.sum)#,ker=ker.local)
				score0.local$scoi = c(score0.local$scoi,svd=svd.local.sum)#,ker=ker.local)
				score0.local$scor[1] = (score0.local$scor[1]+svd.local.sum)#+ker.local)
				score0.local$scoi[1] = (score0.local$scoi[1]+svd.local.sum)#+ker.local)

				#score0.global = score0.local
				#		}
				xyz2t0 = transform_by_rot_matrix(xyz2,rot0)#;print("o")

				#if(1){
				labels = fix_labels(labels=names(dist0i),id1=tid0,id2=sid0)
				names(dist0i) = labels;
				#print(dist0i);readline()
				dist0 = expand_to_nearest_global_nodes(dist0i,xyz1,xyz2t0,tid0,sid0)#;print(dist0);readline()
				#dist0 = NULL
				global0=T
				if(!is.null(dist0)){
					#svd.global.sum = score_svd(dist0i,dist0,xyz1,xyz2t0)#;print(svd.global.sum);readline()
					svd.global.sum = svd.local.sum
					if ( (!is.null(hot1)) & (!is.null(hot2)) ){
						h1 = get_hot_score(hot1,dim(a1))#;print(h1)
						h2 = get_hot_score(hot2,dim(a2))#;print(h2);readline()
					}else{
						h1 = NULL
						h2 = NULL
					}					
					score0.global = make_score_graph_superimposition(dist0,g1,g2,a1,a2,xyz1,xyz2t0,h1,h2,w=w)									
					score0.global$scor = c(score0.global$scor,svd=svd.global.sum)#,ker=ker.global)
					score0.global$scoi = c(score0.global$scoi,svd=svd.global.sum)#,ker=ker.global)
					score0.global$scor[1] = (score0.global$scor[1]+svd.global.sum)#+ker.global)
					score0.global$scoi[1] = (score0.global$scoi[1]+svd.global.sum)#+ker.global)
					if (score0.global$scor[1]>score0.local$scor[1]){
						base0 = score0.global$scor[1]
					}else{
						base0 = score0.local$scor[1]
					}
				}else{
					global0 = F
					#score0.global = score0.local
					base0 = score0.local$scor[1]
					
				}
				#	}
				#}else{
				#	score0.global = score0.local#;print(score0.global);readline()
					#if (svd.global.sum>svd.local.sum){
					#	score0.global$scor["svd"] = svd.global.sum
					#	score0.global$scoi["svd"] = svd.global.sum
					#	score0.global$scor[1] = score0.global$scor[1]+svd.global.sum-svd.local.sum
					#	score0.global$scoi[1] = score0.global$scoi[1]+svd.global.sum-svd.local.sum
					#}
				#}

				#print(score0.global);readline()
					
				#v = auto1.svd$d
				#w = auto2.svd$d
				#cos_t=vector_matching_pair(v,w)		


				#score0.local = make_score_graph_superimposition(xyz1,xyz2,g1,g2,a1,a2,tid0,sid0,rot0,type="local",w=w)
				#score0.global = make_score_graph_superimposition(xyz1,xyz2,g1,g2,a1,a2,tid0,sid0,rot0,type="global",w=w)
				#print(tid0);print(sid0);print(score0.local);readline()
				if(0){
				print(score0.local)
				print(sum(score0.local$distn))
				print(sum(score0.local$edsv))
				#print(sum(score0.local$edwi))
				#print(sum(score0.local$degi))	
				print(sum(score0.local$nodv))
				#print(sum(score0.local$noci))
				}
				#base0 = score0.local$scor[1]+score0.global$scor[1]
				base0 = 2*base0
				if (base0 > base){
					imax = c(i,j)
					rot = rot0
					#dist = dist0
					score.local = score0.local
					score.local$tid = tid0
					score.local$sid = sid0
					global = global0
					if (global) score.global = score0.global
					else score.global = score0.local
					#xyz2t = xyz2t0
					#tid = tid0
					#sid = sid0
					base = base0
				}
			}#;print("o")
			if (i == -105){
			#if (i > -105){
			#if ((i == 9)|(i==22)){
				#print(i);print(tid0);print(sid0)#;print(score0.local$scor);print(score0.global$scor);print(base0)
				#print(score0.local$dist);print(score0.local$distn);print(score0.global$dist);print(score0.global$distn)
				#print(score0.local$scor);print(score0.global$scor);print(base0)
				#print(score0.local);print(score0.global);print(base0)
				#print(score0.local$scoi);print(score0.local$scor)
				#print(score0.global$scoi);print(score0.global$scor);print(base0)
				##labels = fix_labels(labels=names(score.local$dist),id1=tid0,id2=sid0)
				##print(score0.local$scor);print(score0.global$scor);print(score0.local$dist);print(labels)
			#	print(rot0);print(xyz2t0)
				if(0){
				view3d_with_base_axis(m=xyz1[tid0,],centralize=T,new=T,text=T,also2d=F)
				xyz2t0 = transform_by_rot_matrix(xyz2,rot0)
				view3d_with_base_axis(m=xyz2t0[sid0,],col="blue",k=0.5,add=T,centralize=T,text=T,also2d=F)#;readline()
				}
				if(1){
				par(mfrow=c(1,2))
				aux = list()
				aux$n = "1o"
				#aux$a = a1[tid0,tid0]
				aux$a = a1
				auxd = diag(aux$a)
				auxg1 = graph.adjacency(aux$a,weighted=TRUE,mode="undirected",diag=FALSE)
				#plot(auxg);readline()
				aux$g = auxg1 #g1
				coord = xyz1
				coord = transform_by_rot_matrix(coord,rotref)#;print(coord)
				#xyz1.svd = svd(scale(coord,scale=F))
				#n180.1 = rot180(xyz1.svd$v[,1]);print(coord)
				#rotref1 = matrix(rotref,ncol=4,byrow=F)				
				#n180.1[4,] = c(-colMeans(coord),1)
				#n180.1[,4] = c(colMeans(coord),1);print(n180.1)#;readline()
				#n180[1:3,1:3] = rotref1[1:3,1:3] %*% n180.1[1:3,1:3];print(n180)
				#coord = transform_by_rot_matrix(coord,n180.1);print(coord)
				vlabel = 1:length(auxd)
				#vcolor = "white"
				vcolor = rep("white",length(auxd))
				vcolor[tid0]="skyblue4"
				elabel = E(aux$g)$weight
				elabel = NULL
				sub=""
				if (!is.null(dist0i)) sub=paste(round(base0,2),round(score0.local$scor[1],2),round(score0.global$scor[1],2))
				xlim=c(0,60)
				ylim=c(40,80)
				ylim=c(0,40)
				#xlim=c(-15,30)
				#ylim=xlim
				sizef=10
				graph_plot(aux,labcex=0.7,size=sizef*sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,coordi=coord, ewidth=0.2,ecolor="gray75",lcolor="yellow",xlim=xlim,ylim=ylim,sub=sub)#;readline()
				aux$n = "2o"
	#			aux$a = a2[sid0,sid0]
				aux$a = a2
				auxd = diag(aux$a)
				auxg2 = graph.adjacency(aux$a,weighted=TRUE,mode="undirected",diag=FALSE)
				aux$g = auxg2 #g2
				coord = xyz2t0
				#coord = xyz2
				coord = transform_by_rot_matrix(coord,rotref)#;print(coord)
				vlabel = 1:length(auxd)
				#vcolor = "white"
				vcolor = rep("white",length(auxd))
				vcolor[sid0]="skyblue4"
				elabel = E(aux$g)$weight
				elabel = NULL
				sub=""
				if (!is.null(dist0i)) sub=paste(round(base0,2),round(score0.local$scor[1],2),round(score0.global$scor[1],2))
				#xlim=c(0,60)
				#ylim=c(0,85)
				sizef=10
				graph_plot(aux,labcex=0.7,size=sizef*sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,coordi=coord, ewidth=0.2,ecolor="gray75",lcolor="yellow",xlim=xlim,ylim=ylim,sub=sub)#;readline()
				#gl = list()
				#gl[[1]]= auxg1
				#gl[[2]]= auxg2
				#x=shortest_path_kv_kernel(gl)
				#print(x)
				#x1 = CalculateShortestPathKernel(gl)
				#print(x1)
				#x2=CalculateGeometricRandomWalkKernel(gl,0.1)
				#print(x2)
				#x3=CalculateWLKernel(gl,5)
				#print(x3)
				#;readline()
				}
			}
			#i = i + 1
		}
	}

	#a1k = a1[tid,tid]
	#g1k = graph.adjacency(a1k,weighted=TRUE,mode="undirected",diag=FALSE)
	#a2k = a2[sid,sid]
	#g2k = graph.adjacency(a2k,weighted=TRUE,mode="undirected",diag=FALSE)
	#gl = list()
	#gl[[1]] = g1
	#gl[[2]] = g2
	#ker = shortest_path_kv_kernel(gl)
	#ker = CalculateWLKernel(gl,5)
	#ker.diag = sum(diag(ker))
	#if (ker.diag){
	#	score.kernel = (2*ker[1,2])/ker.diag
	#}else{
	#	print(paste("WARNING: zero kernel found"));readline()
	#	score.kernel = 0
	#}
	#print(ker);print(score.kernel)#;readline()

	#print(paste("tid",tid),coll="");print(paste("sid",sid))
	#print(tam1);print(tam2)
	#print(score.local$distn)
	#print(score.global$distn)
	#if (tam1<=tam2){
	if (!is.null(score.local)){
		#labels = fix_labels(labels=names(score.local$dist),id1=tid,id2=sid)
		labels = fix_labels(labels=names(score.local$dist),id1=score.local$tid,id2=score.local$sid)

		#}else{
		#	labels = fix_labels(labels=names(score.local$dist),id1=sid,id2=tid)
		#}
		names(score.local$dist)=labels
		names(score.local$distn)=labels
		#names(score.local$edsi)=labels
		#names(score.local$edsv)=labels
		names(score.local$nodi)=labels
		names(score.local$nodv)=labels
		names(score.local$cosi)=labels
		names(score.local$cosv)=labels
		names(score.local$scov)=labels
		if (global){
			d = score.global$distn
			dd = matrix(as.numeric(unlist(strsplit(names(d),sep))),ncol=2,byrow=T)
			v1 = dd[,2]!=0
			v2 = d > zero
			v = v1 & v2
			score.global$tid = dd[v,1]
			score.global$sid = dd[v,2]		
		}else{
			score.global = score.local
		}
		#print(score.local$distn)
		#print(score.global$distn)
		#print(score.local);print(score.global)

		#score.kernel = 1
		res$imax = imax
		res$score = list()
		res$score$local = score.local
		res$score$global = score.global
		#res$score$kernel = score.kernel
		#res$score$scot = base+score.kernel*base
		res$score$scot = base*(1+ker.global)
		#res$score$tid = tid
		#res$score$sid = sid
		res$score$ker = ker.global
		 
		res$rot = round(rot,r)
		xyz2t = transform_by_rot_matrix(xyz2,rot)#;print(rotref)
		if (!is.null(rotref)){
			res$supxyz = transform_by_rot_matrix(xyz2t,rotref)
		}else{
			res$supxyz = xyz2t
		}
	}else{
		print(paste("WARNING: it seems it was not possible to align something!"))
		res$supxyz = xyz2
		res$score$local = list()
		res$score$local$scor = 0
		res$score$global = list()
		res$score$global$scor = 0 
		res$score$scot = 0
		res$score$ker = ker.global
	}
	res$names = c(parg1$auxname,parg2$auxname)
	#print(tid);print(sid)
	#print(res$score$local$distn)
	#print(res$score$global$distn)
	#print(res$score$local$scoi);print(res$score$local$scor)
	#print(res$score$global$scoi);print(res$score$global$scor)
	#print(res$imax)
	#print(res);readline()
	#view3d_with_base_axis(m=xyz1)
	#view3d_with_base_axis(m=xyz2t,col="lightblue",k=0.5,add=T)#;readline()
	if (0){
		par(mfrow=c(1,2))
		aux = list()
		aux$n ="1o"
		#aux$a = a1[tid,tid]
		aux$a = a1
		auxd = diag(aux$a)
		auxg1a = graph.adjacency(aux$a,weighted=TRUE,mode="undirected",diag=FALSE)
		aux$g = auxg1a #g1
		#coord = xyz1[tid,]
		coord = xyz1
		coord = transform_by_rot_matrix(coord,rotref)
		vlabel = 1:length(auxd)
		vcolor = rep("white",length(auxd))
		vcolor[tid]="skyblue4"
		#vcolor = "white"
		elabel = E(aux$g)$weight
		elabel = NULL
		sub=paste(round(base,2),round(score.local$scor[1],2),round(score.global$scor[1],2))
		xlim=c(0,60)
		ylim=c(40,80)
		sizef=10
		graph_plot(aux,labcex=0.7,size=sizef*sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,coordi=coord, ewidth=0.2,ecolor="gray75",lcolor="yellow",xlim=xlim,ylim=ylim,sub=sub)#;readline()
	}
	if(0){
		aux$n ="2o"
		#aux$a = a2[sid,sid]
		aux$a = a2
		auxd = diag(aux$a)
		auxg2a = graph.adjacency(aux$a,weighted=TRUE,mode="undirected",diag=FALSE)
		aux$g = auxg2a#g2
		#coord = xyz2t[sid,]
		coord = xyz2t
		coord = transform_by_rot_matrix(coord,rotref)
		vlabel = 1:length(auxd)
		vcolor = rep("white",length(auxd))
		vcolor[sid]="skyblue4"
		elabel = E(aux$g)$weight
		elabel = NULL
		sub=paste(round(base,2),round(score.local$scor[1],2),round(score.global$scor[1],2))
		xlim=c(0,60)
		ylim=c(40,80)
		sizef=10
		graph_plot(aux,labcex=0.7,size=sizef*sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,coordi=coord, ewidth=0.2,ecolor="gray75",lcolor="yellow",xlim=xlim,ylim=ylim,sub=sub)
		#gl = list()
		#gl[[1]]= auxg1a
		#gl[[2]]= auxg2a
		#x=shortest_path_kv_kernel(gl)
		#print(x)
		#x1 = CalculateShortestPathKernel(gl)
		#print(x1)
		#x2=CalculateGeometricRandomWalkKernel(gl,0.1)
		#print(x2)
		#x3=CalculateWLKernel(gl,5)
		#print(x3)
		#readline();
	}

	return(res)
}


#align_graph_pdb = function(name1,name2,command="python",script="script-python-v1.py",workdir="Pymol/",minwin=3,sufix=".pdb"){
best_align_graph_pdb = function(name1,name2,xyz1,xyz2,g1,g2,a1,a2,w=c(1,1,1,1,1),crit="dist",command="pymol",script="script-align-v2",workdir="Pymol/",minwin=3,sufix=".pdb",r=5,guide=0){

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
	auxr = list()

	i=minwin		
	#i = 13
	flags = c("-cq")
	#args = c("-cp",script,name1,name2,i,guide) #ultimo se guide=0 (todos os atomos), guide=1 so CA
	args = c(name1,name2,i,guide)
	#auxlist = align_graph_pdb(command,args)
	#print(command);print(flags);print(script);print(args);readline()
	auxlist = align_graph_pdb(command,flags,script,args)
	#print(length(auxlist));print(auxlist);readline()
	auxrmsd = auxlist[[1]]$rmsd
	auxwin = i
	auxrot = auxlist[[1]]$rot
	auxsym = auxlist[[1]]$sym
	auxdist = 0
	#print(name1);print(name2)
	if (crit=="rmsd") auxmin = auxrmsd
	if (crit=="dist"){
		#print(xyz2);print(auxrot);readline()
		xyz2t = transform_by_rot_matrix(xyz2,matrix(auxrot,ncol=4,byrow=T))#;print(i)
		auxdist = find_nearest_nodes(xyz1,xyz2t)
		#auxscore = sum(auxdist)
		#print(xyz1);print(xyz2)
		#if(rgl.cur()!=0) rgl.close()
		#plot3d(xyz1,box=F,type="s",col="blue",add=T)
		#plot3d(xyz2t,box=F,type="s",col="red",add=T)#;readline()
		auxscore = make_score_graph_superimposition(auxdist,g1,g2,a1,a2,w)
		#print(i);print(auxrmsd);print(auxscore);readline()
		#auxmin = auxscore$scor
		auxbase = auxscore$scor
	}
	j=1	
	while(auxlist[[j]]$rmsd !=0 ){
		i=i+1
		#args = c("-cq",script,name1,name2,i)
		#args[5] = i
		#print(args);
		args[3] = i
		#print(args);readline()
		#auxlist = align_graph_pdb(command,args)
		#print(i)
		auxlist = c(auxlist,align_graph_pdb(command,flags,script,args))
		#print(length(auxlist));readline()#;print(auxlist);readline()
		j=j+1
		#print(auxlist[[j]]);readline()
		if (auxlist[[j]]$rmsd == 0) break;
		if (crit=="rmsd"){
			if (auxlist[[j]]$rmsd < auxmin){
				auxrmsd = auxlist[[j]]$rmsd
				auxwin = i
				auxrot = auxlist[[j]]$rot
				auxmin = auxrmsd
				auxsym = auxlist[[j]]$sym
			}
		}
		if (crit=="dist"){
			auxrmsd0 = auxlist[[j]]$rmsd
			auxrot0 = auxlist[[j]]$rot
			auxsym0 = auxlist[[j]]$sym
			xyz2t0 = transform_by_rot_matrix(xyz2,matrix(auxrot0,ncol=4,byrow=T))#;print(i)
			auxdist0 = find_nearest_nodes(xyz1,xyz2t0)
			#auxscore0 = sum(auxdist0)
			#auxscore0 = median(auxdist0)
			auxscore0 = make_score_graph_superimposition(auxdist0,g1,g2,a1,a2,w)#;print(auxscore0);readline()
			#if(rgl.cur()!=0) rgl.close()
			#plot3d(xyz1,box=F,type="s",col="blue",add=T)
			#plot3d(xyz2t,box=F,type="s",col="red",add=T)#;readline()
			#print(j);print(auxrmsd0);print(auxscore0);print(auxmin);readline()
			#if (auxscore0$scor < auxmin){
			if (auxscore0$scor > auxbase){
			#print(xyz2);print(auxrot);readline()
				auxrmsd = auxrmsd0
				auxwin = i
				auxrot = auxrot0
				auxdist = auxdist0
				auxscore = auxscore0
				#auxmin = auxscore0$scor
				auxbase = auxscore0$scor
				auxsym = auxsym0
				xyz2t = xyz2t0
				#xyz2t = transform_by_rot_matrix(xyz2,matrix(auxrot0,ncol=4,byrow=T))
				#plot3d(xyz1,box=F,type="s",col="blue",add=F)
				#plot3d(xyz2t,box=F,type="s",col="red",add=T);readline()
			}
		}
	}
	#auxlist$dist = round(auxdist,r)
	#auxlist$score = round(auxscore,r)
	#auxlist$score = mean(auxlist$dist)
	### ERROR ###	
	###auxlist$score = auxscore
	###auxlist$rmsd = round(auxrmsd,r)
	###auxlist$win = auxwin
	###auxlist$rot = round(auxrot,r)
	###auxlist$sym = auxsym
	###auxlist$score$pain = verify_penalty(auxscore$distn,a1,a2)
	###auxlist$score$scot = auxlist$score$scor-sum(auxlist$score$pain)
	###auxlist$supxyz = xyz2t 
	### ERROR ###
	auxr$score = auxscore
	auxr$rmsd = round(auxrmsd,r)
	auxr$win = auxwin
	auxr$rot = round(auxrot,r)
	auxr$sym = auxsym
	auxr$score$pain = verify_penalty(auxscore$distn,a1,a2)
	auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)
	auxr$supxyz = xyz2t 	
	
	#print(x);readline()
	#print(auxlist);readline()
	return(auxr)
}

test = function(x){
	
	y = list()
	y$x = x
	return(y)

}

# ANTIGO: metodo de similiridade usuando cealign pre-computado
#pre_make_similarity_assessment = function(db,dball,id,pdbname,outdir="Datalign/",sufix="align.csv",sep="_"){
#auxr[[j]] = pre_best_align_graph_pdb(auxname1,auxname2,auxlist,xyz1,xyz2,g1,g2,a1,a2,w)
pre_best_align_graph_pdb = function(auxlist,xyz1,xyz2,g1,g2,a1,a2,rotref,w=c(1,1,1,1,1),crit="dist",minwin=3,r=5,comp=F){

	if (!is.matrix(auxlist)){
		auxlist = t(as.matrix(auxlist))
	}
	#print(auxlist[1,])
	
	##script = paste(workdir,script,sep="")
	#command = paste(workdir,command,sep="")
	##name1 = paste(workdir,name1,sufix,sep="")
	##name2 = paste(workdir,name2,sufix,sep="")
	auxrmsd = c()
	auxwin = c()
	auxrot = c()
	auxsym = c()
	auxr = c()
	colrot = 5:20
	#auxrotlist = list()
	##auxlist = list()

	i=minwin		
	#i = 13
	##flags = c("-cq")
	#args = c("-cp",script,name1,name2,i,guide) #ultimo se guide=0 (todos os atomos), guide=1 so CA
	##args = c(name1,name2,i,guide)
	#auxlist = align_graph_pdb(command,args)
	##print(command);print(flags);print(script);print(args);readline()
	##auxlist = align_graph_pdb(command,flags,script,args)
	#print(length(auxlist));print(auxlist);readline()
	#auxrmsd = auxlist[[1]]$rmsd
	#print(auxlist[1,])
	#print(colnames(auxlist))
	auxrmsd = auxlist[1,"rmsd"]
	#print(auxrmsd);readline()
	auxwin = auxlist[1,"win"]
	auxrot = auxlist[1,colrot]
	auxsym = auxlist[1,"sym"]
	auxdist = 0
	#print(name1);print(name2)
	if (crit=="rmsd") auxmin = auxrmsd
	if (crit=="dist"){
		#print(xyz2);print(auxrot);print(matrix(auxrot,ncol=4,byrow=T))#;readline()
		#x = matrix(auxrot,ncol=4,byrow=T)
		#x = data.matrix(x)
		#print(x[4,]);readline()
		xyz2t = transform_by_rot_matrix(xyz2,matrix(auxrot,ncol=4,byrow=T))#;print(i)
		#print(xyz2t);readline()
		auxdist = find_nearest_nodes(xyz1,xyz2t)
		#auxscore = sum(auxdist)
		#print(xyz1);print(xyz2)
		#if(rgl.cur()!=0) rgl.close()
		#plot3d(xyz1,box=F,type="s",col="blue",add=T)
		#plot3d(xyz2t,box=F,type="s",col="red",add=T)#;readline()
		auxscore = make_score_graph_superimposition(auxdist,g1,g2,a1,a2,w)
		#print(i);print(auxrmsd);print(auxscore);readline()
		#auxmin = auxscore$scor
		auxbase = auxscore$scor
	}
	j=1#;print("oi")	
	while(auxlist[j,"rmsd"] !=0 ){
		i=i+1
		#args = c("-cq",script,name1,name2,i)
		#args[5] = i
		#print(args);
		##args[3] = i
		#print(args);readline()
		#auxlist = align_graph_pdb(command,args)
		#print(i)
		##auxlist = c(auxlist,align_graph_pdb(command,flags,script,args))
		#print(length(auxlist));readline()#;print(auxlist);readline()
		j=j+1
		#print(auxlist[[j]]);readline()

		if (auxlist[j,"rmsd"] == 0) break;
		if (crit=="rmsd"){
			if (auxlist[j,"rmsd"] < auxmin){
				auxrmsd = auxlist[j,"rmsd"]
				auxwin = auxlist[j,"win"]
				auxrot = auxlist[j,colrot]
				auxmin = auxlist[j,"rmsdj"]
				auxsym = auxlist[j,"sym"]
			}
		}
		if (crit=="dist"){
			auxrmsd0 = auxlist[j,"rmsd"]
			auxrot0 = auxlist[j,colrot]
			auxsym0 = auxlist[j,"sym"]
			xyz2t0 = transform_by_rot_matrix(xyz2,matrix(auxrot0,ncol=4,byrow=T))#;print(i)
			auxdist0 = find_nearest_nodes(xyz1,xyz2t0)#;print("oi")
			#auxscore0 = sum(auxdist0)
			#auxscore0 = median(auxdist0)
			auxscore0 = make_score_graph_superimposition(auxdist0,g1,g2,a1,a2,w)#;print(auxscore0);readline()
			#if(rgl.cur()!=0) rgl.close()
			#plot3d(xyz1,box=F,type="s",col="blue",add=T)
			#plot3d(xyz2t,box=F,type="s",col="red",add=T)#;readline()
			#print(j);print(auxrmsd0);print(auxscore0);print(auxmin);readline()
			#if (auxscore0$scor < auxmin){
			if (auxscore0$scor > auxbase){
			#print(xyz2);print(auxrot);readline()
				auxrmsd = auxrmsd0
				auxwin = auxlist[j,"win"]
				auxrot = auxrot0
				auxdist = auxdist0
				auxscore = auxscore0
				#auxmin = auxscore0$scor
				auxbase = auxscore0$scor
				auxsym = auxsym0
				xyz2t = xyz2t0
				#xyz2t = transform_by_rot_matrix(xyz2,matrix(auxrot0,ncol=4,byrow=T))
				#plot3d(xyz1,box=F,type="s",col="blue",add=F)
				#plot3d(xyz2t,box=F,type="s",col="red",add=T);readline()
			}
		}
	}
	#auxlist$dist = round(auxdist,r)
	#auxlist$score = round(auxscore,r)
	#auxlist$score = mean(auxlist$dist)
	auxr$score = auxscore
	auxr$rmsd = round(auxrmsd,r)
	auxr$win = auxwin
	auxr$rot = round(auxrot,r)
	#auxr$rotref = round(rotref,r)
	auxr$sym = auxsym
	auxr$score$pain = verify_penalty(auxscore$distn,a1,a2)
	if (comp) auxr$score$gain = sum(auxscore$nocv) - sum(auxscore$nodv)
	#auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)
	if (comp) auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)+auxr$score$gain
	else auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)
	#auxr$score$sco1 = auxr$score$scot/(length(auxr$score$scov)*4)
	#auxzz = transform_by_rot_matrix(auxyz2,auxrr)
	if (!is.null(rotref)){
		auxr$supxyz = transform_by_rot_matrix(xyz2t,rotref)
	}else{
		auxr$supxyz = xyz2t
	} 
	#print(x);readline()
	#print(auxlist);readline()
	return(auxr)
}


exclude_sum_lines_matrix = function(v,a){
	
	auxa = a[v,]
	#print(!z);print(auxat);print(is.matrix(auxat));readline()
	if (!is.matrix(auxa)){
		auxa = t(as.matrix(auxa))
	}
	auxs = apply(auxa,1,sum)
	return(auxs)

}
verify_penalty = function(d,a1,a2,w=c(1,4),sep="-",k=1,noalign=0){

	auxn = names(d)

	auxn=as.numeric(unlist(strsplit(auxn,sep)))
	m = matrix(auxn,ncol=2,byrow=T)

	tam1 = dim(a1)[1]
	tam2 = dim(a2)[2]

	if (tam1<tam2){
		auxa1 = a1
		auxa2 = a2
	}else{
		auxa1 = a2
		auxa2 = a1
	}
	tam1 = dim(auxa1)[1]
	tam2 = dim(auxa2)[1]

	d = round(d,k)
	z = !(d>0)
	auxs1 = exclude_sum_lines_matrix(z,auxa1)
	auxa1s = sum(auxa1[upper.tri(auxa1,diag=T)])
	auxt1 = sum(auxs1)/auxa1s
	auxt2 = 0
	#print(auxs1);print(auxa1s);print(auxt1);readline()
	if (noalign){
		x = m[,2]
		y = 1:tam2
		z = !(y %in% x)
		auxs2 = exclude_sum_lines_matrix(z,auxa2)
		#z = c(F,F,T,T)
		#mt = upper.tri(auxa2,diag=T)
		#auxa2[!mt]=0
		#auxat = auxa2[!z,]
		#print(!z);print(auxat);print(is.matrix(auxat));readline()
		#if (!is.matrix(auxat)){
		#	auxat = t(as.matrix(auxat))
		#}
		#auxs1 = apply(auxat,1,sum)
	
		auxa2s = sum(auxa2[upper.tri(auxa2,diag=T)])
		auxt2 = sum(auxs2)/auxa2s
	}
	#print(auxn);print(auxa1);print(auxa2);print(z);print(sum(auxs2));print(auxs2);print(auxa2s);print(auxs);readline()

	return(c(auxt1,auxt2)*w)
}

resume_score = function(score,value=4,sep="_",rel=T){
#	colnames(auxtab)=c("n1","n2","rmsd","scorep","pain","scoret","sarea")

	#print(score);readline()
	tam = dim(score)[1]

	auxsc = score
	auxr = score[1,]
	#auxp1 = unlist(strsplit(as.character(auxsc[,"n1"]),sep))
	#auxp2 = unlist(strsplit(as.character(auxsc[,"n2"]),sep))
	###auxc2 = as.numeric(pdbsplit(auxsc[,"n2"])[,6])
	#print(auxc2)
	###auxv = auxc2[1]<=auxc2
	###auxc2[auxv]=auxc2[1]
	###auxc2 = auxc2*value
	#print(auxc2);readline()
	auxp = unlist(strsplit(as.character(auxsc[1,"n2"]),sep))[1]
	#print(auxc);print(auxp);readline()
	#auxw1 = score[1,3]
	#if (rel) auxr[1,3] = auxr[1,3]/auxw1
	#if (rel) auxr[1,3] = auxr[1,3]/auxw1
	#print(auxp)
	j=2
	for (i in 2:tam){
		auxi = unlist(strsplit(as.character(auxsc[i,"n2"]),sep))[1]
		if (!(auxi %in% auxp)){
			auxp = c(auxp,auxi)
			auxr = rbind(auxr,auxsc[i,])
			#if (rel) auxr[j,3] = auxr[j,3]/auxw1
			j = j + 1
		}	
		#print(auxr);readline()
	}
	#auxr$relocal = auxr[,3]/auxr[1,3]
	###auxc2 = as.numeric(pdbsplit(auxr[,"n2"])[,6])
	###auxv = auxc2[1]<=auxc2
	###auxc2[auxv]=auxc2[1]
	###auxc2 = auxc2*value
	aux = auxr[,"sarea"]
	auxc2 = rep(1,length(aux))
	#print(auxr);print(auxc2);readline()
	#auxr[,"sarea"] = auxr[,"scoret"]/auxr[1,"scoret"]
	#### TESTAR ELIMINACAO DO SCORE RELATIVO
	auxr[,"sarea"] = auxr[,"scoret"]/auxc2
	auxr$scglobal = aux
	#colnames(auxr)[3:4]=c("score","relocal")
	colnames(auxr)[c(2,7)]=c("pdbid","sclocal")
	auxr$pain = NULL
	auxr$n1 = NULL
	#print(auxr);readline()
	#dball$scorefinal = auxr
	#print(dball$scorefinal);readline()
	return(auxr)
}

remove_extra_spaces= function(text,split=""){

	#print(pdb$title)#;readline()
	extra = text
	#extra = "titleaaaaaaaaaaa   aa   title aaa"
	extra = toupper(extra)
	#print(extra)
	extra = paste(unlist(strsplit(extra,split)),collapse="")
	extra = gsub("\\s+", " ", str_trim(extra))
	#print(extra);readline()
	return(extra)

}

resume_similarity_assessment = function(dball,parg,score,pdblist,id=1){

	auxall = list()
	idtab = t(sapply(parg, "[[",3))#;print(idtab);readline()

	ids = unique(idtab[,2])
	tami = length(ids)
	
	for (i in 1:tami){
	#print(length(aux));readline()
		auxtab = data.frame()
		tamj=which(idtab[,2]==ids[i])
		for (j in tamj){
			stat = score[[j]]$score$local$scor[1]
			stat = c(stat,score[[j]]$score$global$scor[1])
			stat = c(stat,score[[j]]$score$scot)
			stat = c(stat,score[[j]]$score$ker)
			id0=idtab[j,1]#;print(id0)
			id1=idtab[j,2]#;print(id1)
			id2=idtab[j,3]#;print(id2);readline()
			auxtab0 = data.frame(id0,id1,id2,score[[j]]$names[1],score[[j]]$names[2],stat[1],stat[2],stat[3],stat[4])
			auxtab = rbind(auxtab,auxtab0)
		}
		colnames(auxtab)=c("id0","id1","id2","n1","n2","local","global","scoret","kernel")
		auxo = order(auxtab[,8],decreasing=T)
		auxscore = auxtab[auxo,]#;print(auxscore);readline()
		
		auxall[[i]] = list()
		auxall[[i]]$aligndata = score[tamj]
		auxall[[i]]$prescore = auxscore
		auxall[[i]]$prescore$scorea = 0
		
	}
	dball$res[[1]]=list()
	dball$res[[1]]$alignment = auxall
	dball$res[[1]]$pdbname = toupper(pdblist[id])
	dball$res[[1]] = new_consolidate_pdb_similarity_assessment(db,dball$res[[1]],pdblist)
	return(dball)

}

make_queue_for_par = function(db,dball,dbhot,groupid){

#	auxres = list()
	#pdbname = toupper(pdblist[idb])
	#auxres$pdbname = pdbname
	#auxgrp = which(dball$group[,"pdb"]==pdbname)
	#auxres$groupid = auxgrp
	#auxres$alignment = list()
	id = groupid
	parg12 = list()
	k = 1
	n = 1
	for (i in id){
		id0 = i
		id1 = n
		pdbid = as.numeric(as.vector(dball$group[i,"id"]))
		clusid = as.numeric(as.vector(dball$group[i,"nclus"]))#;print(pdbid);print(clusid)
		parg1 = list()
		parg1$xyz = db[[pdbid]]$exp$superclus[[clusid]]$geomc#;print(xyz1)
		parg1$g = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print("o");print(g1);print("o")
		parg1$a = db[[pdbid]]$exp$superclus[[clusid]]$a#;print(a1)
		parg1$auxname = as.character(dball$group[i,"fakename"])
		parg1$pdbname = as.character(dball$group[i,"pdb"])
		if (!is.null(dbhot)) parg1$hot = dbhot[[pdbid]]
		else parg1$hot = NULL

		#auxr = list()
		tamj = 1:dim(dball$group)[1]
		#tamj = 1:6

		for (j in tamj){
			parg2 = list()
			id2 = j
			pdbid = as.numeric(as.vector(dball$group[j,"id"]))
			clusid = as.numeric(as.vector(dball$group[j,"nclus"]))
			parg2$xyz = db[[pdbid]]$exp$superclus[[clusid]]$geomc
			parg2$g = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print("o");print(g2);print("o")
			parg2$a = db[[pdbid]]$exp$superclus[[clusid]]$a
			#sa2 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
			#sa = sa2/sa1#;print("o")
			parg2$auxname = as.character(dball$group[j,"fakename"])#;print(auxname1);print(auxname2)
			parg2$pdbname = as.character(dball$group[j,"pdb"])
			if (!is.null(dbhot)) parg2$hot = dbhot[[pdbid]]
			else parg2$hot = NULL
			parg12[[k]]=list()
			parg12[[k]]$parg1 = parg1
			parg12[[k]]$parg2 = parg2
			parg12[[k]]$ids = c(id0,id1,id2)
			k = k + 1
		}
		n=n+1
	}
	return(parg12)
	
}


# VERSAO NOVA QUE FICOU VELHA: usando funcao propria!!! ;-) de alinhamento para construir similaridade in memory
new_make_similarity_assessment = function(db,dball,id0,pdblist,dbhot=NULL,rotref=diag(4),self=F,w=c(1,1,1,1,1),r=3,verbose=T,outdir="Datalign/", sufix="align.csv",sep="_", comp=F){

	auxres = list()
	pdbname = toupper(pdblist[id0])
	auxres$pdbname = pdbname
	auxgrp = which(dball$group[,"pdb"]==pdbname)#;print(pdbname);print(id0);print(pdblist);print(dball$group)#;readline()
	auxres$groupid = auxgrp
	auxres$alignment = list()
	id = auxgrp#;print(id)

	auxall=list()
	#id = 5
	id.len = length(id)
	#print(id.len);readline()
	if (self){
#		tami = id[1]
#		tamj = id
		if (id.len>1){
			tami = id[1]:id[(id.len-1)]
		}
		else{
			print(paste("WARNING: id length for self similarty is insufficient in",pdbname))
			tami = id[1]
			#return(auxres)
		}
	}else{
		tami = id
		tamj = 1:dim(dball$group)[1]#;tam=3
	}
	k = 1
	#id = 1 # test 1ppf 5
	#id=1
	#tami = 5
	#tamj = 3
	#tami = 2
	#tamj = 8
	empty=F
	for (i in tami){
	#for (i in 15){
		#print(paste("i:",i))
		if (dball$group$empty[i]){
			#print("fudeu")
			empty=T
			break
		}
		parg1 = list()
		pdbid = as.numeric(as.vector(dball$group[i,"id"]))
		clusid = as.numeric(as.vector(dball$group[i,"nclus"]))#;print(pdbid);print(clusid)
		#if (is.null(db[[pdbid]]$A)) break;
		parg1$xyz = db[[pdbid]]$exp$superclus[[clusid]]$geomc#;print(xyz1)
		parg1$g = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print("o");print(g1);print("o")
		parg1$a = db[[pdbid]]$exp$superclus[[clusid]]$a#;print(a1)
		#sa1 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
		parg1$auxname = as.character(dball$group[i,"fakename"])
		parg1$pdbname = as.character(dball$group[i,"pdb"])
		if (!is.null(dbhot)) parg1$hot = dbhot[[pdbid]]
		else parg1$hot = NULL

		auxr = list()
		auxtab = data.frame()
		id0 = i
		id1 = k
		if (self) {
			if (id.len>1){
				tamj = c(i,i+1)
			}else{
				#tamj = 1 ### correcao de bug quando id.len == 1
				tamj = i	
			}
		}
		for (j in tamj){
		#for (j in 16){ #test 1oyv 7
		#for (j in 5){ #test 1acb 5
		#for (j in 2){ #test 1ppf 3
			#print(paste("j:",j));readline()
			parg2 = list()
			id2 = j
			pdbid = as.numeric(as.vector(dball$group[j,"id"]))
			clusid = as.numeric(as.vector(dball$group[j,"nclus"]))
			parg2$xyz = db[[pdbid]]$exp$superclus[[clusid]]$geomc
			parg2$g = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print("o");print(g2);print("o")
			parg2$a = db[[pdbid]]$exp$superclus[[clusid]]$a
			#sa2 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
			#sa = sa2/sa1#;print("o")
			parg2$auxname = as.character(dball$group[j,"fakename"])#;print(auxname1);print(auxname2)
			parg2$pdbname = as.character(dball$group[j,"pdb"])
			if (!is.null(dbhot)) parg2$hot = dbhot[[pdbid]]
			else parg2$hot = NULL
			#if (pre){
			#auxpdbname = dball$group$pdb[id]
			#auxv = (auxlist[,"id1"]==i)&(auxlist[,"id2"]==j)#;print("oi")
			#print(sum(auxv));readline()
			#return(auxlist)
			#auxr[[j]] = pre_best_align_graph_pdb(auxlist[auxv,],xyz1,xyz2,g1,g2,a1,a2,rotref,w)
			#if (verbose) print(sprintf("%s - %s ",auxname1,auxname2))
			#print("o")
			### auxr[[j]] = new_best_align_graph_pdb(xyz1,xyz2,g1,g2,a1,a2,pdbname1,pdbname2,hot1,hot2,rotref,w)#;print("o")
			#readline()
			if (self){
				#print(parg1);print(parg2);readline()
				auxr[[j]] = new_best_self_graph_pdb(parg1,parg2)
			}else{	
				auxr[[j]] = new_best_align_graph_pdb(parg1,parg2)#;print("o")	
			}

			##auxr[[k]] = new_best_align_graph_pdb(xyz1,xyz2,g1,g2,a1,a2,pdbname1,pdbname2,rotref,w)#;print("o")			
			#print(auxr);readline()
			#}else{
			#	auxr[[i]] = best_align_graph_pdb(auxname1,auxname2,xyz1,xyz2,g1,g2,a1,a2,w)
			#}
			#auxr[[j]]$names = c(auxname1,auxname2)
			##auxr[[k]]$names = c(auxname1,auxname2)
			stat=c()
			#print(auxr[[j]]$score$global$scor);readline()
			##stat = auxr[[k]]$score$local$scor[1]
			##stat = c(stat,auxr[[k]]$score$global$scor[1])
			##stat = c(stat,auxr[[k]]$score$scot)
			##stat = c(stat,auxr[[k]]$score$ker)
			stat = auxr[[j]]$score$local$scor[1]
			stat = c(stat,auxr[[j]]$score$global$scor[1])
			stat = c(stat,auxr[[j]]$score$scot)
			stat = c(stat,auxr[[j]]$score$ker)
			if (verbose) print(sprintf("%s - %s - %.3f - %.3f - %.3f - %.3f",parg1$auxname,parg2$auxname,stat[1],stat[2],stat[3],stat[4]))#;readline()
			auxtab0 = data.frame(id0,id1,id2,parg1$auxname,parg2$auxname,stat[1],stat[2],stat[3],stat[4])
			if (stat[3]==0) print(paste("WARNING: it was not possible to align",parg1$auxname,"-",parg2$auxname))
			auxtab = rbind(auxtab,auxtab0)
			#stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$scot,sa)#;print(stat);readline()
			#if (comp) {
			#	stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$gain,auxr[[j]]$score$scot,sa)
			#	if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[j]]$rmsd,stat[4]))
			#	auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,auxr[[j]]$rmsd,stat[1],stat[2],stat[3],stat[4],stat[5])
			#	auxtab = rbind(auxtab,auxtab0)
			#}else{
			#	stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$scot,sa)
			#	if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[j]]$rmsd,stat[3]))
			#	auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,auxr[[j]]$rmsd,stat[1],stat[2],stat[3],stat[4])
			#	auxtab = rbind(auxtab,auxtab0)
			#}
			#k = k+1
		}
	#}
		#dball$aligndata = auxr
		colnames(auxtab)=c("id0","id1","id2","n1","n2","local","global","scoret","kernel")
		rownames(auxtab)=NULL#;print(auxtab);readline()
		if (!self){
			auxo = order(auxtab[,8],decreasing=T)
			auxscore = auxtab[auxo,]#;prbyint(auxscore);readline()
		}else{
			auxscore = auxtab
		}
		#auxall[[1]] = list()
		#auxall[[1]]$aligndata = auxr
		#auxall[[1]]$prescore = auxscore
		#auxall[[1]]$prescore$scorea = 0
		auxall[[k]] = list()
		auxall[[k]]$aligndata = auxr
		auxall[[k]]$prescore = auxscore
		auxall[[k]]$prescore$scorea = 0
		
		#if (comp){
		#	colnames(auxtab)=c("id0","id1","id2","n1","n2","rmsd","scorep","pain","gain","scoret","sarea")
		#}else{
		#	colnames(auxtab)=c("id0","id1","id2","n1","n2","rmsd","scorep","pain","scoret","sarea")
		#}
		#rownames(auxtab)=NULL
		#auxo = auxtab
		#auxscore = auxo[order(auxo[,9]),c(4:5,9)] #by mean
		#auxscore = auxo[order(auxo[,9],decreasing=T),c(4:5,9)] #by mean
		#auxscore = auxo[order(auxo[,11],decreasing=T),c(4:5,11:12)] #by sum
		#auxscore = auxo[order(auxo[,"scoret"],decreasing=T),] #by sum
		#print(auxscore);readline()
		#dball$tab = auxtab
		#dball$tab = auxtab
		#dball$prescore = auxscore
		#dball$score = resume_score(score=dball$prescore)
		#auxall[[k]]$prescore = auxscore
		#auxall[[k]]$score = resume_score(score=auxall[[k]]$prescore)
		k=k+1
	}
	if (!empty) auxres$alignment = auxall
	#return(auxall)
	return(auxres)

}



# VERSAO ANTIGA: usava dados pre-computados em arquivos para construir similaridade
#pre_make_similarity_assessment(db=db,dball=dball,id=auxgrp,pdbname=pdbname)
pre_make_similarity_assessment = function(db,dball,id,pdbname,rotref,w=c(1,1,1,1,1),r=3,verbose=T,outdir="Datalign/",sufix="align.csv",sep="_", comp=F){


	filename = paste(outdir,toupper(pdbname),sep,sufix,sep="")
	print(paste("Reading pre-computation file:",filename))
	auxlist = read.table(filename,header=T,sep=",")
	#auxcoln = colnames(auxlist)
	auxlist = as.matrix(auxlist)
	#colnames(auxlist) = auxcoln
	#print(auxlist[1,]);readline()
	#print(colnames(auxlist));readline()
	#print(auxlist[1,])
	#print(dim(auxlist));readline()

	auxall=list()
	tam = dim(dball$group)[1]#;tam=3
	k = 1
	for (i in id){

		pdbid = as.numeric(as.vector(dball$group[i,"id"]))
		clusid = as.numeric(as.vector(dball$group[i,"nclus"]))#;print(pdbid);print(clusid)
		xyz1 = db[[pdbid]]$exp$superclus[[clusid]]$geomc#;print(xyz1)
		g1 = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print(g1)
		a1 = db[[pdbid]]$exp$superclus[[clusid]]$a#;print(a1)
		sa1 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
		auxname1 = as.character(dball$group[i,"fakename"])

		auxr = list()
		auxtab = data.frame()
		id0 = i
		id1 = k
		for (j in 1:tam){
			id2 = j
			pdbid = as.numeric(as.vector(dball$group[j,"id"]))
			clusid = as.numeric(as.vector(dball$group[j,"nclus"]))
			xyz2 = db[[pdbid]]$exp$superclus[[clusid]]$geomc
			g2 = db[[pdbid]]$exp$superclus[[clusid]]$ga
			a2 = db[[pdbid]]$exp$superclus[[clusid]]$a
			sa2 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
			sa = sa2/sa1
			auxname2 = as.character(dball$group[j,"fakename"])#;print(auxname1);print(auxname2)
			#if (pre){
			#auxpdbname = dball$group$pdb[id]
			auxv = (auxlist[,"id1"]==i)&(auxlist[,"id2"]==j)#;print("oi")
			#print(sum(auxv));readline()
			#return(auxlist)
			auxr[[j]] = pre_best_align_graph_pdb(auxlist[auxv,],xyz1,xyz2,g1,g2,a1,a2,rotref,w)
			#print(auxr);readline()
			#}else{
			#	auxr[[i]] = best_align_graph_pdb(auxname1,auxname2,xyz1,xyz2,g1,g2,a1,a2,w)
			#}
			auxr[[j]]$names = c(auxname1,auxname2)
			#stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$scot,sa)#;print(stat);readline()
			if (comp) {
				stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$gain,auxr[[j]]$score$scot,sa)
				if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[j]]$rmsd,stat[4]))
				auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,auxr[[j]]$rmsd,stat[1],stat[2],stat[3],stat[4],stat[5])
				auxtab = rbind(auxtab,auxtab0)
			}else{
				stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$scot,sa)
				if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[j]]$rmsd,stat[3]))
				auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,auxr[[j]]$rmsd,stat[1],stat[2],stat[3],stat[4])
				auxtab = rbind(auxtab,auxtab0)
			}
		}
	#}
		#dball$aligndata = auxr
		auxall[[k]] = list()
		auxall[[k]]$aligndata = auxr
		if (comp){
			colnames(auxtab)=c("id0","id1","id2","n1","n2","rmsd","scorep","pain","gain","scoret","sarea")
		}else{
			colnames(auxtab)=c("id0","id1","id2","n1","n2","rmsd","scorep","pain","scoret","sarea")
		}
		rownames(auxtab)=NULL
		auxo = auxtab
		#auxscore = auxo[order(auxo[,9]),c(4:5,9)] #by mean
		#auxscore = auxo[order(auxo[,9],decreasing=T),c(4:5,9)] #by mean
		#auxscore = auxo[order(auxo[,11],decreasing=T),c(4:5,11:12)] #by sum
		auxscore = auxo[order(auxo[,"scoret"],decreasing=T),] #by sum
		#print(auxscore);readline()
		#dball$tab = auxtab
		#dball$tab = auxtab
		#dball$prescore = auxscore
		#dball$score = resume_score(score=dball$prescore)
		auxall[[k]]$prescore = auxscore
		#auxall[[k]]$score = resume_score(score=auxall[[k]]$prescore)
		k=k+1
	}
	return(auxall)

}

# VERSAO ANTIGA: computava a similaridade in memory
#make_general_similarity_assessment = function(db,pdbid,clusid,minclusid=3,workdir="Pymol"){
make_similarity_assessment = function(db,dball,id,w=c(1,1,1,1,1),r=3,verbose=T,pre=T){


	auxall=list()
	tam = dim(dball$group)[1]#;tam=3

	pdbid = as.numeric(as.vector(dball$group[id,"id"]))
	clusid = as.numeric(as.vector(dball$group[id,"nclus"]))#;print(pdbid);print(clusid)
	xyz1 = db[[pdbid]]$exp$superclus[[clusid]]$geomc#;print(xyz1)
	g1 = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print(g1)
	a1 = db[[pdbid]]$exp$superclus[[clusid]]$a#;print(a1)
	sa1 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
	auxname1 = as.character(dball$group[id,"fakename"])

	auxr = list()
	auxtab = data.frame()


	for (i in 1:tam){

		pdbid = as.numeric(as.vector(dball$group[i,"id"]))
		clusid = as.numeric(as.vector(dball$group[i,"nclus"]))
		xyz2 = db[[pdbid]]$exp$superclus[[clusid]]$geomc
		g2 = db[[pdbid]]$exp$superclus[[clusid]]$ga
		a2 = db[[pdbid]]$exp$superclus[[clusid]]$a
		sa2 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
		sa = sa2/sa1
		auxname2 = as.character(dball$group[i,"fakename"])#;print(auxname1);print(auxname2)
		auxr[[i]] = best_align_graph_pdb(auxname1,auxname2,xyz1,xyz2,g1,g2,a1,a2,w)
		auxr[[i]]$names = c(auxname1,auxname2)
		stat = c(auxr[[i]]$score$scor,sum(auxr[[i]]$score$pain),auxr[[i]]$score$scot,sa)#;print(stat);readline()
		if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[i]]$rmsd,stat[3]))
		auxtab0 = data.frame(auxname1,auxname2,auxr[[i]]$rmsd,stat[1],stat[2],stat[3],stat[4])
		auxtab = rbind(auxtab,auxtab0)
		
	}
	#dball$aligndata = auxr
	auxall$aligndata = auxr
	colnames(auxtab)=c("n1","n2","rmsd","scorep","pain","scoret","sarea")
	rownames(auxtab)=NULL
	auxo = auxtab
	#auxscore = auxo[order(auxo[,9]),c(4:5,9)] #by mean
	#auxscore = auxo[order(auxo[,9],decreasing=T),c(4:5,9)] #by mean
	#auxscore = auxo[order(auxo[,11],decreasing=T),c(4:5,11:12)] #by sum
	auxscore = auxo[order(auxo[,"scoret"],decreasing=T),] #by sum
	#print(auxscore);readline()
	#dball$tab = auxtab
	#dball$tab = auxtab
	#dball$prescore = auxscore
	#dball$score = resume_score(score=dball$prescore)
	auxall$prescore = auxscore
	auxall$score = resume_score(score=auxall$prescore)
	return(auxall)

}


pdbsplit = function(v,sep="_",join=T){

	#print(v)
	if (!is.vector(v)){
		v = as.character(v)
	}	
	#v=v[v!=non]
	#print(v)
	vs = strsplit(v,sep)
	tam = length(vs[[1]])
	aux=matrix(unlist(vs),ncol=tam,byrow=T)#;print(dim(aux))
	if (join){
		auxpdb = apply(tomatrix(aux[,1:3]),1,paste,collapse=sep)
		aux[,1] = auxpdb
	}
	#;print(aux);readline()
	if (tam!=7){
		print(paste("WARNING: non-standard name for",vs[[1]]))
	#	aux = cbind(aux[,1:4],aux[,6:7])
	}
	#print(aux);readline()
	colnames(aux)=c("pdb","chain1","chain2","method","polarity","id","nclus")
	return(aux)
}


make_pdb_similarity_assessment_old = function(db,dball,pdbname,rotref,pre=T){

	auxr = list()
	pdbname = toupper(pdbname)
	auxr$pdbname = pdbname

	auxgrp = which(dball$group[,"pdb"]==pdbname)
	#print(auxgrp);readline()
	auxr$groupid = auxgrp

	auxr$alignment = list()
	k=1

	if (pre){
		print(paste("Making similarity alignment with pre-computation for",pdbname))
		auxr$alignment = pre_make_similarity_assessment(db=db,dball=dball,id=auxgrp,pdbname=pdbname,rotref=rotref)
	} else {
		for (i in auxgrp){
			print(paste("Making similarity alignment for",pdbname," - cluster",dball$group[i,"nclus"]))
			auxr$alignment[[k]] = make_similarity_assessment(db=db,dball=dball,id=i)
			k=k+1
		}
	}
			

	return(auxr)

}

make_index_tab_of_assessment = function(dbali,sep="_"){

	tami = length(dbali)
	auxpdbr = c()
	auxclun = c()
	for (i in 1:tami){
		#auxn1 = c(auxn1,dbali[[i]]$score[1,"pdbid"])
		#print(dbali[[i]]$score[,"pdbid"])
		auxpli = pdbsplit(dbali[[i]]$score[,"pdbid"])
		#auxpdb = auxpli[,1]
		auxpdb = auxpli[,"pdb"]  
		#auxclu = auxpli[,6]
		#auxr[[i]]=cbind(auxpdb,dbali[[i]]$score[,c("rmsd","sclocal","scglobal")])
		auxpdbr = cbind(auxpdbr,auxpdb)
		#auxclur = cbind(auxclur,auxclu)
		#print(auxr[[i]]);readline()
	}
	
	#print(auxpdbr);readline()
	tamr = dim(auxpdbr)[1]
	tamc = dim(auxpdbr)[2]
	auxin = c()
	auxn = c()
	for (i in 1:tamr){
		auxpdb = auxpdbr[i,1]
		#auxclu = auxclur[i,1]
		#auxn = c(auxn,auxpdb)
		ind = which(auxpdbr==auxpdb)%%tamr
		ind[ind==0]=tamr
		auxin = rbind(auxin,ind)
	}
	#print(auxin);print(pdbsplit(dbali[[1]]$score[1,"pdbid"]))
	for (i in 1:tami){
		auxpli = pdbsplit(dbali[[i]]$score[1,"pdbid"])
		#auxid = auxpli[5]
		#auxclu = auxpli[6]
		auxid = auxpli["id"]
		auxclu = auxpli["nclus"]
		auxclun = c(auxclun,paste(auxid,sep,auxclu,sep=""))
	}
	auxpdbn = auxpdbr[,1]
	#print(auxpdbr[,1])
	colnames(auxin) = auxclun
	rownames(auxin) = auxpdbn
	#auxa = list()
	#auxa$indn = auxr
	#auxa$indi = auxin
	return(auxin)

}


evaluate_best_consensus = function(db,dbres,tab,value=4,k=2,xo=3,w=c(1,1,0),edgw=T){

#g2 = db[[pdbid]]$exp$superclus[[clusid]]$ga
#a2 = db[[pdbid]]$exp$superclus[[clusid]]$a
#auxdist0 = find_nearest_nodes(xyz1,xyz2t0)#;print("oi")
#auxscore0 = make_score_graph_superimposition(auxdist0,g1,g2,a1,a2,w)#;print(auxscore0);readline()
#auxr$score$pain = verify_penalty(auxscore$distn,a1,a2)
#auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)
#verify_penalty = function(d,a1,a2,w=c(1,4),sep="-",k=1,noalign=0){

	tami = dim(tab)[1]
	aux=list()
	auxr = data.frame()
	auxa = data.frame()
	auxbase = 0
	auxonly=0
	#print(tab)
	if (tami>1){
		auxonly = length(unique(tab$n2))#;print(auxonly);readline()
		for (i in 1:(tami-1)){
		#for (i in 6){
			#par(mfrow=c(1,2))
			auxp1 = pdbsplit(tab[i,"n2"])#;print(auxp1)
			#pdbid1 = as.numeric(as.vector(auxp1[5]))
			#clusid1 = as.numeric(as.vector(auxp1[6]))
			pdbid1 = as.numeric(as.vector(auxp1["id"]))
			clusid1 = as.numeric(as.vector(auxp1["nclus"]))
			g1 = db[[pdbid1]]$exp$superclus[[clusid1]]$ga#;plot(g1)
			a1 = db[[pdbid1]]$exp$superclus[[clusid1]]$a
			id1 = as.numeric(as.vector(tab[i,c("id1","id2")]))#;print(id1)
			xyz1 = dbres$alignment[[id1[1]]]$aligndata[[id1[2]]]$supxyz

			#print(a1);readline()			
			for (j in (i+1):(tami)){
			#for (j in 10){
				#if (tami==1) j=j-1
				auxp2 = pdbsplit(tab[j,"n2"])#;print(auxp1);print(auxp2)
				#pdbid2 = as.numeric(as.vector(auxp2[5]))
				#clusid2 = as.numeric(as.vector(auxp2[6]))
				pdbid2 = as.numeric(as.vector(auxp2["id"]))
				clusid2 = as.numeric(as.vector(auxp2["nclus"]))
				if ((clusid1!=clusid2)|(auxonly==1)){
					g2 = db[[pdbid2]]$exp$superclus[[clusid2]]$ga#;plot(g2)
					a2 = db[[pdbid2]]$exp$superclus[[clusid2]]$a
					id2 = as.numeric(as.vector(tab[j,c("id1","id2")]))#;print(id2)
					xyz2 = dbres$alignment[[id2[1]]]$aligndata[[id2[2]]]$supxyz
					auxdist = find_nearest_nodes(xyz1,xyz2)#;print(auxdist);print(a1);print(a2)
					auxscore = make_score_graph_superimposition(dist=auxdist,g1=g1,g2=g2,a1=a1,a2=a2,xyz1,xyz2,k=k,xo=xo)#;print(auxscore);readline()
					auxscore$pain = verify_penalty(auxscore$distn,a1,a2,noalign=1)
					if (edgw) {
						auxscore$scot = auxscore$scor-sum(auxscore$pain)+sum(auxscore$edwv)
					}else{
						auxscore$scot = auxscore$scor-sum(auxscore$pain)
					}
					#auxri = data.frame(tab[i,],scorei=auxscore$scot)
					#auxrj = data.frame(tab[j,],scorei=auxscore$scot)
					#auxscorea = tab[i,"scoret"]+tab[j,"scoret"]+auxscore$scot
					#auxscorea = tab[i,"scoret"]+tab[j,"scoret"]+tab[i,"scoreo"]+tab[j,"scoreo"]+auxscore$scot
					#auxscorea = tab[i,"scoret"]+tab[j,"scoret"]+tab[i,"scoreo"]+tab[j,"scoreo"]#+auxscore$scot
					#### ATE i=16 para dball[[i]]					
					#if (auxscorea>auxbase){
					####	auxri = data.frame(tab[i,],scorei=auxscore$scot)
					####	auxrj = data.frame(tab[j,],scorei=auxscore$scot)
					####	auxri = data.frame(auxri,scorea=auxscorea)
					####	auxrj = data.frame(auxrj,scorea=auxscorea)
					####	auxr = rbind(auxri,auxrj)
					####	auxa = rbind(auxa,auxr)
						#auxbase=auxscorea
					#}
						auxri = data.frame(tab[i,],scorei=auxscore$scot)
						auxscorea = w[1]*tab[i,"scoret"]+w[2]*tab[i,"scoreo"]+w[3]*auxscore$scot
						#print(auxri);print(auxscorea);readline()
						auxri = data.frame(auxri,scorea=auxscorea)
						auxrj = data.frame(tab[j,],scorei=auxscore$scot)
						auxscorea = w[1]*tab[j,"scoret"]+w[2]*tab[j,"scoreo"]+w[3]*auxscore$scot
						auxrj = data.frame(auxrj,scorea=auxscorea)
						auxr = rbind(auxri,auxrj)
						auxa = rbind(auxa,auxr)
					#;print(auxa);readline()
					#plot(g1,layout=xyz1);plot(g2,layout=xyz2)
					#print(auxscore$dist);print(auxscore$scot);readline()
					#print(a2);readline()
					
				}
				
			}
		}
	}else if (tami==1){
		### TO DO HERE ###
		print(paste("WARNING: it was found a consensus evaluating of a single input in",tab[1,"n1"],"Please, check it!"))
		auxscorea = tab[1,"scoret"]
		auxri = data.frame(tab[1,],scorei=auxscorea)
		auxrj = data.frame(tab[1,],scorei=auxscorea)
		auxri = data.frame(auxri,scorea=auxscorea)
		auxrj = data.frame(auxrj,scorea=auxscorea)
		auxr = rbind(auxri,auxrj)
	}else{
		print(paste("WARNING: something is wrong with consensus evaluating of",tab[1,"n1"]))
	}

	auxo = order(auxa[,"scorea"],decreasing=T)
	auxa = auxa[auxo,]
	#print(auxa);readline()
	auxr=auxa
	return(auxr)
}

found_best_pair_consensus = function(db,dbres,tab,id,value=4,k=2,xo=3,w=c(1,1,0),edgw=F){
#auxr=found_best_pair_consensus(db,dbres,auxr,1)

#g2 = db[[pdbid]]$exp$superclus[[clusid]]$ga
#a2 = db[[pdbid]]$exp$superclus[[clusid]]$a
#auxdist0 = find_nearest_nodes(xyz1,xyz2t0)#;print("oi")
#auxscore0 = make_score_graph_superimposition(auxdist0,g1,g2,a1,a2,w)#;print(auxscore0);readline()
#auxr$score$pain = verify_penalty(auxscore$distn,a1,a2)
#auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)
#verify_penalty = function(d,a1,a2,w=c(1,4),sep="-",k=1,noalign=0){

	tami = dim(tab)[1]
	aux=list()
	auxr = data.frame()
	auxa = data.frame()
	auxbase = 0
	auxonly=0
	#print(tab)
	if (tami>1){
		auxonly = length(unique(tab$n2))#;print(auxonly);readline()
		#for (i in 1:(tami-1)){
		for (i in id){
		#for (i in 6){
			#par(mfrow=c(1,2))
			auxp1 = pdbsplit(tab[i,"n2"])#;print(auxp1)
			#pdbid1 = as.numeric(as.vector(auxp1[5]))
			#clusid1 = as.numeric(as.vector(auxp1[6]))
			pdbid1 = as.numeric(as.vector(auxp1["id"]))
			clusid1 = as.numeric(as.vector(auxp1["nclus"]))
			g1 = db[[pdbid1]]$exp$superclus[[clusid1]]$ga#;plot(g1)
			a1 = db[[pdbid1]]$exp$superclus[[clusid1]]$a
			id1 = as.numeric(as.vector(tab[i,c("id1","id2")]))#;print(id1)
			xyz1 = dbres$alignment[[id1[1]]]$aligndata[[id1[2]]]$supxyz

			#print(a1);readline()			
			for (j in (i+1):(tami)){
			#for (j in 10){
				#if (tami==1) j=j-1
				auxp2 = pdbsplit(tab[j,"n2"])#;print(auxp1);print(auxp2)
				#pdbid2 = as.numeric(as.vector(auxp2[5]))
				#clusid2 = as.numeric(as.vector(auxp2[6]))
				pdbid2 = as.numeric(as.vector(auxp2["id"]))
				clusid2 = as.numeric(as.vector(auxp2["nclus"]))
				if ((clusid1!=clusid2)|(auxonly==1)){
					g2 = db[[pdbid2]]$exp$superclus[[clusid2]]$ga#;plot(g2)
					a2 = db[[pdbid2]]$exp$superclus[[clusid2]]$a
					id2 = as.numeric(as.vector(tab[j,c("id1","id2")]))#;print(id2)
					xyz2 = dbres$alignment[[id2[1]]]$aligndata[[id2[2]]]$supxyz
					auxdist = find_nearest_nodes(xyz1,xyz2)#;print(auxdist);print(a1);print(a2)
					auxscore = make_score_graph_superimposition(dist=auxdist,g1=g1,g2=g2,a1=a1,a2=a2,xyz1,xyz2,k=k,xo=xo)#;print(auxscore);readline()
					#auxscore$pain = verify_penalty(auxscore$distn,a1,a2,noalign=1)
					auxscore$pain = 0
					if (edgw) {
						auxscore$scot = auxscore$scor-sum(auxscore$pain)+sum(auxscore$edwv)
					}else{
						auxscore$scot = auxscore$scor-sum(auxscore$pain)
					}
					#auxri = data.frame(tab[i,],scorei=auxscore$scot)
					#auxrj = data.frame(tab[j,],scorei=auxscore$scot)
					#auxscorea = tab[i,"scoret"]+tab[j,"scoret"]+auxscore$scot
					#auxscorea = tab[i,"scoret"]+tab[j,"scoret"]+tab[i,"scoreo"]+tab[j,"scoreo"]+auxscore$scot
					#auxscorea = tab[i,"scoret"]+tab[j,"scoret"]+tab[i,"scoreo"]+tab[j,"scoreo"]#+auxscore$scot
					#### ATE i=16 para dball[[i]]					
					if (auxscore$scot[1]>auxbase){
						#auxri = data.frame(tab[i,],scorei=auxscore$scot)
						#auxrj = data.frame(tab[j,],scorei=auxscore$scot)
						tab[i,"scorei"]=auxscore$scot[1]
						auxri = tab[i,]
						tab[j,"scorei"]=auxscore$scot[1]
						auxrj = tab[j,]
					####	auxri = data.frame(auxri,scorea=auxscorea)
					####	auxrj = data.frame(auxrj,scorea=auxscorea)
						auxr = rbind(auxri,auxrj)
					####	auxa = rbind(auxa,auxr)
						auxbase=auxscore$scot[1]
					}
						#auxri = data.frame(tab[i,],scorei=auxscore$scot)
						#auxscorea = w[1]*tab[i,"scoret"]+w[2]*tab[i,"scoreo"]+w[3]*auxscore$scot
						#print(auxri);print(auxscorea);readline()
						#auxri = data.frame(auxri,scorea=auxscorea)
						#auxrj = data.frame(tab[j,],scorei=auxscore$scot)
						#auxscorea = w[1]*tab[j,"scoret"]+w[2]*tab[j,"scoreo"]+w[3]*auxscore$scot
						#auxrj = data.frame(auxrj,scorea=auxscorea)
						#auxr = rbind(auxri,auxrj)
						#auxa = rbind(auxa,auxr)
					#;print(auxa);readline()
					#plot(g1,layout=xyz1);plot(g2,layout=xyz2)
					#print(auxscore$dist);print(auxscore$scot);readline()
					#print(a2);readline()
					
				}
				
			}
		}
	}else if (tami==1){
		### TO DO HERE ###
		print(paste("WARNING: it was found a consensus evaluating of a single input in",tab[1,"n1"],"Please, check it!"))
		auxscorea = tab[1,"scoret"]
		tab[1,"scorei"]=auxscorea
		auxri = tab[1,]
		tab[1,"scorei"]=auxscorea
		auxrj = tab[1,]
		#auxri = data.frame(auxri,scorea=auxscorea)
		#auxrj = data.frame(auxrj,scorea=auxscorea)
		auxr = rbind(auxri,auxrj)
		#auxa = rbind(auxri,auxrj)
	}else{
		print(paste("WARNING: something is wrong with consensus evaluating of",tab[1,"n1"]))
	}

	#auxo = order(auxa[,"scorea"],decreasing=T)
	#auxa = auxa[auxo,]
	#print(auxa);readline()
	#auxr=auxa
	return(auxr)
}


evaluate_origin_consensus = function(db,dbres,tab,value=4,k=2,xo=3,comp=F){

	tami = dim(tab)[1]
	aux=list()
	#auxr = data.frame()
	auxr = c()
	auxa = data.frame()
	#auxscoro = c()
	auxbase = 0
	auxonly=0
	#print(tab)
	auxn1 = unique(tab$n1)
	tamj = length(auxn1)
	#print(auxn1)#;readline()
	if (tami>1){
		#auxonly = length(unique(tab$n2))#;print(auxonly);readline()
		#for (i in 1:(tami-1)){
		for (i in 1:tami){
		#for (i in 6){
			#par(mfrow=c(1,2))
			auxp1 = pdbsplit(tab[i,"n2"])#;print(auxp1)
			#pdbid1 = as.numeric(as.vector(auxp1[5]))
			#clusid1 = as.numeric(as.vector(auxp1[6]))
			pdbid1 = as.numeric(as.vector(auxp1["id"]))
			clusid1 = as.numeric(as.vector(auxp1["nclus"]))
			g1 = db[[pdbid1]]$exp$superclus[[clusid1]]$ga#;plot(g1)
			a1 = db[[pdbid1]]$exp$superclus[[clusid1]]$a
			id1 = as.numeric(as.vector(tab[i,c("id1","id2")]))#;print(id1)
			xyz1 = dbres$alignment[[id1[1]]]$aligndata[[id1[2]]]$supxyz
			auxr = c()
			#print(a1);readline()			
			for (j in 1:tamj){
			#for (j in 10){
				#if (tami==1) j=j-1
				auxp2 = pdbsplit(auxn1[j])#;print(auxp1);print(auxp2)
				#pdbid2 = as.numeric(as.vector(auxp2[5]))
				#clusid2 = as.numeric(as.vector(auxp2[6]))
				pdbid2 = as.numeric(as.vector(auxp2["id"]))
				clusid2 = as.numeric(as.vector(auxp2["nclus"]))
				#if ((clusid1!=clusid2)){#|(auxonly==1)){
				#print("---");print(auxn1[j]);print(tab[i,"n1"]);print("---")
				if (auxn1[j]!=tab[i,"n1"]){
					#print("---");print(as.character(tab[i,"n1"]));print(as.character(auxn1[j]))
					#print("entrei");print(auxp1);print(auxp2);print("---")
					g2 = db[[pdbid2]]$exp$superclus[[clusid2]]$ga#;plot(g2)
					a2 = db[[pdbid2]]$exp$superclus[[clusid2]]$a
					idn2 = which(auxn1[j]==tab$n1)#;print(idn2);readline()
					#id2 = as.numeric(as.vector(tab[idn2[1],c("id1","id2")]))#;print(id2)
					id2 = as.numeric(as.vector(tab[idn2[1],c("id1","id0")]))#;print(id2)
					#if (id2[1]!=id2[2]) {
					#	print(paste("WARNING: check idn2 in origin consensus for",auxp2));readline()
					#}
					#id2 = as.numeric(as.vector(tab[j,c("id1","id1")]));print(id1);print(id2);readline()
					xyz2 = dbres$alignment[[id2[1]]]$aligndata[[id2[2]]]$supxyz
					auxdist = find_nearest_nodes(xyz1,xyz2);#print(auxdist);print(a1);print(a2);readline()
					#make_score_graph_superimposition = function(dist,g1,g2,a1,a2,xyz1,xyz2,w=c(1,1,1,1,1),consol=5,sep="-",k=2,xo=5,comp=F,edgw=F){
					auxscore = make_score_graph_superimposition(dist=auxdist,g1=g1,g2=g2,a1=a1,a2=a2,xyz1,xyz2,k=k,xo=xo)#;print(auxscore);readline()
					#auxscore$pain = verify_penalty(auxscore$distn,a1,a2)
					auxscore$pain = 0
					if (comp) auxscore$gain = sum(auxscore$nocv) - sum(auxscore$nodv)
					#auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)
					if (comp) auxscore$scot = auxscore$scor-sum(auxscore$pain)+auxscore$gain
					else auxscore$scot = auxscore$scor-sum(auxscore$pain)
					auxr = c(auxr,auxscore$scot)#;print(auxr);readline()
					#auxri = data.frame(tab[i,],scorei=auxscore$scot)
					#auxrj = data.frame(tab[j,],scorei=auxscore$scot)
					#auxscoro = tab[i,"scoret"]+tab[j,"scoret"]+auxscore$scot
					#if (auxscorea>auxbase){
					#	auxri = data.frame(tab[i,],scorei=auxscore$scot)
					#	auxrj = data.frame(tab[j,],scorei=auxscore$scot)
					#	auxri = data.frame(auxri,scorea=auxscorea)
					#	auxrj = data.frame(auxrj,scorea=auxscorea)
					#	auxr = rbind(auxri,auxrj)
					#	auxa = rbind(auxa,auxr)
					#	#auxbase=auxscorea
					#}
					#;print(auxa);readline()
					#plot(g1,layout=xyz1);plot(g2,layout=xyz2)
					#print(auxscore$dist);print(auxscore$scot);readline()
					#print(a2);readline()
					
				}
				
			}
			#print(auxr);readline()
			auxri = data.frame(tab[i,],scoreo = sum(auxr))
			auxa = rbind(auxa,auxri)
			#	print(auxa);readline()
		}
	}else if (tami==1){
		print(paste("WARNING: it was found a origin consensus evaluating of a single input in",tab[1,"n1"],"Please, check it!"))
		#auxscoreo = tab[1,"scoret"]
		auxa = data.frame(tab[1,],scoreo=0)
		#auxrj = data.frame(tab[1,],scorei=auxscorea)
		#auxri = data.frame(auxri,scorea=auxscorea)
		#auxrj = data.frame(auxrj,scorea=auxscorea)
		#auxr = rbind(auxri,auxrj)
	}else{
		print(paste("WARNING: something is wrong with origin consensus evaluating of",tab[1,"n1"]))
	}

	#auxo = order(auxa[,"scorea"],decreasing=T)
	#auxa = auxa[auxo,]
	#print(auxa);readline()
	#auxr=auxa
	return(auxa)

}

making_ranking_consensus = function(tab,rto=T,rankto=F,levelto=T){

	auxsct = tab$scoret
	auxsco = tab$scoreo
	auxsci = tab$scorei		
	

	if (rto){
		auxsct = roundto(auxsct)
		auxsco = roundto(auxsco)
		auxsci = roundto(auxsci)		

	}

	if (levelto){
		auxsct = factor(auxsct)
		auxsco = factor(auxsco)
		auxsci = factor(auxsci)
		levels(auxsct)=length(levels(auxsct)):1
		levels(auxsco)=length(levels(auxsco)):1
		levels(auxsci)=length(levels(auxsci)):1
		auxsct = as.numeric(as.vector(auxsct))
		auxsco = as.numeric(as.vector(auxsco))
		auxsci = as.numeric(as.vector(auxsci))			
	}	

	#else{
	#	auxsct = rank(tab$scoret)
	#	auxsco = rank(tab$scoreo)
	#	auxsci = rank(tab$scorei)
	#}
	if (rankto){
		auxsct = rank(auxsct)
		auxsco = rank(auxsco)
		auxsci = rank(auxsci)		
	}

	#auxm = cbind(rscoret=auxsct,rscoreo=auxsco,rscorei=auxsci)
	auxm = cbind(rscoret=auxsct,rscoreo=auxsco)#,rscorei=auxsci)

	auxscr = apply(auxm,1,sum)
	
	tab = cbind(tab,auxm,rscorea=auxscr)

	if (levelto){
		auxo = order(tab[,"rscorea"])
		tab = tab[auxo,]
	}else{
		auxo = order(tab[,"rscorea"],decreasing=T)
		tab = tab[auxo,]
	}

	#print(tab);readline()	

	return(tab)

}

reduce_for_one_pdb = function(tab){

	#pdblist = unique(pdbsplit(tab$n2)[,1])
	pdblist = unique(pdbsplit(tab$n2)[,"pdb"])
	auxr = data.frame()


	for (pdbname in pdblist){
		id = grep(pdbname,tab$n2)
		if (length(id)<1) print(paste("WARNING: something is wrong in reduce pdb tab for",pdbname))
		auxr = rbind(auxr,tab[id[1],])
	}
	return(auxr)

}

join_prescore_by_pdb = function(dbres,pdbname){


	tami = length(dbres$alignment)#;print(dbres$alignment)
	auxa = data.frame()
	#print(tami);readline()
	for (i in 1:tami) {
		auxid=grep(pdbname,dbres$alignment[[i]]$prescore$n2,fixed=T)
		#print(pdbname);print(dbres$alignment[[i]]$prescore);print(auxid);readline()
		if (length(auxid)<1){
			print(paste("WARNING: something is wrong with join prescore for",pdbname))
		}
		#auxn2 = rownames(dbres$alignment[[i]]$prescore)[auxid]### PAREI AQUI
		
		auxr = cbind(dbres$alignment[[i]]$prescore[auxid,],auxid)
		auxa = rbind(auxa,auxr)
	}
	#print(auxa);readline()
	return(auxa)
}

new_consolidate_pdb_self_similarity_assessment = function(db,dball,pdblist,debug=F){

	dbsel = dball$res
	tamj = length(dbsel)
	tami = length(pdblist)
	auxr = data.frame()
	auxa = data.frame()
	auxl = list()

	if (tamj){
		i=1
		#for (i in 1){
		for (j in 1:tamj){#;print(j)#;print(tamj)
		#for (j in 3){
			pdbnamei = pdblist[j]#;print(pdbnamei)#;readline()
			#print(paste("NEW Consolidating similarity between",toupper(dbres$pdbname),"-",toupper(pdbnamei)))
			#print(head(dbsel[[j]]$alignment));readline()
			if (length(dbsel[[j]]$alignment)!=0) {
				aux1 = join_prescore_by_pdb(dbsel[[j]],pdbnamei)#;print(aux1);readline()
				aux1.len = dim(aux1)[1] ########## REVER AQUI ###########
				id = seq(1,aux1.len,as.integer(sqrt(aux1.len)))#;print(id);readline()
				aux1 = aux1[id,];
				dball$res[[j]]$scorend = aux1
				dball$res[[j]]$scorend$auxid = NULL
				#auxo = order(auxr$scoret,decreasing=T)[1:id]
				#auxa = rbind(auxa,auxr[auxo,]);print(auxa);readline()
				auxa = rbind(auxa,aux1)
			}
		}
		auxa$auxid = NULL
		#print(auxa);readline()
		dball$scoreself = auxa
	}
	#return(auxa)
	return(dball)

}


new_consolidate_pdb_similarity_assessment = function(db,dbres,pdblist,debug=F){

	
	tami = length(pdblist)
	auxr = data.frame()
	auxa = data.frame()
	auxl = list()
	
	for (i in 1:tami){
	#for (i in 2){
		pdbnamei = pdblist[i]
		print(paste("NEW Consolidating similarity between",toupper(dbres$pdbname),"-",toupper(pdbnamei)))
		auxr = join_prescore_by_pdb(dbres,pdbnamei)#;print(auxr);readline()
		auxo = order(auxr$scoret,decreasing=T)[1]
		auxa = rbind(auxa,auxr[auxo,])#;print(auxa);readline()
	}
	auxo = order(auxa$scoret,decreasing=T)
	auxa = auxa[auxo,]#;print(auxa);readline()
	auxa$auxid = NULL
	dbres$scorend = auxa
	return(dbres)
}


#		auxr = join_prescore_by_pdb(dbres,pdbnamei)#;print(auxr);readline()
		#print(pdbnamei);print(auxr);readline()
#		auxr = evaluate_origin_consensus(db,dbres,auxr)
#		auxr$scorea = (auxr$scoret+auxr$scoreo)/length(dbres$groupid)
#		print(auxr);readline()
#		if (debug) {
#			auxr0 = evaluate_best_consensus(db,dbres,auxr)
#			auxr0[,"scorea"] = (auxr0$scoret+auxr0$scoreo)/length(dbres$groupid)
#			auxl[[i]] = auxr0
#		}
#		auxr = data.frame(auxr,scorei=0,scorea=(auxr$scoret+auxr$scoreo)/length(dbres$groupid))
#		auxo = order(auxr$scorea,decreasing=T)
#		auxr = auxr[auxo,]
#		print(auxr);readline()
		#if (debug) auxr0 = evaluate_best_consensus(db,dbres,auxr)
		#print(pdbnamei);print(x);readline()
		#auxr = making_ranking_consensus(auxr)
		#print(auxr);readline()
		#if (debug) auxl[[i]] = auxr0
#		auxr=found_best_pair_consensus(db,dbres,auxr,1)
		#print(auxr);readline()
		#auxa = rbind(auxa,auxr[1:2,])
#		auxa = rbind(auxa,auxr)
#	} 
	#print(auxa);readline()
	#auxo = order(auxa$scorea,decreasing=T)
#	auxo = order(auxa$scoret,decreasing=T)
#	auxa = auxa[auxo,]
	#print(auxa);readline()
#	auxone=reduce_for_one_pdb(auxa)
	#print(auxone);readline()
	#dbres$scorend = auxa[auxo,]
#	dbres$scorend = auxone
#	dbres$allend = auxl
#	return(dbres)
#}


consolidate_pdb_similarity_assessment = function(db,dbres,pdblist,debug=F){

	
	tami = length(pdblist)
	auxr = data.frame()
	auxa = data.frame()
	auxl = list()
	
	for (i in 1:tami){
	#for (i in 2){
		pdbnamei = pdblist[i]
		print(paste("Consolidating similarity between",toupper(dbres$pdbname),"-",toupper(pdbnamei)))
		auxr = join_prescore_by_pdb(dbres,pdbnamei)#;print(auxr);readline()
		#print(pdbnamei);print(auxr);readline()
		auxr = evaluate_origin_consensus(db,dbres,auxr)
		auxr$scorea = (auxr$scoret+auxr$scoreo)/length(dbres$groupid)
		#print(auxr);readline()
		if (debug) {
			auxr0 = evaluate_best_consensus(db,dbres,auxr)
			auxr0[,"scorea"] = (auxr0$scoret+auxr0$scoreo)/length(dbres$groupid)
			auxl[[i]] = auxr0
		}
#		auxr = data.frame(auxr,scorei=0,scorea=(auxr$scoret+auxr$scoreo)/length(dbres$groupid))
#		auxo = order(auxr$scorea,decreasing=T)
#		auxr = auxr[auxo,]
#		print(auxr);readline()
		#if (debug) auxr0 = evaluate_best_consensus(db,dbres,auxr)
		#print(pdbnamei);print(x);readline()
		#auxr = making_ranking_consensus(auxr)
		#print(auxr);readline()
		#if (debug) auxl[[i]] = auxr0
#		auxr=found_best_pair_consensus(db,dbres,auxr,1)
		#print(auxr);readline()
		#auxa = rbind(auxa,auxr[1:2,])
		auxa = rbind(auxa,auxr)
	} 
	#print(auxa);readline()
	#auxo = order(auxa$scorea,decreasing=T)
	auxo = order(auxa$scoret,decreasing=T)
	auxa = auxa[auxo,]
	#print(auxa);readline()
	auxone=reduce_for_one_pdb(auxa)
	#print(auxone);readline()
	#dbres$scorend = auxa[auxo,]
	dbres$scorend = auxone
	dbres$allend = auxl
	return(dbres)
}


resume_for_web_output_old3 = function (db,dbsel,dbgrp,id,k=0){

	dbhot = list()

	tami = length(dbsel$hot)

	#for (i in 1:tami){
		tamj = mapply(dim,dbsel$hot[[id]])[1,]
		#dbhot[[]] = list()
		for (j in tamj){
			dbhot[[j]] = list()
			m = dbgrp[[i]]$connections_group_bsr[[j]]
			dbhot[[j]]$matrix = filter(m,source!=target)
			#;print(m_f);readline()
		}
	#}

	return(dbhot)
		



}


resume_for_web_output_old = function (db,dball,id,k=0){

#> connetions_group_bsr[3]
#[[1]]
#  source source_volume target target_volume edge_sum      s_x      s_y      s_z      t_x      t_y      t_z
#1      1      5553.502      2     19821.187 1646.727 30.46871 69.08758 11.42233 27.80587 54.77508 12.79056
#2      1      5553.502      3      2085.385  147.141 30.46871 69.08758 11.42233 38.81296 61.25319 19.91191
#3      2     19821.187      3      2085.385  164.266 27.80587 54.77508 12.79056 38.81296 61.25319 19.91191

#auxp2 = pdbsplit(auxn1[j])#;print(auxp1);print(auxp2)
#pdbid2 = as.numeric(as.vector(auxp2[5]))
#clusid2 = as.numeric(as.vector(auxp2[6]))
#auxt1 = ends(auxg1,E(auxg1))
#dball[[35]]$res[[1]]$alignment[[3]]$aligndata[[3]]$supxyz
#matrix(dball[[1]]$align[[51]]$rot,ncol=4,byrow=T)
#find_nearest_nodes
#c(5,x[!(x %in% 5)])

	auxl = list()
	if (!is.null(dball$res[[id]]$scorend)){
		cols = brewer.pal(n = 9, name = 'YlOrRd')
		tami = dim(dball$res[[id]]$scorend)[1]
		#auxr = data.frame()
		#auxa = data.frame()
		#auxl = list()
		print(paste("Resuming for web output..."))
		pdbname0 = dball$res[[id]]$pdbname
		#i0all = grep(pdbname0,dball$res[[id]]$scorend$n2)[1]
		i0 = grep(pdbname0,dball$res[[id]]$scorend$n2)[1]
		#auxall = list()
		#i0k = 1
		#for (i0 in i0all){
		#auxl = list()
		#tami = dim(dball$res[[id]]$scorend)[1]
		al01 = dball$res[[id]]$scorend$id1[i0]
		al02 = dball$res[[id]]$scorend$id2[i0]
		xyz0 = dball$res[[id]]$alignment[[al01]]$aligndata[[al02]]$supxyz
		#xyz0 = transform_by_rot_matrix(xyz0,rotref)
		#i0 = 3
		ilist = 1:tami
		ilist = c(i0,ilist[!(ilist %in% i0)])#;print(ilist)
		#print(i0);print(al01);print(al02);print(xyz0);print(ilist);readline()
		#for (i in 1:tami){#print(dball$res[[id]]$scorend$n2[1])
		for (i in ilist){
			auxp = pdbsplit(dball$res[[id]]$scorend$n2[i])#;print(auxp);readline()		
			#pdbid = as.numeric(as.vector(auxp[5]))
			#clusid = as.numeric(as.vector(auxp[6]))
			pdbid = as.numeric(as.vector(auxp[,"id"]))#;print(pdbid)
			clusid = as.numeric(as.vector(auxp[,"nclus"]))#;print(clusid)
			a = round(db[[pdbid]]$exp$superclus[[clusid]]$a,k)
			g = db[[pdbid]]$exp$superclus[[clusid]]$ga
			e = ends(g,E(g))
			e = cbind(e,round(E(g)$weight,k))
			al1 = dball$res[[id]]$scorend$id1[i]
			al2 = dball$res[[id]]$scorend$id2[i]#;print(al1);print(al2)#;readline()
			xyz = dball$res[[id]]$alignment[[al1]]$aligndata[[al2]]$supxyz
			ngroup0 = dim(xyz)[1]
			if (!is.null(db[[pdbid]]$exp$preclus)){
				ngroup = ngroup0 + length(db[[pdbid]]$exp$preclus$ids)
			}#;print(i)
			#xyz = transform_by_rot_matrix(xyz,rotref)
			if (as.character(dball$res[[id]]$scorend$n1[i])==as.character(dball$res[[id]]$scorend$n2[i])){
				rot = db[[pdbid]]$exp$rotref
			}else{
				rot = dball$res[[id]]$alignment[[al1]]$aligndata[[al2]]$rot
			}
			#plot(g,layout=xyz);print(a);print(e);print(xyz);readline()
			tamj = dim(e)[1]
			auxa = data.frame()
			for (j in 1:tamj){
				auxr = data.frame(source=e[j,1],source_volume=a[e[j,1],e[j,1]],target=e[j,2],target_volume=a[e[j,2],e[j,2]],edge_sum=e[j,3])
				xyzj = c(xyz[e[j,1],],xyz[e[j,2],])
				auxr = data.frame(auxr,s_x=xyzj[1],s_y=xyzj[2],s_z=xyzj[3],t_x=xyzj[4],t_y=xyzj[5],t_z=xyzj[6])
				auxa = rbind(auxa,auxr)
				#print(auxa);readline()
			}
			#x = round(find_nearest_nodes(xyz0,xyz),k)
			#print(as.character(dball$res[[id]]$scorend[i,5]));print(x);readline()

			auxl[[i]]=list()#; print(i);print(j)
			auxl[[i]]$matrix = auxa#;print(auxa)
			auxl[[i]]$scorend = dball$res[[id]]$scorend[i,]
			auxl[[i]]$rot = matrix(rot,ncol=4,byrow=F) #ATENCAO!!! mudou byrow de T para F
			auxl[[i]]$ngroup = ngroup
			nhot = sapply(dball$hot[[id]],dim)#;print(nhot)#;readline()
			idhot=which(ngroup0==nhot[1,])#;print(ngroup0);print(id);print(idhot);readline()
			auxl[[i]]$hot = dball$hot[[id]][[idhot]]$score#;print(auxl);readline()
			#auxl[[i]]$rotref = dball$rotref#;print(id);readline() 
			#auxl[[i]]$corr = names(x)
			vcol = auxl[[i]]$hot
			idcol = round(mapply(renorm_interval,vcol,inv=F),0)
			hotcols = cols[idcol]
			auxl[[i]]$hotcols = hotcols
			##vcol = dbhoti[[id]]$score
			##cols = gray.colors(100,start=0.3,end=1)#;print(vcol)
			##cols = brewer.pal(n = 9, name = 'YlOrRd')
			##idcol = round(mapply(renorm_interval,vcol,inv=F),0)#;print(idcol);readline()
			##cols=cols[idcol]

		}
		#auxall[[i0k]] = auxl
		#i0k=i0k+1
		#}
		#print(auxall);readline()
	}
	return(auxl)
	#return(auxall)
}

verify_hot_nodes = function(db,dist,pdbname1,pdbname2,limit=c(1.4,0.85),sep="-",zero=0.1,k=3){

	
	dist.m=matrix(as.numeric(unlist(strsplit(names(dist),sep))),ncol=2,byrow=T)
	dist.m = as.data.frame(cbind(dist.m,dist))
	dist.m$a1 = 0 #;print(dist.m);readline()
	dist.m$a2 = 0
	dist.m$a = 0#;print(dist.m);readline()
	dist.id = which(dist.m[,3]<limit[1])
	tam = length(dist.id)
	#if (tam>0){
	#dist.hot = tomatrix(dist.m[v,]);print(dist.hot)
	#tam = dim(dist.hot)[1]
	#print(pdbsplit(pdbname1));readline()
	#id1 = as.numeric(pdbsplit(pdbname1)[,5:6])#;print(id1)
	#id2 = as.numeric(pdbsplit(pdbname2)[,5:6])#;print(id2)
	id1 = as.numeric(pdbsplit(pdbname1)[,c("id","nclus")])#;print(id1)
	id2 = as.numeric(pdbsplit(pdbname2)[,c("id","nclus")])#;print(id2)
	d1 = diag(db[[id1[1]]]$exp$superclus[[id1[2]]]$a)#;print(d1)
	d2 = diag(db[[id2[1]]]$exp$superclus[[id2[2]]]$a)#;print(d2);readline()
	did1 = dist.m[,1]
	did2 = dist.m[,2]

	f1 = did1!=0
	f2 = did2!=0
	fd1 = did1[f1]
	fd2 = did2[f2]

	dist.m$a1[f1] = d1[fd1]
	dist.m$a2[f2] = d2[fd2]

	#dist.m$a1 = d1[dist.m[,1]];
	#dist.m$a2 = d2[dist.m[,2]]
	#dist.m[!f1,"a1"] = d1[d.m.1];
	#dist.m[!f2,"a2"] = d2[d.m.2];
	#print(dist.m);readline()
	if (tam>0){
		for (i in dist.id){
			i1 = dist.m[i,1]
			i2 = dist.m[i,2]
			#d12 = abs(d1[i1]-d2[i2]);print(d12)
			d12 = round(min(d1[i1],d2[i2])/max(d1[i1],d2[i2]),k)
			if (d12>limit[2]){
				#dist.m[i,"a1"]=d1[i1]
				#dist.m[i,"a2"]=d2[i2]
				dist.m[i,"a"]=d12
			}else{
				print(paste("WARNING? there is doubt if these parameters (",i,dist.m[i,3],d1[i1],d2[i2],d12,") indicate a hot node in",pdbname1))
			}
		}
		#print(dist.m);readline()
		#return(dist.m)
	}else{
		print(paste("WARNING: it was not found a hot nodes for",pdbname1))
		#print(dist.m);readline()
		#return(dist.m)
	}
	return(dist.m)
}


make_hot_index = function(hot,r=3){


	tam = length(hot)
	newhot = hot

	for (i in 1:tam){
		newhot[[i]]$score = 0
	}
	for (i in 1:tam){
		#print(hot[[i]]);print(newhot[[i]])#;readline()
		repeat{
			id = which(hot[[i]]$a>0)#;print(id);readline()
			s = c()
			pilha = list()
			n = 1
			if (length(id)){
				id = id[1]#;print(id)
				#k1 = hot[[k]][id,"V1"]
				V2 = hot[[i]][id,"V2"]
				a = hot[[i]][id,"a"]
				s = c(s,a)
				pilha[[n]] = c(i,id)#;print(pilha)
				n = n + 1 
				hot[[i]][id,"a"] = 0			
				#;print(hot[[i]]);readline()
				k = i + 1
				repeat{
					#if (is.null(hot[[k]])) break
					if (k>tam) break
					else{
						id = which(hot[[k]]$V1 == V2)#;print("o");print(id)
						V2 = hot[[k]][id,"V2"]
						a = hot[[k]][id,"a"]
						if (a>0) {
							hot[[k]][id,"a"] = 0
							s = c(s,a)
							pilha[[n]] = c(k,id)#;print(pilha)
							n = n + 1
						} else {
							id0 = hot[[k-1]]$V2 == hot[[k]][id,"V1"]
							a = newhot[[k-1]][id0,"a"]
							if (a>0){
								s = c(s,a)
								pilha[[n]] = c(k,id)
								n = n + 1
							} 
							break
						}
						#print(hot[[k]]);readline()
						k = k + 1
					}
				}
			}else{
				break;
			}			
			#print(s);print(sum(s)/tam);print(pilha);readline()
			score = round(sum(s)/tam,r)
			while(length(pilha)){
				p.i = pilha[[1]][1]#;print(p.i)
				p.id = pilha[[1]][2]#;print(p.id)
				newhot[[p.i]][p.id,"score"] = score
				pilha[[1]] = NULL
				#print(pilha);readline()
			}
		
		}
		#print(hot);print(newhot);readline()
	}
	#print(hot);print(newhot);readline()
	return(newhot)

}

### ATENCAO: para desativar hotspot, go=F
set_hotspot_candidates = function(db,dbres,go=T){

	tam = length(dbres$alignment)#;print(tam)
	#print(dbres$pdbname)#;readline()

	hot = list()
	#print(tam);readline()
	if (tam>1){
		for (i in 1:tam){#print(i)
		#for (i in 2){
			if (go){
				print(dbres$alignment[[i]]$prescore);
				id1 = dbres$alignment[[i]]$prescore$id1[2]#;print(id1)
				id2 = dbres$alignment[[i]]$prescore$id2[2]#;print(id2)
				pdbname1 = as.character(dbres$alignment[[i]]$prescore$n1[2])#;print(pdbname1);readline()		
				pdbname2 = as.character(dbres$alignment[[i]]$prescore$n2[2])
				dist = dbres$alignment[[id1]]$aligndata[[id2]]$score$global$dist#;print(dist)#;readline()
				hot[[i]] = verify_hot_nodes(db,dist,pdbname1,pdbname2)#;print(hot[[i]]);readline()
			}else{
				hot[[i]] = data.frame(V1=0,V2=0,dist=0,a1=0,a2=0,a=0,score=0)
			}
			#hot[[i]] = dist.m
			#print(i);print(hot);readline()
		}
		if (go){
			hot=make_hot_index(hot)
		}
	}else{
		#print("EEEEEEEEEEEEEEEEEEEEEEEEENTREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEI")
		hot[[1]]=data.frame(V1=0,V2=0,dist=0,a1=0,a2=0,a=0,score=0)
		#hot[[1]][[1]]=list(score=0)
	}
	#print(hot);readline()
	return(hot)

}

set_transform_reference_rot = function(dball,id=1){

	tam = dim(dball$res[[id]]$scorend)[1]

	#auxv = which(dball$group$id==id)
	#auxn = as.character(dball$group$fakename[auxv])

	#print(auxn);readline()

	auxr = list()

	for (i in 1:tam){
		auxr[[i]]=list()
		#auxend1 = as.character(dball$res[[id]]$scorend[i,"pdb1"])
		#auxend2 = as.character(dball$res[[id]]$scorend[i,"pdb2"])
		auxend1 = as.character(dball$res[[id]]$scorend[i,"n1"])
		auxend2 = as.character(dball$res[[id]]$scorend[i,"n2"])
		#print(auxend1);print(auxend2)
		#clusid = which(auxend1==auxn)
		#pdbid = which(auxend2==dball$group$fakename)
		clusid = as.numeric(as.vector(dball$res[[id]]$scorend[i,"id1"]))
		pdbid = as.numeric(as.vector(dball$res[[id]]$scorend[i,"id2"]))
		pdbname = dball$group$pdb[pdbid]
		auxr[[i]]$pdb = as.character(pdbname)
		auxr[[i]]$pdbpair = c(auxend1,auxend2)
		#print(clusid);print(pdbid);readline()
		auxr[[i]]$rot = dball$res[[id]]$alignment[[clusid]]$aligndata[[pdbid]]$rot
		#print(auxr);readline()
	}
	dball$rotref = auxr

	return(dball)
}


make_group_table = function(dball,nogroup="none",filename=NULL,outdir="Datalign/",sufix=".csv"){

#dball[[i]]$groups = pdbsplit(dball[[i]]$info[,6])
#          gacc lacc  pdb chain1 chain2 method id nclus            fakename
#1_1ppf_4  0.90 0.85 1PPF      E      I  BSRAA  1     4  1PPF_E_I_BSRAA_1_4
#1_1ppf_5  0.85 0.70 1PPF      E      I  BSRAA  1     5  1PPF_E_I_BSRAA_1_5
#1_1ppf_6  0.82 0.56 1PPF      E      I  BSRAA  1     6  1PPF_E_I_BSRAA_1_6

#          gacc lacc  pdb chain1 chain2 method id nclus            fakename
#1_1ppf_4  0.90 0.85 1PPF      E      I  BSRAA  1     4  1PPF_E_I_BSRAA_1_4
#1_1ppf_5  0.85 0.70 1PPF      E      I  BSRAA  1     5  1PPF_E_I_BSRAA_1_5
#1_1ppf_6  0.82 0.56 1PPF      E      I  BSRAA  1     6  1PPF_E_I_BSRAA_1_6

	auxv = dball$info$fakename!=nogroup
	auxr = dball$info[auxv,]
	auxp = data.frame(pdbsplit(auxr$fakename))
	colnames(auxp)=c("pdb","chain1","chain2","method","polarity","id","nclus")
	#print(auxp);readline()
	auxr = data.frame(auxr[,c("gacc","lacc")],auxp,fakename=auxr$fakename)
	#print(auxr);readline()
	dball$groups = auxr
	if (!is.null(filename)){
		filename = paste(outdir,toupper(filename),sufix,sep="")
		write.csv(auxr,filename)
	}
	return(dball)

}

sim_corbi = function(a,b,k=3){

	aux = 1 - abs(a-b)/(max(a,b))

	return(round(aux,k))

}

build_similarity_to_corbi = function(d1,d2,s1="v",s2="u"){

	print(d1);print(d2)

	auxr = c()

	for (i in 1:length(d1)){
		n1 = paste(s1,i,sep="")
		for (j in 1:length(d2)){
			n2 = paste(s2,j,sep="")
			auxs = sim_corbi(d1[i],d2[j])
			auxl = data.frame(n1,n2,auxs)
			auxr = rbind(auxr,auxl)
			#print(auxr);readline()
		}
	}
	#print(auxr);readline()
	return(auxr)
}


#BIHARCK: Retorna o índice de oclusão OD como a razao da distancia d
# da distancia da esfera cn interveniente ao eixo 
get_occlusion_degree = function(c1, c2, cn, rc2, rcn, is_verbose=F,dd=3){
  od = 0
  
  v = cn-c1
  w = c2-c1
  
  u = round(w*(dot(v,w)/(Norm(w)^2)),dd)
  
  #D=||u-w||
  #d = round(Norm(u-w),dd)
  d = round(Norm(v-u),dd)
  duw = round(Norm(w)-Norm(u),dd)
  drc2n = round((rc2+rcn),dd)
  
  cos_t = round(dot(w,v)/(Norm(w)*Norm(v)),dd)
  
  if(is_verbose){
    print(paste0('v = [', v[1], ', ', v[2], ', ',v[3], ']'))
    print(paste0('w = [', w[1], ', ', w[2], ', ',w[3], ']'))
    print(paste0('u = [', u[1], ', ', u[2], ', ',u[3], ']'))
    print(paste0('d = ', round(d,dd), ' | cos_t = ', cos_t))
    print(paste0('rc2+rcn= ',drc2n,' ||u||<||w|| = ',duw))	
  }
  
  od = d/drc2n

  # primeiro teste
  if((cos_t > 0) && (cos_t <1) && (duw>0)){
    if (is_verbose) print('válido')
    #od = d/drc2nd
  }else{
    if (is_verbose) print('Inválido')
    #od=0
  }
  
  # segundo teste
  #if(d <(rc2 + rcn)){
  #  od = (rc2+rcn)/d
    #print('Existe oclusão')
  #}
  #od = (rc2+rcn)/d
  return(od)
}


#SABRINA: monta colunas conforme ids dos grupos e subgrupos
make_columns = function(db){

	tam = unique(db$group)
	auxr = c()

	for (i in tam){

		auxv = db$group == i
		auxc = unique(db[auxv,"subgroup"])
		#print(as.character(i))
		auxr=c(auxr,strcat(as.character(i),paste(auxc)))
		#print(auxr);readline()

	}
	return(auxr)

}

#SABRINA: 
ids_by_groups = function(db){

	auxcols = make_columns(db)
	auxrows = unique(db$id)

	auxm = matrix(0,nrow=length(auxrows),ncol=length(auxcols),byrow=T)
	rownames(auxm) = auxrows
	colnames(auxm) = auxcols
	#print(auxm);readline()
	tam = dim(db)[1]

	for (i in 1:tam){
		x = as.character(db[i,"id"])
		y = strcat(as.character(db[i,"group"]),as.character(db[i,"subgroup"]))
		auxm[x,y]=1
		#print(x);print(y);readline()		
	}

	#print(auxm);readline()
	return(auxm)
}

#SABRINA: imprime vetor/matriz de ocorrencia de atributos
print_attributes_from_fingers = function(db,expid,g,subg){

	auxv=rep(0,length(db$colnames))

	auxatt=paste("att:e",expid,"g",g,sep="")
	auxs = data.frame(att=db$colnames)
	#print(auxs);readline()
	auxn = c()


	for (i in subg){

		auxt = paste("s",i,sep="")
		auxf = db$exp[[expid]]$groups[[g]]$subs[[i]]$finger
		auxs = cbind(auxs,auxf)
		auxn = c(auxn,auxt)
		#print(auxs);readline()
		auxv=auxv+auxf
	}
	#print(auxv)
	auxs = cbind(auxs,auxv)
	colnames(auxs)=c(auxatt,auxn,"total")
	print(auxs)#;readline()

	#return(auxs)
	#tam = dim(auxs)[2]
	filter = auxs[,"total"]!=0
	print(auxs[filter,c(auxatt,"total")])#;readline()



	#auxv = 

}

area_under_curves = function(x,m,abs=T){

	if (is.vector(m)) m = as.matrix(m)
	tam = 1:dim(m)[2]

	auxa = c()

	for (i in tam){
		#print(i)
		#auxa = c(auxa,auc(x,m[,i],type="spline"))
		y = m[,i]
		if (abs) {
			y = abs(y)
		}
		auxa = c(auxa,auc(x,y,type="linear"))
		#print(i)
	}
	return(auxa)

}


#ERNESTO...
#	a = db[[1]]$exp1[[10]]$mr
#	b = db[[1]]$exp1[[9]]$mr
#	x= a - b
stat_res_by_decomp = function(svd,k){

	auxmc = svd$c
	aux = list()
	aux[[1]] = list()
	aux[[1]]$mr = auxmc
	aux[[1]]$dr = 0
	aux[[1]]$cp = 0
	aux[[1]]$acp = 0
	if (k<2) print(paste("WARNING: k must be > 2"))
	else {
		for (j in 2:k){
			aux[[j]]=list()
			auxmr = rebuild_matrix(svd,1:j)
			aux[[j]]$mr = auxmr
			nc = 2:dim(auxmr)[2]
			auxdr = auxmr[,nc]-auxmc[,nc]
			auxdr = cbind(auxmr[,1],auxdr)
			aux[[j]]$dr = auxdr
			auxdiff = aux[[j]]$mr - aux[[j-1]]$mr
			aux[[j]]$cp = auxdiff
			aux[[j]]$acp = area_under_curves(auxmc[,1],aux[[j]]$cp)
		}
	}	
	return(aux)

}

plot_ernest_residues_old = function(m,mc,type,ylim){

	i = 2:dim(mr)[2]
	mi = m[,i]-mc[,i]
	mr = cbind(m[,1],mi)
	plot_ernest(mr,type=type,ylim=ylim)

}

plot_ernest = function(m,tam=NULL,add=F,lty=1,type="l",ylim=c(-1,1),cex=0.001,pch=".",main="",xlab="",ylab=""){
	if (is.null(tam)) {
		tam = 2:dim(m)[2]
		cols = blue2red_v1(length(tam))
	}else{
		auxtam = 2:dim(m)[2]
		cols = blue2red_v1(length(auxtam))
		#cols = cols[tam]
	}
	j = 1
	for (i in tam){
		x = m[,1]
		y = m[,i]
		#print(cols[i])
		if (!add){
			if (j==1){
				plot(x,y,type=type,col=cols[i-1],ylim=ylim,cex=cex,pch=pch,main=main,xlab=xlab,ylab=ylab,lty=lty)
				j=j+1
			}else{
				lines(x,y,type=type,col=cols[i-1],cex=cex,pch=pch,lty=lty)
			}
		}else{
			lines(x,y,type=type,col=cols[i-1],cex=cex,pch=pch,lty=lty)			
		}
		#readline()	
	}

}

plot_ernest_component_analysis = function(db,comp,tam,inter,mfrow,ylim=c(-0.1,1),lty=3){

	yl1 = ylim[1]
	yl2 = ylim[2]

	if (!is.null(mfrow)) par(mfrow=mfrow)

	b=1	
	for (a in inter){
		x = rebuild_component(db[[a]]$A$svd$c,db[[a]]$A$svd$r,db[[a]]$A$a[,1],all=1)
		y = rebuild_component(db[[a]]$exp1[[comp-1]]$mr,db[[a]]$A$svd$r,db[[a]]$A$a[,1],all=1)
		z = rebuild_component(db[[a]]$exp1[[comp]]$cp,db[[a]]$A$svd$r,db[[a]]$A$a[,1],all=0)
		#print(dim(x));print(dim(y));print(dim(z))
		plot_ernest(x,tam=tam[b],ylim=c(yl1,yl2))
		plot_ernest(y,tam=tam[b],ylim=c(yl1,yl2),lty=lty,add=T)
		plot_ernest(z,tam=tam[b],ylim=c(yl1,yl2),lty=lty,add=T)
		b=b+1
		#readline()
	}
}

rebuild_matrix = function(m,inter){

	auxm = m$us[,inter] %*% t(m$v[,inter])
	return(auxm)

}

as_bipartite = function(m,matrixtype="normal"){


	tami = dim(m)[1]
	tamj = dim(m)[2]

	if (matrixtype=="normal"){
		m11 = matrix(rep(0,tami*tami),nrow=tami)
		m22 = matrix(rep(0,tamj*tamj),nrow=tamj)
	}else{
		#m11 = Matrix(rep(0,tami*tami),nrow=tami)
		#m22 = Matrix(rep(0,tamj*tamj),nrow=tamj)
		m11 = Diagonal(tami,x=0)
		m22 = Diagonal(tamj,x=0)
	}
	m12 = m
	m21 = t(m)
	#print(dim(m11));print(dim(m22));print(dim(m12));print(dim(m21));readline()
	m1 = cbind(m11,m12)
	m2 = cbind(m21,m22)

	mr = rbind(m1,m2)

	#print(dim(mr));readline()
	return(mr)
	
	#return(db)

}


# GENERAL: (dbv[[1]],expid=1,g=5,krange=1:5,mfrow=c(2,3))
view_3D_decomp_pam_general = function(db,expid,g,axes,xlim,abcline=T,pal=1,size=1,us=F){

	#g=5
	#axes=c(1,1,1)
	#axes=4:6
	#axes=7:9
	#axes=1:3
	#xlim=c(-0.5,0.5)
	ylim=xlim
	zlim=ylim
	m.svd = db$A$svd
	m.pam = db$exp[[expid]]$pam[[g]]
	clus = m.pam$clustering
	if (pal==1){
		bcol=clus
		clus=1:length(clus)
	}
	if (pal==2){
		bcol = rainbow(g)
	}
	if (us){
		#print(bcol[clus])
		#print(m.svd$us[,axes[1]]);print(m.svd$u[,axes[1]]);readline()
		if (is.null(xlim)){
			plot3d(m.svd$us[,axes[1]],m.svd$us[,axes[2]],m.svd$us[,axes[3]],size=size,aspect="iso",box=FALSE,type="s",col=bcol[clus],xlab=axes[1],ylab=axes[2],zlab=axes[3])
		} else {
			plot3d(m.svd$us[,axes[1]],m.svd$us[,axes[2]],m.svd$us[,axes[3]],size=size,aspect="iso",box=FALSE,type="s",col=bcol[clus],xlim=xlim,ylim=ylim,zlim=zlim,xlab=axes[1],ylab=axes[2],zlab=axes[3])
		}
	}else{
		plot3d(m.svd$u[,axes[1]],m.svd$u[,axes[2]],m.svd$u[,axes[3]],size=size,aspect="iso",box=FALSE,type="s",col=bcol[clus],xlim=xlim,ylim=ylim,zlim=zlim,xlab=axes[1],ylab=axes[2],zlab=axes[3])
	}
	#a=t(m.svd$v[1:3,1:3])
	#if(abcline) abclines3d(0, 0, 0, a = a, col = "blue")
	if(abcline) abclines3d(0, 0, 0, a = diag(3), col = "gray")

}

red2blue_v1 = function (Ncolors){
	
	col = c()

	if (Ncolors==1){
		return("#FF0000")
	}else{

		Ncolors = Ncolors - 1

		for (i in Ncolors:0)
			col = c(col,rgb(i,0,(Ncolors-i),maxColorValue=Ncolors))

		return (col)
	}
}

blue2red_v1 = function(Ncolors){

	auxv = red2blue_v1(Ncolors)
	return(rev(auxv))

}

red2blue_v2 = function (Ncolors){
	
	col = c()

	#Ncolors = Ncolors - 1
	half = round(Ncolors/2,0)

	for (i in half:0)
		#col = c(col,rgb(Ncolors-i,Ncolors-i,Ncolors,maxColorValue=Ncolors))
		#col = c(col,rgb(half,half-i,half-i,maxColorValue=half))
		col = c(col,rgb(half-i,half-i,half,maxColorValue=half))
	half = half
	for (i in 1:half)
		#col = c(col,rgb(Ncolors-i,Ncolors-i,Ncolors,maxColorValue=Ncolors))
		#col = c(col,rgb(half-i,half-i,half,maxColorValue=half))
		col = c(col,rgb(half,half-i,half-i,maxColorValue=half))

	return (col)

}



select_sub = function (m,mlist){

	len = length(mlist)
	mret = c()

	for (i in 1:len){
		aux = m[which(m[,2]==mlist[i]),]
		if (dim(aux)[1]!=0){
			mret = rbind(mret,aux)
		} 
	}
	return(mret)

}

filter_asa = function(m,error){

	aux = m[which(m[,6]>error),]
	return(aux)

}


#pre_process_name_file =function(dirnamein,dirnameout,pdbids,prestr,posstr){

#	tam = length(pdbids)

#}

###################################################################################################
# converter uma lista numa string
preprocess01 = function(n){

	tamn = length(n)
	x = c()

	for (i in 1:tamn){
		x = c(x,paste(n[i][[1]]))
	}

	return(x)

}

# Transforma codigo de aminoacido de 3 letras para 1 letra
res = function(a){

	aux.r = c("GLY","ALA","VAL","LEU","ILE","CYS","MET","PRO","PHE","TYR","TRP","LYS","ARG","HIS","GLU","ASP","SER","THR","ASN","GLN")
	aux.l = c("G","A","V","L","I","C","M","P","F","Y","W","K","R","H","E","D","S","T","N","Q")


	return(aux.l[which(aux.r==a)])



}

# Adjusta os nomes chain-res-atom
preprocess02 = function(n){

	tamv = length(n)
	#aux = "e"
	#enz = "e"
	#i = 1
	aux = c()
	for (i in 1:tamv) {
		aux.i = unlist(strsplit(n[i],"[.]"))
		aux.c = aux.i[1]
		aux.r3 = substr(aux.i[2],1,3)
		aux.r =res(aux.r3)
		aux.n =gsub(aux.r3,"",aux.i[2])
		aux.a = aux.i[3]
		#aux = c(aux,paste(aux.c,aux.r,aux.n,aux.a,sep=""))
		aux = c(aux,paste(aux.c,aux.r,aux.n,".",aux.a,sep=""))
	}

	#print(aux)
	return(aux)

}

preprocess03 = function(m){


  tamr = dim(m)[1]
  tamc = dim(m)[2]
  
  for (i in 1:tamr){
    for (j in 1:tamc){
      if (is.na(m[i,j])) m[i,j] = 0
    }
  
  }
  return(m)
}

preprocess04 = function(m){


  tamr = dim(m)[1]
  tamc = dim(m)[2]
  
  for (i in 1:tamr){
    for (j in 1:tamc){
      if (is.infinite(m[i,j])) m[i,j] = 0
    }
  
  }
  return(m)
}

#Faz split de "eR146.C" num vetor v=c("e" "R" "146" "C")
splitatom = function(str,sep="[.]",typename=1){

	if (typename==1){
		aux1 = strsplit(str,sep)	
		aux1n = nchar(aux1[[1]][1])
		auxresi = substr(aux1[[1]][1],3,aux1n)
		auxchain = substr(aux1[[1]][1],1,1)
		auxresn = substr(aux1[[1]][1],2,2)
		auxatom = aux1[[1]][2]
		aux = c(auxchain,auxresn,auxresi,auxatom)
	}else
	if (typename==2){
		aux1 = unlist(strsplit(str,sep))
		auxchain = aux1[1]
		auxlen = nchar(aux1[2])
		auxresn = substr(aux1[2],1,3)
		auxresi = substr(aux1[2],4,auxlen)
		auxatom = aux1[3]
		#print(aux1)#;readline()
		aux = c(auxchain,aa321(auxresn),auxresi,auxatom)
		#print(aux);readline()
	}else
	if(typename==3){
		aux1 = unlist(strsplit(str,sep))
		auxchain = aux1[1]
		auxlen = nchar(aux1[2])
		auxresn = substr(aux1[2],1,3)
		auxresi = substr(aux1[2],4,auxlen)
		auxatom = unlist(strsplit(aux1[3],"[-]"))[1]
		aux = c(auxchain,aa321(auxresn),auxresi,auxatom)
	}else
	if(typename==4){
		aux1 = unlist(strsplit(str,sep))
		auxchain = aux1[1]
		auxlen = nchar(aux1[2])
		auxresn = substr(aux1[2],1,3)
		auxresi = substr(aux1[2],4,auxlen)
		#auxatom = unlist(strsplit(aux1[3],"[-]"))[1]#;print(unlist(strsplit(aux1[3],"[-]")));readline()
		auxatom = unlist(strsplit(aux1[3],"[-]"))
		auxatom = c(auxatom[1],auxatom[3])
		aux = c(auxchain,aa321(auxresn),auxresi,auxatom)
	}else
	if(typename==5){
		aux1 = unlist(strsplit(str,sep))
		auxchain = aux1[1]
		auxlen = nchar(aux1[2])
		auxresn = substr(aux1[2],1,3)
		auxresi = substr(aux1[2],4,auxlen)
		#auxatom = unlist(strsplit(aux1[3],"[-]"))[1]#;print(unlist(strsplit(aux1[3],"[-]")));readline()
		auxatom = unlist(strsplit(aux1[3],"[-]"))#;print(aux1[3]);print(auxatom);readline()
		auxatom = c(auxatom[1],auxatom[3],auxatom[4])
		aux = c(auxchain,aa321(auxresn),auxresi,auxatom)

	}
	return(aux)

}

# Retorna vetor com quantidade de arestas
edges = function(m){
	
	tamr = dim(m)[1]
	aux = c()

	for (i in 1:tamr){
		aux = c(aux,sum(m[i,]>0))		
	}

	return(aux)

}

#Colore conforme e=enzima e i=inibidor

colv = function(n,cole="RED",coli="BLUE"){

	tam = length(n)
	aux = c()
	
	for (i in 1:tam){
		if (substr(n[i],1,1)=="e") aux = c(aux,cole)
		if (substr(n[i],1,1)=="i") aux = c(aux,coli)
	}
	return(aux)

}

colvap = function(n,cola="BLACK",colp=NA,colc="GRAY100"){

	tam = length(n)
	aux = c()
	
	for (i in 1:tam){
		if (n[i]=="a") aux = c(aux,cola)
		if (n[i]=="p") aux = c(aux,colp)
		if (n[i]=="c") aux = c(aux,colc)
	}
	return(aux)

}


#make_atom_polarity_code = function(){

#	pol = c()

#	aux = list()
#	aux$apolar = list()
#	aux$apolar$code = "a"
#	aux$apolar$list = c("C.CB","C.SG","V.CB","V.CG1","V.CG2","P.CB","P.CG","P.CD","A.CB","I.CB","I.CG1","I.CG2","I.CD1","Q.CB","Q.CG")

#	aux$polar = list()
#	aux$polar$code = "p"
#	aux$polar$list = c("OE1","NE2")

#	aux$meso = list()
#	aux$meso$code = "m"
#	aux$meso$list = c("Q.CD","S.CB","N.CG","E.CD","D.CG","K.CE","T.CB","R.CZ")

#	aux$backbone = list()
#	aux$backbone$code = "b"
#	aux$backbone$list = c("N","C","O")

#	aux$calpha = list()
#	aux$calpha$code = "c"
#	aux$calpha$list = c("CA")

#}


# Verifica se um carbono eh meso (anfipático)
is.meso = function(res,atom,meso){

	auxmeso = strsplit(meso,"[.]")
	tammeso = length(auxmeso)
	#print(auxmeso);readline()
	for (j in 1:tammeso){
		auxres = auxmeso[[j]][1]
		auxatom = auxmeso[[j]][2]
		if ((auxres == res)&(auxatom == atom)) return(TRUE)
	}
	return(FALSE)
}


new_codify_polarity7 = function(n,meso,aacode,vdw,bin=F){

	#"E.ILE16.N-1-N-A",...,"E.ILE16.CB-5-C-A"
	tam = length(n)
	aux = c()
	aa_canonical = c("G","P","A","V","L","I","H","F","W","Y","C","M","S","T","N","Q","D","E","K","R")

	for (i in 1:tam){
		atom.v = splitatom(n[i],type=5)#;print(head(n));print(head(atom.v));readline()
		res = atom.v[2]
		atom = atom.v[4]
		type = atom.v[5]
		tag = atom.v[6] #A - ATOM, H - HETATM
		if (is.meso(res,atom,meso)) aux = c(aux,aacode[2])
		else if ((atom == "CA")&(type == "C")) aux = c(aux,aacode[5])
		else if (atom == "C") aux = c(aux,aacode[6])
		else if (atom %in% c("N","O")) aux = c(aux,aacode[4])
		else if (atom %in% c("SG","SD")) aux = c(aux,aacode[1])
		else if (type == "C") aux = c(aux,aacode[1])
		else if ((type == "S")&(tag == "A")) aux = c(aux,aacode[1])
		else if ((type == "S")&(tag == "H")) aux = c(aux,aacode[3])
		else if (type == "N") aux = c(aux,aacode[3])
		else if (type == "O") aux = c(aux,aacode[3])
		else if (type == "P") aux = c(aux,aacode[3])
		else if (type == "H") aux = c(aux,aacode[8])
		else aux = c(aux,aacode[7])
		if ( (res %in% aa_canonical) & (aux[i]=="o")){
			print(paste("WARNING: Atom",n[i],"is canonical and it was considered as another"))
		}
	}
	if (length(aux)!=tam){
		print(paste("WARNING: something wrong in codify polarity 6"))
	}
	return(aux)
}


new_codify_polarity6 = function(n,meso,aacode,vdw,bin=F){

	tam = length(n)
	aux = c()
	aa_canonical = c("G","P","A","V","L","I","H","F","W","Y","C","M","S","T","N","Q","D","E","K","R")

	for (i in 1:tam){
		atom.v = splitatom(n[i],type=5)#;print(atom.v);readline()
		res = atom.v[2]
		atom = atom.v[4]
		type = atom.v[5]
		tag = atom.v[6]
		if (is.meso(res,atom,meso)) aux = c(aux,aacode[2])
		else if ((atom == "CA")&(type == "C")) aux = c(aux,aacode[5])
		else if (atom == "C") aux = c(aux,aacode[6])
		else if (atom %in% c("N","O")) aux = c(aux,aacode[4])
		else if (atom %in% c("SG","SD")) aux = c(aux,aacode[1])
		else if (type == "C") aux = c(aux,aacode[1])
		else if ((type == "S")&(tag == "A")) aux = c(aux,aacode[1])
		else if ((type == "S")&(tag == "H")) aux = c(aux,aacode[3])
		else if (type == "N") aux = c(aux,aacode[3])
		else if (type == "O") aux = c(aux,aacode[3])
		else if (type == "P") aux = c(aux,aacode[3])
		else aux = c(aux,aacode[7])
		if ( (res %in% aa_canonical) & (aux[i]=="o")){
			print(paste("WARNING: Atom",n[i],"is canonical and it was considered as another"))
		}
	}
	if (length(aux)!=tam){
		print(paste("WARNING: something wrong in codify polarity 6"))
	}
	return(aux)
}

new_codify_polarity5 = function(n,meso,aacode,vdw,bin=F){

	tam = length(n)
	aux = c()
	aa_canonical = c("G","P","A","V","L","I","H","F","W","Y","C","M","S","T","N","Q","D","E","K","R")

	for (i in 1:tam){
		atom.v = splitatom(n[i],type=4)#;print(atom.v);readline()
		res = atom.v[2]
		atom = atom.v[4]
		type = atom.v[5]
		if (is.meso(res,atom,meso)) aux = c(aux,aacode[2])
		else if ((atom == "CA")&(type == "C")) aux = c(aux,aacode[5])
		else if (atom == "C") aux = c(aux,aacode[6])
		else if (atom %in% c("N","O")) aux = c(aux,aacode[4])
		else if (atom %in% c("SG","SD")) aux = c(aux,aacode[1])
		else if (type == "C") aux = c(aux,aacode[1])
		else if (type == "S") aux = c(aux,aacode[1])
		else if (type == "N") aux = c(aux,aacode[3])
		else if (type == "O") aux = c(aux,aacode[3])
		else if (type == "P") aux = c(aux,aacode[3])
		else aux = c(aux,aacode[7])
		if ( (res %in% aa_canonical) & (aux[i]=="o")){
			print(paste("WARNING: Atom",n[i],"is canonical and it was considered as another"))
		}
	}
	if (length(aux)!=tam){
		print(paste("WARNING: something wrong in codify polarity 5"))
	}
	return(aux)
}

#code=c("a","m","p","b","c","i","o")
new_codify_polarity4 = function(n,meso,aacode,vdw,bin=F){

	tam = length(n)
	aux = c()
	aa_canonical = c("G","P","A","V","L","I","H","F","W","Y","C","M","S","T","N","Q","D","E","K","R")
	carbon = c()

	for (i in 1:tam){

		atom.v = splitatom(n[i],type=3)#;print(atom.v)#;readline()
		res = atom.v[2]
		atom = atom.v[4]

		if (is.meso(res,atom,meso)) aux = c(aux,aacode[2])
		else if (atom %in% c("CA")) aux = c(aux,aacode[5])
	  else if (atom %in% c("CB","CD","CE","CG","CZ","CH","CD1","CD2","CE1","CE2","CE3","CZ1","CZ2","CZ3","CH1","CH2","CH3")) aux = c(aux,aacode[1])
    else if (atom %in% c("SG","SD")) aux = c(aux,aacode[1])
		else if (atom %in% c("N","O")) aux = c(aux,aacode[4])
		else aux = c(aux,aacode[3])

	}
	if (length(aux)!=tam){
		print(paste("WARNING: something wrong in codify polarity 4"))
	}

	return(aux)
}

#code=c("a","m","p","b","c","i","o")
new_codify_polarity3 = function(n,meso,aacode,bin=F){

	tam = length(n)
	aux = c()
	aa_canonical = c("G","P","A","V","L","I","H","F","W","Y","C","M","S","T","N","Q","D","E","K","R")

	for (i in 1:tam){

		atom.v = splitatom(n[i],type=3)#;print(atom.v)#;readline()
		res = atom.v[2]
		atom = atom.v[4]
		atom.sep = unlist(strsplit(atom,""))
		atom.sep.tam = length(atom.sep)
		if (grepl("^C",atom)){
			if (atom.sep.tam==1) aux = c(aux,aacode[6]) #C carbonila - backbone
			else if (atom.sep.tam==2){
				if(atom.sep[2]=="A") aux = c(aux,aacode[5]) #CA
				else if (is.meso(res,atom,meso)) aux = c(aux,aacode[2]) # C anfipatico/meso
				else if (atom.sep[2] %in% c("B","D","E","G","Z","H")) aux = c(aux,aacode[1])
				else aux = c(aux,aacode[7])
			} else if (atom.sep.tam==3){
				if (atom.sep[2] %in% c("B","D","E","G","Z","H")){
					if (atom.sep[3] %in% 1:3){
						aux = c(aux,aacode[1])
					} else aux = c(aux,aacode[1])
				} else aux = c(aux,aacode[7])		
			} else aux = c(aux,aacode[7])
 
		} else if (grepl("^S",atom)){
				if (atom.sep.tam==2){
					if (atom.sep[2] %in% c("D","G")) aux = c(aux,aacode[1])
					else aux = c(aux,aacode[7])
				} else	aux = c(aux,aacode[7])
		} else if (grepl("^N",atom)){
			if (atom.sep.tam==1) aux = c(aux,aacode[4])
			else aux = c(aux,aacode[3])
		} else if (grepl("^O",atom)){
			if (atom.sep.tam==1) aux = c(aux,aacode[4])
			else aux = c(aux,aacode[3])
		} else aux = c(aux,aacode[7])

		if ( (res %in% aa_canonical) & (aux[i]=="o")){
			print(paste("WARNING: Atom",n[i],"is canonical and it was considered as other"))
		}

		#print(aux);readline()
	}
	if (length(aux)!=tam){
		print(paste("WARNING: something wrong in codify polarity 3"))
	}

	return(aux)
}


# Experimental: discrimina aromáticos
new_codify_polarity2 = function(n,meso,aacode,aro){

#meso=c("Q.CD","S.CB","N.CG","E.CD","D.CG","K.CE","T.CB","R.CZ")
#code=c("a","m","p","b","c","i","o")

	tam = length(n)	
	#aux.c = grepl(c,n)
	aux = c()


	for (i in 1:tam){
		auxres = splitatom(n[i])[2]
		auxatom = splitatom(n[i])[4]
		aux.sep = unlist(strsplit(auxatom,""))
		#if (is.aromatic(auxres,auxatom,aro)) aux=c(aux,aacode[7])
		#else
		if (grepl("C",auxatom)){
			#aux.sep = unlist(strsplit(auxatom,""))
			if (aux.sep[length(aux.sep)]=="A") aux = c(aux,aacode[5]) # CA
			else if (aux.sep[length(aux.sep)]=="C") aux = c(aux,aacode[6]) #C carbonila - backbone
			else if (is.meso(auxres,auxatom,meso)) aux = c(aux,aacode[2]) # C anfipatico 
			else aux = c(aux,aacode[1]) # C apolar
		} else if (grepl("S",auxatom)) aux = c(aux,aacode[1]) # S apolar
		  else if ((aux.sep[length(aux.sep)]=="N")|(aux.sep[length(aux.sep)]=="O")) aux = c(aux,aacode[4]) # N, O backbone
		  else aux = c(aux,aacode[3]) # atom polar
	}
	if (length(aux)!=tam) print("WARNING: possibly some of the atom names is without coding")
	return(aux)
}

# Nova codificação em 6 classes/tipos de aa
new_codify_polarity = function(n,meso,aacode){

#meso=c("Q.CD","S.CB","N.CG","E.CD","D.CG","K.CE","T.CB","R.CZ")
#code=c("a","m","p","b","c","i")
	
	tam = length(n)	
	#aux.c = grepl(c,n)
	aux = c()


	for (i in 1:tam){
		auxres = splitatom(n[i])[2]
		auxatom = splitatom(n[i])[4]
		aux.sep = unlist(strsplit(auxatom,""))
		if (grepl("C",auxatom)){
			#aux.sep = unlist(strsplit(auxatom,""))
			if (aux.sep[length(aux.sep)]=="A") aux = c(aux,aacode[5]) # CA
			else if (aux.sep[length(aux.sep)]=="C") aux = c(aux,aacode[6]) #C carbonila - backbone
			else if (is.meso(auxres,auxatom,meso)) aux = c(aux,aacode[2]) # C anfipatico 
			else aux = c(aux,aacode[1]) # C apolar
		} else if (grepl("S",auxatom)) aux = c(aux,aacode[1]) # S apolar
		  else if ((aux.sep[length(aux.sep)]=="N")|(aux.sep[length(aux.sep)]=="O")) aux = c(aux,aacode[4]) # N, O backbone
		  else aux = c(aux,aacode[3]) # atom polar
	}
	if (length(aux)!=tam) print("WARNING: possibly some of the atom names is without coding")
	return(aux)
}


# Versão mais aprimorada do charv para codificação do tipo de aminoácido.
codify_polarity = function(n,c,s,b=c("a","p","c")){

	tam = length(n)	
	#aux.c = grepl(c,n)
	aux = c()
	for (i in 1:tam){
		auxatom = splitatom(n[i])[4]
		if (grepl(c,auxatom)){
			aux.sep = unlist(strsplit(auxatom,""))
			if (aux.sep[length(aux.sep)]=="A") aux = c(aux,b[3])
			else if (aux.sep[length(aux.sep)]!="C") aux = c(aux,b[1]) 
			else aux = c(aux,b[2])
		} else if (grepl(s,auxatom)) {		
			aux = c(aux,b[1])
		} else {
			aux = c(aux,b[2])
		}
	}
	return(aux)

}

# OBSOLETED: Classifica colname como apolar, polar ou CA
# Corrige falha das cisteínas.
charv = function(n,c,b=c("a","p","c")){

	tam = length(n)	
	#aux.c = grepl(c,n)
	aux = c()
	for (i in 1:tam){
		auxatom = splitatom(n[i])[4]
		if (grepl(c,auxatom)){
			aux.sep = unlist(strsplit(auxatom,""))
			if (aux.sep[length(aux.sep)]=="A") aux = c(aux,b[3])
			else if (aux.sep[length(aux.sep)]!="C") aux = c(aux,b[1]) 
			else aux = c(aux,b[2])
		} else {		
			aux = c(aux,b[2])
		}
	}
	return(aux)

}

# OBSOLETED: Classifica colname como apolar, polar ou CA 
# Erra para cisteinas...
charv_old = function(n,c,b=c("a","p","c")){

	tam = length(n)	
	aux.c = grepl(c,n)
	aux = c()
	for (i in 1:tam){
		if (aux.c[i]==1){
			aux.sep = unlist(strsplit(n[i],""))
			if (aux.sep[length(aux.sep)]=="A") aux = c(aux,b[3])
			else if (aux.sep[length(aux.sep)]!="C") aux = c(aux,b[1]) 
			else aux = c(aux,b[2])
			
			#else aux = c(aux,b[2])
		}else{
			aux = c(aux,b[2])
		}

	}
	return(aux)

}


find_max_position = function(v,lim=1,p=3){

	if (!is.na(lim)){
		aux = which(round(v,p)==lim)
		if (length(aux)) return (aux)
		else return (0)
	} else {
		auxv = round(v,p)
		aux = which(auxv==max(auxv))
		return(aux[1])
	}

}

# Encontra o id mais numeroso no vetor indicador v
find_max_occurrence = function(v){

	auxf = as.factor(v)
	auxl = levels(auxf)
	
	auxmax=0
	
	#tam = 1:length(auxl)
	
	for (i in auxl){
		auxsum = sum(v==i)
		if (auxsum>auxmax) {
			auxmax = auxsum
			auxi = i
		}
	}
	return(as.numeric(auxi))
}


# "desbipartiza" a matriz
desbipatite = function(m,n){

	tam = dim(m)[1]
	tamv = length(n)
	aux = "e"
	enz = "e"
	i = 1
	while (aux[1]==enz){
		aux = unlist(strsplit(n[i],"[.]"))
		i = i + 1
	}
	i = i-2
	#j = tam
	#print(paste(i-2))
	mr = m[1:i,(i+1):tam]
	return(mr)

}



# imprime pares de atomos dadas distancias
atompardist = function(m,n,v){

	tam = dim(m)[1]
	tamv = length(v)
	for (i in 1:tam){
		for (j in 1:tam){
			for ( z in 1:tamv){
				if (m[i,j]==v[z]) print(paste(n[i]," ",n[j]," : ",m[i,j]))
			}
		}	
	}
}


# Retorna num vetor valores não-zero de m
getnozero = function(m){

	taml = dim(m)[1]
	tamc = dim(m)[2]

	v=c()
	if (symmetric_check(m,"","",FALSE)){
		for (i in 1:taml){
			for (j in i:tamc){
				if (m[i,j]!=0) v = c(v,m[i,j])
			}
		}
	} else {
		print("WARNING: weights are not symmetric in matrix")
		for (i in 1:taml){
			for (j in 1:tamc){
				if (m[i,j]!=0) v = c(v,m[i,j])
			}
		}
	}
	return(v)

}

binarize = function(m){

	taml = dim(m)[1]
	tamc = dim(m)[2]

	aux = m
	
	for (i in 1:taml){
		for (j in 1:tamc){
			if (m[i,j]!=0) aux[i,j] = 1
		}
	}
	return(aux)
}

# gera uma matriz completa binária de n bits
ntobin = function(n){

	tam = 2^n-1
	auxm = c()
	tamaux = length(as.integer(intToBits(0)))
	for (i in 0:tam){
	
		aux = rev(as.integer(intToBits(i)))
		aux = aux[(tamaux-n+1):tamaux]
		auxm = rbind(auxm,aux)
	}
	rownames(auxm)=0:tam
	return (auxm)

}

# Gera a decomposição de um hipercubo dado pela matriz m
make_svd_cube = function(m){

	aux = svd(m)
	aux$us = aux$u %*% diag(aux$d)
	
	return(aux)


}

# Plota a decomposição m$svd de um hipercubo m$a
plot_svd_cube = function(m,main,type=c(1,1),cex=0.6){

	tam = dim(m$a)[1]

	somaum = apply(m$a,1,sum)
	auxcol = as.factor(somaum)
	#print(auxcol)
	taml = length(levels(auxcol))
	levels(auxcol)=red2blue_v1(taml)
	#print(auxcol)
	#readline()
	#aux = svd(m)
	#aux$us = aux$u %*% diag(aux$d)
	if (type[1]==1){
		plot(m$svd$us[,1],m$svd$us[,2],type="n",main=main);
		text(m$svd$us[,1],m$svd$us[,2],log(tam,2)-somaum,cex=cex,col=paste(auxcol))
	}
	if (type[2]==1){
		plot(m$svd$us[,1],m$svd$us[,2],type="n",main=main)
		text(m$svd$us[,1],m$svd$us[,2],0:(tam-1),cex=cex)
	}
	if (type[1]==2){
		plot(m$svd$us[,1],m$svd$us[,2],type="n",main=main)
		text(m$svd$us[,1],m$svd$us[,2],somaum,cex=cex,col=paste(auxcol))
	}
}

gen_svd_test = function(xseq,a,b,sd=0.2){

	y=c()
	#y = 0.75*x+2
	for (x in xseq){
		y = c(y,a*(x+rnorm(1,mean=0,sd=sd))+b)
		#z = a*(x+rnorm(1,mean=0,sd=sd))+b
		#y = c(y,z+rnorm(1,mean=0,sd=sd))
	}
	#return(cbind(x,y))
	return(y)
}

# TO DO: gera matrix de distancias Euclideanas entre os vetores linhas de db
# usa propriedade de que dist = Norm(v1-v2)
#make_matrix_dist_euclidean = function(db,tami=c(),tamj=c()){
#	if (is.null(tami)){
#		tamx = 1:dim(db)[1]
#		#tamy = 1:dim(db)[2]
#		aux = db[tamx,tamy]
#	} else {
#		aux = db[tami,tamj]
#		tamx = 1:dim(aux)[1]
#		#tamy = 1:dim(aux)[2]
#	}
#	
#	for (i in tamx){
#		for (j in tamj){
#			
#		}
#
#	}
#}

# Retorna um vetor de ids contíguos a partir de "i". Aceita uma quebra da contiguidade de diff
ids_in_serie = function(v,i,diff=2){

	tam = length(v)
	auxv = as.numeric(v)
	aux=c(v[i])
	
	if ((tam-i)<0) return(0)
	else {
		while ((tam-i)>0){
			auxdiff = auxv[i+1]-auxv[i]
			if (auxdiff <= diff){
				aux = c(aux,v[i+1])
				i = i+1
			} else {
				break
			}
		}
	}
	return(aux)
}


# Encontra a maior série de ids contiguos em v
find_bigger_ids_in_serie = function(v){

	tam = length(v)
	auxbig = 0
	auxi = 0
	auxlen = 0
	for (i in 1:tam){
		aux=ids_in_serie(v,i)	
		auxlen = length(aux)
		if (auxlen>auxbig){
			auxids = aux
			auxbig = auxlen
		}
	}
	return(auxids)
	
}

#Elimina ids repetidos na lista vl
no_repeat_id = function(vl){

	tam = length(vl$rn)
	#print(tam)
	aux = list()
	aux$rn=c()
	aux$id=c()
	auxdiff = vl$rn[1]
	aux$rn = c(aux$rn,auxdiff)
	aux$id = c(aux$id,vl$id[1])

	for (i in 2:tam){
		if (auxdiff != vl$rn[i]){
			#aux = c(aux,auxdiff)
			aux$rn = c(aux$rn,vl$rn[i])
			aux$id = c(aux$id,vl$id[i])
			auxdiff = vl$rn[i]
		}
	}
	return(aux)

}

#OBSOLETED: Elimina ids repetidos no vetor v
no_repeat_id_OLD = function(v){

	tam = length(v)
	#print(tam)
	aux = c()
	auxdiff = v[1]
	aux = c(aux,auxdiff)

	for (i in 2:tam){
		if (auxdiff != v[i]){
			#aux = c(aux,auxdiff)
			aux = c(aux,v[i])
			auxdiff = v[i]
		}
	}
	return(aux)

}

# OBSOLETED: Extrai o resid do name. Ex: iT17.CG2, extrai 17, separando em lista conforme a cadeia
extract_resi_OLD = function(name,chain_e="e",chain_i="i"){

	tam = length(name)
	aux = list()
	aux$e = c()
	aux$i = c()

	for (i in 1:tam){
		auxatom = splitatom(name[i])[3]
		auxtype = splitatom(name[i])[1]
		if (auxtype == chain_e){
			aux$e = c(aux$e,auxatom)
		} else {
			aux$i = c(aux$i,auxatom)
		}
	}
	return(aux)

}

#Extrai o resid do name. Ex: iT17.CG2, extrai 17, separando em lista conforme a cadeia
extract_resi = function(name){

	tam = length(name)
	aux = c()
	
	for (i in 1:tam){
		auxatom = splitatom(name[i])[3]
		#auxtype = splitatom(name[i])[1]
		aux=c(aux,auxatom)
	}
	return(aux)

}

#Extrai o res name do name. Ex: iT17.CG2, extrai T, separando em lista conforme a cadeia
extract_resn = function(name){

	tam = length(name)
	aux = c()
	
	for (i in 1:tam){
		auxatom = splitatom(name[i])[2]
		#auxtype = splitatom(name[i])[1]
		aux=c(aux,auxatom)
	}
	return(aux)

}

#Extrai o resid do name.
extract_resid = function(name){

	tam = length(name)
	aux = c()
	
	for (i in 1:tam){
		auxatom = i
		#auxatom = splitatom(name[i])[2]
		#auxtype = splitatom(name[i])[1]
		aux=c(aux,auxatom)
	}
	return(aux)

}

# Tenta discriminar átomos e resíduos do inibidor como da alça (loop), não alça (noloop) e da enzima (base)
find_inibitory_loop = function(name,eivector,chain_1="E",chain_2="I"){

	aux = list()
	aux$loop = c()
	aux$noloop = c(0)
	aux$base = c()
	aux$baseid = c()
	aux$loopatom = c()
	aux$loopres = c()
	aux$baseatom = c()
	aux$baseres = c()
	auxresids=list()
	
	ne = name[eivector==chain_1]
	ni = name[eivector==chain_2]
	auxresids$e = list()
	auxresids$i = list()
	auxresids$e$rn = extract_resi(ne)
	auxresids$i$rn = extract_resi(ni)
	auxresids$e$id = extract_resid(ne)
	auxresids$i$id = extract_resid(ni)
	auxresids$i$id = auxresids$i$id + length(ne)
	#print(auxresids$i)
	#auxresids$eres = extract_resn(ne)
	#auxresids$ires = extract_resn(ni)
	#auxnorep.i = as.numeric(no_repeat_id(auxresids$i))
	#auxnorep.e = as.numeric(no_repeat_id(auxresids$e))
	auxnorep.i = list()
	auxnorep.i = no_repeat_id(auxresids$i)
	#Sprint(auxnorep.i)
	auxnorep.e = list()
	auxnorep.e = no_repeat_id(auxresids$e)
	#print(auxnorep.i)
	#print(auxnorep.e)
	#readline()
	aux$base = auxnorep.e$rn
	aux$baseid = auxnorep.e$id
	#print(auxnorep)
	auxbigger = find_bigger_ids_in_serie(auxnorep.i$rn)
	filter = auxnorep.i$rn %in% auxbigger
	#print(filter)
	#print(auxbigger)
	#print(auxnorep.i$rn)
	#readline()
	tamnorep = length(auxnorep.i$rn)
	tambigger = length(auxbigger)

	aux$loop = paste(as.numeric(auxbigger[1]):as.numeric(auxbigger[tambigger]))
	aux$loopid = auxnorep.i$id[filter]
	aux$noloop = auxnorep.i$rn[!filter]
	aux$noloopid = auxnorep.i$id[!filter]
	tamnoloop = length(aux$noloop)
	if (tamnoloop==0){
		aux$noloop=0
		aux$noloopid = 0
	}
	return(aux)
}

# OBSOLETED: Tenta classificar os resíduos do inibidor como da alça (loop) e não alça (noloop)
find_inibitory_loop_OLD2 = function(name,eivector,chain_1="E",chain_2="I"){

	aux = list()
	aux$loop = c()
	aux$noloop = c(0)
	aux$base = c()
	auxresids=list()
	
	auxresids$e = extract_resi(name[eivector==chain_1])
	auxresids$i = extract_resi(name[eivector==chain_2])
	#print(auxresids)
	#readline()
	#auxnorep.i = as.numeric(no_repeat_id(auxresids$i))
	#auxnorep.e = as.numeric(no_repeat_id(auxresids$e))
	auxnorep.i = no_repeat_id(auxresids$i)
	#print(auxnorep.i)
	auxnorep.e = no_repeat_id(auxresids$e)
	aux$base = auxnorep.e
	#print(auxnorep)
	auxbigger = find_bigger_ids_in_serie(auxnorep.i)
	
	tamnorep = length(auxnorep.i)
	tambigger = length(auxbigger)

	aux$loop = paste(as.numeric(auxbigger[1]):as.numeric(auxbigger[tambigger]))

	for (i in 1:tamnorep){
		auxsum = sum(as.numeric(auxnorep.i[i]==auxbigger))
		if (auxsum == 0){
			aux$noloop = c(aux$noloop,auxnorep.i[i])
		}
	}

	tamnoloop = length(aux$noloop)
	if (tamnoloop>1){
		aux$noloop = aux$noloop[2:tamnoloop]
	} 
	return(aux)
}


# OBSOLETED: Tenta classificar os resíduos do inibidor como da alça (loop) e não alça (noloop)
find_inibitory_loop_OLD = function(name){

	aux = list()
	aux$loop = c()
	aux$noloop = c(0)
	aux$base = c()

	auxresids = extract_resi(name)
	#print(auxresids)
	#readline()
	#auxnorep.i = as.numeric(no_repeat_id(auxresids$i))
	#auxnorep.e = as.numeric(no_repeat_id(auxresids$e))
	auxnorep.i = no_repeat_id(auxresids$i)
	#print(auxnorep.i)
	auxnorep.e = no_repeat_id(auxresids$e)
	aux$base = auxnorep.e
	#print(auxnorep)
	auxbigger = find_bigger_ids_in_serie(auxnorep.i)
	
	tamnorep = length(auxnorep.i)
	tambigger = length(auxbigger)

	aux$loop = paste(as.numeric(auxbigger[1]):as.numeric(auxbigger[tambigger]))

	for (i in 1:tamnorep){
		auxsum = sum(as.numeric(auxnorep.i[i]==auxbigger))
		if (auxsum == 0){
			aux$noloop = c(aux$noloop,auxnorep.i[i])
		}
	}

	tamnoloop = length(aux$noloop)
	if (tamnoloop>1){
		aux$noloop = aux$noloop[2:tamnoloop]
	} 
	return(aux)
}

#Estima melhor numero de clusters por PDBid dado o numero de elementos conexos extras (numelems) e um deslocamento d
find_better_k_clustering_per_id = function(m,numelems,d){

	tam = dim(m)[1]
	aux = c()
	
	if (d==0) print("WARNING: displacement can not be zero")
	
	for (i in 1:tam){
	
		aux = c(aux,m[i,numelems[i]+d])

	}
	return(aux)
}

#retorna vetor de atom pertencente a um res conforme name
res_to_atom_ids = function(resid,name){

	tam = 1:length(name)
	aux=c()
	for (i in tam){
		auxresid = splitatom(name[i])[3]
		if (auxresid==resid) aux=c(aux,i)
	}
	return(aux)
}

# Retorna vetor contendo os átomos pertencentes a vetor res
res_to_atom_vector = function(v,name){

	tam = 1:length(v)
	aux = c()
	
	for (i in tam){
		aux = c(aux,res_to_atom_ids(v[i],name))
	}
	return(aux)
}

# Acrescenta informações atômicas do loop em db
add_atom_loop_list = function(db,chain_1="E",chain_2="I"){

	aux = db
	ne = aux$ei==chain_1
	ni = aux$ei==chain_2
	tame = length(aux$name[ne])
	#tam = 1:length(aux$inter$loop)
	#print(aux$inter$loop)
	#print(ni)
	#readline()
	aux$inter$loopatomid =  res_to_atom_vector(aux$inter$loop,aux$name[ni])
	aux$inter$loopatomid = aux$inter$loopatomid + tame
	aux$inter$baseatomid =  res_to_atom_vector(aux$inter$base,aux$name[ne])
	if (aux$inter$noloopid[1]==0) aux$inter$noloopatomid = 0
	else {
		aux$inter$noloopatomid =  res_to_atom_vector(aux$inter$noloop,aux$name[ni])
		aux$inter$noloopatomid = aux$inter$noloopatomid + tame
	}
	return(aux)
}

# OBSOLETED: Retorna um vetor indicado das posicões do vetor v contendo o escalar a
# Faça v %in% a
select_scale_in_vector = function(a,v){

	auxfilter = as.factor(v)
	auxlevels = levels(auxfilter)
	auxprefilter = auxlevels==a
	levels(auxfilter) = auxprefilter
	auxfilter = as.logical(auxfilter)
	return(auxfilter)

}

# Associa atomos do looping com k cluster numbers definidos em inter
add_loopid_per_clustering = function(inter,kclus){

	aux = list()

	used = c()
	using = kclus$loopatomclus[1]

	#print(using)
	#readline()
	
	#auxlevels = levels(inter$loopatomid)
	
	tam = length(levels(as.factor(kclus$loopatomclus)))
	#print(tam)
	#readline()
	i=1
	while (!is.na(using)){
		aux[[i]]=list()
		#auxfilter = select_scale_in_vector(using,inter$loopatomclus)
		auxfilter = kclus$loopatomclus %in% using
		#if (length(auxfilter)==0) aux[[i]]$loop = 0
		if (sum(auxfilter)==0) aux[[i]]$loop = 0
		else aux[[i]]$loop=inter$loopatomid[auxfilter]
		
		#auxfilter = select_scale_in_vector(using,inter$baseatomclus)
		auxfilter = kclus$baseatomclus %in% using
		if (sum(auxfilter)==0) aux[[i]]$base = 0
		else aux[[i]]$base=inter$baseatomid[auxfilter]
		
		#auxfilter = select_scale_in_vector(using,inter$noatomclus)
		auxfilter = kclus$noloopatomclus %in% using
		if (sum(auxfilter)==0) aux[[i]]$noloop = 0
		else aux[[i]]$noloop=inter$noloopatomid[auxfilter]
		
		used = c(used,using)
		if (length(used)<tam){
			auxfilter = !(kclus$loopatomclus %in% used)
			using = kclus$loopatomclus[auxfilter]
			using = using[1]
		} else if (kclus$noloopatomclus[1] != 0){
			#print(inter$noloopatomclus)
			auxfilter = !(kclus$noloopatomclus %in% used)
			#print(auxfilter)
			using = kclus$noloopatomclus[auxfilter]
			#print(using)
			using = using[1]
			#print(using)
			#readline()
		} else using = NA
		#print(using)
		#print(aux[[1]])
		#readline()
		i=i+1
	}
	return(aux)
}


# OBSOLETED: Associa atomos do looping com k cluster numbers definidos em inter
add_loopid_per_clustering_OLD = function(inter){

	aux = list()

	used = c()
	using = inter$loopatomclus[1]
	
	#auxlevels = levels(inter$loopatomid)
	
	tam = length(levels(as.factor(inter$loopatomclus)))
	i=1
	while (!is.na(using)){
		aux[[i]]=list()
		#auxfilter = select_scale_in_vector(using,inter$loopatomclus)
		auxfilter = inter$loopatomclus %in% using
		#if (length(auxfilter)==0) aux[[i]]$loop = 0
		if (sum(auxfilter)==0) aux[[i]]$loop = 0
		else aux[[i]]$loop=inter$loopatomid[auxfilter]
		
		#auxfilter = select_scale_in_vector(using,inter$baseatomclus)
		auxfilter = inter$baseatomclus %in% using
		if (sum(auxfilter)==0) aux[[i]]$base = 0
		else aux[[i]]$base=inter$baseatomid[auxfilter]
		
		#auxfilter = select_scale_in_vector(using,inter$noatomclus)
		auxfilter = inter$noloopatomclus %in% using
		if (sum(auxfilter)==0) aux[[i]]$noloop = 0
		else aux[[i]]$noloop=inter$noloopatomid[auxfilter]
		
		used = c(used,using)
		if (length(used)<tam){
			auxfilter = !(inter$loopatomclus %in% used)
			using = inter$loopatomclus[auxfilter]
			using = using[1]
		} else if (inter$noloopatomclus[1] != 0){
			#print(inter$noloopatomclus)
			auxfilter = !(inter$noloopatomclus %in% used)
			#print(auxfilter)
			using = inter$noloopatomclus[auxfilter]
			#print(using)
			using = using[1]
			#print(using)
			#readline()
		} else using = NA
		#print(using)
		#print(aux[[1]])
		#readline()
		i=i+1
	}
	return(aux)
}

add_v_betweenness_old = function(g,mean=T,onlynoz=T,norm=F,k=2){

	aux=betweenness(g)
	if (onlynoz){
		aux = aux[aux!=0]
	}
	#print(aux)
	if (length(aux)==0){
		aux = 0
	}else if (norm){
		#aux = mean(aux/max(aux))
		aux = aux/max(aux)
	}
	if (mean){
		aux = round(mean(aux),k)
	}
	#print(aux);readline()
	return (aux)

}

add_betweenness = function(g,type,mean=T,onlynoz=T,norm=F,k=2){

	if (type=="e"){
		aux=edge_betweenness(g,direct=F)
	}
	if (type=="v"){
		aux=betweenness(g,direct=F)
	}
	if (onlynoz){
		aux = aux[aux!=0]
	}
	if (length(aux)==0){
		aux = 0
	} else if (norm){
		#print(aux)#;readline()
		aux = aux/max(aux)
		#print(aux);readline()
	}
	if (mean){
		aux = round(mean(aux),k)
	}
	#print(aux);readline()
	return (aux)

}

add_clus_net_info = function(db,net){

	tam = 1:length(db$exp$pam)
	db$exp$clus = list()
	print(db$pdbname)

	for (i in tam){
		auxv = db$exp$pam[[i]]$clustering
		auxf = as.numeric(levels(factor(auxv)))
		#print(auxf);readline()
		
		db$exp$clus[[i]] = list()
		#db$clus[[i]]$subclus = list()
		for (j in auxf){
			db$exp$clus[[i]]$subclus[[j]] = list()
			#print(i);print(j)#;readline()
			auxfj = auxv==j
			auxm = db$A$a[auxfj,auxfj]
			if (is.matrix(auxm)){# print("ops")
				#print(dim(auxm));readline()
				auxg = graph.adjacency(auxm,weighted=TRUE,mode="undirected",diag=F)
				db$exp$clus[[i]]$subclus[[j]]$g = auxg
				if (net[1]){
					#db$exp$clus[[i]]$subclus[[j]]$vbet = add_betweenness(auxg,type="v")
					#auxgi = auxg
					#E(auxgi)$weight = 1/E(auxgi)$weight
					auxcbet = centr_betw(auxg,directed=F)$centralization
					if (is.na(auxcbet)|is.nan(auxcbet)){
						#print(db$pdbname);print(i);print(j);readline()
						auxcbet = 0
					}
					db$exp$clus[[i]]$subclus[[j]]$cbet = auxcbet
				}
				#db$exp$clus[[i]]$subclus[[j]]$cbet = auxcbet
				if (net[2]){
					#db$exp$clus[[i]]$subclus[[j]]$vbet = add_betweenness(auxg,type="v")
					#auxgi = auxg
					#E(auxgi)$weight = 1/E(auxgi)$weight
					auxcclo = centr_clo(auxg,mode="all")$centralization
					if (is.na(auxcclo)|is.nan(auxcclo)){
						#print(db$pdbname);print(i);print(j);readline()
						auxcclo = 0
					}
					db$exp$clus[[i]]$subclus[[j]]$cclo = auxcclo	
				}
				#if (net[2]){
					#db$exp$clus[[i]]$subclus[[j]]$ebet = add_betweenness(auxg,type="e")
					#print(i);print(j)
					#db$exp$clus[[i]]$subclus[[j]]$cbet = centr_betw(auxg)	
				#}
			}else{
				auxm = as.matrix(auxm)
				auxg = graph.adjacency(auxm,weighted=TRUE,mode="undirected",diag=T)
				db$exp$clus[[i]]$subclus[[j]]$g = auxg
				if (net[1]){
					#db$exp$clus[[i]]$subclus[[j]]$vbet = 0
					#db$exp$clus[[i]]$subclus[[j]]$ebet = 0
					db$exp$clus[[i]]$subclus[[j]]$cbet = 0
					#db$exp$clus[[i]]$subclus[[j]]$cclo = 0
					#print(i);print(j)	
				}
				if (net[2]){
					#db$exp$clus[[i]]$subclus[[j]]$vbet = 0
					#db$exp$clus[[i]]$subclus[[j]]$ebet = 0
					#db$exp$clus[[i]]$subclus[[j]]$cbet = 0
					db$exp$clus[[i]]$subclus[[j]]$cclo = 0
					#print(i);print(j)	
				}
			}	
		}
	}
	return(db)
}

# Associa a cada atomid um cluster number conforme kclus
map_atomic_loop_clustering = function(db,kclus){

	#> f=1:88 %in% dbaAA[[1]]$inter$loopid
	#> dbaAA[[1]]$exp$pam[[6]]$clustering[f]

	aux = list()
	aux$loopatomclus = c()
	aux$baseatomclus = c()
	aux$noloopatomclus = c()
	
	tamv = 1:length(db$name)

	#print(tamv)
	#print(db$inter$loopid)
	auxf = tamv %in% db$inter$loopatomid
	#print(auxf)
	#readline()
	aux$loopatomclus = db$exp$pam[[kclus]]$clustering[auxf]
	
	auxf = tamv %in% db$inter$baseatomid
	aux$baseatomclus = db$exp$pam[[kclus]]$clustering[auxf]
	
	if (db$inter$noloopid[1]==0) aux$noloopatomclus = 0
	else {
		auxf = tamv %in% db$inter$noloopatomid
		aux$noloopatomclus = db$exp$pam[[kclus]]$clustering[auxf]
	}
	return(aux)
}

# Calcula o centro geométrico da matrix m 
geometric_center = function(m){

	#auxyz = pdb[atomids,]
	#print(auxyz)
	if (sum(is.na(m))!=0){
		print("WARNING: NAs found in geometric center calculation")
	}
	geomc = c(mean(m[,1]),mean(m[,2]),mean(m[,3]))
	#print(geomc)
	#readline()
	return (geomc)
}


add_dists_matrix_per_super_clustering = function(db,geom_dist=T,shorted_path=T,gd=F,all_path=T){

	tam = 1:length(db$exp$superclus)

	for (i in tam){
		#tamclus = 1:length(db$kclus[[i]]$clus)
		tamclus = 1:i
		auxgeomc = c()
		for (j in tamclus){
			#print("ok");readline()
			if (geom_dist){
				geomc = db$exp$superclus[[j]]$geomc
				auxd = dist(geomc)
				if (length(auxd)==0){
					auxd = as.matrix(0)
				}
				db$exp$superclus[[j]]$d = auxd
				if (gd){
					db$exp$superclus[[j]]$gd = graph.adjacency(as.matrix(auxd),weighted=TRUE,mode="undirected")
				}
				#print(db$exp$superclus[[j]]$gd);readline()
			}
			if (shorted_path){
				#round(distances(dba2[[5]]$exp$superclus[[4]]$g),0)
				if (all_path){
					spath = distances(db$exp$superclus[[j]]$ga)
					auxd = as.dist(spath)
					if (length(auxd)==0){
						auxd = as.matrix(0)
					}
				}else{
					spath = E(db$exp$superclus[[j]]$ga)$weight
					auxd = spath
				}
				db$exp$superclus[[j]]$p = auxd
			}
			#auxv = unlist(db$kclus[[i]]$clus[[j]])
			#filter = which(db$exp$pam[[i]]$clustering==j)
			#filter = which(auxv!=0)
			#if (length(filter)==0) print(paste("WARNING: no id found in",db$pdbname,"superclus",i,"subclus",j,"for geometric center"))
			#atomids = auxv[filter]
			#atomids = filter
			#if (length(atomids)==1){
			#	auxyz = db$pdb$xyz[atomids,]
			#	auxyz = t(auxyz)
			#}
			#else auxyz = db$pdb$xyz[atomids,]
			#print(i);print(j);print(auxyz);readline()
			#auxgeomc = rbind(auxgeomc,geometric_center(auxyz))
			#aux$kclus[[i]]$clus[[j]]$geomc = auxgeomc
			#db$exp$superclus[[i]]$
		}
		#db$exp$superclus[[i]]$geomc = auxgeomc
	}
	return(db)

}

#add_center_node = function(db){
#
#}

# Adiciona o centro geométrico dos átomos de cada cluster em clus
add_geometric_center_per_super_clustering = function(db,k=3){

	#tam = 2:length(db$kclus)
	#tam = 2:length(db$exp$superclus)
	tam = 1:length(db$exp$superclus)
	#aux = db
	
	if (is.list(db$pdb)){
		for (i in tam){#;print(paste("i:",i))
			#tamclus = 1:length(db$kclus[[i]]$clus)
			tamclus = 1:i
			auxgeomc = c()
			for (j in tamclus){#;print(paste("j:",j))
				#auxv = unlist(db$kclus[[i]]$clus[[j]])
				filter = which(db$exp$pam[[i]]$clustering==j)
				#filter = which(auxv!=0)
				if (length(filter)==0) print(paste("WARNING: no id found in",db$pdbname,"superclus",i,"subclus",j,"for geometric center"))
				#atomids = auxv[filter]
				atomids = filter#;print(c("atomids: ",atomids))
				if (length(atomids)==1){
					auxyz = tomatrix(db$pdb$xyz[atomids,])#;print(auxyz)
					if (dim(auxyz)[1]>1)	auxyz = t(auxyz)
				}
				else auxyz = db$pdb$xyz[atomids,]
				#print(i);print(j);print(auxyz);readline()
				auxgeomc = rbind(auxgeomc,geometric_center(auxyz))
				#aux$kclus[[i]]$clus[[j]]$geomc = auxgeomc
				#db$exp$superclus[[i]]$
			}
			#print(auxgeomc)
			db$exp$superclus[[i]]$geomc = round(auxgeomc,k)
		}
	} else {
		print(paste("WARNING: no PDB coordinates found for",db$pdbname))
	}
	return(db)
}


# Adiciona o centro geométrico dos átomos de cada cluster em clus
add_geometric_center_per_clustering = function(db){

	tam = 2:length(db$kclus)
	aux = db
	
	if (is.list(db$pdb)){
		for (i in tam){
			tamclus = 1:length(db$kclus[[i]]$clus)
			for (j in tamclus){
				auxv = unlist(db$kclus[[i]]$clus[[j]])
				filter = which(auxv!=0)
				if (length(filter)==0) print(paste("WARNING: no id found in",db$pdbname,"klus",i,"clus",j,"for geometric center"))
				atomids = auxv[filter]
				if (length(atomids)==1){
					auxyz = db$pdb$xyz[atomids,]
					auxyz = t(auxyz)
				}
				else auxyz = db$pdb$xyz[atomids,]
				#print(i);print(j);print(auxyz);readline()
				auxgeomc = geometric_center(auxyz)
				aux$kclus[[i]]$clus[[j]]$geomc = auxgeomc
			}
		}
	} else {
		print(paste("WARNING: no PDB coordinates found for",db$pdbname))
	}
	return(aux)
}

make_cluster_filter = function(db,cutpar){

	tam = length(db$exp$superclus)

	auxvol = diag(db$exp$superclus[[tam]]$a)
	auxid = order(auxvol)
	#print(auxid);readline()

	if (cutpar<1){
		auxvolcut = cutpar*sum(auxvol)
	}else{
		auxvolcut = cutpar
	}
	#print(auxvolcut)

	auxs = 0
	auxids = c()

	for (i in 1:tam){
		auxi = auxs + auxvol[auxid[i]]
		#print(auxi);print(auxvolcut);realine()
		if (auxi>auxvolcut){
			break;
		}		
		auxs = auxi
		auxids = c(auxids,auxid[i])
	}
	#print(auxs);print(auxids);readline()

	#for (i in auxids){
	#}
	auxpam = db$exp$pam[[tam]]$clustering
	auxf = !(auxpam %in% auxids)
	auxr = list()
	auxr$atomfil = auxf
	auxr$mfil = auxids	

	#return(auxf)
	return(auxr)

	#aux = db
	#aux$A$a = db$A$a[auxf,auxf]
	#aux$name = db$name[auxf]
	#aux$id = db$id[auxf]
	#aux$ape = db$ape[auxf]
	#aux$ap = db$ap[auxf]
	#aux$pdb$atid = db$pdb$atid[auxf]
	#aux$pdb$xyz = db$pdb$xyz[auxf]
	#aux$exp=list()
	#aux$g = graph.adjacency(aux$A$a,weighted=TRUE,mode="undirected")
	#aux$L = NULL
	#aux$Lw = NULL
	return(auxf)
}


make_centrality = function(db){

	if(!is.null(db$A)){
		g = db$g
		tnet = data.frame(cbind(get.edgelist(g),E(g)$weight))
		tnet=filter(tnet,X1!=X2)
		tnet=symmetrise_w(tnet)
		db$centrality = closeness_w(tnet,gconly=F)
	}
	return(db)

}

make_all_cluster_filter = function(db,cutpar){


	if (!is.null(db$A)){
		aux = list()
		auxr = list()

		#auxf = make_cluster_filter(db[[i]],cutpar)
		auxr = make_cluster_filter(db,cutpar)
		#print(auxf);readline()
		#aux[[i]] = apply_filter(db[[i]],auxf)	
		aux = apply_filter(db,auxr$atomfil)
		#aux[[i]]$mv = auxr$mfil	
		aux$exp$preclus = list()
		auxlen = length(db$exp$superclus)
		if (!is.null(auxr$mfil)){
			aux$exp$preclus$ids = auxr$mfil
			a = db$exp$superclus[[auxlen]]$a[auxr$mfil,auxr$mfil]
			aux$exp$preclus$a = tomatrix(a)
			geomc = db$exp$superclus[[auxlen]]$geomc[auxr$mfil,]
			aux$exp$preclus$geomc = tomatrix(geomc)
		}
		return(aux)
	}else{
		return(db)
	}	
}

# OBSOLETED: Adiciona o centro geométrico dos átomos de cada cluster em clus
add_geometric_center_per_clustering_OLD = function(db){

	tam = 1:length(db$clus)
	aux = db
	
	for (i in tam){
		auxv = unlist(db$clus[[i]])
		filter = which(auxv!=0)
		atomids = auxv[filter]
		auxyz = db$pdb$xyz[atomids,]
		auxgeomc = geometric_center(auxyz)
		aux$clus[[i]]$geomc = auxgeomc
	}
	return(aux)
}

contact_area_vector = function(v,r1,r2,probe){

	mapply(contact_area,v,MoreArgs = list(r1=r1,r2=r2,probe=probe))

}

# Calcula área de contato conforme Silveira-Romanelli equation
contact_area = function(d,r1,r2,probe){

	#print(d);readline()

	R1 = r1+probe
	R2 = r2+probe

	if (d>0) {
		aux1 = 2*pi*(R1^2+R2^2)
		aux2 = pi*(R1+R2)*d
		aux3 = (1+((R1-R2)/d)^2)
		aux = aux1-aux2*aux3
		#print(aux1);print(aux2);print(aux3);
		if (aux>0) return(aux)
		else return (0)
	}
	else {
		#print("WARNING: division by zero or negative d em contact_area")
		return(0)
	}

}


# Retorna c(id,x,y,z) do atom em name 
return_xyz_from_pdb = function(name,pdb,verbose=FALSE,typename=1){

	auxm = c()
	aux = splitatom(name,typename=typename)
	print(aux)#;readline()
	#auxch = aux[1] #chain
	#auxrn = aux[2] #resn
	#auxri = aux[3] #resid
	#auxan = aux[4] #atomn

	f1 = pdb$atom$chain == aux[1]
	f2 = aa321(pdb$atom$resid) == aux[2]
	f3 = pdb$atom$resno == aux[3]
	f4 = pdb$atom$elety == aux[4]
	
	f = as.logical(f1 * f2 * f3 * f4)
	
	if (sum(f)!=1) print(paste("WARNING: no xyz atom or more than one xyz atom found for ",name,sep="")) 
	auxm = c(which(f),pdb$atom$x[f],pdb$atom$y[f],pdb$atom$z[f])
	
	if (verbose){
		print(name)
		print(which(f))
		print(pdb$atom$chain[f])
		print(aa321(pdb$atom$resid[f]))
		print(pdb$atom$resno[f])
		print(pdb$atom$elety[f])
		print(auxm)#;readline()
	}
	#readline()
	#names(auxm)=NULL
	return(auxm)
}


# Acrescenta coordenadas XYZ do PDB conforme atomos em db$name
#insert_pdb_xyz = function(db,dirpathin,prefix="pdb_align_all/",sufix="a.pdb"){
insert_pdb_xyz = function(db,pdbxyz,prefix="",sufix=".pdb",typename=1){

	#auxfile = paste(dirpathin,prefix,db$pdbname,sufix,sep="")
	#print(auxfile)
	#if (!file.exists(auxfile)){
	#	auxfile = paste(dirpathin,prefix,toupper(db$pdbname),sufix,sep="")
	#}
	#if (!file.exists(auxfile)){
	#	auxfile = paste(dirpathin,prefix,tolower(db$pdbname),sufix,sep="")
	#}
	#if (file.exists(auxfile)){
	#	auxpdb = read.pdb(auxfile,ATOM.only=TRUE)
	if (length(pdbxyz)){
		auxpdb = pdbxyz#;print(auxpdb);readline()
		auxname = db$element_name
		auxids = db$atom_ids
		auxyz=list()
		#auxyz$atid = c()
		mxyz = matrix(pdbxyz$xyz,ncol=3,byrow=T)
		f = pdbxyz$atom$eleno %in% db$atom_ids
		auxyz$xyz = mxyz[f,]
		auxyz$atid = pdbxyz$atom$eleno[f]
		
		#tam = 1:length(auxname)
		#print("ok")
		#for (i in tam){
		#	tempxyz = return_xyz_from_pdb(auxname[i],auxpdb,typename=typename)#;print(tempxyz);readline()
		#	auxyz$atid = c(auxyz$atid,tempxyz[1])
		#	auxyz$xyz = rbind(auxyz$xyz,tempxyz[2:4])
		#}
		rownames(auxyz)=NULL
	} else {
		print(paste("Some problem with insert xyz for",db$pdbname))
		auxyz=c()
	}
	db$pdb = auxyz
	#return(auxyz)
	return(db)
}

# Acrescenta coordenadas XYZ do PDB conforme atomos em db$name
#insert_pdb_xyz = function(db,dirpathin,prefix="pdb_align_all/",sufix="a.pdb"){
insert_pdb_xyz_old = function(db,dirpathin,prefix="",sufix=".pdb",typename=1){

	auxfile = paste(dirpathin,prefix,db$pdbname,sufix,sep="")
	#print(auxfile)
	if (!file.exists(auxfile)){
		auxfile = paste(dirpathin,prefix,toupper(db$pdbname),sufix,sep="")
	}
	if (!file.exists(auxfile)){
		auxfile = paste(dirpathin,prefix,tolower(db$pdbname),sufix,sep="")
	}
	if (file.exists(auxfile)){
		auxpdb = read.pdb(auxfile,ATOM.only=TRUE)
		auxname = db$element_name
		auxyz=list()
		auxyz$atid = c()
		auxyz$xyz = c()
		
		tam = 1:length(auxname)
		#print("ok")
		for (i in tam){
			tempxyz = return_xyz_from_pdb(auxname[i],auxpdb,typename=typename)#;print(tempxyz);readline()
			auxyz$atid = c(auxyz$atid,tempxyz[1])
			auxyz$xyz = rbind(auxyz$xyz,tempxyz[2:4])
		}
		rownames(auxyz)=NULL
	} else {
		print(paste("WARNING: file",auxfile,"do not exist"))
		auxyz=c()
	}
	db$pdb = auxyz
	#return(auxyz)
	return(db)
}


# Testa e retorna um vetor de atom names C,O,N,S para cálculo de Ac
# Detecta algumas ambiguidades...
atom_name_test_ok = function(vatom,vdw){

	tamv = length(vatom)
	tamvdw = length(vdw)

	auxv=c()

	aux = rep(0,tamvdw)

	for (i in 1:tamv){		
		auxatom = splitatom(vatom[i])[4]
		aux = rep(0,tamvdw)
		for (j in 1:tamvdw){
			auxname = names(vdw)[j]
			auxfilter = grepl(auxname,auxatom)
			#print(auxfilter);print(names(vdw)[j])
			if (auxfilter) {
				aux[j]=aux[j]+1
				auxv = c(auxv,auxname)
			}
		}
		auxsum = sum(aux)
		if (auxsum!=1) {
			print("WARNING: atom name ambiguation detected")
			return(FALSE)
		}
	}
	return(auxv)

}

set_atom_numbers = function(label,symbol,radius,num=T){

	#label = c()
	#label = c(symbol,paste0(symbol,toupper(letters)))
	if (num) label = c(label,paste0(symbol,0:9))
	vdw = rep(radius,length(label))
	names(vdw)=label
	return(vdw)

}

#http://www.cgl.ucsf.edu/chimera/1.1700/docs/UsersGuide/midas/vdwtables.html
#http://periodictable.com/Properties/A/VanDerWaalsRadius.v.html
set_atom_names = function(type=2){

	if (type==1){
		vdw = c()
		vdw0 = c("C","CA","CB","CD","CE","CG","CZ","CH")
		vdw = c(vdw,set_atom_numbers(vdw0,"C",1.700))
		vdw0 = c("O","OD","OE","OH","OG","OA","OX")
		vdw = c(vdw,set_atom_numbers(vdw0,"O",1.520))
		vdw0 = c("N","ND","NE","NH","NZ")
		vdw = c(vdw,set_atom_numbers(vdw0,"N",1.550))
		vdw0 = c("S","SD","SG")
		vdw = c(vdw,set_atom_numbers(vdw0,"S",1.800))
		vdw = c(vdw,MEAN=1.657) #media dos mais frequentes
		vdw0 = c("P","PA","PB")
		vdw = c(vdw,set_atom_numbers(vdw0,"P",1.800))
		vdw0 = c("ZN")
		vdw = c(vdw,set_atom_numbers(vdw0,"ZN",1.390,F))
		vdw0 = c("CU")
		vdw = c(vdw,set_atom_numbers(vdw0,"CU",1.400,F))
	}else{
		vdw = c()
		vdw0 = c("C")
		vdw = c(vdw,set_atom_numbers(vdw0,"C",1.700,F))
		vdw0 = c("O")
		vdw = c(vdw,set_atom_numbers(vdw0,"O",1.520,F))
		vdw0 = c("N")
		vdw = c(vdw,set_atom_numbers(vdw0,"N",1.550,F))
		vdw0 = c("S")
		vdw = c(vdw,set_atom_numbers(vdw0,"S",1.800,F))
		vdw0 = c("P")
		vdw = c(vdw,set_atom_numbers(vdw0,"P",1.800,F))
		vdw0 = c("NA")
		vdw = c(vdw,set_atom_numbers(vdw0,"NA",1.020,F))
		vdw0 = c("CA")
		vdw = c(vdw,set_atom_numbers(vdw0,"CA",1.000,F))
		vdw0 = c("ZN")
		vdw = c(vdw,set_atom_numbers(vdw0,"ZN",0.73,F))
		vdw0 = c("CU")
		vdw = c(vdw,set_atom_numbers(vdw0,"CU",0.73,F))
		vdw0 = c("FE")
		vdw = c(vdw,set_atom_numbers(vdw0,"FE",0.61,F))
		vdw0 = c("MG")
		vdw = c(vdw,set_atom_numbers(vdw0,"MG",0.72,F))
		vdw0 = c("MN")
		vdw = c(vdw,set_atom_numbers(vdw0,"MN",0.83,F))
		vdw0 = c("NI")
		vdw = c(vdw,set_atom_numbers(vdw0,"NI",0.69,F))
		vdw0 = c("CO")
		vdw = c(vdw,set_atom_numbers(vdw0,"CO",0.65,F))
		vdw0 = c("CL")
		vdw = c(vdw,set_atom_numbers(vdw0,"CL",1.81,F))
		vdw0 = c("BR")
		vdw = c(vdw,set_atom_numbers(vdw0,"BR",1.96,F))
		vdw0 = c("I")
		vdw = c(vdw,set_atom_numbers(vdw0,"I",2.20,F))
		vdw0 = c("H")
		vdw = c(vdw,set_atom_numbers(vdw0,"H",1.000,F))
		vdw0 = c("SE")
		vdw = c(vdw,set_atom_numbers(vdw0,"SE",1.900,F))
		vdw = c(vdw,MEAN=1.00)#media dos nao C-O-N-S-P ate aqui

	}
	return(vdw)

}

#vdw = c(1.700,1.520,1.550,1.800,1.200,1.470,1.800,1.750, 1.400)
#names(vdw)=c("C","O","N","S","H","F","P","Cl","CU")

get_vdw = function(atom_name,atom_list,vdw,df=vdw["MEAN"],nf="X"){

	x = vdw[atom_list]
	na.id = which(is.na(x))#;print(na.id)
	if (length(na.id)>0){
		print(paste("WARNING: unspecified atom type found for",atom_name[na.id],"Associating a mean radius to it"))
		x[na.id] = df
		names(x)[na.id]=nf
	}
	return(x)

}


typify_atom_names = function(a,vdw,sep=c("[.]","[-]"),type=2){

	#print("ok");readline()
	#print(a[1:10])
	x = data.frame(strsplit(a,sep[1]))[3,]

	if (type==2){#print(data.frame(strsplit(as.character(unlist(x)),sep[2])));readline()
		x = data.frame(strsplit(as.character(unlist(x)),sep[2]))[3,]#;print(x);readline()
		x = as.character(unlist(x))
		x = get_vdw(a,x,vdw)
	}
	if (type==1){
		x = data.frame(strsplit(as.character(unlist(x)),sep[2]))[1,]#;print(x);readline()
		x = as.character(unlist(x))#;print(x);readline()
		#x = strtrim(x,1)
		x = strtrim(x,2)
	}
	#x = a
	#print(x);readline()
	#x = vdw[x]
	#na.id = which(is.na(x))#;print(na.id)
	#if (length(na.id)>0){
	#	print(paste("WARNING: unspecified atom type found for",a[na.id],"Associating a mean radius to it"))
	#	x[na.id] = df
	#	names(x)[na.id]=nf
	#}
	#print(x[1:10]);readline()
	return(x)

}

dist2area = function(v,vdw,probe,liminf){

	#print(ind);print(vdw.v[ind[1]]);print(vdw.v[ind[2]]);readline()
	area = contact_area(v[3],r1=vdw[v[1]],r2=vdw[v[2]],probe=probe)
	if (area>liminf){
		return(area)
	}else{
		return(0)
	}
	#print(x);readline()

}
#Transforma matriz de distância em matriz de área de contato através da equação de Silveira-Romanelli
silveira_romanelli_transformation = function(db,vdw,probe,liminf=5,p=3,maxdist=7){

#"EL35.CB"

	if (typeof(db)=="list"){
		m = db$A$a
		auxv = atom_name_test_ok(db$name,vdw)
		tami = dim(m)[1]
		tamj = dim(m)[2]
	
		auxm = matrix(rep(0,tami*tamj),nrow=tami,ncol=tamj)

		for (i in 1:tami){
			ri = vdw[auxv[i]]
			for (j in 1:tamj){
				rj = vdw[auxv[j]]
				auxca = contact_area(m[i,j],ri,rj,probe)
				if (auxca>=liminf) auxm[i,j] = round(auxca,p)
				else auxm[i,j] = 0
				#if (m[i,j]!=0) {print(auxm[i,j]);print(ri);print(rj);print(m[i,j]);readline()}
			}
		}
		db$A$a = auxm
		db$typename = "Area"
		return(db)
	}else{
		if (dim(db)[1]>1){
			m = db
			atom_names = colnames(m)
			vdw.v = typify_atom_names(atom_names,vdw)
			tam = dim(m)[1]
			auxm = Diagonal(tam,x=0)
			colnames(auxm)=atom_names
			ind = which((m>0&m<maxdist), arr.ind=T)
			indd = cbind(ind,m[ind])#;return(indd)#;print(indd);print(dim(indd));readline();return(indd)
			auxm[ind]=apply(indd,1,dist2area,vdw=vdw.v,probe=probe,liminf=liminf)
			#print(auxm[1:5,1:5]);readline()
			#print(dim(auxm));print(symmetric_check(auxm))#;readline()
			f = apply(auxm,1,Norm)!=0
			auxm = auxm[f,f]#;print(auxm[1:5,1:5]);readline()
			#print(dim(auxm));print(symmetric_check(auxm));readline()
			#print(x[1:10]);readline()
		}else{
			auxm=Matrix()
		}
		return(auxm)
	}
	


}

#Transforma matriz de distância em matriz de área de contato através da equação de Silveira-Romanelli
silveira_romanelli_transformation_old = function(db,vdw,probe,liminf=0,p=3){

	m = db$A$a
	tami = dim(m)[1]
	tamj = dim(m)[2]
	auxv = atom_name_test_ok(db$name,vdw)
	auxm = matrix(rep(0,tami*tamj),nrow=tami,ncol=tamj)

	for (i in 1:tami){
		ri = vdw[auxv[i]]
		for (j in 1:tamj){
			rj = vdw[auxv[j]]
			auxca = contact_area(m[i,j],ri,rj,probe)
			if (auxca>=liminf) auxm[i,j] = round(auxca,p)
			else auxm[i,j] = 0
			#if (m[i,j]!=0) {print(auxm[i,j]);print(ri);print(rj);print(m[i,j]);readline()}
		}
	}
	#return(auxm)
	db$A$a = auxm
	db$typename = "Area"
	return(db)

}


# Aplica filter à lista m, compondo outro conjunto de dados
apply_filter = function(m,filter){


	#onlyc = as.logical(diag(m$Lw$a))
	#print(sum(filter));readline()

	aux = list()
	if (!is.null(m$matrixname)) aux$matrixname = m$matrixname
	if (!is.null(m$matrixtype)) aux$matrixtype = m$matrixtype
	if (!is.null(m$pdbname)) aux$pdbname = m$pdbname
	if (!is.null(m$typename)) aux$typename = m$typename
	if (!is.null(m$enzcode)) aux$enzcode = m$enzcode
	if (!is.null(m$inbcode)) aux$inbcode = m$inbcode
	if (!is.null(m$chains)) aux$chains = m$chains
	if (!is.null(m$title)) aux$title = m$title

	aux$A = list()
	aux$A$a = m$A$a[filter,filter]
	if (!is.null(m$A$b)){
		aux$A$b = m$A$b[filter,filter]
		aux$bin = m$bin[filter]
	}
	if (!is.null(m$element_name)) aux$element_name = m$element_name[filter]
	if (!is.null(m$element_polarity)) aux$element_polarity = m$element_polarity[filter]
	if (!is.null(m$vdw)) aux$vdw = m$vdw[filter]
	if (!is.null(m$atom_ids)) aux$atom_ids = m$atom_ids[filter]
	#aux$weight = getnozero(aux$A$a)
	aux$id = m$id[filter]
	if (!is.null(m$ap)) aux$ap = m$ap[filter]
	if (!is.null(m$apg)) aux$apg = m$apg[filter]
	if (!is.null(m$ape)) aux$ape = m$ape[filter]
	if (!is.null(m$idn)) aux$idn = m$idn[filter]
	if (!is.null(m$ide)) aux$ide = m$ide[filter]
	if (!is.null(m$ei)) aux$ei = m$ei[filter]
	if (!is.null(m$rownames)) aux$rownames = m$rownames[filter]
	if (!is.null(m$colnames)) aux$colnames = m$colnames[filter]

	if (is.list(m$pdb)){
		#print(length(m$pdb))
		#print(m$pdb)
		#freadline()
		aux$pdb = list()
		aux$pdb$atid = m$pdb$atid[filter]
		aux$pdb$xyz = m$pdb$xyz[filter,]
	}
	aux$exp = list()
	#aux$exp$vlist = list()
	aux$exp$pam = list()
	#aux$exp$pam = c("")
	
	return(aux)
}

conn_to_adj_matrix = function(conn){

	
	dij = as.matrix(filter(conn,source==target)[c("source","target","source_volume")])
	tam = dim(dij)[1]
	a = matrix(0,nrow=tam,ncol=tam)#;print(conn)
	aij = as.matrix(conn[c("source","target","edge_sum")])#;print(aij);print(dij);readline()
	a[aij[,1:2]] = aij[,3]
	a[aij[,2:1]] = aij[,3] 
	a[dij[,1:2]] = dij[,3]
	colnames(a) = NULL	
	rownames(a) = NULL
	#print(a);readline()
	return(a)

}


conn_to_xyz_matrix = function(conn){

	tam = dim(conn)[1]
	#a = matrix(0,nrow=tam,ncol=3)#;print(conn) 
	aij = as.matrix(filter(conn,source==target)[c("source","target","s_x","s_y","s_z")])#;print(aij)#;readline()
	a = tomatrix(aij[,3:5])
	colnames(a) = NULL	
	rownames(a) = NULL
	#;print(a);readline()
	#if (!is.null(rot)){
	#	a = round(transform_by_rot_matrix(a,rot),k)
	#}

	return(a)

}


get_info_connections_group = function(conn,local,global,spot,dense,center){

	tam = length(conn)
	res=list()

	for (i in 1:tam){
		#i=5
		res[[i]] = list()#;print(conn[[i]]);readline()
		res[[i]]$a = conn_to_adj_matrix(conn[[i]])
		res[[i]]$xyz = conn_to_xyz_matrix(conn[[i]])#;print(res);readline()
		res[[i]]$sensitivity = local[[i]]
		res[[i]]$accuracy = global[[i]]
		res[[i]]$spot = list()
		res[[i]]$spot$score = spot$score[[i]]
		res[[i]]$spot$colors = spot$colors[[i]]
		res[[i]]$dense = dense[[i]]
		res[[i]]$center = center[[i]]
		#res[[i]]$ga = graph.adjacency(res[[i]]$a,weighted=TRUE,mode="undirected",diag=FALSE)
		#vertex_attr(res[[i]]$ga) = list(label=rep(1,i))
		#print(res);readline()
	}

	return(res)

}

add_superclus_from_align = function(align,n=5){

	tam = length(align)
	super = list()
	#matrix(rot,ncol=4,byrow=F)
	k = 1
	for (i in 1:tam){#print(paste("i:",i))
		super[[k]] = list()
		super[[k]]$a = align[[i]]$p1$a
		super[[k]]$geomc = align[[i]]$p1$xyz
		super[[k]]$rot = diag(4)
		super[[k]]$score_align = NA
		super[[k]]$score_group = align[[i]]$t1[1,6:7]
		super[[k]]$spot = align[[i]]$p1$spot
		super[[k]]$dense = align[[i]]$p1$dense
		super[[k]]$qualy = align[[i]]$p1$qualy
		o = order(align[[i]]$base_vector,decreasing=T)
		score_group = align[[i]]$t2[1,6:7]
		spot = align[[i]]$p2$spot
		dense = align[[i]]$p2$dense
		qualy = align[[i]]$p2$qualy
		k = k + 1
		for (j in 1:n){#print(paste("j:",j))
		#for (j in 5){
			#print(o)
			if (is.null(align[[i]]$score_list[[o[j]]])){
				j = lastj
			}
			super[[k]] = list()
			super[[k]]$a = align[[i]]$p2$a
			super[[k]]$geomc = align[[i]]$score_list[[o[j]]]$xyz2t
			super[[k]]$rot = matrix(align[[i]]$score_list[[o[j]]]$rot,ncol=4,byrow=F)
			super[[k]]$score_align = align[[i]]$score_list[[o[j]]]$score$scor[1]
			super[[k]]$score_group = score_group
			super[[k]]$spot = spot
			super[[k]]$dense = dense
			super[[k]]$qualy = qualy
			lastj = j
			k = k + 1
		}
	}
	return(super)
}

prune_superclus = function(super,maxcut){

	tamclus = length(super)
	repeat{
		if (tamclus>maxcut){
			super[tamclus] = NULL
			tamclus = length(super)
		}else{
			return(super)
		}
	}

}

get_connections_group_bsr = function(super,original=T,k=1){

	group = list()
	group$conn = list()
	group$vol = list()
	group$local = list()
	group$global = list()
	group$dense_t = list()
	group$dense_b = list()
	group$rot = list()
	group$spot = list()
	group$center = list()
	#group$spot = list()
	tam = length(super)

	for (i in 1:tam){
		ids = which(upper.tri(super[[i]]$a,diag=T),arr.ind=T)#;print(ids)
		idz = which(super[[i]]$a[ids]>0)#;print(idz)
		ids = tomatrix(ids[idz,])#;print(ids);readline()
		edge_sum = round(super[[i]]$a[ids],k)
		f = ids[,1] == ids[,2]
		edge_sum[f] = 0
		a_diag = diag(super[[i]]$a)
		tam = length(a_diag)
		source = ids[,1]
		target = ids[,2]
		source_volume = a_diag[ids[,1]]
		target_volume = a_diag[ids[,2]]
		#gscore = rep(super[[i]]$accall,tam)
		#source_gscore = gscore[ids[,1]]
		if (!is.null(super[[i]]$accbyclass)) {
			source_score = super[[i]]$accbyclass[,1][ids[,1]]
			#target_gscore = gscore[ids[,2]]
			target_score = super[[i]]$accbyclass[,1][ids[,2]]
		}else{
			source_score = rep(NA,dim(ids)[1])
			target_score = rep(NA,dim(ids)[1])
		}
		#source_density = super[[i]]$density$bipar[ids[,1]]
		#target_density = super[[i]]$density$bipar[ids[,2]]
		if (!is.null(super[[i]]$density)) {
			source_density_t = super[[i]]$density$total[ids[,1]]
			target_density_t = super[[i]]$density$total[ids[,2]]
			source_density_b = super[[i]]$density$bipar[ids[,1]]
			target_density_b = super[[i]]$density$bipar[ids[,2]]
		}else{
			source_density_t = rep(NA,dim(ids)[1])
			target_density_t = rep(NA,dim(ids)[1])
			source_density_b = rep(NA,dim(ids)[1])
			target_density_b = rep(NA,dim(ids)[1])
		}
		if(original){
			source_xyz = tomatrix(tomatrix(super[[i]]$geomc0)[ids[,1],])
			target_xyz = tomatrix(tomatrix(super[[i]]$geomc0)[ids[,2],])#;print(source_xyz);readline()
		}else{
			source_xyz = tomatrix(tomatrix(super[[i]]$geomc)[ids[,1],])
			target_xyz = tomatrix(tomatrix(super[[i]]$geomc)[ids[,2],])#;print(source_xyz);readline()
		}
		colnames(source_xyz) = c("s_x","s_y","s_z")
		colnames(target_xyz) = c("t_x","t_y","t_z")
		conn_tab = data.frame(source,source_volume,target,target_volume,edge_sum,source_xyz,target_xyz,source_score,target_score, source_density_t,target_density_t,source_density_b,target_density_b)
		rownames(conn_tab) = NULL
		vol_tab = filter(conn_tab,source==target)[c("source","source_volume","s_x","s_y","s_z")]
		colnames(vol_tab) = c("group","sum","x","y","z")
		#spot_tab = filter(conn_tab,source!=target)
		#spot_tab$source_score = NULL
		#spot_tab$target_score = NULL
		group$conn[[i]] = conn_tab
		group$vol[[i]] = vol_tab
		if (!is.null(super[[i]]$accbyclass)) {
			group$local[[i]] = as.vector(t(filter(group$conn[[i]],source==target)["source_score"]))
			group$global[[i]] = as.vector(super[[i]]$accall)
		}else{
			group$local = NULL
			group$global = NULL
		}
		if (!is.null(super[[i]]$density)) {
			group$dense_t[[i]] = as.vector(t(filter(group$conn[[i]],source==target)["source_density_t"]))
			group$dense_b[[i]] = as.vector(t(filter(group$conn[[i]],source==target)["source_density_b"]))
		}else{
			group$dense_t = NULL
			group$dense_b = NULL
		}
		if (!is.null(super[[i]]$rot)) {
			group$rot[[i]] = super[[i]]$rot
		}else{
			group$rot = NULL
		}
		if (!is.null(super[[i]]$score_align)) {
			group$score_align[[i]] = super[[i]]$score_align
		}else{
			group$score_align = NULL
		}
		if (!is.null(super[[i]]$score_group)) {
			group$score_group[[i]] = super[[i]]$score_group
		}else{
			group$score_group = NULL
		}
		if (!is.null(super[[i]]$spot)) {
			group$spot[[i]] = super[[i]]$spot
		}else{
			group$spot = NULL
		}
		if (!is.null(super[[i]]$dense)) {
			group$dense[[i]] = super[[i]]$dense
		}else{
			group$dense = NULL
		}
		if (!is.null(super[[i]]$qualy)) {
			group$qualy[[i]] = super[[i]]$qualy
		}else{
			group$qualy = NULL
		}
		if (!is.null(super[[i]]$centrality)) {
			group$center[[i]] = super[[i]]$centrality
		}else{
			group$center = NULL
		}
		group$conn[[i]]$source_score = NULL
		group$conn[[i]]$target_score = NULL
		group$conn[[i]]$source_density_t = NULL
		group$conn[[i]]$target_density_t = NULL
		group$conn[[i]]$source_density_b = NULL
		group$conn[[i]]$target_density_b = NULL
		#group$spot[[i]] = spot_tab
		#print(conn_tab);readline()
	
	}
	

	return(group)

}

re_list = function(db,n=2){

	res = list()
	tam = length(db)
	k = 0
	for (i in 1:tam){
		
		if (i%%n==1){
			k = k + 1
			res[[k]] = list()
			j = 1
		}
		res[[k]][[j]] = list()
		res[[k]][[j]] = db[[i]]
		j = j + 1

	}
	return(res)

}

glance = function(v){

	print(head(v));
	print(tail(v));

}

get_bsr_info_interactions = function(bsr,atom_names,atom_codes,k=1){

	#tami = dim(bsr)[1]
	#tamj = dim(bsr)[2]
	#tam = length(atom_names)
	#;print(atom_names);readline()
	#atom_codes = (atom_names,meso,aacode)#;print(pol);readline()

	#ids = which(bsr>0,arr.ind=T)#;print(dim(ids))
	bsr_tab = data.frame()
	ids = which(upper.tri(bsr),arr.ind=T) # indices triagulares superiores
	idz = which(bsr[ids]>0) # valores bsr > 0 
	ids = tomatrix(ids[idz,])#;print(dim(ids));readline()
	

	#ids = which(bsr>0,arr.ind=T);print(dim(ids))
	#ids = unique(ids);print(dim(ids));readline()
	bsr_value = round(bsr[ids],k)
	tam = length(bsr_value)
	source = atom_names[ids[,1]]
	target = atom_names[ids[,2]]
	source_atom_codes = atom_codes[ids[,1]]
	target_atom_codes = atom_codes[ids[,2]]
	f = source_atom_codes == "a"
	source_polarity = rep("POLAR",tam)
	source_polarity[f] = "NONPOLAR"
	f = target_atom_codes == "a"
	target_polarity = rep("POLAR",tam)
	target_polarity[f] = "NONPOLAR"
	bsr_tab = data.frame(source=source,source_polarity=source_polarity,target=target,target_polarity=target_polarity,bsr_value=bsr_value)

	return(bsr_tab)

}

### CENTRALITY ### verificar possibilidade de acrescenta-la aqui depois
get_bsr_info_elements = function(atom_name,polarity,value,xyz){

	tam = length(polarity)

	pol = rep("POLAR",tam)

	f = polarity=="a"
	pol[f] = "NONPOLAR"

	bsr_tab = data.frame(atom_name,pol,value,xyz)
	colnames(bsr_tab)=c("name","polarity","value","x","y","z")
	
	return(bsr_tab)

}

get_bsr_info_groups =function(atom_name,pam,maxcut){

	tam = length(pam)
	if (tam>maxcut) tam = maxcut
	bsr_group = list()	
	bsr_group[[1]] = data.frame(element=atom_name,group=pam[[1]]$clustering)

	if (tam>1){
		for (i in 2:tam){
			bsr_tab = data.frame(element=atom_name,group=pam[[i]]$clustering)
			o = order(pam[[i]]$clustering)
			bsr_group[[i]] = bsr_tab[o,]
		}
	}
	
	return(bsr_group)
}

# OBSOLETED: Elimina linhas zeros de m retornando outra matriz aux
only_connected_old = function(m,onlyc){


	#onlyc = as.logical(diag(m$Lw$a))

	aux = list()
	aux$pdbname = m$pdbname
	aux$typename = m$typename
	aux$enzcode = m$enzcode
	aux$inbcode = m$inbcode
	#aux$d = m$d[onlyc,onlyc]
	aux$A = list()
	aux$A$a = m$A$a[onlyc,onlyc]
	aux$name = m$name[onlyc]
	aux$weight = getnozero(aux$A$a)
	aux$id = m$id[onlyc]
	aux$ap = m$ap[onlyc]
	aux$apg = m$apg[onlyc]
	aux$ape = m$ape[onlyc] #new codification
	aux$idn = m$idn[onlyc]
	aux$ide = m$ide[onlyc]
	aux$exp = list()
	aux$exp$vlist = list()
	aux$exp$pam = list()
	#aux$exp$pam = c("")
	
	return(aux)
}

make_gaussians_for_ernest = function(){

	x = seq(400,900,by=0.2)
	g1 = 1000*dnorm(x,650,65)
	g2 = 80*dnorm(x,750,20)
	g12 = g1+g2
	auxa = c()

	auxa = cbind(x,g12)

	r = seq(0.9,0.1,by=-0.1)

	for (i in 1:length(r)){
		auxa = cbind(auxa,r[i]*g12)
	}
	return(auxa)

}

#GENERAL: Entrada geral de dados
input_general_data = function(filename,matrixname,sufix,typename="Values",row.names=0,col1onull=T,as.matrix=T,bin=F,weight=F,header=T,rnd=2,matrixtype="normal"){

	data = list()
	data$matrixname = matrixname
	data$matrixtype = matrixtype
	#data$typename = typename
	data$A = list()

	if (row.names>0){
		if (sufix==".csv"){
			m = read.csv(filename,header=header,sep=",",row.names=row.names)
			if (matrixtype=="sparse"){
				m = Matrix(m)
			}
		}
		if (sufix==".mtx"){
			m = readMM(filename)#,header=header,sep=",",row.names=row.names)
		}
		#data$colnames = 1:dim(m)[2]
		data$colnames = colnames(m)
		data$rownames = rownames(m)
	}
	else {
		if (sufix==".csv"){
			m = read.csv(filename,header=header,sep=",")
			if (matrixtype=="sparse"){
				m = Matrix(m)
			}
		}
		if (sufix==".mtx"){
			m = readMM(filename)#,header=header,sep=",",row.names=row.names)
		}
		data$colnames = colnames(m)
		data$rownames = rownames(m)
	}
	
	
	#print(data$rownames);print(data$colnames);readline()
	
	if (col1onull){
		tami = dim(m)[1]
		tamj = dim(m)[2]
		m = m[1:tami,2:tamj] #mudar no futuro...
		#if (as.matrix) m = as.matrix(m)
		#print(m);readline()
		#data$colnames = colnames(m)
		#data$rownames = rownames(m)
		if (row.names==0){
			data$colnames = data$colnames[2:tamj]
			data$rownames = data$rownames[1:tami]
		}
	}
	#print(data$colnames);print(data$rownames);readline()
	if (as.matrix){
		if (matrixtype=="normal") m = as.matrix(m)
		else m = Matrix(m)
	}
	colnames(m)=NULL
	rownames(m)=NULL
	data$A$a = round(m,rnd)
	if (bin) {
		data$A$b = binarize(data$A$a)
		data$bin = apply(data$A$b,1,sum)
	}
	if (weight){
		data$weight = getnozero(data$A$a)
	}
	data$exp = list()
	#data$exp$vlist = list()
	#data$exp$pam = list()
	return(data)
}

# GENERAL: Aplica filter à lista m, compondo outro conjunto de dados
apply_general_filter = function(m,rowfilter,colfilter=rowfilter,typename="Revalues"){


	aux = list()
	aux$matrixname = m$matrixname
	if(!is.null(typename)){
		aux$typename = typename
	} else {
		aux$typename = m$typename
	}
	aux$A = list()
	#m = read.csv(filename,header=T,sep=",")
	#print(dim(m$A$a))
	#print(rowfilter)
	aux$A$a = m$A$a[rowfilter,colfilter]
	#aux$colnames = colnames(m)[colfilter]
	#aux$rownames = rownames(m)[rowfilter]
	aux$colnames = m$colnames[colfilter]
	aux$rownames = m$rownames[rowfilter]	
	#print(aux$colnames);readline()
	#if (col1onull){
	#	tami = dim(m)[1]
	#	tamj = dim(m)[2]
	#	m = m[1:tami,2:tamj] #mudar no futuro...
	#	data$colnames = colnames(m)
	#	data$rownames = rownames(m)
	#}
	#if (as.matrix){
	#	m = as.matrix(m)
	#}
	#colnames(m)=NULL
	#rownames(m)=NULL
	#data$A$a = m
	if (!is.null(m$A$b)) {
		aux$A$b = m$A$b[rowfilter,colfilter]
		aux$bin = apply(aux$A$b,1,sum)
	}
	if (!is.null(m$weight)){
		aux$weight = getnozero(aux$A$a)
	}
	aux$exp = list()
	#data$exp$vlist = list()
	#data$exp$pam = list()
	return(aux)

	#onlyc = as.logical(diag(m$Lw$a))

}

threshold_filter = function(db,k,diagz=T){

	auxv = db$A$a<k
	db$A$a[auxv]=0
	if (diagz) diag(db$A$a)=0
	return(db)

}


make_super_matrix_areacut_old = function(db,areacutpar,id,minpar=4,pre=F,k=2){

	auxsuperclus = db$exp$superclus
	#minsensi = auxsuperclus[[i]]$accbyclass[,1]
	auxn = db$pdbname
	tamsuper = length(auxsuperclus)
	#tam = 2:tamsuper
	tam = 1:tamsuper
	
	auxr = c()
	auxnr = c()
	goahead = T	

	for (i in tam){
		#minsensi = round(min(auxsuperclus[[i]]$accbyclass[,1]),k)
		minsensi = round(min(auxsuperclus[[i]][[minpar]][,1]),k)
		#print(minsensi);readline()
		#if (minsensi==1){
		#	goahead = F
		#	if ((i+1)<=tamsuper){ 
		#		minsensiprox = round(min(auxsuperclus[[i+1]]$accbyclass[,1]),k)
		#		if (minsensiprox<1){
		#			goahead=T
		#		}
		#	}else{
		#		print(paste("WARNING: problem in find next min sensitivity less than 1"))
		#	}	
		#}
		#if (goahead){
		if (pre){
			if (minsensi<areacutpar[4]){
				break;
			}

		}
		if (goahead_superclus(superclus=auxsuperclus,i=i,minsensi=minsensi,tamsuper=tamsuper,minpar=minpar)){
			auxd = diag(auxsuperclus[[i]]$a)		
			#auxd = dist_vector(auxsuperclus[[i]]$geomc)
			auxcs = cutoff_scanning(auxd,init=areacutpar[1],step=areacutpar[2],imax=areacutpar[3])
			#print(auxcs);readline()
			auxr = rbind(auxr,auxcs)
			auxnr = c(auxnr,paste(id,auxn,i,sep="_"))
			#print(auxr);print(auxnr);readline()
			#minsensi = round(min(auxsuperclus[[i]]$accbyclass[,1]),k)
		}
		if (minsensi<areacutpar[4]){
			break;
		}
	}
	rownames(auxr)=auxnr
	return(auxr)
}

make_super_matrix_metric_old = function(db,metricpar,id,minpar=4,pre=F,k=2,gotest=T){

	auxsuperclus = db$exp$superclus
	#minsensi = auxsuperclus[[i]]$accbyclass[,1]
	auxn = db$pdbname
	#tam = 2:length(auxsuperclus)
	tamsuper = length(auxsuperclus)
	#tam = 2:tamsuper
	tam = 1:tamsuper

	auxr = c()
	auxnr = c()
	#auxd = c()
	#goahead = T	
	go = T

	for (i in tam){
		#minsensi = round(min(auxsuperclus[[i]]$accbyclass[,1]),k)
		minsensi = round(min(auxsuperclus[[i]][[minpar]][,1]),k)
		minaccu  = round(min(auxsuperclus[[i]][[minpar-1]][,1]),k)
		if (minaccu<metricpar[5]){
			break;
		}
		#if (pre){
		#	if (minsensi<metricpar[4]){
		#		break;
		#	}
	
		#}
		if (gotest){
			go = goahead_superclus(superclus=auxsuperclus,i=i,minsensi=minsensi,tamsuper=tamsuper,minpar=minpar)
		}
		if (go){

			auxd = c()
			if (metricpar[1]){
				auxd = c(auxd,round(auxsuperclus[[i]]$accall[1],k))
			}
			if (metricpar[2]){
				auxd = c(auxd,round(min(auxsuperclus[[i]]$accbyclass[,1]),k))
			}
			auxr = rbind(auxr,c(auxd,i,length(E(auxsuperclus[[i]]$ga)$weight)))
			auxnr = c(auxnr,paste(id,auxn,i,sep="_"))
		}
		if (minsensi<metricpar[4]){
			break;
		}
	}
	rownames(auxr)=auxnr
	return(auxr)
}

goahead_superclus = function(superclus,i,minsensi,tamsuper,minpar=4,k=2){

	goahead = T
	if (i==1){
		return(goahead)
	}
	if (minsensi==1){
		goahead = F
		if ((i+1)<=tamsuper){ 
			#minsensiprox = round(min(superclus[[i+1]]$accbyclass[,1]),k)
			minsensiprox = round(min(superclus[[i+1]][[minpar]][,1]),k)#;print(
			if (minsensiprox<1){
				goahead = T
			}
		}else{
			#goahead = F
			goahead = T
			print(paste("WARNING: problem in find next min accuracy or sensitivity less than 1"))
		}	
	}
	return(goahead)

}

c_vector = function(clus,n,type){

	auxr = c()

	for (i in 1:n){
		if (type=="bet"){
			auxi = clus$subclus[[i]]$cbet
		}

		if (type=="clo"){
			auxi = clus$subclus[[i]]$cclo
		}
		if (is.na(auxi)|is.nan(auxi)) print ("ops")
		auxr = c(auxr, auxi)

	}

	return(auxr)	

}

diff_zero = function(v) { 

	return(sum(v!=0)) 

}

make_all_rollback_test_super_clustering = function(superclus,i,par){

		auxa = superclus[[i]]$a
		auxga = superclus[[i]]$ga
		auxaccby = superclus[[i]]$accbyclass

		if (sum(auxaccby[,1]==0)|sum(auxaccby[,2]==0)) return(T)

		auxec = edge_connectivity(auxga)
		#print(superm)#;readline()
		if (auxec == 1){
			auxdz = apply(auxa,2,diff_zero)
			#print(auxdz)
			tam = which(auxdz==2)
			#print(tam)#;readline()
			for (j in tam){
				auxnb = as.vector(neighbors(auxga,j))
				auxnb = c(auxnb,j)
				#print(auxnb)
				auxmean = mean(auxaccby[auxnb,1])
				#print(auxmean)
				if (auxmean<par){
					#print(auxmean<par);readline()
					return(T)
				}
			}			
	
		}
		return(F)

}

make_all_rollback_test_super_clustering_old1 = function(superclus,i,par){

		auxa = superclus[[i]]$a
		auxga = superclus[[i]]$ga
		auxaccby = superclus[[i]]$accbyclass
		auxec = edge_connectivity(auxga)
		#print(superm)#;readline()
		if (auxec == 1){
			auxdz = apply(auxa,2,diff_zero)
			#print(auxdz)
			tam = which(auxdz==2)
			#print(tam)#;readline()
			for (j in tam){
				auxnb = as.vector(neighbors(auxga,j))
				auxnb = c(auxnb,j)
				#print(auxnb)
				auxmean = mean(auxaccby[auxnb,1])
				#print(auxmean)
				if (auxmean<par){
					#print(auxmean<par);readline()
					return(T)
				}
			}			
	
		}
		return(F)

}

sub_rownames =function(m,row_id,pos,value,sep="_"){

	rown = rownames(m)[row_id]
	rown = unlist(strsplit(rown,sep))
	if (pos==0){
		rown[length(rown)]=as.character(value)
	}else{
		rown[pos]=as.character(value)
	}
	rown = paste(rown,collapse=sep)
	rownames(m)[row_id]=rown
	return(m)

}

find_good_super_clus = function(db,i,pdbname,lowcut,ids,pol=c("ALL","PP","AA"),sep="_"){

	tam = length(db)
	good_table = data.frame()
	#print(tam)
	if (tam){
		for (j in 1:tam){
			minsen = min(db[[j]]$sensitivity)
			if (minsen<lowcut[2]) break
			acc = db[[j]]$accuracy
			if (acc<lowcut[1]) break
			pdbnamej = paste(c(pdbname,pol[ids[3]],j),collapse=sep)#;print(pdbname);readline()
			good_line = data.frame(pdbnamej,i,ids[2],ids[3],j,acc,minsen)
			good_table = rbind(good_table,good_line)

		}
	}else{
		good_table = data.frame(paste(c(pdbname,pol[ids[3]],0),collapse=sep),i,ids[2],ids[3],0,0,0)
	}
	colnames(good_table)=c("pdbname","pdbid","inter","pol","nclus","acc","minsen")
	return(good_table)

}


make_pre_superclus_table = function(db,ids,lowcut,sep="_"){

	tam = length(db)
	n_table = data.frame()	

	for (i in 1:tam){
		pdbname = db[[i]][[ids[2]]]$pdb
		#i_table = find_good_super_clus(db[[i]][[ids[2]]]$clus[[ids[3]]],i,pdbname,lowcut,ids)
		i_table = find_good_super_clus(db[[i]][[ids[2]]]$clus[[ids[3]]]$nclus,i,pdbname,lowcut,ids)
		n_table = rbind(n_table, i_table)
		#print(res_table);readline()
	}
	#print(n_table);readline()
	return(n_table)

}

find_nearest = function(v,a){

	if (!is.vector(v)){
		v = as.vector(unlist(v))
	}
	#print(v);print(a);
	id = which(v-a<0)[1]#;print(id);readline()
	return(id)


}

trim_tab = function(tab){

	grp = unique(tab$group)
	tr = data.frame()
	for (i in grp){
		z = filter(tab,group==i)#;print(median(z$nclus))
		f = ceiling(median(z$nclus))==z$nclus
		tr = rbind(tr, z[f,])#;print(res1);readline()
	}
	return(tr)	

}

try_best_pair_tabs = function(pdbids,t1,t2,n=3){

	tam1 = dim(t1)[1]
	tam2 = dim(t2)[1]
	x = c(t1$minsen,t2$minsen)#;print(x)
	y = pam(x,n)#;print(y$clustering)

	#print(y$clustering[(tam1+1):(tam1+tam2)]);readline()

	p1 = data.frame(t1,group=y$clustering[1:tam1])
	p2 = data.frame(t2,group=y$clustering[(tam1+1):(tam1+tam2)])

	t1 = p1[round(tam1/2):tam1,]
	t2 = p2[round(tam2/2):tam2,]

	if (t1$nclus[1]==1) print(paste("WARNING: failure to find a best pair table for",pdbids[1]))
	if (t2$nclus[1]==1) print(paste("WARNING: failure to find a best pair table for",pdbids[2]))

	t1 = trim_tab(t1)
	t2 = trim_tab(t2)

	res = list()
	res$p1 = p1
	res$p2 = p2
	res$t1 = t1
	res$t2 = t2

	return(res)

	#print(t1);print(t2);readline()

	#pdbn = c(rep(pdbids[1],dim(t1)[1]),rep(pdbids[2],dim(t2)[1]))
	#z = data.frame(pdbn,x,y$clustering)
	#z1 = z[z[,1]==pdbids[1],]
	#z2 = z[z[,1]==pdbids[2],]
	#print(z);print(z1);print(z2);readline()


}


filter_pre_table = function(tab,cluspar){

	#print(tab)
	tamtab = dim(tab)[1]
	if (tamtab==1){
		print(paste("WARNING: impossible to prefilter align table"))
		return(tab)
	}
	tseq = seq(cluspar[3],by=-0.1,length.out=cluspar[2])
	tamseq = length(tseq)	
	res = c()
	#print(tseq)
	for (i in 1:tamseq){
		#print(tseq[i])
		id=find_nearest(tab["minsen"],tseq[i])#;print(id)
		if (!is.na(id)){
			if (tab$nclus[id]<cluspar[1]){
				id1 = id + 1
				if (id1 < tam){
					if (!(id1 %in% res)){
						res=c(res,id1)
					}
				}else{
					if (!(id %in% res)){
						res=c(res,id)
					}
				}
			}else{
				if (!(id %in% res)){
					res=c(res,id)
				}
			}
		}else{
			print(paste("WARNING: it was not possible to find a nearest minsen"))
			break;
		}	
	}
	if (length(res)>1){
		med = round(median(tab[res,"nclus"]))#;print(med)
		if (!(med %in% tab[res,"nclus"])){
			res = c(res,which(med==tab[,"nclus"]))
			res = sort(res)
		}
		return(tab[res,])
	}else{
		return(tab[(tamtab-1):tamtab,])
	}
	#print(res);
	#print(tab[res,]);readline()
}


make_all_metric_super_clustering = function(db,id,metricpar,lowcut,minclusid,maxclusid,join=F,only_last=0,minpar=4,pre=F,sep="_", continuum=F,dcut=c(0.2,0.01)){

	if (!is.null(db$A)){
		auxmet = make_all_super_clustering(db=db,type="met",par=c(metricpar,lowcut),id=id,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
		#print(auxmet);readline()
		f = ((auxmet[,3]>=minclusid)&(auxmet[,3]<=maxclusid))
		tamf = sum(f)
		if (tamf>0){
			if (tamf>1){
				auxmet = auxmet[f,]
			}else{
				aux = auxmet[f,]
				aux = t(data.frame(aux))
				rownames(aux) = rownames(auxmet)[f]
				auxmet = aux
			}
		}else{
			print(paste("WARNING: it was not possible to reach the limits between minclusid ",minclusid,"and",maxclusid,"for",db$pdbname))
			auxmet = t(data.frame(c(rep(0,4),id)))
			rownames(auxmet) = paste(id,db$pdbname,"0",sep=sep)
		}
	}else{
		auxmet = t(data.frame(c(rep(0,4),id)))
		rownames(auxmet) = paste(id,db$pdbname,"0",sep=sep)
		#print(auxmet);readline()
	}
	#print(auxmet);readline()
	return(auxmet)

}


# Antiga make_super_matrix = function(db,type,par,id,minpar=4,pre=F,k=2,gotest=T,only_last=T,join=T){
make_all_super_clustering = function(db,type,par,id,minpar=4,pre=F,k=2,gotest=F,join=T,only_last=1,continuum=F){
	auxsuperclus = db$exp$superclus
	#minsensi = auxsuperclus[[i]]$accbyclass[,1]
	auxn = db$pdbname
	tamsuper = length(auxsuperclus)
	#tam = 2:tamsuper
	tam = 1:tamsuper;

	go = T
	
	auxr = c()
	auxnr = c()	
	#goahead = T
	#print(auxn);readline()
	for (i in tam){
		#print(i)
		minsensi = round(min(auxsuperclus[[i]][[minpar]][,1]),k)
		minaccu  = round(min(auxsuperclus[[i]][[minpar-1]][,1]),k)
		#print(minaccu);print(par[5]);readline()
		###if (minaccu<par[5]){ ### ATENCAO: retirei min accuracy fora. Verificar possiveis efeitos !
			###i=i-1
			###break
		###}
		#if (pre){
			#if (minsensi<par[4]){
		#	if ((minsensi<par[4])|(minaccu<par[5]){
		#		break;
		#	}

		#}
		if (gotest){
			go = goahead_superclus(superclus=auxsuperclus,i=i,minsensi=minsensi,tamsuper=tamsuper,minpar=minpar)
		}
		#print(auxn)
		#print(type)
		#print(go);readline()
		if (go){
			if (type=="geo"){
				#auxd = dist_vector(auxsuperclus[[i]]$geomc)
				auxd = as.vector(auxsuperclus[[i]]$d)
				see_inf=F
				scanning = T
				#print(auxd);readline()
			}
			if (type=="area"){
				auxd = diag(auxsuperclus[[i]]$a)
				see_inf=F
				scanning = T
			}
			if (type=="bet"){
				auxclus = db$exp$clus[[i]]
				auxd = c_vector(auxclus,i,type)
				see_inf=F
				scanning = T
				#print(auxd);readline()
			}
			if (type=="clo"){
				auxclus = db$exp$clus[[i]]
				auxd = c_vector(auxclus,i,type)
				see_inf=F
				scanning = T
				#print(auxd);readline()
			}
			if (type=="sho"){
				#auxd = dist_vector(auxsuperclus[[i]]$geomc)
				auxd = as.vector(auxsuperclus[[i]]$p)
				see_inf=F
				scanning = T
				#print(auxd);readline()
			}
			if (type=="met"){
				scanning = F
				auxd = c()
				if (par[1]){
					auxd = c(auxd,round(auxsuperclus[[i]]$accall[1],k))
				}
				if (par[2]){
					auxd = c(auxd,round(min(auxsuperclus[[i]]$accbyclass[,1]),k))
				}
				auxr = rbind(auxr,c(auxd,i,length(E(auxsuperclus[[i]]$ga)$weight),id))
				#auxnr = c(auxnr,paste(id,auxn,i,sep="_"))
			}
			if (length(auxd)==0){
				auxd = 0
			}
			if (scanning){
				auxcs = cutoff_scanning(auxd,init=par[1],step=par[2],imax=par[3],see_inf=see_inf,continuum=continuum)
				auxr = rbind(auxr,auxcs)
				#auxnr = c(auxnr,paste(id,auxn,i,sep="_"))
			}
			#print(auxd);print(auxcs);readline()
			#auxr = rbind(auxr,auxcs)
			auxnr = c(auxnr,paste(id,auxn,i,sep="_"))
		}
		if (minsensi<par[4]){
			#i = i - 1
		#if ((minsensi<par[4])|(minaccu<par[5])
			break
		}
		#print(auxr);print(auxnr);readline()
	}
	#print(auxnr);readline()
	rownames(auxr)=auxnr
	#print(auxr);print(i);readline()
	if (i>0){
		auxr = auxr[1:i,]
		auxnr = auxnr[1:i]
		if (!is.matrix(auxr)){
			auxr = t(as.matrix(auxr))
			rownames(auxr)=auxnr
		}
	}else{
		print(paste("WARNING: roll back fail for ",auxn,"in make_all_super_clustering"))
	}

	###print(db$pdbname);print(auxr);print(i);print(par[4])
	##roll_test = make_all_rollback_test_super_clustering(auxsuperclus,i,par[4])
	###print(roll_test)
	##if (roll_test){
	##	if ((i-1)>0){
	##		auxr = auxr[1:(i-1),]
	##		auxnr = auxnr[1:(i-1)]
	##		if (!is.matrix(auxr)){
	##			auxr = t(as.matrix(auxr))
	##			rownames(auxr)=auxnr
	##		}
	##	}else{
	##		print(paste("WARNING: roll back fail for ",auxn,"in make_all_super_clustering"))
	##	}
	##}
	#print(auxnr);print(auxr);readline()
	if (only_last){
		auxdim = dim(auxr)[1]
		#auxr = auxr[auxdim,]
		auxdim0 = auxdim-only_last+1
		if (auxdim0>0){
			auxr = auxr[auxdim0:auxdim,]
			auxnr = auxnr[auxdim0:auxdim]
		}else{
			auxr = auxr[1:auxdim,]
			auxnr = auxnr[1:auxdim]
		}
		#print(auxr);print(auxnr);readline()
		if (!is.matrix(auxr)){
			auxr = t(as.matrix(auxr))
			#rownames(auxr)=auxnr
		}
		rownames(auxr)=auxnr
		#auxgeor = rbind(auxgeor,auxgeo)
	}
	#print(auxr);readline() 
	if (join){
		auxdim = dim(auxr)[1]
		auxrni = rownames(auxr)[auxdim]
		if (type=="met"){
			#auxdim = dim(auxr)[1]
			#auxrni = rownames(auxr)[auxdim]
			auxr = auxr[auxdim,]
			#if (!is.matrix(auxr)){
			#	auxr = t(as.matrix(auxr))
			#	rownames(auxr)=auxrni
			#}
		}else{
			auxr=apply(auxr,2,sum)
			#print(auxr)
		}
		if (!is.matrix(auxr)){
			auxr = t(as.matrix(auxr))
			rownames(auxr)=auxrni
		}
		#rownames(auxr)=auxrni
		#auxrn = c(auxrn,auxrni)
		#print(auxr);readline()
		#auxr=apply(auxr,2,mean)
		#auxgeo=c(auxgeo,nc=auxgeodim)
		#if (!is.matrix(auxr)){
		#	auxr = t(as.matrix(auxr))
		#	rownames(auxr)=auxnr
		#}
	}
	#print("oi");print(auxr);readline()
	return(auxr)
}

#make_all_super_clustering(db=db[[i]],type="geo",par=c(geopar,lowcut),id=i,minpar=minpar,pre=pre,only_last=only_last,join=join)
make_all_super_clustering_old = function(db,type,par,id,minpar,pre,only_last,join){

	aux = make_super_matrix(db=db,type=type,par=par,id=id,minpar=minpar,pre=pre)
	#print(auxgeo);readline()
	if (only_last){
		auxdim = dim(aux)[1]
		aux = aux[auxdim,]
		#auxgeor = rbind(auxgeor,auxgeo)
	} else if (join){
		#auxdim = dim(auxgeo)[1]
		aux=apply(aux,2,sum)
		#auxgeo=c(auxgeo,nc=auxgeodim)
	}
	return(aux)
	#auxgeor = rbind(auxgeor,auxgeo)
}

make_super_matrix_geocut_old = function(db,par,id,minpar=4,pre=F,k=2){

	auxsuperclus = db$exp$superclus
	#minsensi = auxsuperclus[[i]]$accbyclass[,1]
	auxn = db$pdbname
	tamsuper = length(auxsuperclus)
	#tam = 2:tamsuper
	tam = 1:tamsuper
	
	auxr = c()
	auxnr = c()	
	#goahead = T

	for (i in tam){
		#minsensi = round(min(auxsuperclus[[i]]$accbyclass[,1]),k)
		minsensi = round(min(auxsuperclus[[i]][[minpar]][,1]),k)
		if (pre){
			if (minsensi<geocutpar[4]){
				break;
			}

		}
		if (goahead_superclus(superclus=auxsuperclus,i=i,minsensi=minsensi,tamsuper=tamsuper,minpar=minpar)){
			auxd = dist_vector(auxsuperclus[[i]]$geomc)
			if (length(auxd)==0){
				auxd = 0
			}
			auxcs = cutoff_scanning(auxd,init=par[1],step=par[2],imax=par[3])
			#print(auxd);print(auxcs);readline()
			auxr = rbind(auxr,auxcs)
			auxnr = c(auxnr,paste(id,auxn,i,sep="_"))
		}
		if (minsensi<geocutpar[4]){
			break;
		}
	}
	rownames(auxr)=auxnr
	return(auxr)
}

dist_vector = function(m){

	auxd = dist(m)
	auxd = as.vector(auxd)

	return(auxd)

}

cutoff_scanning = function(auxv,init=0,step=0.01,imax=0.5,acc=F,bycol=F,k=0,see_inf=T,zero_to_init=T,continuum=T,scale=""){

	#auxrx = c(0)
	auxrx = c()
	auxry = c()
	auxrz = c()
	auxr = c()

	auxryss = 0

	#auxv = as.vector(db$A$a)

	i1 = init
	if (zero_to_init){
		i1 = 0
		i2 = init
		auxrx = c(auxrx,i2)
		if (continuum){
			auxryi = which((auxv>=i1)&(auxv<i2))
			if (scale=="sqrt"){
				auxrys = round(sum(sqrt(auxv[auxryi])),k)
			}else{
				auxrys = round(sum(auxv[auxryi]),k)	
			}
		}else{
			auxrys = sum((auxv>=i1)&(auxv<i2))
		}
		auxry = c(auxry,auxrys)
		i1 = i2
	}
	
	while(i1<imax){

		i2 = i1+step
		auxrx = c(auxrx,i2)
		if (continuum){
			auxryi = which((auxv>=i1)&(auxv<i2))
			if (scale=="sqrt"){
				auxrys = round(sum(sqrt(auxv[auxryi])),k)
			}else{
				auxrys = round(sum(auxv[auxryi]),k)
			}
		}else{
			auxrys = sum((auxv>=i1)&(auxv<i2))
		}
		auxry = c(auxry,auxrys)
		if (acc){
			auxryss = auxryss+auxrys
			auxrz = c(auxrz,auxryss)
		}
		i1 = i2
	}
	if (acc){
		#db$A$cs = cbind(auxrx,auxry,auxrz)
		auxr = cbind(auxrx,auxry,auxrz)
	}else{
		#db$A$cs = cbind(auxrx,auxry)
		if (bycol){
			auxr = cbind(auxrx,auxry)
		} else {
			#auxr = rbind(auxrx,auxry)
			auxr = auxry
			names(auxr)=round(auxrx,k)
		}
	}
	#db$A$st = auxrx
	#db$A$cs = auxry
	if (see_inf){
		auxinf = sum(auxv==Inf)
		auxr = c(auxr,inf=auxinf)
	}
	return(auxr)
}

cutoff_scanning_old = function(db,init=0,step=0.01,imax=0.5,acc=F){

	auxrx = c(0)
	auxry = c(0)
	auxrz = c(0)

	auxryss = 0

	auxv = as.vector(db$A$a)

	i1 = init

	while(i1<imax){

		i2 = i1+step
		auxrx = c(auxrx,i2)
		auxrys = sum((auxv>=i1)&(auxv<i2))
		auxry = c(auxry,auxrys)
		if (acc){
			auxryss = auxryss+auxrys
			auxrz = c(auxrz,auxryss)
		}
		i1 = i2
	}
	if (acc){
		db$A$cs = cbind(auxrx,auxry,auxrz)
	}else{
		db$A$cs = cbind(auxrx,auxry)
	}
	#db$A$st = auxrx
	#db$A$cs = auxry
	return(db)
}

# Megaentrada de dados - OLD
input_data = function(filename,pdbname,typename,meso,aacode,sufix,matrixtype="normal",bin=FALSE,sep=",",sepn="[-]",codnew=0){

	data = list()
	data$matrixname = pdbname
	data$matrixtype = matrixtype
	data$pdbname = tolower(substr(pdbname,1,4))
	#data$typename = typename
	#data$d = as.matrix(read.csv(filename,header=F,sep=","))
	data$A = list()
	print(paste("Reading ",filename))
	if (sufix == ".csv"){
		m = read.csv(filename,header=T,sep=sep)
		tami = dim(m)[1]
		tamj = dim(m)[2]
		elen = as.character(m[1:tami,1])#;print(elen[1])
		m = m[1:tami,2:tamj]
		#elen = colnames(m);print(elen[1:10]) #obscuramente, o read.csv transforma "-" em "." para column
		elen1 = unlist(strsplit(elen[1],sepn))#;print(elen1);readline()
		if (length(elen1)==1){
			data$element_name = colnames(m)
			print(paste("WARNING: no atomid found for",pdbname))
		}else{
			elen = matrix(unlist(mapply(strsplit,elen,MoreArgs = list(split=sepn))),ncol=2,byrow=T)
			data$element_name = elen[,1]
			data$atom_ids = as.numeric(elen[,2])
		}
		#print(data$element_name[1:10])
		#print(data$atom_ids[1:10]);readline()
		names(data$element_name) = NULL
		names(data$atom_ids) = NULL
		m = as.matrix(m)
		colnames(m) = NULL
		rownames(m) = NULL
		data$A$a = m
		colnames(data$A$a) = NULL
		rownames(data$A$a) = NULL
		if (matrixtype=="sparse"){	
			data$A$a = Matrix(data$A$a)		
		}
	}
	if (sufix == ".mtx"){
		readMM(filename)
	}
	
	#print(m[1:10,1:10]);readline()
	#print(colnames(m[1]))
	#auxcodes=unlist(strsplit(colnames(m)[1],"[.]"))
	#print(colnames(m));readline()
	#print(substr(auxcodes[1],2,4))
	#readline()
	#data$enzcode = auxcodes[2]
	#data$inbcode = auxcodes[3]
	#print(auxcodes)
	#readline()
	#tami = dim(m)[1]
	#tamj = dim(m)[2]
	#m1 = m[1:tami,2:tamj]
	#print(colnames(m1));readline()
	#m1 = as.matrix(m1)
	#colnames(m1) = NULL
	#rownames(m1) = NULL
	#data$A$a = m1
	#colnames(data$A$a) = NULL
	#rownames(data$A$a) = NULL
	if (bin) {
		data$A$b = binarize(data$A$a)
		data$bin = apply(data$A$b,1,sum)
	}
	#data$A$a = as.matrix(read.csv(filename,header=F,sep=","))
	#data$name = preprocess02(preprocess01((read.csv(rotname,header=F,sep=";"))))
	#print("ok")
	#data$weight = getnozero(data$A$a) # vetor com os pesos não zeros 
	#print(data$weight)
	#readline()
	#data$id = 1:length(data$name) # numeração de cada átomo conforme entrada
	if (codnew==1){
		data$id = 1:length(data$element_name)
		auxcodes=unlist(strsplit(colnames(m)[1],"[.]"))
		data$element_name = preprocess02(colnames(m1))
		data$enzcode = auxcodes[2]
		data$inbcode = auxcodes[3]
		data$ape = new_codify_polarity(data$element_name,meso,aacode) #Novo aprimoramento com 6 definições de aacode
		auxcod = as.factor(data$ape)
		levels(auxcod) = c("a",rep("p",5))
		data$ap = as.character(auxcod)
		data$idn = paste(data$id,substr(data$element_name,2,2),sep="") #id e codigo 1 letra aa
		data$ide = substr(data$element_name,1,1) #"e" = enzima, "i" = inibidor
		names(data$weight) = NULL
		names(data$id) = NULL
		names(data$idn) = NULL
		names(data$ide) = NULL
		names(data$ap) = NULL
		names(data$element_name) = NULL
	}
	#OBSOLETED
	if (codnew==2){
		data$ap = charv(data$element_name,"C") #Apolar, Polar e C (alfa)
		data$apg = codify_polarity(data$element_name,"C","S") # Aprimoramento de charv, considera S apolar
		data$ape = new_codify_polarity(data$element_name,meso,aacode) #Novo aprimoramento com 6 definições de aacode
		data$idn = paste(data$id,substr(data$element_name,2,2),sep="") #id e codigo 1 letra aa
		data$ide = substr(data$element_name,1,1) #"e" = enzima, "i" = inibidor
		data$ei = rep("-",length(data$ide))
	}
	data$exp = list()

	#data$exp$vlist = list()
	#data$exp$pam = list()
	#data$exp$pam = c("")
	#data$exp$kmeans = list()
	#colnames(data$d) = NULL
	#rownames(data$d) = NULL
	#colnames(data$A$a) = NULL
	#rownames(data$A$a) = NULL
	#names(data$weight) = NULL
	#names(data$id) = NULL
	#names(data$idn) = NULL
	#names(data$ide) = NULL
	#names(data$ap) = NULL

	return(data)

}


# OBSOLETED: Megaentrada de dados - OLD

input_data_old = function(filename,rotname,pdbname,typename){

	data = list()
	data$pdbname = pdbname
	data$typename = typename
	#data$d = as.matrix(read.csv(filename,header=F,sep=","))
	data$A = list()
	data$A$a = as.matrix(read.csv(filename,header=F,sep=","))
	data$name = preprocess02(preprocess01((read.csv(rotname,header=F,sep=";"))))
	data$weight = getnozero(data$A$a) # vetor com os pesos não zeros 
	data$id = 1:length(data$name) # numeração de cada átomo conforme entrada
	data$ap = charv(data$name,"C") #Apolar, Polar e C (alfa)
	data$idn = paste(data$id,substr(data$name,2,2),sep="") #id e codigo 1 letra aa
	data$ide = substr(data$name,1,1) #"e" = enzima, "i" = inibidor
	data$exp = list()
	data$exp$vlist = list()
	data$exp$pam = list()
	#data$exp$pam = c("")
	#data$exp$kmeans = list()
	#colnames(data$d) = NULL
	#rownames(data$d) = NULL
	colnames(data$A$a) = NULL
	rownames(data$A$a) = NULL
	names(data$weight) = NULL
	names(data$id) = NULL
	names(data$idn) = NULL
	names(data$ide) = NULL
	names(data$ap) = NULL
	names(data$name) = NULL

	return(data)

}

# Monta as Laplacianas L, Lw, Ls
laps = function(m,anorm=FALSE,r=10,only=c(T,T)){

	if (!is.null(m$A)){
		aux = m
		aux = make_laplacian(aux,anorm=anorm,r=r,only=only)
		return(aux)
	}else{
		return(m)
	}
	#aux$Lw = make_decomp(aux$Lw)
	#aux$L = make_decomp(aux$L)
	

}

#GENERAL: decomposições EIG e SVD
decomps_general = function(m,A,L,Lw,EIG,SVD,centered=1,weight=c()){

	aux = m
	if (A) aux$A = make_decomp_general(aux$A,EIG,SVD,aux$pdbname,centered=centered,weight=weight)
	if (L) aux$L = make_decomp_general(aux$L,EIG,SVD,aux$pdbname)
	if (Lw) aux$Lw = make_decomp_general(aux$Lw,EIG,SVD,aux$pdbname)
	return(aux)

}


# decomposições EIG e SVD
decomps = function(m,L,Lw,EIG,SVD,karpack=0,r=10){

	aux = m
	if (!is.null(m$A)){
		if (!karpack){
			if (L) aux$L = make_decomp(m=aux$L,EIG=EIG,SVD=SVD,pdbname=aux$pdbname,r=r)
			if (Lw) aux$Lw = make_decomp(m=aux$Lw,EIG=EIG,SVD=SVD,pdbname=aux$pdbname,r=r)
		} else {
			if (L) aux$L = make_decomp_arpack(m=aux$L,EIG=EIG,SVD=SVD,pdbname=aux$pdbname,r=r,karpack=karpack)
			if (Lw) aux$Lw = make_decomp_arpack(m=aux$Lw,EIG=EIG,SVD=SVD,pdbname=aux$pdbname,r=r,karpack=karpack)
		}
		return(aux)
	}else{
		return(m)
	}

}


# OBSOLETED: Monta as Laplacianas L e Lw e suas decomposições EIG e SVD
laps_and_decomps = function(m){

	aux = m
	aux = make_laplacian(aux)
	aux$Lw = make_decomp(aux$Lw,TRUE,TRUE,"")
	aux$L = make_decomp(aux$L,TRUE,TRUE,"")
	return(aux)

}

# Monta a matriz Walk Laplacia a partir de m (Laplaciana) 
### QUE HORROR !!! REFAZER !!!
walk_laplacian = function(m){

	tamr = dim(m)[1]
	tamc = dim(m)[2]
	
	for (i in 1:tamr){
		d = m[i,i]
		if (d!=0){
			#print(d)
			for (j in 1:tamc){
				m[i,j] = m[i,j]/d
			}
		}
	}
	return(m)

}

#Somente sparse
walk.laplacian = function(A,L){

	d = apply(A,1,sum)
	D1 = Diagonal(x=1/d)
	Lw = D1 %*% L

	return(Lw)

}

# Monta as matrizes Laplacians L e Lw
make_laplacian = function(m,anorm=FALSE,r=10,only=c(T,T),verbose=T){

	m$g = graph.adjacency(m$A$a,weighted=TRUE,mode="undirected")
	#plot(ppf5$g)
	if (verbose) print(paste("Making Non-Normalized Laplacian L"))
	m$L = list()
	if (m$matrixtype=="normal"){
		m$L$a = as.matrix(graph.laplacian(m$g,norm=FALSE))
	}else{
		m$L$a = Matrix(graph.laplacian(m$g,norm=FALSE))
	}
	m$L$a = round(m$L$a,r)

	#m$Ls = list()
	#m$Ls$a = as.matrix(graph.laplacian(m$g,norm=TRUE))
	#m$Ls$a = round(m$Ls$a,r)
	if (only[2]){
		if (verbose) print(paste("Making Normalized Laplacian Lw"))
		if (anorm){
			m$Lw = list()
			if (m$matrixtype=="normal"){
				m$Lw$a = walk_laplacian(m$L$a)
			}else{
				m$Lw$a = walk.laplacian(m$A$a,m$L$a)			
			}
			m$Lw$a = round(m$Lw$a,r)
		}else{
			m$Lw = list()
			if (m$matrixtype=="normal"){
				m$Lw$a = as.matrix(graph.laplacian(m$g,norm=TRUE))
			}else{
				m$Lw$a = Matrix(graph.laplacian(m$g,norm=TRUE))
			}
			m$Lw$a = round(m$Lw$a,r)	
		}
	}
	if (!only[1]){
		if (verbose) print(paste("Flushing Non-Normalized Laplacian L"))
		m$L = NULL
	}		
	return(m)

}

subva = function(v,a){

	return(v-a)

}

#center_column = function(nc){
#}

eliminate_all_equal_column = function(db,bkp=T){

	if (bkp) {
		db$A$b = db$A$a
		db$colnamez = db$colnames
	}
	auxa = db$A$a
	#auxv = as.logical(apply(auxa,2,sum))
	auxv = apply(auxa,2,sd)!=0
	auxa = auxa[,auxv]
	db$A$a = auxa
	db$colnames = db$colnames[auxv]
	#print(dim(db$A$b))
	#print(dim(auxa))
	return(db)
}


eliminate_zero_column = function(db,bkp=T){

	if (bkp) {
		db$A$b = db$A$a
		db$colnamez = db$colnames
	}
	auxa = db$A$a
	auxv = as.logical(apply(auxa,2,sum))
	auxa = auxa[,auxv]
	db$A$a = auxa
	db$colnames = db$colnames[auxv]
	#print(dim(db$A$b))
	#print(dim(auxa))
	return(db)
}

eliminate_zero_column_old = function(db,bkp=T){

	if (bkp) {
		db$A$b = db$A$a
	}
	auxa = db$A$a
	auxv = as.logical(apply(auxa,2,sum))
	auxa = auxa[,auxv]
	db$A$a = auxa
	#print(dim(db$A$b))
	#print(dim(auxa))
	return(db)
}


# GENERAL: Faz decomposição EIG e SVD sobre m
make_decomp_general = function(m,EIG,SVD,pdbname="",r=6,centered=0,diag=T,weight=c()){

	warn = "WARNING: complex numbers in eigendecomposition of:"
	#print("ok")

	if (centered) {
		if (centered==1){
			aux = m$a-mean(m$a)
		}
		if (centered==2){
			aux = scale(m$a)
			if(sum(attr(aux,"scaled:scale")==0)>0){
				print("WARNING: in centering data, scale:rms may be zero generating NaN")
				#print(attr(aux,"scaled:scale")==0)
			}
			attr(aux,"scaled:center")=NULL
			attr(aux,"scaled:scale")=NULL

		}
		if (centered==3){
			aux = scale(m$a,scale=F)
			#if(sum(attr(mts,"scaled:scale")==0))>0){
			#	print("WARNING: scale:rms may be zero generating NaN")
			#}
			attr(aux,"scaled:center")=NULL
			attr(aux,"scaled:scale")=NULL
		}
		if (centered==4){
			
			aux1 = scale(m$a[,1])
			auxa = m$a[,2:dim(m$a)[2]]
			auxmin = min(auxa)
			auxa = auxa - auxmin
			auxmax = max(auxa)
			auxa = auxa/auxmax
			auxmean = mean(auxa)
			aux2 = auxa - auxmean
			aux = cbind(aux1,aux2)

		}
		if (centered==5){
			
			aux1 = scale(m$a[,1])
			auxa = m$a[,2:dim(m$a)[2]]
			auxmin = apply(auxa,2,min)
			#print(auxmin);readline()
			auxa = sweep(auxa,2,auxmin)
			#plot(auxa[,1],type="l");readline()
			#auxmin = min(auxa)
			#auxa = auxa - auxmin
			auxmax = apply(auxa,2,max)
			#print(auxmax);readline()
			auxa = sweep(auxa,2,auxmax,FUN="/")
			#plot(auxa[,1],type="l");readline()
			#auxmean = mean(auxa)
			#aux2 = auxa - auxmean
			aux2 = scale(auxa,scale=F)
			#plot(aux2[,1],type="l");readline()
			aux = cbind(aux1,aux2)

		}
		#if (centered==4){
		#	aux = a
		#	for (i in 1:dim(a)[2]){
		#		aux[,i] = center_column(m$a[,i])
		#	}
		#
		#}
	}
	else aux = m$a
	if (SVD){
		if (length(weight)>0){
			aux = aux*weight
		}
		#print(weight);readline()
		svd = list()
		svd$a = svd(aux)
		svd$c = aux
		svd$u = round(svd$a$u,r)
		auxs = round(diag(svd$a$d),r)
		svd$v = round(svd$a$v,r)
		svd$us = svd$u %*% auxs
		if (diag){
			svd$s = diag(auxs)
		}else{
			svd$s = auxs
		}
		if (centered == 4) svd$r = c(auxmin,auxmax,auxmean)
		#if (centered == 4) svd$r = c(auxmin,auxmax,auxmean)
		svd$a = NULL
		m$svd = svd
	}

	if (EIG){
		#print(aux[1:10,1:10])
		eig = list()
		eig$a = eigen(aux)
		#print("ok")
		eig$c = aux
		if (is.complex(eig$a$values)) print(paste(warn,pdbname))
		eig$v = round(eig$a$vectors,r)
		auxs = round(diag(eig$a$values),r)
		eig$vs = eig$v %*% auxs
		if (diag){
			eig$s = diag(auxs)
		}else{
			eig$s = auxs
		}
		eig$a = NULL
		m$eig = eig
	}
	return(m)
		#aux[,19:20] = 1.5*aux[,19:20]
		#aux[,1:9] = aux[,1:9]*c(0,seq(1,by=1,length.out=8))
		#aux[,1:9] = aux[,1:9]*c(0,seq(600,by=100,length.out=8)/600)
		#aux[,10:18] = aux[,10:18]*c(0,seq(1,by=0.1,length.out=8))
		#aux[,10:18] = aux[,10:18]*c(0,seq(60,by=10,length.out=8)/60)
		#aux[,19] = aux[,19]*8
}

re_scale = function(v,vmin=1,vmax=10){

	auxmax = max(v)
	auxmin = min(v)

	auxr = vmin + (vmax-vmin)*((v-auxmin)/(auxmax-auxmin))
	return(auxr)

}

re_scale_weight = function(v,part=10,vmin=1,vmax=2){

	#auxwei = as.numeric(colnames(dball[[1]]$A$a))
	#if (clear_na){
	auxv = is.na(v)
	v = v[!auxv]
	#}
	auxlen = length(v)
	if (part){
		auxpart = which(v==part)
		#auxlen = length(v)
		#print(auxweip);readline()
		auxwei1 = v[1:(auxpart-1)]/v[1]
		auxwei1 = re_scale(auxwei1,vmin=vmin,vmax=vmax)
		auxwei2 = v[auxpart:auxlen]/v[auxpart]
		#print(auxwei2);readline()
		auxwei2 = re_scale(auxwei2,vmin=vmin,vmax=vmax)
		auxwei3 = rep(1,sum(auxv))
		auxwei = c(auxwei1,auxwei2,auxwei3)
	} else {
		auxwei = v[1:auxlen]/v[1]
		#print(auxwei);readline()
		auxwei = re_scale(auxwei,vmin=vmin,vmax=vmax)
	}
	return(auxwei)

}

#Faz decomposicao ARPACK computando apenas k valores menores
#Funcionando apenas para EIG
make_decomp_arpack = function(m,EIG,SVD,pdbname="",centered=F,r=10,imgp=0.0001,diag=T,vs=F,karpack=20,which="SR",verbose=T){

	warne = "WARNING: complex numbers in eigen decomposition of:"
	warns = "WARNING: complex numbers in svd decomposition of:"

	tam = dim(m$a)[1]
	if (karpack>tam){
		print(paste("WARNING: karpack=",karpack,"is larger than dim matrix",tam,"Adjusting karpack to dim..."))
		karpack = tam
	}

	#print(karpack);readline()

	if (SVD){
		print(paste("ARPACK SVD decomposition has not been implemented yet!"))
	}
	if (EIG){
		eig = list()
		if (verbose) print(paste("Making ARPACK eigen decomposition..."))
		#eig$a = eigen(m$a)
		#print(karpack);print(which);
		eig$a = eigs(m$a,k=karpack,which=which)
		#print(eig$a);readline()
		complex = is.complex(eig$a$values)
		if (complex) {
			if (verbose) {
				print(paste(warne,pdbname))
				print(paste("trying round to near real.."))
			}
			sum.es = sum(Im(eig$a$values))
			if (sum.es<imgp){
				if (verbose) print(paste("SUCCESS!!!"))
				eig$complex2real = T
				eig$v = round(Re(eig$a$vectors),r)
				if (vs){
					auxs = abs(round(diag(Re(eig$a$values)),r))
					eig$vs = eig$v %*% auxs
					if (diag){
						eig$s = diag(auxs)
					}else{
						eig$s = auxs
					}
				}else{
					eig$s = abs(round(Re(eig$a$values),r))

				}
			}else{
				print(paste("Sorry... fail :-("))
				eig$complex2real = F
			}
		}else{
			if (verbose) print(paste("No complex results to deal..."))
			eig$v = round(eig$a$vectors,r)
			if (vs){
				auxs = round(diag(eig$a$values),r)
				eig$vs = eig$v %*% auxs	
				if (diag){
					eig$s = diag(auxs)
				}else{
					eig$s = auxs
				}
			}else{
				eig$s = abs(round(eig$a$values,r))
			}
		}
		eig$a = NULL
		m$eig = eig
	}
	#print(round(m$eig$s,2));readline()
	return(m)
}

# Faz decomposição EIG e SVD sobre m
make_decomp = function(m,EIG,SVD,pdbname="",centered=F,r=10,imgp=0.0001,diag=T,vs=F,verbose=T){

	warne = "WARNING: complex numbers in eigen decomposition of:"
	warns = "WARNING: complex numbers in svd decomposition of:"

	if (SVD){
		svd = list()
		svd$a = svd(m$a)
		complex = is.complex(svd$a$d)
		if (complex) {
			print(paste(warns,pdbname))
		}
		svd$u = round(svd$a$u,r)
		auxs = round(diag(svd$a$d),r)
		svd$v = round(svd$a$v,r)
		svd$us = svd$u %*% auxs
		if (diag){
			svd$s = diag(auxs)
		}else{
			svd$s = auxs
		}
		svd$a = NULL
		m$svd = svd
	}

	if (EIG){
		eig = list()
		if (verbose) print(paste("Making eigen decomposition..."))
		eig$a = eigen(m$a)
		if (verbose) print(paste("Verifying complex results..."))
		complex = is.complex(eig$a$values)
		if (complex) {
			print(paste(warne,pdbname))
			print(paste("trying round to near real.."))
			sum.es = sum(Im(eig$a$values))
			if (sum.es<imgp){
				print(paste("SUCCESS!!!"))
				eig$complex2real = T
				eig$v = round(Re(eig$a$vectors),r)
				if (vs){
					auxs = round(diag(Re(eig$a$values)),r)
					eig$vs = eig$v %*% auxs
					if (diag){
						eig$s = diag(auxs)
					}else{
						eig$s = auxs
					}
				}else{
					eig$s = round(eig$a$values,r)

				}
			}else{
				print(paste("Sorry... fail :-("))
				eig$complex2real = F
			}
		}else{
			if (verbose) print(paste("No complex result to deal..."))
			eig$v = round(eig$a$vectors,r)
			if (vs){
				auxs = round(diag(eig$a$values),r)
				eig$vs = eig$v %*% auxs	
				if (diag){
					eig$s = diag(auxs)
				}else{
					eig$s = auxs
				}
			}else{
				eig$s = round(eig$a$values,r)
			}
		}
		eig$a = NULL
		m$eig = eig
	}
	return(m)

}


rebuild_component = function(m,r,c1,all=1){

	auxm = m[,2:dim(m)[2]]
	auxmin = r[1]
	auxmax = r[2]
	auxmean = r[3]

	if (all==1){
		auxm = auxm + auxmean
	}
	if (all>1){
		auxm = auxm*auxmax
		auxm = auxm + auxmin
	}

	auxm = cbind(c1,auxm)

	return(auxm)

}


make_superclus_general = function(auxids,i,auxa){

	#auxids = aux$exp$pam[[i]]$clustering
	#tamids = levels(auxids)
	#print(tamids);readline()
	#auxm = matrix(rep(0,i*i),i,i)
	auxwclus[[i]] = list()
	#auxwclus[[i]]$gi = 
	#for (j in 1:i){
	#auxj = which(auxids == j)
	#if (!is.null(ma)){
	#	auxsupa = ma
	#}else{		
	auxsupa = make_matrix_superclus(auxids,auxa)
	#}
	#print(auxsupa);readline()
	auxwclus[[i]]$a = auxsupa
	if (!only_a){
		#diag(auxsupa) = 0 E(aux$g)$weight
		auxg = graph.adjacency(auxsupa,weighted=TRUE,mode="undirected",diag=FALSE)
		#print(i);print(E(auxg)$weight);readline()
		if (is.null(E(auxg)$weight)){
			diag(auxsupa)=1
			#print(auxsupa)
			auxg = graph.adjacency(auxsupa,weighted=NULL,mode="undirected")		
		}# else {
		#print("ok")
		auxwclus[[i]]$ga = auxg
		#}
		if (length(levels(factor(auxsupa)))>1){
			auxcm = confusionMatrix(as.table(auxwclus[[i]]$a))
			auxwclus[[i]]$accall = round(t(as.matrix(auxcm$overal[1:2])),k)
			#print(as.matrix(auxcm$byClass));readline()
			if (!is.matrix(auxcm$byClass)){
				auxcm$byClass = round(t(as.matrix(auxcm$byClass[1:2])),k)
				#print(auxcm$byClass);readline()
				auxcm$byClass = rbind(auxcm$byClass,auxcm$byClass[1,]) 
				auxwclus[[i]]$accbyclass = auxcm$byClass
				#print(auxwclus[[i]]$accbyclass);readline()
			}else{
			#print(auxcm$byClass);readline()
				auxwclus[[i]]$accbyclass = round(auxcm$byClass[,1:2],k)
			}
		}else{
			auxwclus[[i]]$accall = t(as.matrix(c(Accuracy=1,Kappa=1)))
			auxwclus[[i]]$accbyclass = t(as.matrix(c(Sensitivity=1,Specificity=1)))
		}
	}
	return(auxwclus[[i]])

}


make_superclus_reorientation = function(db,rotref,refid,k=3){

	if (!is.null(db$A)){
		#tamclus = length(db$exp$superclus)
		tam = length(db$exp$superclus)
		if (refid>tam) refid = tam
		m = db$exp$superclus[[refid]]$geomc
		target = diag_x(3)
		target[1,2] = 1
		target[2,1] = -1
		target[3,3] = 1 ## TESTE ##
		target = diag(3)#;print(m);print(target);print(refid);readline()
		if (dim(m)[1]>2){
			rotref = get_rot_matrix(source=m,target=target,svd=c(T,F),center=c(T,F))#;print(rotref);readline()
		}else{
			rotref = diag(4)
		}
		#db$pdb$xyz0 = db$pdb$xyz
		#db$pdb$xyz =  round(transform_by_rot_matrix(db$pdb$xyz0,rotref),k)

		for (i in 1:tam){
			geomc = db$exp$superclus[[i]]$geomc
			db$exp$superclus[[i]]$geomc0 = geomc
			db$exp$superclus[[i]]$geomc = round(transform_by_rot_matrix(geomc,rotref),k)
		}
		db$exp$rotref = rotref
	}
	return(db)

}	


make_superclus_for = function(ids,a,e_name,centrality=NULL,only_a,density,k){

	#tamids = levels(auxids)
	#print(tamids);readline()
	#auxm = matrix(rep(0,i*i),i,i)
	auxwclus = list()
	#auxwclus[[i]]$gi = 
	#for (j in 1:i){
	#auxj = which(auxids == j)
	#if (!is.null(ma)){
	#	auxsupa = ma
	#}else{		
	auxsupa = make_matrix_superclus(ids,a)#;print(auxsupa)#;readline()
	#}
	#print(auxsupa);readline()
	auxwclus$a = auxsupa
	if (!only_a){
		#diag(auxsupa) = 0 E(aux$g)$weight
		auxg = graph.adjacency(auxsupa,weighted=TRUE,mode="undirected",diag=FALSE)
		#vertex_attr(auxg) = list(label=1:dim(auxsupa)[1])
		vertex_attr(auxg) = list(label=rep(1,dim(auxsupa)[1]))
		#print("0k")
		#print(i);print(E(auxg)$weight);readline()
		if (is.null(E(auxg)$weight)){
			diag(auxsupa)=1
			#print(auxsupa)
			auxg = graph.adjacency(auxsupa,weighted=NULL,mode="undirected")
			#vertex_attr(auxg) = list(label=1:dim(auxsupa)[1])
			vertex_attr(auxg) = list(label=rep(1,dim(auxsupa)[1]))		
		}# else {
		#print("ok")
		auxwclus$ga = auxg
		#}
		if (length(levels(factor(auxsupa)))>1){
			auxcm = confusionMatrix(as.table(auxwclus$a))#;print(auxcm);readline()
			auxwclus$accall = round(t(as.matrix(auxcm$overal[1:2])),k)
			#print(as.matrix(auxcm$byClass));readline()
			if (!is.matrix(auxcm$byClass)){
				auxcm$byClass = round(t(as.matrix(auxcm$byClass[1:2])),k)
				#print(auxcm$byClass);readline()
				auxcm$byClass = rbind(auxcm$byClass,auxcm$byClass[1,]) 
				auxwclus$accbyclass = auxcm$byClass
				#print(auxwclus[[i]]$accbyclass);readline()
			}else{
			#print(auxcm$byClass);readline()
				auxwclus$accbyclass = round(auxcm$byClass[,1:2],k)
			}
		}else{
			auxwclus$accall = t(as.matrix(c(Accuracy=1,Kappa=1)))
			auxwclus$accbyclass = t(as.matrix(c(Sensitivity=1,Specificity=1)))
		}
		#print("ok")
		#auxwclus[[i]]$hot = rep(0,dim(auxsupa)[1])
		if (density){
			#print(a);print(e_name);print(ids);
			auxwclus$density = calculate_superclus_density(a,e_name,ids)#;print("ok")
			#readline()
		}
		#print("ok")
		if (!is.null(centrality)){
			auxwclus$centrality = calculate_superclus_centrality(centrality,ids)
		}
		#print("ok")
	}
	#print(auxwclus);readline()
	return(auxwclus)

}

calculate_superclus_centrality = function(center,ids,factor=100){

	#center = data.frame(center)#;print(center)
	res = c()
	id = unique(ids)#;print(ids)
	for (i in id){
		f = ids == i
		cm = mean(center[f,"n.closeness"])
		res = c(res,cm)#;print(res);readline()
	}
	res = round(factor*res)
	names(res) = id
	#print(res);readline()
	return(res)

}


mult_v = function(v,w){

	w[v[1]]*w[v[2]]

}

density_bipartite_graph = function(g,p,weilim=60,zero="0"){

	p0 = NULL
	p.n = names(p)#;print(p.n)
	f=p.n==zero
	if (sum(f)) p0 = p[f]
	#print(p);print(p0)#;readline()
	if (weilim){
		w = E(g)$weight
		w[w>weilim]=weilim
		w = w/weilim
		n = sum(w) 
		#;print(n);print(ecount(g));print(length(w));readline()
	}else{
		n = ecount(g)
	}
	s = 0
	if (n){
		#t1 = choose(sum(p),2)-sum(mapply(choose,p,MoreArgs=list(2)));print(t1)
		#s1=n/t1;print(s1);readline()
		#;print(p)
		tamp = length(p)
		t=0
		if (tamp>1){
			tc = combn(tamp,2)#;print(tc);readline()
			t = sum(apply(tc,2,mult_v,p))#;print(t)
		}
		if (!is.null(p0)){
			t0 = choose(p0,2)#;print(t0)
			t = t + t0
		}
		
		#print(t)
		s=n/t#;print(s);readline()
	}
	return(s)

}

density_weigth_graph = function(g,weilim=60){

	w = E(g)$weight
	w[w>weilim]=weilim
	w = w/weilim
	n = sum(w) 
	vc = vcount(g)
	t = choose(vc,2)
	s = n/t
	if (is.nan(s)){
		return(0)
	}else{
		return(s)
	}
}

calculate_superclus_density = function(m,e_name,ids,sep="[.]",ncol=3,weilim=60,factor=100){

	#print(table(ids));readline()
	res = list()#;print(ids)
	id = unique(ids)#;print(id)
	rex = c()
	rey = c()
	
	for (i in id){
	#	print(i)
	#for (i in 4){
		f = ids == i#;print(ids);print(i);print(f);print(m)
		mi = tomatrix(m[f,f])#;print(mi);readline()
		gi = graph.adjacency(mi,weighted=TRUE,mode="undirected",diag=FALSE)#;print(sep_atom_name_id(e_name[f],sep=sep,ncol=ncol));readline()
		e_frame = sep_atom_name_id(e_name[f],sep=sep,ncol=ncol)#;print(e_frame)
		f = e_frame[,1]==1
		if (sum(f)){
			e_frame[f,1] = e_frame[f,2]
		}
		e_table = table(e_frame[,1])#;print(e_table)#;readline()
		x = density_bipartite_graph(gi,e_table,weilim)#;print(x)#;readline()
		#x1 = density_bipartite_graph(gi,e_table,weilim=0);print(x1);readline()
		#if (length(e_table)>1){
			#parts = #;print(parts);readline()
		#	x = density_bipartite_graph(gi,e_table)
		#}else{
		#	e_factor = factor(e_frame[,2])
		#	levels(e_factor) = 1:length(levels(e_factor))
		#	e_v = as.character(e_factor)
		#	e_table = table(e_v);print(e_table);readline()
		#	x = density_bipartite_graph(gi,e_table);print(x);readline()
		#}
		if (weilim){
			y = density_weigth_graph(gi,weilim)
			#print(y);print(edge_density(gi));readline()
		}else{
			y = edge_density(gi)
		}
		#print(x);print(y);readline()
		rex = c(rex,x)
		rey = c(rey,y)
	}
	names(rex)=id
	names(rey)=id
	res$total = round(rey*factor)
	res$bipar = round(rex*factor)
	#print(res);readline()
	return(res)
}

#Monta mantriz de adjacencias do supercluster
make_matrix_superclus = function(idlist,a,zero=0.01,k=2){

	tam = levels(factor(idlist))

	#auxid = which(idlist == id)

	tamlen = length(tam)	

	auxm = matrix(rep(0,tamlen*tamlen),tamlen,tamlen)
	#print(tam);print(tamlen)#;readline()
	#print(sum(a));print(idlist)
	for (i in 1:tamlen){
		auxidlin = which(idlist == i)
		#auxidlin = which(idlist == i)
		#auxm2 = a[auxid,auxidcol]
		for (j in i:tamlen){
			auxidcol = which(idlist == j)
			#print(auxid);print(auxidcol);readline()
			auxm2 = a[auxidlin,auxidcol]
			auxs = sum(as.vector(auxm2))
			if (i == j){
				auxs = auxs/2
				if (auxs<zero) auxs = zero ### DANGER! ###
			}
			auxm[i,j] = round(auxs,k)
			auxm[j,i] = auxm[i,j]
			#print(i);print(j);print(auxs);print(auxm);readline()
		}
	}
	return(auxm)
}

make_pam_k_scanning_test = function(db,maxrange,init=1,type=1,oneg=T,conex=F,verbose=T,general=F,pamin=T,time=F){

	print("ok")
	return("ok")

}

make_pam = function(k,m,pamin=T,pdbname=""){

	ne = length(m$s)
	no = ne-k+1
	if (k>=ne){
		#print(paste("WARNING:k=",k,"must be less than n=",ne,"in",db$pdbname))
		print(paste("WARNING:k=",k,"must be less than n=",ne,"in",pdbname))
		return(NULL)
	}
	aux = pam(m$v[,no:ne],k)
	if (pamin){
		auxr = list()
		auxr$clustering = aux$clustering#;print("ok")
		auxr$silinfo = list()				
		auxr$silinfo$clus.avg.widths = aux$silinfo$clus.avg.widths
		auxr$silinfo$avg.width = aux$silinfo$avg.width
	}else{
		auxr = aux
	}
	return(auxr)

}


#Monta lista superclus onde cada cluster vira um nó. O rotulo do nó eh o volume das arestas internas
#e o rotulo das arestas do superclus eh o volume das arestas externas. Com isso, monta-se uma 
#matriz do supergrafo, que tambem pode ser vista como uma matriz de confusao. 
make_superclus = function(db,only_last=NULL,only_a=F,ma=NULL,k=3,verbose=T,geo=T,density=T){

	#aux = db
	#;print("o") 
	if (!is.null(db$A)){
		auxwclus = list()
		auxlen = length(db$exp$pam)#;print(auxlen);readline()
		if (!is.null(only_last)){
			auxids = db$exp$pam[[only_last]]$clustering#;print(auxids)#;readline()
			if (verbose) print(paste("Making superclus from pam for last",only_last,"clusters"))
			auxwclus[[1]] = make_superclus_for(auxids,db$A$a,db$element_name,db$centrality,only_a,density,k)#;print("ok")
			db$exp$superclus = c(db$exp$superclus,auxwclus)#;print("ok");readline()
		}else{
			tam = 1:auxlen
			for (i in tam){#;print(i)
			#for (i in 8){
				if (verbose) print(paste("Making superclus from pam with",i,"clusters"))
				auxids = db$exp$pam[[i]]$clustering
				auxwclus[[i]] = make_superclus_for(auxids,db$A$a,db$element_name,db$centrality,only_a,density,k)
			}
			db$exp$superclus = auxwclus
		}
		if (geo){
			db = add_geometric_center_per_super_clustering(db)
		}
		#print(db$exp$superclus);readline()
	}
	return(db)

}

make_pam_par = function(db,k0,maxrange,pamin,type,ncores){

	pamk=list()
	if (ncores){
		registerDoMC(ncores)
		pamk = foreach(k=k0:maxrange) %dopar% {
			if (type==1){
				make_pam(k=k,m=db$Lw$eig,pamin=pamin,pdbname=db$pdbname)
			}else{
				make_pam(k=k,m=db$Lw$svd,pamin=pamin,pdbname=db$pdbname)
			}
		}
	}else{
		i = 1
		for (k in k0:maxrange){			
			if (type==1){
				pamk[[i]] = make_pam(k=k,m=db$Lw$eig,pamin=pamin,pdbname=db$pdbname)
			}else{
				pamk[[i]] = make_pam(k=k,m=db$Lw$svd,pamin=pamin,pdbname=db$pdbname)
			}
			i = i + 1
		}
	}
	return(pamk)

}

#make_pam_superclus(db[[i]],lowcut,type=type,oneg=oneg,conex=conex,general=general,init=init,ncores=ncores)
make_pam_superclus = function(db,lowcut,init=1,type=1,oneg=T,conex=F,general=F,pamin=T,ncores=0){

	#ncores = 10
	if (!is.null(db$A)){
		ranges = data.frame(cut=c(0.9,0.7,0.5),n=c(10,5,5))
		n = dim(db$A$a)[1]
		k = init#;print("o")
		if (k==1){
			#oneg=F
			if (oneg){
				db$exp$pam[[1]] = list()#;print("o")
				db$exp$pam[[1]]$clustering = rep(1,n)#;print("o")
				db = make_superclus(db,only_last=T,geo=F)#;print("o")
			}else{
				db$exp$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
			}
		}
		#print(db$exp);readline()	
		if (k==1) k0 = 2
		maxrange=10
		stop = F
		tameigs = length(db$Lw$eig$s)-1
		if (maxrange>tameigs){
				print(paste("WARNING: maxrange before repeat does not reached; length of pam is smaller"))
			  print(paste("Probabily all zeros in eigenvalues for",db$pdbname))
				maxrange = tameigs
		}
		#print(k0);print(maxrange);print(lowcut[3]);readline()
		repeat{
			if (k0>maxrange) break;
			pamk = list()
			pamk = make_pam_par(db,k0,maxrange,pamin,type,ncores)
			#print(length(pamk))#;readline()
			db$exp$pam = c(db$exp$pam,pamk)
			tampam = length(db$exp$pam)#;print(tampam);readline()
			#print(mapply(table,lapply(db$exp$pam,"[[",1)));readline()
			#### PARALELIZAR AQUI NO FUTURO ####
			for (i in k0:maxrange){#;print(i)
			#for (i in 16){
				#if (i <= tampam){
					db = make_superclus(db,only_last=i,geo=F)#;print(db$exp$superclus[[i]])
					minsen = min(db$exp$superclus[[i]]$accbyclass[,1])
					if (minsen < lowcut[2]){
						stop = T
						break;
					}
					accall = db$exp$superclus[[i]]$accall[,1]
					if (accall<lowcut[1]){
						stop = T
						break;
					}
				#}else{
				#	break;
				#}
				#print(db$exp$pam[[k]]);readline()
				#print(table(db$exp$pam[[k]]$clustering));
				#k=k+1
				#if (k==lowcut[3]){
				#	break;
				#}
			}
			#if (maxrange>tampam){
			#	print(paste("WARNING: maxrange inside repeat does not reached; length of pam is smaller"))
			#  print(paste("Probabily all zeros in eigenvalues for",db$pdbname))
			#	break;
			#}
			#print(mapply(table,lapply(db$exp$superclus,"[[",1)));readline()
			if (stop){
				repeat{
					#print(i)
					tam = length(lapply(db$exp$pam,"[[",1))
					#print(tam);readline()
					if (tam>i){
						db$exp$pam[[tam]] = NULL
					}else{
						break;
					}
				}
				break;
			}#;print("ok")
			k0=maxrange+1
			if (k0>lowcut[3]){
				break;
			}
			dr = ranges[minsen>=ranges[,1],2][1]
			maxrange = maxrange+dr
			if (maxrange>lowcut[3]){
				#maxrange = lowcut[3]-1
				maxrange = lowcut[3]
			}
			if (maxrange>tameigs){
				print(paste("WARNING: maxrange after repeat does not reached; length of pam is smaller"))
			  print(paste("Probabily all zeros in eigenvalues for",db$pdbname))
				maxrange = tameigs
			}
			#print(k0);print(maxrange);print(lowcut[3]);readline()
		}
	}#;print("o")
	db = add_geometric_center_per_super_clustering(db)#;print("o")
	#print(length(lapply(db$exp$pam,"[[",1)))
	#print(length(lapply(db$exp$superclus,"[[",1)))
	#readline();
	return(db)
}


make_pam_k_scanning_par = function(db,maxrange,init=1,type=1,oneg=T,conex=F,verbose=T,general=F,pamin=T,time=F,ncores=2){

	if (!is.null(db$A)){
		#if (time) t = c()
		k = init
		#print(dim(db$A$a));readline()
		if (verbose){
			#paste(print("Over range:",k))
			print("Over range:")
			print(k:maxrange)
		}
		#n = dim(db$Lw$eig$v)[1]
		if (!general){
			if (!is.null(db$element_name)){
				n = length(db$element_name)
			}else{
				print(paste("WARNING: no size length found for",db$pdbname))
			}
		}else{
			if (!is.null(db$A$a)){
				n = dim(db$A$a)[1]
			}else{
				print(paste("WARNING: no size length found for",db$matrixname))
			}
		}

		if (k==1){
			if (oneg){
				db$exp$pam[[1]] = list()
				db$exp$pam[[1]]$clustering = rep(1,n) 
			}else{
				db$exp$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
			}
		}

		if (conex){
			auxsz = sum(db$Lw$eig$s==0)#;print(auxsz)
			maxrange = auxsz
		}
		#print(maxrange);readline()
		if (k==1) k = 2
		#for (k in 2:maxrange){
		registerDoMC(ncor)
		k0 = k
		res = foreach(k=k0:maxrange) %dopar% {
			if (type==1){
				make_pam(k=k,m=db$Lw$eig,pamin=pamin,pdbname=db$pdbname)
			}else{
				make_pam(k=k,m=db$Lw$svd,pamin=pamin,pdbname=db$pdbname)
			}
		}
		db$exp$pam = c(db$exp$pam,res) 
		##while (k<=maxrange){
			##ne = length(db$Lw$eig$s)
			##if (verbose) print(k)
			#expam = list()
			#db$exp$pam[[k]] = aux
			#n = dim(db[[i]]$Lw$svd$us)[1]
			##no = ne-k+1
			#print(paste(db[[i]]$pdbname,n,k))
			##if (k>=ne){
			##	print(paste("WARNING:k=",k,"must be less than n=",ne,"in",db$pdbname))
			##	break;
				#db[[i]]$exp$pam[[expid]] = 0
			##}else{
				#if (time) t = c(proc.time()[3],t)
				##if (type==1){
					#print(no);print(n);print(k);readline()
					##aux = pam(db$Lw$eig$v[,no:ne],k)
				##} else {
				##	aux = pam(db$Lw$svd$v[,no:ne],k)
				##}
				#aux = pam(db[[i]]$Lw$svd$us[,no:n],k)
				#exp1 = db[[i]]$exp$vlist[[1]]
				#exp2 = make_vector_indicator_list(aux$clustering)	
				#exp = make_logical_vector_product(exp1,exp2)
				#db[[i]]$exp$vlist[[expid]]=exp
				##if (pamin){
					#print(sum(db$exp$pam[[1]]$clustering==1));print(k);readline()
				##	db$exp$pam[[k]] = list()
				##	db$exp$pam[[k]]$clustering = aux$clustering#;print("ok")
				##	db$exp$pam[[k]]$silinfo = list()				
				##	db$exp$pam[[k]]$silinfo$clus.avg.widths = aux$silinfo$clus.avg.widths
				##	db$exp$pam[[k]]$silinfo$avg.width = aux$silinfo$avg.width
				##}else{
				##	db$exp$pam[[k]] = aux
				##}
				#if (time) t = c(proc.time()[3],t)
				#if (time) print(paste("Time PAM k:",k,round((t[1]-t[2])/60,1)))
				#exp = list()
				#exp = make_vector_indicator_list(db[[i]]$ide)
				#db[[i]]$exp[[expid]]=exp
			##}
			##k = k + 1
		##}
	}
	return(db)
}


make_pam_k_scanning = function(db,maxrange,init=1,type=1,oneg=T,conex=F,verbose=T,general=F,pamin=T,time=F,partam=0,ncores=10){


	if (!is.null(db$A)){
		if ((partam)&( dim(db$Lw$eig$v)[1]>partam )){#NAO HABILITAR PARTAM... com problemas nao identificado
			db = make_pam_k_scanning_par(db,maxrange=maxrange,type=type,oneg=oneg,conex=conex,general=general,init=init,ncores=ncores)
		}else{
			#if (time) t = c()
			k = init
			#print(dim(db$A$a));readline()
			if (verbose){
				#paste(print("Over range:",k))
				print("Over range:")
				print(k)
			}
			#n = dim(db$Lw$eig$v)[1]
			if (!general){
				if (!is.null(db$element_name)){
					n = length(db$element_name)
				}else{
					print(paste("WARNING: no size length found for",db$pdbname))
				}
			}else{
				if (!is.null(db$A$a)){
					n = dim(db$A$a)[1]
				}else{
					print(paste("WARNING: no size length found for",db$matrixname))
				}
			}

			if (k==1){
				if (oneg){
					db$exp$pam[[1]] = list()
					db$exp$pam[[1]]$clustering = rep(1,n) 
				}else{
					db$exp$pam[[1]] = list(t="Experiment with 1 cluster not exist for PAM")
				}
			}

			if (conex){
				auxsz = sum(db$Lw$eig$s==0)#;print(auxsz)
				maxrange = auxsz
			}
			#print(maxrange);readline()
			if (k==1) k = 2
			#for (k in 2:maxrange){
			while (k<=maxrange){
				ne = length(db$Lw$eig$s)
				if (verbose) print(k)
				#expam = list()
				#db$exp$pam[[k]] = aux
				#n = dim(db[[i]]$Lw$svd$us)[1]
				no = ne-k+1
				#print(paste(db[[i]]$pdbname,n,k))
				if (k>=ne){
					print(paste("WARNING:k=",k,"must be less than n=",ne,"in",db$pdbname))
					break;
					#db[[i]]$exp$pam[[expid]] = 0
				}else{
					if (time) t = c(proc.time()[3],t)
					if (type==1){
						#print(no);print(n);print(k);readline()
						aux = pam(db$Lw$eig$v[,no:ne],k)
					} else {
						aux = pam(db$Lw$svd$v[,no:ne],k)
					}
					#aux = pam(db[[i]]$Lw$svd$us[,no:n],k)
					#exp1 = db[[i]]$exp$vlist[[1]]
					#exp2 = make_vector_indicator_list(aux$clustering)	
					#exp = make_logical_vector_product(exp1,exp2)
					#db[[i]]$exp$vlist[[expid]]=exp
					if (pamin){
						#print(sum(db$exp$pam[[1]]$clustering==1));print(k);readline()
						db$exp$pam[[k]] = list()
						db$exp$pam[[k]]$clustering = aux$clustering#;print("ok")
						db$exp$pam[[k]]$silinfo = list()				
						db$exp$pam[[k]]$silinfo$clus.avg.widths = aux$silinfo$clus.avg.widths
						db$exp$pam[[k]]$silinfo$avg.width = aux$silinfo$avg.width
					}else{
						db$exp$pam[[k]] = aux
					}
					#if (time) t = c(proc.time()[3],t)
					#if (time) print(paste("Time PAM k:",k,round((t[1]-t[2])/60,1)))
					#exp = list()
					#exp = make_vector_indicator_list(db[[i]]$ide)
					#db[[i]]$exp[[expid]]=exp
				}
				k = k + 1
			}
		}
	}
	return(db)
}

# OBSOLETED: transforma uma matrix estatistica (ex: ks.test()) numa matrix centralizada com diagonal zero(se o caso)
make_stat_matrix = function(db,inter,p=3,ttype=c("No Centered","Centered","Diag=0")){
	aux=list()
	aux$A = list()
	if (inter==1){
		aux$type=ttype[1]
	} else if (inter==2){
		aux$type=paste(ttype[1],ttype[3],sep=" - ")
		diag(db)=0
	} else if (inter==3){
		aux$type=ttype[2]
		db = db - mean(db)
	} else {
		aux$type=paste(ttype[2],ttype[3],sep=" - ")
		db = db - mean(db)
		diag(db)=0
	}
	#if (inter==3)
	#if (center) aux$type = "Centered"
	#else aux$type = "No Centered"
	#if (inter) aux$type = paste(aux$type,ttype,sep="-")
	#diag(db)=0
	#if (inter==0) 
	#	if (center) db = db - mean(db)
	#while (inter){
	#	if (center) db = db - mean(db) #centrifica
	#	diag(db)=0
	#	inter = inter-1
	#}
	aux$A$a = round(db,p)
	return(aux)
}


make_general_stats_distribution = function(db,tam,names,verbose=TRUE,p=3){

	tamlen = tam
	auxm=matrix(rep(0,tamlen*tamlen),nrow=tamlen,ncol=tamlen)
	#print(auxm)
	#readline()
	for (i in 1:tamlen){
		for (j in 1:tamlen){
			x1 = db[,i]
			x2 = db[,j]
			aux = ks.test(x1,x2)
			auxm[i,j] = round(aux$p.value,p)
		}
	}
	#colnames(auxm) = names
	#rownames(auxm) = names
	if (verbose) print(auxm)
	return(auxm)
}

# Cria matriz de p-values ks.test modo continuum com opção de diag = 0
make_stat_input = function(db,diagz=1,type=1,name="dba",verbose=TRUE,p=3,subtext="Continuum - Diag="){

	index=1:length(db)
	aux=list()
	aux$A=list()
	if (is.list(db)){
		index=1:length(db)
		aux$A$a=make_mass_stats_weight_distribution(db,index,type=type,name,verbose,p)
		}
	else {
		index = dim(db)[2]
		aux$A$a=make_general_stats_distribution(db,index,"")
	}
	if (diagz) {
		diag(aux$A$a)=0
		aux$type=paste(subtext,as.numeric(!diagz),sep="")
	} else {
		aux$type=paste(subtext,as.numeric(!diagz),sep="")
	}
	return(aux)
}

shortest_path_kv_kernel = function(g,type="dk",bin=T,k=0){

	tam = length(g)#;print("0k")
	auxr = matrix(rep(0,tam*tam),nrow=tam)

	#print(g[[1]])
	#i=5;j=5
	for (i in 1:tam){
		for (j in i:tam){
			#print(paste(i,j))#;readline()
			if (type=="dk") auxr[i,j] = sp_kv_dk_kernel(g[[i]],g[[j]],bin=bin,k=k)
			if (type=="lk") auxr[i,j] = sp_kv_lk_kernel(g[[i]],g[[j]],bin=bin,k=k)
			#print(auxr[i,j]);readline()
			auxr[j,i] = auxr[i,j]

		}
	}
	#print("ok")
	return(auxr)


}

sp_kv_dk_kernel =  function(g1,g2,bin=T,k=0,verbose=F){

	if (bin){
		E(g1)$weight = rep(1,length(E(g1)$weight))
		E(g2)$weight = rep(1,length(E(g2)$weight))

	}
	#if (verbose) {
		#par(mfrow=c(1,2));plot(g1);plot(g2)
	#}
	#auxd1 = factor(round(distances(g1),k))
	#auxd2 = factor(round(distances(g2),k))
	#auxmax1 = max(levels(auxd1))
	#auxmax2 = max(levels(auxd2))
	#auxd1 = as.vector(round(distances(g1),k))
	#auxd2 = as.vector(round(distances(g2),k))
	auxd1 = round(distance_table(g1)$res,k)
	auxd2 = round(distance_table(g2)$res,k)
	if (verbose) {
		print(auxd1);print(auxd2)#;readline()
	}
	auxmax1 = length(auxd1)
	auxmax2 = length(auxd2)
	#print(paste("max1",auxmax1))
	#print(paste("max2",auxmax2))
	if (auxmax1<2){
		if (auxmax1==0){
			#auxv1 = c(0,0)
			return(0)
		} else {
			auxd1 = c(auxd1,0)
			auxmax1 = length(auxd1)
		}
	}
	if (auxmax2<2){
		if (auxmax2==0){
			#auxv2 = c(0,0)
			return(0)
		} else {
			auxd2 = c(auxd2,0)
			auxmax2 = length(auxd2)
		}
	}
	#print("a");print(auxd1);print(auxd2);readline()
	if (auxmax1>auxmax2){
		#auxmax = as.numeric(auxmax1)
		auxdiff = auxmax1-auxmax2
		auxd2 = c(auxd2,rep(0,auxdiff))
		auxr = dot(auxd1,auxd2)
		return(auxr)
	}#else{
	if (auxmax2>auxmax1){
		#auxmax = as.numeric(auxmax2)
		auxdiff = auxmax2-auxmax1
		auxd1 = c(auxd1,rep(0,auxdiff))
		auxr = dot(auxd1,auxd2)
		return(auxr)
	}
	#print("b");print(auxd1);print(auxd2);readline()
	auxr = dot(auxd1,auxd2)
	return(auxr)
	#if (verbose){
	#	print(auxmax);print(auxd1);print(auxd2)
	#}
	#auxv1 = tabulate(auxd1,nbins=auxmax+1)
	#auxv2 = tabulate(auxd2,nbins=auxmax+1)
	#if (no_zero){
	#	auxv1[1]=0
	#	auxv2[1]=0
	#}
	#if (verbose){
	#	print(auxd1);print(auxd2);readline()
	#}#;print("ok")
	#if (length(auxd1)>0){
	#	auxv1 = auxd1
	#}else{
	#	auxv1 = c(0,0)
	#}
	#if (length(auxd1)>0){
	#if (length(auxd2)>0){
	#	auxv2 = auxd2
	#}else{
	#	auxv2 = c(0,0)
	#}
	#print(paste("auxv1",auxv1))
	#print(paste("auxv2",auxv2))
	#print(auxv1);print(auxv2);readline()
	#auxr = dot(auxv1,auxv2)
	#print(paste("posdot",auxr));readline()
	#if (verbose){
#		print(auxd1);print(auxd2);print(auxr);readline()
	#}

	#return(auxr)

}
 

sp_kv_lk_kernel =  function(g1,g2,bin=T,k=0,verbose=F){

	if (bin){
		E(g1)$weight = rep(1,length(E(g1)$weight))
		E(g2)$weight = rep(1,length(E(g2)$weight))

	}
	#auxd1 = sum()

	#auxd1 = factor(round(distances(g1),k))
	#auxd2 = factor(round(distances(g2),k))
	#auxmax1 = max(levels(auxd1))
	#auxmax2 = max(levels(auxd2))
	#auxd1 = as.vector(round(distances(g1),k))
	#auxd2 = as.vector(round(distances(g2),k))
	auxd1 = round(distance_table(g1)$res,k)
	auxd2 = round(distance_table(g2)$res,k)
	auxmax1 = length(auxd1)
	auxmax2 = length(auxd2)
	if (auxmax1>auxmax2){
		#auxmax = as.numeric(auxmax1)
		auxdiff = auxmax1-auxmax2
		auxd2 = c(auxd2,rep(0,auxdiff))
	}#else{
	if (auxmax2>auxmax1){
		#auxmax = as.numeric(auxmax2)
		auxdiff = auxmax2-auxmax1
		auxd1 = c(auxd1,rep(0,auxdiff))
	}
	
	#if (verbose){
	#	print(auxmax);print(auxd1);print(auxd2)
	#}
	#auxv1 = tabulate(auxd1,nbins=auxmax+1)
	#auxv2 = tabulate(auxd2,nbins=auxmax+1)
	#if (no_zero){
	#	auxv1[1]=0
	#	auxv2[1]=0
	#}
	if (verbose){
		print(auxd1);print(auxd2)#;readline()
	}#;print("ok")
	if (length(auxd1)>0){
		auxv1 = auxd1
	}else{
		auxv1 = c(0,0)
	}
	if (length(auxd1)>0){
		auxv2 = auxd2
	}else{
		auxv2 = c(0,0)
	}
	#print(auxv1);print(auxv2);readline()
	auxr = dot(auxv1,auxv2)

	return(auxr)

}



sp_kv_dk_kernel_old_old =  function(g1,g2,bin=T,k=0,verbose=T){

	if (bin){
		E(g1)$weight = rep(1,length(E(g1)$weight))
		E(g2)$weight = rep(1,length(E(g2)$weight))

	}
	if (verbose) par(mfrow=c(1,2));plot(g1);plot(g2)
	#auxd1 = factor(round(distances(g1),k))
	#auxd2 = factor(round(distances(g2),k))
	#auxmax1 = max(levels(auxd1))
	#auxmax2 = max(levels(auxd2))
	#auxd1 = as.vector(round(distances(g1),k))
	#auxd2 = as.vector(round(distances(g2),k))
	auxd1 = round(distance_table(g1)$res,k)
	auxd2 = round(distance_table(g2)$res,k)
	if (verbose) {
		print(auxd1);print(auxd2)#;readline()
	}
	auxmax1 = length(auxd1)
	auxmax2 = length(auxd2)
	print(paste("max1",auxmax1))
	print(paste("max2",auxmax2))
	if (auxmax1>auxmax2){
		#auxmax = as.numeric(auxmax1)
		auxdiff = auxmax1-auxmax2
		auxd2 = c(auxd2,rep(0,auxdiff))
	}#else{
	if (auxmax2>auxmax1){
		#auxmax = as.numeric(auxmax2)
		auxdiff = auxmax2-auxmax1
		auxd1 = c(auxd1,rep(0,auxdiff))
	}
	
	#if (verbose){
	#	print(auxmax);print(auxd1);print(auxd2)
	#}
	#auxv1 = tabulate(auxd1,nbins=auxmax+1)
	#auxv2 = tabulate(auxd2,nbins=auxmax+1)
	#if (no_zero){
	#	auxv1[1]=0
	#	auxv2[1]=0
	#}
	#if (verbose){
		print(auxd1);print(auxd2)#;readline()
	#}#;print("ok")
	if (length(auxd1)>0){
		auxv1 = auxd1
	}else{
		auxv1 = c(0,0)
	}
	#if (length(auxd1)>0){
	if (length(auxd2)>0){
		auxv2 = auxd2
	}else{
		auxv2 = c(0,0)
	}
	print(paste("auxv1",auxv1))
	print(paste("auxv2",auxv2))
	#print(auxv1);print(auxv2);readline()
	auxr = dot(auxv1,auxv2)
	print(paste("posdot",auxr))#;readline()
	if (verbose){
		print(auxd1);print(auxd2);print(auxr)#;readline()
	}

	return(auxr)

}

sp_kv_kernel_old =  function(g1,g2,bin=T,k=0,verbose=T){

	if (bin){
		E(g1)$weight = rep(1,length(E(g1)$weight))
		E(g2)$weight = rep(1,length(E(g2)$weight))

	}
	#auxd1 = factor(round(distances(g1),k))
	#auxd2 = factor(round(distances(g2),k))
	#auxmax1 = max(levels(auxd1))
	#auxmax2 = max(levels(auxd2))
	auxd1 = as.vector(round(distances(g1),k))
	auxd2 = as.vector(round(distances(g2),k))
	auxmax1 = max(auxd1)
	auxmax2 = max(auxd2)
	if (auxmax1>=auxmax2){
		auxmax = as.numeric(auxmax1)
	}else{
		auxmax = as.numeric(auxmax2)
	}
	if (verbose){
		print(auxmax);print(auxd1);print(auxd2)
	}
	auxv1 = tabulate(auxd1,nbins=auxmax+1)
	auxv2 = tabulate(auxd2,nbins=auxmax+1)
	#if (no_zero){
	#	auxv1[1]=0
	#	auxv2[1]=0
	#}
	if (verbose){
		print(auxmax);print(auxd1);print(auxd2);print(auxv1);print(auxv2)#;readline()
	}
	auxr = dot(auxv1,auxv2)

	return(auxr)

}

shortest_path_kv_kernel_old =  function(g1,g2,bin=T,k=0){

	if (bin){
		E(g1)$weight = rep(1,length(E(g1)$weight))
		E(g2)$weight = rep(1,length(E(g2)$weight))

	}
	auxd1 = round(distances(g1),k)
	auxd2 = round(distances(g2),k)
	#print(as.vectorauxd1);print(auxd2);readline()
	tamg1 = dim(auxd1)[1]	
	#tamg1 = 1dim(auxd2)[1]
	auxr = 0
	for (i in 1:tamg1){
		for (j in 1:tamg1){
			aux1 = auxd1[i,j]
			auxs = sum(aux1==auxd2)
			#print(aux1);print(auxs);readline()
			auxr = auxr + auxs
		}
	}
	auxr=auxr-tamg1*tamg1
	return(auxr)
}

# Adiciona matriz de p-values binarizados conforme lim
add_bin_stat_input = function(db,pdbids,subtext=c("Bin="," - Diag="),diagz=1,lim=0.05){

	#lim=0.05
	#dks[[2]]$A=list()
	#dks[[2]]$A$a = apply(dks[[1]]$A$a>lim,1,as.numeric)
	#dbtb5 = apply(dbt>lim,1,as.numeric)
	#rownames(dks[[2]]$A$a)=pdbids
	tam = length(db)
	aux=list()
	aux$A=list()
	aux$A$a = apply(db[[tam]]$A$a>lim,1,as.numeric)
	if (diagz) {
		diag(aux$A$a)=0
		aux$type=paste(subtext[1],lim,subtext[2],as.numeric(!diagz),sep="")
	} else {
		aux$type=paste(subtext[1],lim,subtext[2],as.numeric(!diagz),sep="")
	}
	rownames(aux$A$a)=pdbids
	tam = length(db)
	db[[tam+1]]=aux
	return(db)
}

# Centraliza matriz estatística já presente em db
add_center_stat_input = function(db,subtext=" - Centered"){


	tam = length(db)
	j = tam + 1
	for (i in 1:tam){
		aux = list()
		aux$A = list()
		aux$A$a=db[[i]]$A$a-mean(db[[i]]$A$a)
		aux$type=paste(db[[i]]$type,subtext,sep="")
		db[[j]]=aux
		j = j + 1
	}
	return(db)
}


# transforma uma matrix estatistica (ex: ks.test()) numa matrix centralizada com diagonal zero(se o caso)
# OBSOLETED
make_stat_matrix_OLD = function(db,center,inter,p=3,ttype="Diag=0"){
	aux=list()
	aux$A = list()
	if (center) aux$type = "Centered"
	else aux$type = "No Centered"
	if (inter) aux$type = paste(aux$type,ttype,sep="-")
	#diag(db)=0
	if (inter==0) 
		if (center) db = db - mean(db)
	while (inter){
		if (center) db = db - mean(db) #centrifica
		diag(db)=0
		inter = inter-1
	}
	aux$A$a = round(db,p)
	return(aux)
}

#Checa simetria de m
symmetric_check = function(m,pdbname="",typename="",verbose=FALSE){

	if (sum(apply(m-t(m),1,Norm))){
		if (verbose) print(paste("Warning: adjacency matriz of",pdbname,"by",typename,"is not symmetric"))
		return(FALSE)
	} else {
		if (verbose) print(paste("Good: adjacency matriz of",pdbname,"by",typename,"is symmetric"))
		return(TRUE)
	}
}

# OBSOLOTED: Codifica numa matriz as relações entre o tipo de resíduo
codify_matrix = function(type){


	if (type == 1){
		aux = matrix(c(1,2,3,
				2,4,5,
				3,5,6),nrow=3,ncol=3)
		rownames(aux) = c("a","c","p")
		colnames(aux) = c("a","c","p")
	}
	if (type == 2){
		aux = matrix(c(1,1,2,
				1,1,2,
				2,2,3),nrow=3,ncol=3)
		rownames(aux) = c("a","c","p")
		colnames(aux) = c("a","c","p")

	}
	return(aux)
}

# Idenfica res/atom com dist == 0 e area != 0
verify_dist_area_anomaly = function(dbarea,dbdist,noredun=FALSE){
# Não funciona ainda para noredun=TRUE

	tami = dim(dbarea$A$a)[1]
	tamj = dim(dbarea$A$a)[2]
	aux = c()
	#if (noredun){
	#	x=dbarea$ide=="i"
	#	y=dbarea$ide=="e"
	#} else {
	#	x = rep(TRUE,tamj)
	#	y = rep(TRUE,tami)
	#}
	#ma = dbarea$A$a[x,y]
	#md = dbdist$A$a[x,y]
	#tamai = dim(ma)[1]
	#tamaj = dim(ma)[2]

	for (i in 1:tami){
		for (j in 1:tamj){
			if ((dbdist$A$a[i,j]==0)&(dbarea$A$a[i,j]!=0))
				aux = rbind(aux,c(dbdist$name[i],dbdist$idn[i],"-",dbdist$name[j],dbdist$idn[j],dbdist$A$a[i,j],dbarea$A$a[i,j]))
		}

	}
	return(aux)
}



# OBSOLETED: Retorna um vetor conforme a codificação de resíduos em type
codify_residue_type = function(ve,vi,type){

	tame = length(ve)
	tami = length(vi)
	aux = c()
	codm = codify_matrix(type)

	for (e in 1:tame){
		for (i in 1:tami){
				aux = c(aux,codm[ve[e],vi[i]])
		}
	}
	return(aux)

}

# Imprime o numero de componentes conexos em db
connected_component_number = function(db,p){

	tami = dim(db$Lw$eig$v)[1]
	#tamj = dim(db$Lw$a)[2]

	aux = sum(round(diag(db$Lw$eig$s),p)==0)
	auxtab = table(as.factor(round(db$Lw$eig$v[,tami],p)))
	auxvalues = as.numeric(names(auxtab))
	print(paste(db$pdbname,"weigthed by",db$typename,":",aux))
	print(as.vector(auxtab))
	print(paste("(",paste(auxvalues),")"))
	#print(paste("comp #:",paste(auxtab),"comp values:",paste(auxvalues)))
	#print("")

}


# Particiona o range de pesos entre min e max em n partições 
make_edge_partition = function(minmaxn){

	rangemax = ceiling(minmaxn[2])
	rangemin = floor(minmaxn[1])
	range = rangemax - rangemin
	diff = ceiling(range/minmaxn[3])
	aux = c(rangemax)

	for (i in 1:minmaxn[3]){
		auxdiff = aux[i]-diff
		if ((auxdiff)<=0){
			 aux = c(aux,0.1) # a menor partição não pode ser zero
			 break
		} else {
			aux = c(aux,auxdiff)
		}
	}
	return(aux)

}


# Decompõe vetor indicador vind em vetores lógicos conforme número de categorias (levels)
make_vector_indicator_list = function(vind,verbose=FALSE){

	n = length(levels(as.factor(vind)))
	mlogic = diag(rep(1,n))>0

	tami = dim(mlogic)[1]

	vlist = list()


	for (i in 1:tami){

		vf = as.factor(vind)
		levels(vf) = mlogic[i,]
		vlist[[i]] = as.logical(vf)

	
	}
	return(vlist)
	
}

#Recebe duas listas de vetores logicos e faz o produto cartesiano deles
make_logical_vector_product = function(v1,v2){

	aux = list()
	k = 1	

	tam1 = length(v1)
	tam2 = length(v2)

	for (i in 1:tam2){
		for (j in 1:tam1){

			aux[[k]] = v1[[j]] & v2[[i]];k=k+1
		}
	}
	return(aux)

}

#OBSOLETED: Classifica cadeias em db conforme códigos em chainids
classify_chains = function(db,chainids,pdbids,chain_1="E",chain_2="I"){

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

view_inibitory_loop_as_char = function(clus,name,clusep="|",resep=":"){

	tam = 1:length(clus)
	auxclus=c()
	aux=c()
	auxname = extract_resn(name)
	#print(tam)
	auxclus=c(clusep)
	for (i in tam){
		auxid = clus[[i]]$loop
		auxclus = c(auxclus,auxname[auxid],clusep)
		#auxclus=c("|")
	}
	auxclus=strcat(auxclus)
	return(auxclus)
}

#Mostra os códigos de cadeias registradas em db
show_chains = function(db){

	tam = length(db)
	for (i in 1:tam){
		aux = c(db[[i]]$pdbname,levels(as.factor(db[[i]]$ide)))
		aux = as.vector(aux)
		print(aux)
	}
	print("-----")
}

# Prepara vetor de string para compor nomes dos vértices no Pymol
make_vertex_sufix = function(n,base){

	aux = c()

	for (i in 1:n){
		aux = c(aux,paste(base[1],i,sep=""),paste(base[2],i,sep=""))
	}
	return(aux)

}

#Define c(xlim,ylim) conforme as coordenadas dos centros geometricos
find_minmax_geomxy_in_clus = function(clus){

	tam = 1:length(clus)
	aux = c()
	for (i in tam){
		aux = rbind(aux,clus[[i]]$geomc)
	}
	auxminmax = c(min(aux[,1]),max(aux[,1]),min(aux[,2]),max(aux[,2]))
	return(auxminmax)
}


# transforma uma lista de atomos em lista de residuos (apenas diferentes)
get_res_from_atom_list = function(atoml,name){

	aux = c()

	auxname = name[atoml]
	#print(auxname)
	auxresi = apply(as.matrix(auxname),1,extract_resi)
	auxresn = apply(as.matrix(auxname),1,extract_resn)
	auxresu = unique(auxresi)
	
	for (i in auxresu){
	
		filter = auxresi == i
		#print(unique(auxresn[filter]))
		#readline()
		aux = c(aux,unique(auxresn[filter]))
	}
	
	#print(auxresn)
	return(aux)

}

make_res_table = function(atoms,xyz,area,no_lig_all=T){
  
  #browser()
  atom.tab = tibble(atoms=atoms,xyz,area)
  atom.tab = atom.tab %>% separate(atoms,c("chain","residue","atom"))
  #browser()
  res.tab = atom.tab %>% select(chain,residue,x,y,z,area) %>% group_by(chain,residue) %>%
            summarise(x=mean(x),y=mean(y),z=mean(z),area=sum(area)) %>% ungroup()
  
  if (no_lig_all){
    res.tab = res.tab %>% filter(chain == "A")
  }
  
  res.lig = atom.tab %>% select(chain,residue,x,y,z,area) %>% filter(chain != "A")
  res.tab = bind_rows(res.tab,res.lig)
  #browser()
  # res.tab = res.tab %>% separate(residue,c("resn","resi"),sep="(?<=[A-Za-z])(?=[0-9])",
  #                                remove=F) ### Falha em U5G
  res.tab = res.tab %>% extract(residue,c("resn","resi"),"(^.{0,2}\\D)(\\d*$)",
                                remove=F)
 
  
  
  res.tab = res.tab %>% mutate(across(where(is.double),~round(.x, digits=3)))
  res.tab = res.tab %>% mutate(resc=aa321(resn)) 
  res.tab = res.tab %>% mutate(resC = case_when(
                               resc == "X" ~ resn,
                               TRUE ~ resc)
                               )
  res.tab = res.tab %>% mutate(resCi = paste0(resC,resi))
  res.tab = res.tab %>% relocate(any_of(c("resc","resC","resCi")), .after = resn)
  
  
  
  
  return(res.tab)
  
  
}

add_res_information_to_clusters = function(db){
  
  atoms = db$element_name
  xyz = db$pdb$xyz
  #browser()
  ngroup = length(db$exp$pam)
  #ngroup = 1
  db$exp$residues = list()
  for (i in 1:ngroup){
    db$exp$residues[[i]] = list()
    for (j in 1:i){
      db$exp$residues[[i]][[j]] = list()
      id.j = grep(j,db$exp$pam[[i]]$clustering)
      atoms.j = atoms[id.j]
      xyz.j = xyz[id.j,,drop=F]
      area.j = apply(db$A$a[id.j,],1,sum)
      res.table = make_res_table(atoms.j,xyz.j,area.j)
      #browser()
      db$exp$residues[[i]][[j]]$res.id = id.j
      db$exp$residues[[i]][[j]]$res.table = res.table
      #residues[[i]][[j]]$res.cod = sapply(atoms.j,splitatom,typename=3)[2,]
      #residues[[i]][[j]]$res.id = id.j
      #residues[[i]][[j]]$gc = geometric_center(xyz.i)
    }
  }
  return(db)
  
}


# OBSOLETED: transforma uma lista de atomos em lista de residuos (diferentes)
get_res_from_atom_list_OLD = function(atoml,name){

	auxname = name[atoml]
	#print(auxname)
	#auxresi = apply(as.matrix(auxname),1,extract_resi)
	auxresn = apply(as.matrix(auxname),1,extract_resn)
	#print(auxresn)
	return(unique(auxresn))

}


res_text = function(mid,res,r,text.cex,angt=360){

	if (res[1]!=0){

		tam = length(res)
		#ang = 0
		ang = angt/tam

		r = c(mid[1]+r,mid[2])

		for (i in 1:tam){
			
			text(r[1],r[2],res[i],cex=text.cex)
			r = rotatexy(r,ang,mid=mid)
			
		}

	}

} 

renumerate_res_table = function(db,adj=-24){
  
  
  ma = db %>% filter(chain == "A", resn != "HOH") 
  ma$resi = as.character(as.numeric(ma$resi)+adj)
  ma = ma %>% mutate(resCi = paste0(resC,resi),
                     residue = paste0(resn,resi)
                     )
  mb = db %>% filter(chain == "B" | resn == "HOH")
  mab = bind_rows(ma,mb)
  #browser()
  #View(mab)
  
  return(mab)
  
  #browser()
}

consolidate_res_table = function(db,seq.id,noid=NULL){
  
  m = tibble()
  #browser()
  for (i in seq.id){
    aux = db[[1]][[i]]$exp$residues[[1]][[1]]$res.table
    if (!is.null(noid)){
      f = aux$resn %in% noid
      aux$resCi[f] = aux$resn[f]
      aux$resi[f] = " "
      aux$residue[f] = aux$resn[f]
    }
    m = bind_rows(m,aux)
    
  }
  #res.tab = atom.tab %>% select(chain,residue,x,y,z,area) %>% group_by(chain,residue) %>%
  #  summarise(x=mean(x),y=mean(y),z=mean(z),area=sum(area)) %>% ungroup()
  
  ma = m %>% filter(chain == "A") %>% group_by (chain,residue,resn,resc,resC,resCi,resi) %>%
       summarise(x=mean(x),y=mean(y),z=mean(z),area=mean(area)) %>% ungroup()
  mb = m %>% filter(chain == "B")
  
  mab = bind_rows(ma,mb)
  
  return(mab)
  #browser()
  
  
}

generate_spheres_mesh = function(m,n=24,r=1){
  
  xyz = as_tibble(fibonaccisphere(n=n,r=r))
  m.out = tibble()
  
  for (i in 1:dim(m)[1]){
    
    xyz.i = xyz %>% mutate(x=x+m[i,1]$x,y=y+m[i,2]$y,z=z+m[i,3]$z)
    m.out = bind_rows(m.out,xyz.i)
    
  }
  return(m.out)
  #browser()
  
}

plotly_all_together_ligs = function(p1,p2=NULL,p3=NULL,text="",
                                    range=c(),nticks=4,p.sizes=c(1200,600)){
  #browser()
  axx <- list(
    nticks = nticks,
    #range = c(2, 18)
    range = range[1:2]
  )
  axy <- list(
    nticks = nticks,
    #range = c(90, 104)
    range = range[3:4]
  )
  axz <- list(
    nticks = nticks,
    #range = c(4, 22)
    range = range[5:6]
  )
  f = list(
    family = "Roboto Condensed",
    size = 18,
    color = "black")
  
  t1 = list(
    #text = "Crystal iNOScs Sapropterin (A/B)",
    text = text[1],
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
  
  p1 = p1 %>% layout(
    scene1 = list(xaxis = axx,yaxis = axy,zaxis = axz, aspectmode = "cube"),
    annotations = t1 
  )
  
  if(!is.null(p2)){
    t2 = list(
      #text = "Docked iNOScs Sapropterin poses",
      text = text[2],
      font = f,
      xref = "paper",
      yref = "paper",
      yanchor = "bottom",
      xanchor = "center",
      align = "center",
      x = 0.5,
      y = 0.95,
      showarrow = FALSE
    )
    p2 = p2 %>% layout(
      scene2 = list(xaxis = axx,yaxis = axy,zaxis = axz, aspectmode = "cube"),
      annotations = t2 
    )
    
  }
  if(!is.null(p3)){
    t3 = list(
      #text = "Top docked iNOScs poses",
      text = text[3],
      font = f,
      xref = "paper",
      yref = "paper",
      yanchor = "bottom",
      xanchor = "center",
      align = "center",
      x = 0.5,
      y = 0.95,
      showarrow = FALSE
    )
    p3 = p3 %>% layout(
      scene3 = list(xaxis = axx,yaxis = axy,zaxis = axz, aspectmode = "cube"),
      annotations = t3
    )
  }
  #browser()
  if (is.null(p2)&is.null(p3)){
    p = p1
  } 
  if (!is.null(p2)&is.null(p3)){
    p = subplot(p1,p2) 
  } 
  if (!is.null(p2)&!is.null(p3)){
    p = subplot(p1,p2,p3)  
  }
  p$x$layout <- p$x$layout[grep('NA', names(p$x$layout), invert = TRUE)]
  
  return(p)
  
}

#https://plotly.com/r/3d-axes/
plotly_all_together_INOSas_ligs = function(p1,p2,p3){
  
    axx <- list(
      nticks = 4,
      range = c(-8, 8)
    )
    axy <- list(
      nticks = 4,
      range = c(90, 110)
    )
    axz <- list(
      nticks = 4,
      range = c(10, 30)
    )
    f = list(
      family = "Roboto Condensed",
      size = 18,
      color = "black")
    
    t1 = list(
      text = "iNOSas ARG ITU (A/B)",
      font = f,
      xref = "paper",
      yref = "paper",
      yanchor = "bottom",
      xanchor = "center",
      align = "center",
      x = 0.5,
      y = 0.95,
      showarrow = FALSE
    )
    
    t2 = list(
      text = "iNOSas bis-isothiourea poses",
      font = f,
      xref = "paper",
      yref = "paper",
      yanchor = "bottom",
      xanchor = "center",
      align = "center",
      x = 0.5,
      y = 0.95,
      showarrow = FALSE
    )
    
    t3 = list(
      text = "Top iNOSas poses",
      font = f,
      xref = "paper",
      yref = "paper",
      yanchor = "bottom",
      xanchor = "center",
      align = "center",
      x = 0.5,
      y = 0.95,
      showarrow = FALSE
    )
    
    
    p1 = p1 %>% layout(
      scene1 = list(xaxis = axx,yaxis = axy,zaxis = axz, aspectmode = "data"),
      annotations = t1 
    )
    
    p2 = p2 %>% layout(
      scene2 = list(xaxis = axx,yaxis = axy,zaxis = axz, aspectmode = "data"),
      annotations = t2 
    )
    p3 = p3 %>% layout(
      scene3 = list(xaxis = axx,yaxis = axy,zaxis = axz, aspectmode = "data"),
      annotations = t3
    )
    
  #browser()
  p = subplot(p1,p2,p3)
  #p$x$layout <- p$x$layout[grep('NA', names(p$x$layout), invert = TRUE)]
  
  return(p)
  
}

plotly_all_together_INOSas = function(p1,p2,p3,p4,p5,p6){
  
  if(1){
  axx <- list(
    nticks = 4,
    range = c(-6, 6)
  )
  axy <- list(
    nticks = 4,
    range = c(90, 110)
  )
  axz <- list(
    nticks = 4,
    range = c(10, 30)
  )
  f = list(
    family = "Roboto Condensed",
    size = 18,
    color = "black")
  
  t1 = list(
    text = "iNOS ARG (A/B)",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
  
  t2 = list(
    text = "iNOS ARG ITU (A/B)",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
  
  t3 = list(
    text = "iNOS bis-isothiourea poses",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
  
  t4 = list(
    text = "Top iNOSas COX1",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
  
  t5 = list(
    text = "Top iNOSas iNOScs",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
  
  t6 = list(
    text = "Top just iNOSas",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
  
  p1 = p1 %>% layout(
    scene1 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t1 
  )
  
  p2 = p2 %>% layout(
    scene2 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t2 
  )
  p3 = p3 %>% layout(
    scene3 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t3
  )
  
  p4 = p4 %>% layout(
    scene4 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t4
  )
  
  p5 = p5 %>% layout(
    scene5 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t5
  )
  
  p6 = p6 %>% layout(
    scene6 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t6
  )
  }
  #browser()
  p = subplot(p1,p2,p3,p4,p5,p6,nrows=2)
  #p$x$layout <- p$x$layout[grep('NA', names(p$x$layout), invert = TRUE)]
  
  return(p)
  
}


plotly_all_together_COX1 = function(p0,p1,p2){
  
  axx <- list(
    nticks = 5,
    range = c(240, 260)
  )
  axy <- list(
    nticks = 5,
    range = c(80, 120)
  )
  axz <- list(
    nticks = 5,
    range = c(-48, -28)
  )
  f = list(
    family = "Roboto Condensed",
    size = 18,
    color = "black")

  t0 = list(
    text = "ARACHIDONIC ACID CRYSTAL (COX1)",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
    
  t1 = list(
    text = "ARACHIDONIC ACID POSES (COX1)",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
  
  t2 = list(
    text = "AROEIRA LIGAND POSES (COX1)",
    font = f,
    xref = "paper",
    yref = "paper",
    yanchor = "bottom",
    xanchor = "center",
    align = "center",
    x = 0.5,
    y = 0.95,
    showarrow = FALSE
  )
 
  p0 = p0 %>% layout(
    scene1 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t0 
  )
   
  p1 = p1 %>% layout(
    scene2 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t1 
    )
  p2 = p2 %>% layout(
    scene3 = list(xaxis = axx,yaxis = axy,zaxis = axz),
    annotations = t2
  )
  p = subplot(p0,p1,p2)
  p$x$layout <- p$x$layout[grep('NA', names(p$x$layout), invert = TRUE)]

  return(p)
  
}

plot_lig_res_interactions = function(m,scene,area.lim=500,alpha=0.01,n=240,ntext=40){
  
  # m1 = m %>% filter(area>area.lim | chain == "B")
  # 
  # temp = m1$area
  # 
  # f = m1$chain == "A"
  # tempi = temp
  # tempi[f] = 0
  # 
  # m1 = m1 %>% mutate(resCi = case_when(
  #   chain == "B" ~ " ",
  #   TRUE ~ resCi)
  # )
  
  #browser()
  
  ma = m %>% filter(area>area.lim & chain == "A")
  mb = m %>% filter(chain == "B") #%>% mutate(resCi = " ")
  
  mb.mesh = mb %>% select(x,y,z)
  mb.mesh = generate_spheres_mesh(mb.mesh,n=n)
  
  #browser()
  # pb = mb.mesh %>% plot_ly(x= ~x, y= ~y, z = ~z, type = "mesh3d", 
  #                          alphahull=6, opacity=0.1,flatshading=F) #%>%
  #browser()
  
  p = plot_ly(scene=scene)
  p = p %>% add_trace(
    data = mb.mesh,
    x = ~x,
    y = ~y,
    z = ~z,
    type = "scatter3d", 
    mode = 'markers', 
    marker = list(size = 0.5),
    opacity = 0.1,
    hoverinfo = "skip"
  )
  
  p = p %>% add_text(
    data = ma,
    x = ~x,
    y = ~y,
    z = ~z,
    text= ~resCi,
    textposition = 'middle',
    hoverinfo = "text",
    textfont = list(family = "Roboto Condensed", size = ~area/ntext)#
  ) %>%
    hide_legend()
  
  return(p)
  
  # axy <- list(
  #   nticks = 4,
  #   range = c(100, 130)
  # )
  # p <- p %>% layout(
  #   scene = list(xaxis=axy,yaxis=axy,zaxis=axy)
  # )
  # p <- p %>% layout(
  #   scene = list(xaxis = list(range = c(240, 260)),
  #                yaxis = list(range = c(100, 130)),
  #                zaxis = list(range = c(-30, -50))
  #               )
  #)

  
  #p$x$layout <- p$x$layout[grep('NA', names(p$x$layout), invert = TRUE)]

  
  #browser()
  
  # if(0){
  # pb = mb.mesh %>% plot_ly(x= ~x, y= ~y, z = ~z) %>%
  #   add_trace(
  #     type = "scatter3d", 
  #     mode = 'markers', 
  #     marker = list(size = 0.5),
  #     opacity = 0.1,
  #     hoverinfo = "skip"
  #   )
  # 
  # pa = ma  %>% plot_ly(x= ~x, y= ~y, z = ~z) %>%
  #     add_text(
  #       text= ~resCi,
  #       textposition = 'middle',
  #       hoverinfo = "text",
  #       textfont = list(family = "Roboto Condensed", size = ~area/40)#
  #     ) %>%
  #     hide_legend()
  # 
  # p = subplot(pa,pb)
  # }
  # #browser()
  # 
  # if(0){
  # p = m1 %>% plot_ly(x= ~x, y= ~y, z = ~z) %>%
  #   add_mesh(
  #     type = "mesh3d",
  #     alphanull = 0,
  #     opacity = 0.2
  #   ) %>%
  #   add_markers(
  #     size = tempi,
  #     sizes = c(0,100),
  #     alpha = 0.4,
  #     opacity = 0.8
  #   ) %>% 
  #   add_text(
  #     text= ~resCi,
  #     textposition = 'middle',
  #     hoverinfo = "text",
  #     textfont = list(family = "Roboto Condensed", size = temp/40)#
  #   ) %>%
  #   hide_legend()
  # }
  # #browser()

  
}

#plot_new_model_clus = function(db,clusid,main,xlab,ylab,xlim,ylim,text.marg=NULL,add=0,lcol="BLACK",col="BLACK",dlim,dres,ri,text.cex,mtext.cex){
plot_new_model_clus = function(db,clusid,main,xlab,ylab,xlim,ylim,text.marg=NULL,add=0,lcol="GRAY80",col="BLACK",text.cex=0.5,mtext.cex=0.5,rot=NULL,r=2,dr=0.8){
	#print(clus)
	#readline()
	#if (is.null(xlim)){
	#	auxminmax = find_minmax_geomxy_in_clus(clus)
	#	xlim=c(auxminmax[1]-dlim[1],auxminmax[1]+dlim[1])
	#	ylim=c(auxminmax[3]-dlim[2],auxminmax[4]+dlim[2])
	#}
	if (add==0) plot(0,0,xlim=xlim,ylim=ylim,type="n",asp=1,main=paste(main,db$pdbname,sep=" - "),xlab=xlab,ylab=ylab)
	
	if (!is.null(text.marg)) mtext(text.marg,3,0.1,cex=mtext.cex)
	
	#if (is.null(clusvis)) tam = 1:length(clus)
	#else tam = clusvis
	#lenri = length(ri)

	clus = db$exp$superclus[[clusid]]
	cluslist = db$exp$pam[[clusid]]$clustering

	auxchain = factor(substr(db$name,1,1))

	auxlevel = levels(auxchain)

	if (length(auxlevel)>2) print(paste("WARNING: more than 2 chains found in",db$pdbname))

	auxres = substr(db$name,2,2)
	auxv1 = auxlevel[1]==auxchain
	auxv2 = auxlevel[2]==auxchain

	if (sum(auxv1)<sum(auxv2)){
		aux = auxv2
		auxv2 = auxv1
		auxv1 = aux
	}

	tam = 1:clusid

	auxpdbxyz = db$pdb$xyz
	auxgeomc = clus$geomc

	if (!is.null(rot)){
		auxpdbxyz = transform_by_rot_matrix(auxpdbxyz,matrix(rot,ncol=4,byrow=T))
		auxgeomc = transform_by_rot_matrix(clus$geomc,matrix(rot,ncol=4,byrow=T))
	}

	for (i in tam){

		#xgc = clus[[i]]$geomc[1]
		#ygc = clus[[i]]$geomc[2]
		x = auxgeomc[i,][1]
		y = auxgeomc[i,][2]

		auxvr = cluslist == i

		auxvr1 = auxvr & auxv1
		auxvr2 = auxvr & auxv2
		#lenres = length(clus[[i]]$loopres)
		#print(ri[lenres])
		#print(x);print(y)
		#readline()
		#if (lenres > lenri) lenres = lenri
		#plotcircle(mid=c(x,y),col=col,lcol=lcol,r=ri[lenres])
		plotcircle(mid=c(x,y),lcol=lcol,r=r)
		#auxres = strcat(auxres[auxvr2])
		auxres2 = get_res_from_atom_list(auxvr2,db$name)
		auxres1 = get_res_from_atom_list(auxvr1,db$name)
		#print(auxres2);readline()
		auxyz = auxpdbxyz[auxvr,]

		#print(auxres);readline()
		#text(x,y,strcat(clus[[i]]$loopres),cex=text.cex,font=2)
		text(x,y,strcat(auxres2),cex=text.cex,font=1)#;readline()
		#res_text(c(xgc,ygc),clus[[i]]$baseres,ri[lenres]+dres,text.cex)
		res_text(c(x,y),auxres1,r+dr,text.cex)
		auxellips = ellipsoidhull(auxyz[,1:2])
		lines(predict(auxellips),lty=3,col=lcol,lwd=2)

	}
	
}

# Constroi um modelo visual da clusterização baseado no centro geométrico dos clusteres
plot_model_clus = function(clus,clusvis,pdbname,main,xlab,ylab,xlim,ylim,text.marg,add,lcol,col,dlim,dres,ri,text.cex,mtext.cex){

	#print(clus)
	#readline()
	if (is.null(xlim)){
		auxminmax = find_minmax_geomxy_in_clus(clus)
		xlim=c(auxminmax[1]-dlim[1],auxminmax[1]+dlim[1])
		ylim=c(auxminmax[3]-dlim[2],auxminmax[4]+dlim[2])
	}
	if (add==0) plot(0,0,xlim=xlim,ylim=ylim,type="n",asp=1,main=paste(main,pdbname,sep=" - "),xlab=xlab,ylab=ylab)
	
	if (!is.null(text.marg)) mtext(text.marg,3,0.1,cex=mtext.cex)
	
	if (is.null(clusvis)) tam = 1:length(clus)
	else tam = clusvis
	lenri = length(ri)
	
	for (i in tam){

		xgc = clus[[i]]$geomc[1]
		ygc = clus[[i]]$geomc[2]
		x = xgc
		y = ygc
		lenres = length(clus[[i]]$loopres)
		#print(ri[lenres])
		#print(x);print(y)
		#readline()
		if (lenres > lenri) lenres = lenri
		#plotcircle(mid=c(x,y),col=col,lcol=lcol,r=ri[lenres])
		plotcircle(mid=c(x,y),lcol=lcol,r=ri[lenres])
		text(x,y,strcat(clus[[i]]$loopres),cex=text.cex,font=2)
		res_text(c(xgc,ygc),clus[[i]]$baseres,ri[lenres]+dres,text.cex)
		lines(predict(clus[[i]]$ellips),lty=3,col=lcol,lwd=2)

	}
	
}

# Constroi um modelo visual da clusterização baseado no centro geométrico dos resíduos
plot_model_clus2 = function(clus,clusvis,pdbname,main,xlab,ylab,xlim,ylim,text.marg,add,lcol,col,dlim,dres,ri,text.cex,mtext.cex,restype,lwd,lty,font){

	#print(clus)
	#readline()
	if (is.null(xlim)){
		auxminmax = find_minmax_geomxy_in_clus(clus)
		xlim=c(auxminmax[1]-dlim[1],auxminmax[1]+dlim[1])
		ylim=c(auxminmax[3]-dlim[2],auxminmax[4]+dlim[2])
	}
	if (!is.null(pdbname)) aux.main = paste(main,pdbname,sep=" - ")
	else aux.main = main
	if (add==0) plot(0,0,xlim=xlim,ylim=ylim,type="n",asp=1,main=aux.main,xlab=xlab,ylab=ylab)
	
	if (!is.null(text.marg)) mtext(text.marg,3,0.1,cex=mtext.cex)
	
	if (is.null(clusvis)) tami = 1:length(clus)
	else tami = clusvis
	#lenri = length(ri)
	#print(restype);readline()
	for (i in tami){
		lines(predict(clus[[i]]$ellips),lty=lty,col=lcol,lwd=lwd)
		if (restype[1]){
			tamj = 1:length(clus[[i]]$loopresn)
			for (j in tamj){
				if (!is.null(dim(clus[[i]]$loopresgc))){
					x = clus[[i]]$loopresgc[j,1]
					y = clus[[i]]$loopresgc[j,2]
					text(x,y,clus[[i]]$loopresn[j],cex=text.cex,font=font[1],col=col[1])
				} else {
					print(paste("WARNING: ",pdbname,"has no geometric center for clus",i,"in loopresgc")) 
				
				}
			}
		}
		if (restype[2]){
			tamj = 1:length(clus[[i]]$baseresn)
			for (j in tamj){
				if (!is.null(dim(clus[[i]]$baseresgc))){
					x = clus[[i]]$baseresgc[j,1]
					y = clus[[i]]$baseresgc[j,2]
					text(x,y,clus[[i]]$baseresn[j],cex=text.cex,font=font[2],col=col[2])
				}
			}
		}
		if (restype[3]){
			tamj = 1:length(clus[[i]]$noloopresn)
			for (j in tamj){
				if (!is.null(dim(clus[[i]]$noloopresgc))){
					x = clus[[i]]$noloopresgc[j,1]
					y = clus[[i]]$noloopresgc[j,2]
					text(x,y,clus[[i]]$noloopresn[j],cex=text.cex,font=font[3],col=col[3])
				}
			}
		}
		#lines(predict(clus[[i]]$ellips),lty=3,col=lcol,lwd=2)
	}
	
}

multi_plot_model_clus = function(db,krange,mfrow,main,xlab,ylab,xlim,ylim,type=1,clusvis=c(),restype=c(1,1,1),text.marg=c("Avg.width = ","k = "),add=0,lcol="GRAY60",col=c("WHITE","BLUE","RED","SLATEBLUE"), dlim=c(0,10),dres=0.8,ri=c(1.3,1.3,1.5,1.9,2.02,2.05,2.08,2.11),text.cex=0.6,mtext.cex=0.7,p=3,lwd=2,lty=c("11"),font=c(2,1,1),pdbn=1){

	#if (!is.null(mfrow)) par(mfrow=mfrow)
	if (length(mfrow)==2) par(mfrow=mfrow)
	if (!is.null(pdbn)) auxpdb = db$pdbname
	else auxpdb = pdbn
	
	for (i in krange){
		avgwd = round(db$exp$pam[[i]]$silinfo$avg.width,p)
		#if (!is.null(text.marg)) aux.marg=paste(text.marg[1],avgwd," , ",text.marg[2],i," , ",db$enzcode,"/",db$inbcode,sep="")
		if (length(text.marg)==2) aux.marg=paste(text.marg[1],avgwd," , ",text.marg[2],i," , ",db$enzcode,"/",db$inbcode,sep="")
		else aux.marg = text.marg
		if (type==1) plot_model_clus(db$kclus[[i]]$clus,clusvis,auxpdb,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,text.marg=aux.marg, add=add,lcol=lcol,col=col[1],dlim=dlim,dres=dres,ri=ri,text.cex=text.cex,mtext.cex=mtext.cex)#;print("ok")
		if (type==2) plot_model_clus2(db$kclus[[i]]$clus,clusvis,auxpdb,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,text.marg=aux.marg, add=add,lcol=lcol,col=col[1:3],dlim=dlim,dres=dres,ri=ri,text.cex=text.cex,mtext.cex=mtext.cex,restype=restype,lwd=lwd,lty=lty,font=font)
	}
}

# Compõe matrix de cores
compose_palette = function(n,pal="BLUE2RED"){

	tam = 2*n
	
	if (pal=="RED2BLUE")
		auxc = red2blue_v1(tam)
	if (pal=="BLUE2RED")
		auxc = blue2red_v1(tam)
	auxm = c()
	
	for (i in 1:n){
		auxm=rbind(auxm,c(auxc[i],auxc[tam-i+1]))
	}
	rownames(auxm)=NULL
	return(auxm)
}

# Plota mais elabora de modelos, inclusive com sobreposição
multi_plus_plot_model_clus = function(db,krange,main,xlab,ylab,xlim,ylim,clusvis,lcoli,legxy=c(39,75),type=2,mfrow=c(1,1),restype=c(1,1,1), restext=c("I","E","I"),lwd=1,lty=3, font=c(2,2,1),text.cex=0.9,text.marg=c(),pdbn=c(),leg.pch="X",leg.cex=0.8,legrestype=restype){


	tam = length(db)
	
	colv = c()
	colm = compose_palette(tam)
	colm = cbind(colm,colm[,1])
	#print(colm)
	auxname = c()

	if (tam==1) {
		pdbn="ok"
		text.marg=c("Avg.width = ","k = ")
	}
	
	if (is.character(lcoli)) lcol=lcoli
	else lcol = colm[1,lcoli]
	
	if (is.null(text.marg)) text.marg=strcat("k=",paste("(",strcat(paste(krange),collapse=","),")",sep=""))
	#print(text.marg)
	#print(colm);readline()
	multi_plot_model_clus(db[[1]],krange[1],mfrow,main,xlab,ylab,xlim,ylim,type,col=colm[1,],lcol=lcol,restype=restype,lwd=lwd,lty=lty,font=font,text.cex=text.cex,pdbn=pdbn,text.marg=text.marg,clusvis[[1]])
	
	if (legrestype[1]) {
		auxname = c(auxname,paste(db[[1]]$pdbname,restext[1],sep="-"))
		colv = c(colv,colm[1,1])
	}
	if (legrestype[2]) {
		auxname = c(auxname,paste(db[[1]]$pdbname,restext[2],sep="-"))
		colv = c(colv,colm[1,2])
	}
	if (legrestype[3]) {
			auxname = c(auxname,paste(db[[1]]$pdbname,restext[3],sep="-"))
			colv=c(colv,colm[1,3])
	}
	if (!sum(legrestype)){
			auxname = c(auxname,db[[1]]$pdbname)
			colv=c(colv,lcol)
		}
	
	text.marg=c()
	i=2
	
	while (i<=tam){
		#print("ok")
	#for (i in 2:tam){
		if (is.character(lcoli)) lcol=lcoli
		else lcol = colm[i,lcoli]
		multi_plot_model_clus(db[[i]],krange[i],c(),main,xlab,ylab,xlim,ylim,type,col=colm[i,],lcol=lcol,restype=restype,lwd=lwd,lty=lty,font=font,text.cex=text.cex,pdbn=pdbn,text.marg=text.marg,clusvis=clusvis[[i]],add=1)
		#auxname = c(auxname,db[[i]]$pdbname)
		if (legrestype[1]) {
			auxname = c(auxname,paste(db[[i]]$pdbname,restext[1],sep="-"))
			colv=c(colv,colm[i,1])
		}
		if (legrestype[2]) {
			auxname = c(auxname,paste(db[[i]]$pdbname,restext[2],sep="-"))
			colv=c(colv,colm[i,2])
		}
		if (legrestype[3]) {
			auxname = c(auxname,paste(db[[i]]$pdbname,restext[3],sep="-"))
			colv=c(colv,colm[i,3])
		}
		if (!sum(legrestype)){
			auxname = c(auxname,db[[i]]$pdbname)
			colv=c(colv,lcol)
		}
		i=i+1
	
	}
	if (!sum(legrestype))
		legend(legxy[1],legxy[2],auxname,lty=lty,col=colv,text.col=colv,text.font=font[1],cex=leg.cex)
	else
		legend(legxy[1],legxy[2],auxname,pch=leg.pch,col=colv,text.col=colv,text.font=font[1],cex=leg.cex)

}



#OBSOLETED: Constroi um modelo visual da clusterização
plot_model_clus_OLD2 = function(clus,pdbname,main,xlab,ylab,xlim,ylim,text.marg,add,lcol,col,dlim,dres,ri,text.cex,mtext.cex){

	#print(clus)
	#readline()
	if (is.null(xlim)){
		auxminmax = find_minmax_geomxy_in_clus(clus)
		xlim=c(auxminmax[1]-dlim[1],auxminmax[1]+dlim[1])
		ylim=c(auxminmax[3]-dlim[2],auxminmax[4]+dlim[2])
	}
	if (add==0) plot(0,0,xlim=xlim,ylim=ylim,type="n",asp=1,main=paste(main,pdbname,sep=" - "),xlab=xlab,ylab=ylab)
	
	mtext(text.marg,3,0.1,cex=mtext.cex)
	
	tam = 1:length(clus)
	lenri = length(ri)
	
	for (i in tam){

		xgc = clus[[i]]$geomc[1]
		ygc = clus[[i]]$geomc[2]
		x = xgc
		y = ygc
		lenres = length(clus[[i]]$loopres)
		#print(ri[lenres])
		#print(x);print(y)
		#readline()
		if (lenres > lenri) lenres = lenri
		plotcircle(mid=c(x,y),col=col,lcol=lcol,r=ri[lenres])
		text(x,y,strcat(clus[[i]]$loopres),cex=text.cex,font=2)
		res_text(c(xgc,ygc),clus[[i]]$baseres,ri[lenres]+dres,text.cex)

	}
	
}

#OBSOLETED: Constroi um modelo visual da clusterização
plot_model_clus_OLD = function(db,main,xlab,ylab,xlim,ylim,add=0,lcol="GRAY",col="WHITE",dlim=c(0,10),dres=0.8,ri=1.3,text.cex=0.6){

	if (is.null(xlim)){
		auxminmax = find_minmax_geomxy_in_clus(db$clus)
		xlim=c(auxminmax[1]-dlim[1],auxminmax[1]+dlim[1])
		ylim=c(auxminmax[3]-dlim[2],auxminmax[4]+dlim[2])
	}
	if (add==0) plot(0,0,xlim=xlim,ylim=ylim,type="n",asp=1,main=paste(main,db$pdbname,sep=" - "),xlab=xlab,ylab=ylab)
	
	tam = 1:length(db$clus)
	
	for (i in tam){

		xgc = db$clus[[i]]$geomc[1]
		ygc = db$clus[[i]]$geomc[2]
		x = xgc
		y = ygc
		plotcircle(mid=c(x,y),col=col,lcol=lcol,r=ri)
		text(x,y,strcat(db$clus[[i]]$loopres),cex=text.cex,font=2)
		res_text(c(xgc,ygc),db$clus[[i]]$baseres,ri+dres,text.cex)

	}
	
}

# retorna coordenadas xyz em nivel de res de uma lista de atomos (apenas diferentes)
get_res_geom_center_from_atom_list = function(atoml,name,pdb){

	aux = c()

	auxname = name[atoml]
	#print(auxname)
	auxresi = apply(as.matrix(auxname),1,extract_resi)
	auxresn = apply(as.matrix(auxname),1,extract_resn)
	auxresu = unique(auxresi)
	#print(auxresi);print(auxresn);readline()
	for (i in auxresu){
	
		filter = auxresi == i
		filterid = atoml[filter]
		#auxyz = pdb[filterid,]
		if (length(filterid)==1){
			auxyz = pdb[filterid,]
			auxyz = t(auxyz)
		}
		else auxyz = pdb[filterid,]
		auxgc = geometric_center(auxyz)
		#if (is.na(auxgc)[1]) {
		#	print("OPA");readline()
		#}
		#print(auxyz);print(auxgc);readline()
		aux = rbind(aux,auxgc)
		#print(unique(auxresn[filter]))
		#readline()
		#aux = c(aux,unique(auxresn[filter]))
	}
	
	#print(auxresn)
	return(aux)

}

add_res_geom_center_per_clustering = function(db){

	tami = 2:length(db$kclus)
	aux = db	
	
	for (i in tami){
		tamj = 1:length(db$kclus[[i]]$clus)
		for (j in tamj){

			if (db$kclus[[i]]$clus[[j]]$loop[1]!=0){
				if (is.list(db$pdb)) {
					auxyz = get_res_geom_center_from_atom_list(db$kclus[[i]]$clus[[j]]$loop,db$name,db$pdb$xyz)
					rownames(auxyz)=NULL
					#auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$loop,db$name)
					aux$kclus[[i]]$clus[[j]]$loopresgc = auxyz
				} else {
					aux$kclus[[i]]$clus[[j]]$loopresgc = 0
					print(paste("WARNING: no PDB coordinate for geom center",db$pdbname))
				}
			} else aux$kclus[[i]]$clus[[j]]$loopresgc = 0
		
			if (db$kclus[[i]]$clus[[j]]$base[1]!=0){
				if (is.list(db$pdb)){
					auxyz = get_res_geom_center_from_atom_list(db$kclus[[i]]$clus[[j]]$base,db$name,db$pdb$xyz)
					rownames(auxyz)=NULL
					aux$kclus[[i]]$clus[[j]]$baseresgc = auxyz
					#auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$base,db$name)
					#aux$kclus[[i]]$clus[[j]]$baseresn = auxres
				} else {
					aux$kclus[[i]]$clus[[j]]$loopresgc = 0
					print(paste("WARNING: no PDB coordinate for geom center",db$pdbname))
				}
			} else aux$kclus[[i]]$clus[[j]]$baserescg = 0

			if (db$kclus[[i]]$clus[[j]]$noloop[1]!=0){
				if (is.list(db$pdb)){
					auxyz = get_res_geom_center_from_atom_list(db$kclus[[i]]$clus[[j]]$noloop,db$name,db$pdb$xyz)
					rownames(auxyz)=NULL
					aux$kclus[[i]]$clus[[j]]$noloopresgc = auxyz
					#auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$noloop,db$name)
					#aux$kclus[[i]]$clus[[j]]$noloopresn = auxres
				} else {
					aux$kclus[[i]]$clus[[j]]$loopresgc = 0
					print(paste("WARNING: no PDB coordinate for geom center",db$pdbname))
				}
			} else aux$kclus[[i]]$clus[[j]]$noloopresgc = 0
		}

	}
	return(aux)

}


# Adiciona informações dos res em cada cluster (clus)
add_res_per_clustering = function(db){
	
	tami = 2:length(db$kclus)
	aux = db	

	for (i in tami){
		tamj = 1:length(db$kclus[[i]]$clus)
		for (j in tamj){
			if (db$kclus[[i]]$clus[[j]]$loop[1]!=0){
				auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$loop,db$name)
				aux$kclus[[i]]$clus[[j]]$loopresn = auxres
			} else aux$kclus[[i]]$clus[[j]]$loopresn = 0
		
			if (db$kclus[[i]]$clus[[j]]$base[1]!=0){
				auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$base,db$name)
				aux$kclus[[i]]$clus[[j]]$baseresn = auxres
			} else aux$kclus[[i]]$clus[[j]]$baseresn = 0

			if (db$kclus[[i]]$clus[[j]]$noloop[1]!=0){
				auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$noloop,db$name)
				aux$kclus[[i]]$clus[[j]]$noloopresn = auxres
			} else aux$kclus[[i]]$clus[[j]]$noloopresn = 0
		}

	}
	return(aux)

}


# Adiciona elipsoide mínima aos dados
# ATENÇÃO: não funciona bem para clusteres pequenos...
add_ellipsoid_per_clustering = function(db,tami=2:7){

	#tami = 2:length(db$kclus)
	aux = db
	
	
	for (i in tami){
		tamj = 1:length(db$kclus[[i]]$clus)
		for (j in tamj){
			auxe = c()
			if (db$kclus[[i]]$clus[[j]]$loop[1]!=0){
				auxe = c(auxe,db$kclus[[i]]$clus[[j]]$loop)
				#auxell = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$loop,db$name)
				#aux$kclus[[i]]$clus[[j]]$loopresn = auxell
			} #else aux$kclus[[i]]$clus[[j]]$loopresn = 0
		
			if (db$kclus[[i]]$clus[[j]]$base[1]!=0){
				auxe = c(auxe,db$kclus[[i]]$clus[[j]]$base)
				#auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$base,db$name)
				#aux$kclus[[i]]$clus[[j]]$baseresn = auxres
			} #else aux$kclus[[i]]$clus[[j]]$baseresn = 0

			if (db$kclus[[i]]$clus[[j]]$noloop[1]!=0){
				auxe = c(auxe,db$kclus[[i]]$clus[[j]]$noloop)
				#auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$noloop,db$name)
				#aux$kclus[[i]]$clus[[j]]$noloopresn = auxres
			} #else aux$kclus[[i]]$clus[[j]]$noloopresn = 0
			if (!is.null(auxe)){
				#print(auxe);readline()
				if (is.list(db$pdb)){
					auxy = db$pdb$xyz[auxe,1:2]
					if (is.matrix(auxy)){
						#print(paste(db$pdbname,"kclus",i,"clus",j));readline()
						aux$kclus[[i]]$clus[[j]]$ellips = ellipsoidhull(auxy)
					} else {
						aux$kclus[[i]]$clus[[j]]$ellips = 0
					}
				} else {
					print(paste("WARNING: no PDB coordinates found for ellipsoidhull calculation in",db$pdbname))
				}	
			} else {
				print(paste("WARNING: no atom id for ellipsoid calculation in",db$pdbname))
			}
		}

	}
	return(aux)


}

# OBSOLETED: Adiciona elipsoide mínima aos dados
# ATENÇÃO: não funciona bem para clusteres pequenos...
add_ellipsoid_per_clustering_OLD = function(db,tami=2:7){

	#tami = 2:length(db$kclus)
	aux = db
	
	
	for (i in tami){
		tamj = 1:length(db$kclus[[i]]$clus)
		for (j in tamj){
			auxe = c()
			if (db$kclus[[i]]$clus[[j]]$loop[1]!=0){
				auxe = c(auxe,db$kclus[[i]]$clus[[j]]$loop)
				#auxell = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$loop,db$name)
				#aux$kclus[[i]]$clus[[j]]$loopresn = auxell
			} #else aux$kclus[[i]]$clus[[j]]$loopresn = 0
		
			if (db$kclus[[i]]$clus[[j]]$base[1]!=0){
				auxe = c(auxe,db$kclus[[i]]$clus[[j]]$base)
				#auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$base,db$name)
				#aux$kclus[[i]]$clus[[j]]$baseresn = auxres
			} #else aux$kclus[[i]]$clus[[j]]$baseresn = 0

			if (db$kclus[[i]]$clus[[j]]$noloop[1]!=0){
				auxe = c(auxe,db$kclus[[i]]$clus[[j]]$noloop)
				#auxres = get_res_from_atom_list(db$kclus[[i]]$clus[[j]]$noloop,db$name)
				#aux$kclus[[i]]$clus[[j]]$noloopresn = auxres
			} #else aux$kclus[[i]]$clus[[j]]$noloopresn = 0
			if (!is.null(auxe)){
				#print(auxe);readline()
				auxy = db$pdb$xyz[auxe,1:2]
				aux$kclus[[i]]$clus[[j]]$ellips = ellipsoidhull(auxy)
			} else {
				print("WARNING: no atom id for ellipsoid calculation")
			}
		}

	}
	return(aux)


}


###############################################################################################################
# PRIMEIRAS ANALISES
#versão mais simples da manifold_travel
manifold_travel2 = function(x,y,autov,rot,main.text,colei,zoom=c(min(x),max(x),min(y),max(y)),met="SVD"){

	xlabel=" autovector"
	ylabel=xlabel
	#colei=colv(name)
	plot(x,y,pch=paste(""),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""),xlim=zoom[1:2],ylim=zoom[3:4])
	#points(circx,circy,cex=5,lty=3)
	text(x,y,labels=rot,cex=0.60,col=colei)
	abline(v=0,h=0,lty=3)

}

# visualiza planos dos autovetores
#r=manifold_travel(m,lap,autov=c(v1,v2),con=0,dec,scale=sc,psvd=c(0,0,0,0,1),peig=c(0,0,0,0,1),zoom=c(minus,maxus,minus1,maxus1),main.text)
manifold_travel = function(m,lap,autov,con,dec,scale,psvd,peig,zoom,main.text,hubs=c(100,100)){

	if (dec=="eig") eig=1
	else eig=0
	if (dec=="svd") svd=1
	else svd=0
	#maxus = max(m$Lw$svd$uv)
	#minus = min(m$Lw$svd$uv)
	#zoom = c(minus,maxus,minus,maxus)
	rot = m$colname
	#colei = colv(m$name)
	colei = "black"
	xlabel=" autovector"
	ylabel=xlabel
	#print("0k")
	r = list()
	if (con){
		if (svd){
			if (eig){
				aux=sum(c(psvd,peig))
				if (aux==8) {
					#print("ok1")
					par(mfrow=c(2,4))
				} else {
					#print("ok2")
					par(mfrow=c(1,sum(c(psvd,peig))))
				}
			} else { 
				#print("ok3")
				par(mfrow=c(1,sum(psvd)))
			}
		} else {
			#print("ok4")
			par(mfrow=c(1,sum(peig)))
		}
	}
	 #par(mfrow=c(1,4))
	if (svd){
		met = "SVD"
		if (lap=="A"){
			if (scale){
				x = m$A$svd$us[,autov[1]]
				y = m$A$svd$us[,autov[2]]
			} else {
				x = m$A$svd$u[,autov[1]]
				y = m$A$svd$u[,autov[2]]
			}
		}
		if (lap=="Lw"){
			if (scale){
				x = m$Lw$svd$us[,autov[1]]
				y = m$Lw$svd$us[,autov[2]]#;print("ok")
			} else {
				x = m$Lw$svd$u[,autov[1]]
				y = m$Lw$svd$u[,autov[2]]
			}
		}
		if (lap=="L"){
			if (scale){
				x = m$L$svd$us[,autov[1]]
				y = m$L$svd$us[,autov[2]]
			} else {
				x = m$L$svd$u[,autov[1]]
				y = m$L$svd$u[,autov[2]]
			}
		}
		#r$svd = cbind(k=1:length(x),x,y,W=diag(m$L$a),w=edges(m$d))
		#r$svd = cbind(k=m$id,x,y,W=diag(m$L$a),w=edges(m$d))
		#circx = c(x[hubs])
		#circy = c(y[hubs])
		#plot(x,y,pch=substr(rot,1,1),xlab=paste(auto[1]," autovalue"),ylab=paste(auto[2]," autovalue"),main=paste(met," from Laplacian"),col=colei)
		#points(circx,circy,cex=5,lty=3)
		#print(psvd)
		if (psvd[1]){
			plot(x,y,pch=substr(rot,2,2),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""),col=colei)
			#print(paste(substr(rot,2,2),r$svd[,1]))
			#plot(x,y,pch=paste(substr(rot,2,2),r$svd[,1]),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met," from Laplacian"),col=colei)
			#points(circx,circy,cex=5,lty=3)
		}
		if (psvd[2]){
			plot(x,y,pch=charv(rot,"C"),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""),col=colei)
			#points(circx,circy,cex=5,lty=3)
		}
		if (psvd[3]){
			plot(x,y,pch=paste(""),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""))
			#points(circx,circy,cex=5,lty=3)
			text(x,y,labels=paste(r$svd[,5]),cex=0.90,col=colei)
		}
		if (psvd[4]){
			plot(x,y,pch=paste(""),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""))
			#points(circx,circy,cex=5,lty=3)
			text(x,y,labels=paste(r$svd[,1]),cex=0.80,col=colei)
		}
		if (psvd[5]){
			#print(zoom)
			plot(x,y,pch=paste(""),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""),xlim=zoom[1:2],ylim=zoom[3:4])
			#points(circx,circy,cex=5,lty=3)
			#text(x,y,labels=paste(substr(rot,2,2),r$svd[,1],sep=""),cex=0.60,col=colei)
			text(x,y,labels=m$colnames,cex=0.60,col=colei)
			abline(v=0,h=0,lty=3)
		}

	}

	if (eig) {
		met = "EIGEN"
		if (lap=="A"){
			if (scale){
				x = m$A$eig$vs[,autov[1]]
				y = m$A$eig$vs[,autov[2]]
			} else {
				x = m$A$eig$v[,autov[1]]
				y = m$A$eig$v[,autov[2]]
			}
		}
		if (lap=="Lw"){
			if (scale){
				x = m$Lw$eig$vs[,autov[1]]
				y = m$Lw$eig$vs[,autov[2]]
			} else {
				x = m$Lw$eig$v[,autov[1]]
				y = m$Lw$eig$v[,autov[2]]
			}
		}
		if (lap=="L"){
			if (scale){
				x = m$L$eig$vs[,autov[1]]
				y = m$L$eig$vs[,autov[2]]
			} else {
				x = m$L$eig$v[,autov[1]]
				y = m$L$eig$v[,autov[2]]
			}
		}
		#r$eig = cbind(k=1:length(x),x,y,W=diag(m$L$a),w=edges(m$d))
		r$eig = cbind(k=m$id,x,y,W=diag(m$L$a),w=edges(m$d))
		#circx = c(x[hubs])
		#circy = c(y[hubs])
		if (peig[1]){
			plot(x,y,pch=substr(rot,2,2),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""),col=colei)
			#points(circx,circy,cex=5,lty=3)
		}
		if (peig[2]){
			plot(x,y,pch=charv(rot,"C"),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""),col=colei)
			#points(circx,circy,cex=5,lty=3)
		}
		if (peig[3]){
			plot(x,y,pch=paste(""),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""),col=colei)
			#points(circx,circy,cex=5,lty=3)
			text(x,y,labels=paste(r$eig[,5]),cex=0.90,col=colei)
		}
		if (peig[4]){
			plot(x,y,pch=paste(""),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""))
			#points(circx,circy,cex=5,lty=3)
			text(x,y,labels=paste(r$eig[,1]),cex=0.80,col=colei)
		}
		if (psvd[5]){
			plot(x,y,pch=paste(""),xlab=paste(autov[1],xlabel),ylab=paste(autov[2],ylabel),main=paste(met,": ",main.text,sep=""),xlim=zoom[1:2],ylim=zoom[3:4])
			#points(circx,circy,cex=5,lty=3)
			text(x,y,labels=paste(substr(rot,2,2),r$eig[,1],sep=""),cex=0.60,col=colei)
			abline(v=0,h=0,lty=3)
		}
	}
	return(r)
} 



make_line = function(coord,xo,yo,sx,sy,seg,vet){

	tamv = length(vet)
	tamseg = length(seg)
	x = c()
	y = c()
	xi = xo
	yi = yo

	for (i in 1:tamseg){
		for (j in 1:tamv){
			if (seg[i]==vet[j]){
				#x = c(x,xi)
				#y = c(y,yi)
				#com = rbind(com,c(xi,yi))
				coord[j,1] = xi
				coord[j,2] = yi
				xi = xi + sx
				yi = yi + sy
			}	
		}
	}
	return(coord)
}

view3d_with_base_axis = function(m,b=diag(3),w=c(1,1,1),p0=c(0,0,0),lim=NULL,col="BLACK",lty=3,type="rotation",lab=c("x","y","z"),size=2,k=1,limcex=1.1,centralize=F,svd=T,also2d=T,axis=T,new=F,close=F,add=F,text=F,adj=3,textcex=3,plane=F){

	if (close){
		while(rgl.cur()!=0) rgl.close()
	}

	if ((!add)&(new)){
		#if(rgl.cur()!=0) rgl.close()
		open3d()
	}

	col = toupper(col)

	if (svd){
		rownames(m)=NULL
		colnames(m)=NULL
		m0 = scale(m,scale=F)#;print(m0)
		p0 = attributes(m0)[[2]]#;print(p0);readline()
		m0.svd = svd(m0)
		b = m0.svd$v
		w = m0.svd$d		
	}
	#print(is.matrix(m))

	if (centralize){
		p0 = c(0,0,0)
		if (svd){
			m = m0
		}else{
			m = scale(m,scale=F)
		}
	}

	if (is.null(lim)&(!add)){
		maxx = max(abs(m[,1:2]))#;print("o")
		maxx = maxx*limcex
		lim = c(-maxx,maxx)
		#maxlim = max(m)*limcex
		#minlim = min(m)*limcex
		#lim = c(minlim,maxlim)
		#const = sqrt(2)/2
		#D = Norm(p0)
		#lim = c(D-(w[1]*limcex),D+(w[1]*limcex))
		#lim = c(min(m[,1:2])*(1/limcex),max(m[,1:2])*limcex)
		#maxm = max(m)/2
		#maxc = c(maxm,maxm,maxm)
		#d = po-maxc
		#lim = c(maxc,0)+d[1:2],c(0,maxc)+d[1:2])
		#xlim = c(min(m[,1]),max(m[,1]))
		#ylim = c(min(m[,2]),max(m[,2]))
		#zlim = c(min(m[,3]),max(m[,3]))
		#print(p0);print(w[1]);print(lim)
	}

	#print(lim)
	
	w = w*k
	ids=plot3d(x=m[,1],y=m[,2],z=m[,3],type="s",size=size,box=F,col=col,xlab=lab[1],ylab=lab[2],zlab=lab[3],xlim=lim,ylim=lim,zlim=lim,add=add)
	if (text){
		text3d(x=m[,1],y=m[,2],z=m[,3],text=paste(1:dim(m)[1]),adj=adj,cex=textcex)

	}

	if (plane){
		x = m[,1]
		y = m[,2]
		z = m[,3]
		res = lm(z ~ x + y)
		ap = round(res$coefficients["x"],3)
		bp = round(res$coefficients["y"],3)
		cp = -1
		dp = round(res$coefficients[1],3)
		planes3d(ap, bp, cp, dp, alpha = 0.2,col=col)
		b = t(tomatrix(c(ap,bp,cp)));print(b);print(res)
	}

	if (also2d){
		if (!add){
			plot(x=m[,1],y=m[,2],type="p",col=col,xlim=lim,ylim=lim)
		}else{
			points(x=m[,1],y=m[,2],col=col,xlim=lim,ylim=lim)
		}
	}
	if (axis){
#		for (i in 1:dim(b)[1]){
		for (i in 1:dim(b)[2]){
			if (round(w[i],3)){
				p1 = c(w[i]*b[,i])
				arrow3d(p0=p0,p1=p0+p1,col=col,lty=lty,type=type)
			}
		}
	}
	#plot3d(xx[,1],xx[,2],xx[,3],type="s",size=3,box=F,col=col,xlab=lab[1],ylab=lab[2],zlab=lab[3],add=T)


	return(ids)

}

rot180 = function(v,k=5){


	if (round(Norm(v),k)!=1){
		print(paste("WARNING: v must be unitary for rot180. Making it unitary..."))
		v = v/Norm(v)
	}
	
	m = diag(4)
	m[1,1] = 2*v[1]*v[1]-1
	m[1,2] = 2*v[1]*v[2]
	m[1,3] = 2*v[1]*v[3]

	m[2,1] = 2*v[2]*v[1]
	m[2,2] = 2*v[2]*v[2]-1
	m[2,3] = 2*v[2]*v[3]

	m[3,1] = 2*v[3]*v[1]
	m[3,2] = 2*v[3]*v[2]
	m[3,3] = 2*v[3]*v[3]-1

	return(t(round(m,k)))

}

#https://stackoverflow.com/questions/42520301/find-rotation-matrix-of-one-vector-to-another-using-r
rot_vectors = function(x,y){
  u=x/sqrt(sum(x^2))

  v=y-sum(u*y)*u
  v=v/sqrt(sum(v^2))

  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))

  sint=sqrt(1-cost^2);

  return(diag(length(x)) - u %*% t(u) - v %*% t(v) + cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v)))
}

get_rot_matrix_by_plane = function(m,vz=c(0,0,1),vy=c(0,1,0),k=3,svd=T){

	if (dim(m)[1]>2){
		m0 = scale(m,scale=F)
		p0 = attributes(m0)[[2]]
		rot = diag(4)
		if (svd){
			m0.svd = svd(m0)
			#v1 = t(m0.svd$v)[,1]#;print(m0.svd$v);
			v1 = m0.svd$v[,1]#;print(v1);readline()
			rot0 = rot_vectors(v1,vy)
			if (sum(is.nan(rot0))>0){
				print(paste("WARNING: trying to find a rot matrix with two aligned vectors"))
				#return(NULL)
				rot0 = diag(3)
			}else{
				rot[1:3,1:3] = rot0
				rot[,4] = c(p0,1)
				rot[4,] = c(-p0,1)
				m = transform_by_rot_matrix(m,rot)
				m0 = scale(m,scale=F)
				p0 = attributes(m0)[[2]]
			}
		}

		x = m0[,1]
		y = m0[,2]
		z = m0[,3]
		res = lm(z ~ x + y)
	
		print(summary(res))#;readline()
		print(round(res$coef,k))
		

		a = round(res$coefficients["x"],k)
		b = round(res$coefficients["y"],k)
		c = -1
		d = round(res$coefficients[1],k)

		v1 = normalize.vector(c(a,b,c))

		#rot = diag(4)
		rot1 = rot_vectors(v1,vz)
		if (sum(is.nan(rot1))>0){
			print(paste("WARNING: trying to find a rot matrix with two aligned vectors"))
			rot = diag(4)
			return(NULL)
		}else if (svd){
			rot[1:3,1:3] = rot0 %*% rot1
		}else{
			rot[1:3,1:3] = rot1
		}
		
		rot[,4] = c(p0,1)
		rot[4,] = c(-p0,1)
	
		#print(rot);readline()
		return(rot)


	}else{
		print(paste("WARNING: matrix dimension is less than 3; it is not possible to find a plane"))
	}

}

get_rot_matrix_by_plane_old = function(m,v0=c(0,0,1),k=3,svd=T){

	if (dim(m)[1]>2){
		m0 = scale(m,scale=F)
		p0 = attributes(m0)[[2]]
		x = m0[,1]
		y = m0[,2]
		z = m0[,3]
		res = lm(z ~ x + y)
	#	print(res);readline()

		a = round(res$coefficients["x"],k)
		b = round(res$coefficients["y"],k)
		c = -1
		d = round(res$coefficients[1],k)

		v1 = normalize.vector(c(a,b,c))

		rot = diag(4)
		rot0 = rot_vectors(v1,v0)
		if (sum(is.nan(rot0))>0){
			print(paste("WARNING: trying to find a rot matrix with two aligned vectors"))
			return(NULL)
		}else{
			rot[1:3,1:3] = rot0
		}
		rot[,4] = c(p0,1)
		rot[4,] = c(-p0,1)

		if (svd){
			m = transform_by_rot_matrix(m,rot)
			m0 = scale(m,scale=F)
			p0 = attributes(m0)[[2]]
			m0.svd = svd(m0)
			v0 = c(0,1,0)
			#v1 = t(m0.svd$v)[,1]#;print(m0.svd$v);
			v1 = m0.svd$v[,1]#;print(v1);readline()
			rot1 = rot_vectors(v1,v0)
			if (sum(is.nan(rot1))>0){
				print(paste("WARNING: trying to find a rot matrix with two aligned vectors"))
				return(NULL)
			}else{
				rot[1:3,1:3] = rot0 %*% rot1
			}
			rot[,4] = c(p0,1)
			rot[4,] = c(-p0,1)
		}
		
		#print(rot);readline()
		return(rot)


	}else{
		print(paste("WARNING: matrix dimension is less than 3; it is not possible to find a plane"))
	}

}

diag_x = function(n,v=c(1,2)){

	m = diag(n)
	temp = m[v[1],]
	m[v[1],] = m[v[2],]
	m[v[2],] = temp
	return(m)

}

get_rot_matrix = function(source,target,svd=c(T,T),center=c(T,T),basevector=F,k=5){

	source = round(source,k)
	target = round(target,k)

	if (center[1]) source0 = scale(source,scale=F)
	else source0 = source
	if (center[2]) target0 = scale(target,scale=F)
	else target0 = target

	if (svd[1]) {
		source0.svd = svd(source0)
		vs = t(source0.svd$v)
	}else{
		vs = source0
	}
	if (svd[2]){
		target0.svd = svd(target0)
		vt = t(target0.svd$v)
	}else{
		vt = target0
	}			
	#print(vs);print(vt);readline()	
#	vs = t(source0.svd$v)
#	vt = t(target0.svd$v)

	rot = diag(4)

	rot[1,1] = sum(vs[,1]*vt[,1])
	rot[1,2] = sum(vs[,1]*vt[,2])
	rot[1,3] = sum(vs[,1]*vt[,3])

	rot[2,1] = sum(vs[,2]*vt[,1]) 
	rot[2,2] = sum(vs[,2]*vt[,2])
	rot[2,3] = sum(vs[,2]*vt[,3])

	rot[3,1] = sum(vs[,3]*vt[,1]) 
	rot[3,2] = sum(vs[,3]*vt[,2])
	rot[3,3] = sum(vs[,3]*vt[,3])
	

	if (svd[1]&svd[2]){
		rot[,4] = c(colMeans(target),1)
		rot[4,] = c(-colMeans(source),1)
	}else if (svd[1]){
		cM = colMeans(source)
		rot[,4] = c(cM,1)
		rot[4,] = c(-cM,1)
	}else{
		cM = colMeans(target)
		rot[,4] = c(cM,1)
		rot[4,] = c(-cM,1)
	}

	rot = round(rot,k)
	
	if (basevector){
		res = list()
		res$rot = rot
		res$svd.s = source0.svd
		res$svd.t = target0.svd
		return(res)
	}else{
		return(rot)
	}
	

}


transform_by_rot_matrix_fit = function(xyz,rot){

	if (Norm(rot[,4])!=1){
		print(paste("WARNING1: rot matrix is in old format"))
	}else if (sum(rot[,4])!=1){
		print(paste("WARNING2: rot matrix is in old format"))
	}
	
	xyz0 = scale(xyz,scale=F)
	xyz0 = cbind(xyz0,rep(1,dim(xyz0)[1]))
	xyz0 = xyz0 %*% rot

	return(xyz0[,1:3])

}

transform_by_rot_matrix = function(xyz, rot,self=NULL){

	if (!is.matrix(rot)){
		rot=matrix(rot,ncol=4,byrow=F) ### ATENCAO: byrow mudou de T para F !!!
	}

	xyz = cbind(xyz,rep(1,dim(xyz)[1]))

	if (!is.matrix(xyz)){
		xyz = as.matrix(xyz)
		if (sum(dim(xyz)==c(4,1))==2) xyz = t(xyz)
	}

	# TRANSLATION PREROTATION
	prerot = diag(4)
	prerot[4,] = rot[4,]

	#print(prerot)
	#print(scale(xyz,s=F))
	xyz = xyz %*% prerot

	#print(xyz);readline()

	# INNER ROTATION
	inrot = diag(4)
	inrot[1:3,1:3] = rot[1:3,1:3]
	xyz = xyz %*% inrot

	# TRANSLATION POSROTATION
	posrot = diag(4)
	posrot[4,] = rot[,4]

	#print(posrot)#;readline()
	xyz = xyz %*% posrot
	
	#xyz = scale(xyz,scale=F) 
	#equivale xyz[,i] - colMeans(xyz)[i], ou seja, colocar o centroide em (0,0,0)
	
	#xyz = cbind(xyz,rep(1,dim(xyz)[1]))
	#xyz = xyz %*% rot

	#print(xyz);readline()

	return(xyz[,1:3])

}

transform_by_rot_matrix_try = function(xyz, rot, self=F, k=3){

	if (!is.matrix(rot)){
		rot=matrix(rot,ncol=4,byrow=F) ### ATENCAO: byrow mudou de T para F !!!
	}

	#rot=diag(4)
	#print(xyz)
	#print(rot)

	if (self) colmean = round(colMeans(xyz),k)#;print(colmean);readline()

	xyz = cbind(xyz,rep(1,dim(xyz)[1]))

	# TRANSLATION PREROTATION
	prerot = diag(4)
	if (self){
		prerot[4,] = c(-colmean,1)
	}else{
		prerot[4,] = rot[4,]
	}
	#print(prerot)
	#print(scale(xyz,s=F))
	xyz = xyz %*% prerot

	#print(xyz);readline()

	# INNER ROTATION
	inrot = diag(4)
	inrot[1:3,1:3] = rot[1:3,1:3]
	xyz = xyz %*% inrot

	# TRANSLATION POSROTATION
	posrot = diag(4)
	if (self){
		posrot[4,] = c(colmean,1)
	}else{
		posrot[4,] = rot[,4]
	}
	#print(posrot)#;readline()
	xyz = xyz %*% posrot
	
	#xyz = scale(xyz,scale=F) 
	#equivale xyz[,i] - colMeans(xyz)[i], ou seja, colocar o centroide em (0,0,0)
	
	#xyz = cbind(xyz,rep(1,dim(xyz)[1]))
	#xyz = xyz %*% rot

	#print(xyz);readline()

	return(xyz[,1:3])

}


transform_by_rot_matrix_ce = function(xyz,rot){
	
	#print(xyz)
	if (!is.matrix(rot)){
		rot=matrix(rot,ncol=4,byrow=T)
	}

	#print(xyz)
	auxyz = cbind(xyz,rep(1,dim(xyz)[1]))
	auxyz = t(auxyz)

	#print(auxyz);readline()
	
	#translation pre-rotation
	auxr1 = diag(4)#;print(auxr1)
	auxr1[,4]=rot[4,]#;print(rot[4,])
	#print(auxr1);readline()
	auxr =  auxr1 %*% auxyz
	#print(auxr);readline()
	#rotation
	auxr2 = rot
	auxr2[4,] = c(0,0,0,1)
	auxr2[,4] = c(0,0,0,1)
	#print(auxr2);readline()
	auxr =  auxr2 %*% auxr

	#translation pos-rotation
	auxr3 = diag(4)
	auxr3[,4]=rot[,4]
	#print(auxr3);readline()
	auxr =  auxr3 %*% auxr

	auxr = t(auxr)[,1:3]
	#print(auxr)
	return(auxr)

}

get_xylim_from_model = function(db,dball,graphname,pdbid,sep="_",dxy=NULL,rotref=diag(4),self=F){

	#rotid =  grep(pdbname,dball$rotref)
	aux = list()#;print("ok")
	#pars1 = unlist(strsplit(graphname,sep))
	pars1 = pdbsplit(graphname)
	pdbname = pars1[1,"pdb"]#;print(pdbname)
	rotid =  grep(pdbname,dball$rotref)
	#rotid =  grep(graphname,dball$rotref)
	#print(pdbname);readline()
#	rotref = dball$rotref[[rotid[1]]]$rot
	#rotref = diag(4)
	#pdbid1 = as.numeric(pars1[5])
	#clusid1 = as.numeric(pars1[6])
	pdbid1 = as.numeric(pars1[1,"id"])
	clusid1 = as.numeric(pars1[1,"nclus"])
	#print(pdbid1);print(clusid1);readline()
	
#	auxyz1=transform_by_rot_matrix(db[[pdbid1]]$exp$superclus[[clusid1]]$geomc,rotref)

	i = 1
	repeat{
		if (!is.null(dball$res[[pdbid]]$alignment[[1]]$aligndata[[i]])){
			break;
		}else{
			i = i+1
		}
	}
	auxyz1 = dball$res[[pdbid]]$alignment[[1]]$aligndata[[i]]$supxyz#;print(pdbid)

	#nclus = dim(auxyz1)[1]
	#xyzref = db[[1]]$exp$superclus[[nclus]]$geomc
	#rotref = get_rot_matrix(target=xyzref,source=auxyz1,base=F)
	#auxyz1 = transform_by_rot_matrix(auxyz1,rotref,self=NULL)#;print(auxzz);readline()


	auxyz1 = transform_by_rot_matrix(auxyz1,rotref,self=self)
	#auxyz1 = transform_by_rot_matrix(auxyz1,rotref,self=self)
	xlim = range(auxyz1[,1])#;print(xlim);readline()
	ylim = range(auxyz1[,2])
	if (!is.null(dxy)){
		xlim[1] = xlim[1]-dxy[1]
		xlim[2] = xlim[2]+dxy[1]
		ylim[1] = ylim[1]-dxy[2]
		ylim[2] = ylim[2]+dxy[2]
	}
	aux$xlim = xlim
	aux$ylim = ylim
	return(aux)
}


visualize_align_supercluster_hot = function(db,dball,dbhot,pdbid,cols,dxy,rotref=diag(4),self=T){

	tam = length(dball$res[[pdbid]]$alignment)
	mfrow = c(tam,2)
	i = 1
	visualize_align_super_cluster(db=db,dball=dball,dbhot=dbhot,pdbid=pdbid,clusid=i,ids=1:2,cols=cols,mfrow=mfrow, dxy=dxy,rotref=rotref,self=self)
	mfrow = NULL
	for (i in 2:tam){
		visualize_align_super_cluster(db=db,dball=dball,dbhot=dbhot,pdbid=pdbid,clusid=i,ids=1:2,cols=cols,mfrow=mfrow, dxy=dxy,rotref=rotref,self=self)		
	}
}


visualize_align_supercluster_prescore = function(db,dball,dbhot,pdbnames,cols,dxy,rotref=diag(4),self=F){


	id1 = as.numeric(as.vector(dball$group[which(dball$group$pdb==pdbnames[1]),"id"]))#;print(id1)
	id2 = as.numeric(as.vector(dball$group[which(dball$group$pdb==pdbnames[2]),"id"]))#;print(id2)
	tam1 = length(id1)
	tam2 = length(id2)
	parid=c(id1[1],1,0)
	#print(id1);print(id2);print(parid);readline()

	if (!self) mfrow=c(tam1,tam2+1)
	else mfrow=c(1,tam2+1)
	ids=grep(pdbnames[2],dball$res[[parid[1]]]$alignment[[parid[2]]]$prescore$n2)#;print(ids)#;readline()
	visualize_align_super_cluster(db=db,dball=dball,pdbid=parid[1],clusid=parid[2],ids=c(1,ids),dbhot=dbhot,score=parid[3],cols=cols,mfrow=mfrow,dxy=dxy,rotref=rotref,self=self)#;readline()
	if ((tam1>1)&(!self)){
	#if (tam1>1){
		for (j in 2:tam1){
			parid=c(id1[1],j,0)
			#print(parid);readline()
			ids=grep(pdbnames[2],dball$res[[parid[1]]]$alignment[[parid[2]]]$prescore$n2)#;print(ids)
			visualize_align_super_cluster(db=db,dball,pdbid=parid[1],clusid=parid[2],ids=c(1,ids),dbhot=dbhot,score=parid[3],cols=cols,mfrow=NULL,dxy=dxy,rotref=rotref,self=self)#;readline()
		}
	}

}


visualize_align_super_cluster = function(db,dball,pdbid,ids,clusid,dbhot=NULL,score=1,xlim=NULL,ylim=NULL,cols=NULL,mfrow=NULL,sep="_",dxy=c(4,2),sub="",subid=5,k=3,rotref=diag(4),self=F){


	#graphnames = as.character(dball$score[,1])[1]
	#graphnames = as.character(dball$group$fakename)[1]
	
	#print(graphnames);readline()

	#auxlim = get_xylim_from_model(db=db,dball=dball,graphnames[1],dxy=dxy,sep="_")#;print("ok")
	#xlim = auxlim$xlim
	#ylim = auxlim$ylim

	if (clusid==0){#print("ok")
		graphnames = as.character(dball$res[[pdbid]]$scorend$n2)[ids[1]]
		#print(graphnames);readline()
		auxlim = get_xylim_from_model(db=db,dball,graphname=graphnames[1],pdbid=pdbid,dxy=dxy,sep="_",rotref=rotref,self=self)#;print("ok")
		xlim = auxlim$xlim
		ylim = auxlim$ylim
		#print(auxlim)
		graphnames = as.character(dball$res[[pdbid]]$scorend$n2)[ids]#;print(graphnames)
		#dbhot = dbhot[[ids]]
		auxn1 = as.factor(as.character(dball$res[[pdbid]]$scorend$n1)[ids])
		#print(auxn1);readline()
		#auxc1 = pdbsplit(as.character(auxn1))[,6]
		auxc1 = pdbsplit(as.character(auxn1))[,"nclus"]
		#auxc1 = auxc1[,dim(auxc1)[2]]
		#print(auxc1)#;readline()
		#print(auxn1);readline()
		levels(auxn1) = 1:length(levels(auxn1))
		auxn1 = as.numeric(auxn1)
		#print(pdbid);print(ids);readline()
		subscoret =dball$res[[pdbid]]$scorend$scoret
		subscorei =dball$res[[pdbid]]$scorend$scorei
		subscorea =dball$res[[pdbid]]$scorend$scorea
		sub1=0;
		sub2=0;
		sub3=0;
		if (!is.null(subscoret)) sub1=round(subscoret[ids],k)
		if (!is.null(subscorei)) sub2=round(subscorei[ids],k)
		if (!is.null(subscorea)) sub3=round(subscorea[ids],k)
		sub=paste(sub1,sub2,sub3,"(",auxc1,")")#;print(sub);readline()
	}else{
		#graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$score$pdbid)[1]
		#if (self) graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$prescore$n1)[ids[1]]
		#else 
		graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$prescore$n2)[ids[1]]
		#print(graphnames);readline()
		auxlim = get_xylim_from_model(db=db,dball=dball,graphname=graphnames[1],pdbid=pdbid,dxy=dxy,sep="_",rotref=rotref,self=self)#;print("ok")
		xlim = auxlim$xlim
		ylim = auxlim$ylim
		#print(auxlim)
		#if (score==1) {
		#	graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$score$pdbid)[ids]
		#	tam = length(graphnames)
		#	sub = round(dball$res[[pdbid]]$alignment[[clusid]]$score$sclocal,k)[ids]
		#}
		#else {
			#if (self) graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$prescore$n1)[ids]
			#else 
			graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$prescore$n2)[ids]
			tam = length(graphnames)
			#sub1 = round(dball$res[[pdbid]]$alignment[[clusid]]$prescore$scorep,k)[ids]
			#sub2 = round(dball$res[[pdbid]]$alignment[[clusid]]$prescore$scoret,k)[ids]
			sub1 = round(dball$res[[pdbid]]$alignment[[clusid]]$prescore$scoret,k)[ids]
			sub2 = round(dball$res[[pdbid]]$alignment[[clusid]]$prescore$kernel,k)[ids]
			sub3 = round(dball$res[[pdbid]]$alignment[[clusid]]$prescore$scorea,k)[ids]
			#sub4 = round(dball$res[[pdbid]]$alignment[[clusid]]$prescore$scoret,k)[ids]
			sub = paste(sub1,sub2,sub3,sep=sep)
		#}
	
	}
	

	tam = length(graphnames)#;print(tam);readline()
	
	if (!is.null(mfrow)) par(mfrow=mfrow)
	#print(graphnames)

	for (i in 1:tam){
		#print(graphnames[i])#;readline()
		#pars2 = unlist(strsplit(graphnames[i],sep))#;print(pars2);readline()
		pars2 = pdbsplit(graphnames[i])#;print(pars2);readline()
		pdbname = pars2[1]#;print(pdbname)	
		#pdbid2 = as.numeric(pars2[5])
		#clusid2 = as.numeric(pars2[6])
		pdbid2 = as.numeric(pars2[1,"id"])
		clusid2 = as.numeric(pars2[1,"nclus"])
		#auxyz2=db[[pdbid2]]$exp$superclus[[clusid2]]$geomc
		#idgrp = which(graphnames[i]==dball$group$fakename)
		#idgrp = dball$res[[pdbid]]$scorend$id2[ids[i]]
		if (clusid==0) {
			#idali = auxn1[i]
			idgrp = dball$res[[pdbid]]$scorend$id2[ids[i]]
			idali = dball$res[[pdbid]]$scorend$id1[ids[i]]
		}else{
			idgrp = dball$res[[pdbid]]$alignment[[clusid]]$prescore$id2[ids[i]]
			idali = dball$res[[pdbid]]$alignment[[clusid]]$prescore$id1[ids[i]]
			
		}
		#print(idali);print(idgrp);readline()
		auxzz = dball$res[[pdbid]]$alignment[[idali]]$aligndata[[idgrp]]$supxyz#;print(auxzz);readline()
		nclus = dim(auxzz)[1]
		#xyzref = db[[1]]$exp$superclus[[nclus]]$geomc
		#rotref = get_rot_matrix(target=xyzref,source=auxzz,base=F)
		auxzz = transform_by_rot_matrix(auxzz,rotref,self=self)#;print(auxzz);readline()
		
		#print(dball$aligndata[[id]]$names)
		#print(idgrp);print(auxzz);readline()
		#if (!is.null(cols)) plot3d(auxzz,box=F,type="s",col=cols[i],add=T)#;print(pdbid2);print(clusid2);readline()
		rotid =  grep(pdbname,dball$rotref)#;print(pdbname);print(rotid)
#		rot = dball$rotref[[rotid[1]]]$rot
		rot = dball$res[[pdbid]]$alignment[[idali]]$aligndata[[idgrp]]$rot
		sid = dball$res[[pdbid]]$alignment[[idali]]$aligndata[[idgrp]]$score$global$sid#;print(sid)#;readline()
		vcolor = rep("white",dim(auxzz)[1])
		fcolor = rep("gray",dim(auxzz)[1])
		if (is.null(dbhot)){	
			vcolor[sid]="skyblue4"#;print(vcolor)#;readline()
			cols=vcolor
			visualize_super_cluster(db[[pdbid2]],i=clusid2,typeg="ga",resize=1,coordi=auxzz,cols=cols,xlim=xlim,ylim=ylim,sub=sub[i],rot=rot)
		}else{#print(ids[i]);print(dbhot[[ids[i]]]);readline()
			#if (length(dbhot)==1){
			#	dbhoti = dbhot[[1]]
			#}else{
			#	dbhoti = dbhot[[ids[i]]]
			#}
			#;print(dbhot);readline()
			dbhoti = dbhot[[pdbid2]]
			#nclus = dim(auxzz)[1]#;print(nclus)
			nhot = sapply(dbhoti,dim)#;print(nhot);readline()
			id=which(nclus==nhot[1,])#	;print(id)
			if (length(id)){
				vcol = dbhoti[[id]]$score
				cols = gray.colors(100,start=0.3,end=1)#;print(vcol)
				#cols = brewer.pal(n = 9, name = 'YlOrRd')
				cols = colorRampPalette(c("brown","red", "yellow", "white"))(n = 300)
				#idcol = round(mapply(renorm_interval,vcol,inv=F),0)#;print(idcol);readline()
				idcol = round(mapply(renorm_interval,vcol,MoreArgs=list(b=c(1,300),inv=T)),0)
				cols=cols[idcol]#;print(cols)
			}else{
				cols=vcolor
			}
			fcolor[sid]="black"
			visualize_super_cluster(db[[pdbid2]],i=clusid2,typeg="ga",resize=1,coordi=auxzz,cols=cols,xlim=xlim,ylim=ylim,sub=sub[i],rot=rot,fcolor=fcolor)
			#print(nclus);print(id);print(dbhot[[id]]);readline()
			
		}
		#visualize_super_cluster(db[[pdbid2]],i=clusid2,typeg="ga",resize=1,coordi=auxzz,cols=cols,xlim=xlim,ylim=ylim,sub=sub[i],rot=rot)
		#;readline()
		#i=i+1
		#print("o")
	}


}


set_spot_score = function(tam,hot,col_type=2,n=300){

	#print(tam);
	if (col_type==1) cols = brewer.pal(n = 9, name = 'YlOrRd')
	else cols = colorRampPalette(c("brown","red", "yellow", "white"))(n = n)#;print(col)
	#tam = length(conn)
	res = list()#;print(tam)#;print(hot)
	ids = mapply(dim,hot)[1,]#;print(ids)
	for (i in 1:tam){
		res[[i]]=list()
		j = which(ids %in% i)	
		#if (length(j)){
		if (length(j)==1){
			res[[i]]$score = hot[[j]]$score
		}else{
			res[[i]]$score = rep(0,i)
		}
		#print(res[[i]]$score)
		idcol = round(mapply(renorm_interval,res[[i]]$score,MoreArgs=list(b=c(1,n),inv=T)),0)#;print(idcol)
		res[[i]]$colors = cols[idcol]
		#;print(res);readline()
	}
	return(res)

}


renorm_interval = function(x,b=c(1,9),a=c(0,1),inv=F){

	y = b[2]-( ((b[2]-b[1])*(a[2]-x))/(a[2]-a[1]) )

	if (inv){
		return( (b[1]+b[2])-y )
	}else{
		return(y)
	}

}

visualize_super_cluster = function(db,i,typeg="ga",norm=F,dc=600,sizef=10,k=0,resize=0.6,coordi=NULL,cols=NULL,xlim=NULL,ylim=NULL,vtype=2,sub="",rot=NULL,fcolor="black"){
#print("0k")
	#print(i)
	aux = list()
	aux$a = db$exp$superclus[[i]]$a
	aux$n = db$matrixname
	auxd = diag(aux$a)
	if (length(db$exp$preclus)!=0){
		aux$pre = list()
		prelen = length(db$exp$preclus$ids)
		aux$pre$d = round(sizef*sqrt(diag(db$exp$preclus$a)),k)
		aux$pre$g = graph.adjacency(diag(prelen),weighted=NULL,mode="undirected")
		aux$pre$xyz = tomatrix(transform_by_rot_matrix(db$exp$preclus$geomc,rot))
		if (vtype==1){
			aux$pre$vlabel = diag(db$exp$preclus$a)
		}
		else {
			aux$pre$vlabel = (length(auxd)+1):(length(auxd)+prelen)
		}
		#print(aux$pre);readline() 
	}
	
	if (typeg=="ga"){
		aux$g = db$exp$superclus[[i]]$ga
	}
	if (typeg=="gd"){
		aux$g = db$exp$superclus[[i]]$gd
	}
	#auxd = diag(aux$a)
	#print(coordi)
	if(is.null(coordi)){
		auxyz = db$exp$superclus[[i]]$geomc
	}else{
		auxyz = coordi
	}
	if (norm){
		#auxyz = auxyz/max(auxyz)
		coord = auxyz - mean(auxyz)
		coord = coord/max(coord)
		#print(coord)
		coord = coord*dc
		#print(coord)
	}else{
		coord = auxyz
	}
	#diag(auxm) = 0
	#manifold2graph = function(vx,vy,inv=1,sc=1,tz=1,dv=2,tzlim=c(1200,400),dimx=1200,dimy=800){
	#print(auxyz[,1:2])
	#coord = manifold2graph(coord[,1],coord[,2],dv=1,tzlim=c(1000,68),dimx=1000,dimy=136)#;print("ok")
	#print(coord)
	if(is.null(cols)){
		vcolor = "white"
	}else{
		vcolor = cols
	}
	if (vtype==1) {
		vlabel = round(auxd,k)
	}
	else {
		vlabel = 1:length(auxd)
	}
	#print(auxyz)
	#print(coord)
	#print(E(aux$g)$weight)
	if (typeg=="ga"){
		elabel = E(aux$g)$weight
	}
	if (typeg=="gd"){
		elabel = E(aux$g)$weight
	}
	if (!is.null(elabel)){
		elabel = round(elabel,k)
	}
	#if(is.null(coordi)){
	#	coordi = coord
	#}#else{
	#	vcolor = cols
	#}
	#elabel = round(E(aux$g)$weight,k)
	#elabel = NULL
	#print(coordi)
	#print(vlabel)
	#print(aux)
	#plot(coord[,1:2],type="p",col="blue",xlim=xlim,ylim=ylim);readline()
	#print(vcolor);readline()
	r=graph_plot(aux,labcex=1.0,size=sizef*sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,coordi=coord, ewidth=0.2,ecolor="gray75",lcolor="yellow",xlim=xlim,ylim=ylim,sub=sub,fcolor=fcolor)
	#print("0k")
	#r=graph_plot(aux,labcex=0.8,size=sqrt(auxd),vcolor=vcolor,vlabel=vlabel,elabel=elabel,ewidth=0.2,ecolor="gray75",lcolor="yellow",coordi=coordi)
	#abclines3d(0, 0, 0, a = diag(3), col = "gray")
	#tkplot.setcoords(r,resize*coord[,1:2])
	#tkconfigure(igraph:::.tkplot.get(r)$canvas, "bg"="white")
     #tkplot.center(r)
	#print(tkplot.getcoords(r))
	#print(r)
	#abclines3d(0, 0, 0, a = diag(3), col = "gray")
	#print(aux$a)

}

# Plota grafo
graph_plot = function(m,gdim=c(1200,800),labcex=1,labcole="salmon",labcoli="lightblue",glayout="layout.fruchterman.reingold",size=16, vcolor="white",vlabel=m$colname,elabel=NULL,ewidth=1,ecolor="darkgray",lcolor="white",fcolor="black",coordi=NULL,xlim=NULL,ylim=NULL,sub=""){

	#print(m);readline()
	#id=tkplot(m$g,vertex.label=paste(m$id,m$name,sep=":"),vertex.label.cex=labcex,vertex.color=colv(m$name,cole=labcole,coli=labcoli),canvas.width=gdim[1],canvas.height=gdim[2],vertex.frame.color=colvap(m$ap))
	#id=tkplot(m$g,label.color=lcolor,vertex.label=vlabel,vertex.label.cex=labcex,vertex.color=vcolor,canvas.width=gdim[1],canvas.height=gdim[2],vertex.frame.color="black",vertex.size=size,edge.width=ewidth,edge.color=ecolor,edge.label=elabel)
	if (is.null(coordi)){
		id=tkplot(m$g,edge.label.cex=labcex,vertex.label=vlabel,vertex.label.cex=labcex,vertex.color=vcolor,canvas.width=gdim[1],canvas.height=gdim[2],vertex.frame.color="black",vertex.size=size,edge.width=ewidth,edge.color=ecolor,edge.label=elabel)
		tkplot.center(r)
		return(id)
	}else{#print("ok")
	#id=tkplot(m$g,edge.label.cex=labcex,vertex.label=vlabel,vertex.label.cex=labcex,vertex.color=vcolor,canvas.width=gdim[1],canvas.height=gdim[2],vertex.frame.color="black",vertex.size=size,edge.width=ewidth,edge.color=ecolor,edge.label=elabel)
	#tkplot.setcoords(id,coordi)
		ewidth=2
		#print(size)
		#size = 0.7*size
		#plot(m$pre$g);readline()
		#print(vcolor);print("o")
		plot(m$g,edge.label.cex=labcex,vertex.label=vlabel,vertex.label.cex=labcex,vertex.color=vcolor,vertex.frame.color=fcolor,
vertex.size=size,edge.width=ewidth,edge.color=ecolor,edge.label=elabel,layout=coordi,vertex.label.color="black",main=m$n,xlim=xlim,ylim=ylim,rescale=F,axes=F)#;readline()
		if (!is.null(m$pre)){#print("ok")#;print(coordi)
			vcolor = "white"
			plot(m$pre$g,edge.label.cex=labcex,vertex.label=m$pre$vlabel,vertex.label.cex=labcex,vertex.color=vcolor,vertex.frame.color=fcolor,
vertex.size=m$pre$d,layout=m$pre$xyz,vertex.label.color="black",xlim=xlim,ylim=ylim,rescale=F,axes=T,add=T)
		}
		mtext(sub,line=0,cex=1)
		#readline()
	}	
}

#	ppf5iAT$tk$plot1$lines = make_lines(linxy=c(60,600),linstep=c(30,-30,200,30),linum=c(29,12),vum=diag(ppf5iAT$Lw$a),no=c(1100,700,0,-30))
coord_graph = function(diaglog,tr,tc,segs,xo,yo,sx,sy,sdx,sdy,xf=1100,yf=700,sxf=0,syf=-30){

	coord=matrix(rep(0,tr*tc),ncol=tc,nrow=tr)
	nozero = diaglog
	zero = !nozero
	tamsegs = length(segs)

	for (i in 1:(tamsegs-1)){

		seg = segs[[i]]
		coord = make_line(coord,xo,yo,sx,sy,seg,c(1:tr))
		print(coord)
		#readline()
		xo = xo + sdx
		yo = yo + sdy
	}

	return(coord)
	

}

visualize_super_cluster_space = function(dball,v=1:3,f=NULL,sep="_",glabel=1,coltype=NULL){

	
	if (is.null(f)){
		x = dball$A$svd$us[,v[1]]
		y = dball$A$svd$us[,v[2]]
		z = dball$A$svd$us[,v[3]]
		rown = dball$rownames
	}else{
		x = dball$A$svd$us[f,v[1]]
		y = dball$A$svd$us[f,v[2]]
		z = dball$A$svd$us[f,v[3]]
		rown = dball$rownames[f]
	}
	
	rown = unlist(strsplit(rown,sep))
	auxr = c()
	auxi = c()
	#print(rown);readline()	
	#for (i in 1:tam){
	tam = length(rown)
	i=1
	while(i<=tam){
		if (glabel==1){
			auxr = c(auxr,paste(rown[i],rown[i+2],sep=sep))
		}
		if (glabel==2){
			auxr = c(auxr,paste(rown[i]))
		}
		auxi = c(auxi,rown[i])
		i=i+3
	}
	#print(auxi);readline()
	if (is.null(coltype)){
		auxi = as.factor(auxi)
		col = red2blue_v1(length(levels(auxi)))
		levels(auxi) = col
		col = as.character(auxi)
	}
	if (length(coltype)>0){
		auxi = as.factor(coltype)
		col = red2blue_v1(length(levels(auxi)))
		levels(auxi) = col
		col = as.character(auxi)
	}
	#print(auxi);readline()
	#col = c(rep("blue",6),rep("lightblue",4),rep("red",6))
	#text = c(paste("p",2:7,sep=""),paste("r",2:5,sep=""),paste("a",8:13,sep=""))
	text = auxr
	#plot(x,y,type="p")
	plot3d(x,y,z,size=1,col=col,aspect="iso",box=FALSE,type="p")
	text3d(x,y,z,text=text,col=col)

}






# transforma uma reta coordenada em outra
transcoord = function(v,minv,maxv,maxdim){

	tam = length(v)

	auxv = v

	if ((maxv-minv)==0) {
		print("WARNING: possible division by zero")
		return(v)
	}
	for (i in 1:tam){

		auxv[i] = (maxdim*(v[i]-minv)+(maxv-v[i]))/(maxv-minv)
		auxv[i] = round(auxv[i],0)
		#print(v[i]);print(minv);print(auxv[i]);readline()
	

	}
	return(auxv)

}


transzero = function(vxy,coord,tzlim){
	
	tamr = dim(coord)[1]

	for (i in 1:tamr){

		if ((vxy[i,1]==0)&(vxy[i,2]==0)){
			coord[i,1] = tzlim[1]
			coord[i,2] = tzlim[2]
		}

	}

	return(coord)

}

invert_vector = function(vx,vy){

	auxx = vx
	auxy = vy

	tamr = length(vx)

	for (i in 1:tamr){

		if (auxx[i]>0){

			auxy[i] = -auxy[i]
		}
	}
	return(cbind(auxx,auxy))

}


# Reflete o plano manifold na visualização do grafo
manifold2graph = function(vx,vy,inv=1,sc=1,tz=1,dv=2,tzlim=c(1200,400),dimx=1200,dimy=800){

	#dimx=1200
	#dimy=800
	#dv = 2
	#sc = 1.8
	#sc = 1
	if (inv){
		vxy = invert_vector(vx,vy)
		vx = vxy[,1]
		vy = vxy[,2]
	}
	v1x = transcoord(vx,dv*min(vx),dv*max(vx),dimx)
	v1y = transcoord(vy,dv*min(vy),dv*max(vy),dimy)
	v1y = sc*v1y
	v1y = dimy-v1y
	v1y = round(v1y,0)
	coord = cbind(v1x,v1y)
	if (tz) coord = transzero(cbind(vx,vy),coord,tzlim)
	#coord = cbind(rep(1200,length(vx)),rep(800,length(vx)))
	#print(coord)
	#print(cbind(vx,vy,v1x,v1y))
	#ppf5z$tk$plot1$lines2 = coord
	#tkplot.setcoords(ppf5z$tk$plot1$id,ppf5z$tk$plot1$lines2)
	#tkplot.center(ppf5z$tk$plot1$id)
	return(coord)


}

manifold_view = function(m,v,lap,dec,main.text="",sc,x.lim=c(0,1200),y.lim = c(0,800),tz=0,inv=0){

	coord=c()
	r=list()
	v1 = v[1]
	v2 = v[2]
	#x.lim = c(0,1200)
	#y.lim = c(0,800)
	if (lap=="A"){
		if (dec=="eig"){
			vx = m$A$eig$v[,v1]
			vy = m$A$eig$v[,v2]

		}
		if (dec=="svd"){
			vx = m$A$svd$u[,v1]
			vy = m$A$svd$u[,v2]
		}
	
	}
	if (lap=="L"){
		if (dec=="eig"){
			vx = m$L$eig$v[,v1]
			vy = m$L$eig$v[,v2]

		}
		if (dec=="svd"){
			vx = m$L$svd$u[,v1]
			vy = m$L$svd$u[,v2]
		}		
	}
	if (lap=="Lw"){
		if (dec=="eig"){
			vx = m$Lw$eig$v[,v1]
			vy = m$Lw$eig$v[,v2]
		}
		if (dec=="svd"){
			vx = m$Lw$svd$u[,v1]#;print(vx);readline()
			vy = m$Lw$svd$u[,v2]
		}
	}
	maxus = max(vx)
	minus = min(vx)
	maxus1 = max(vy)
	minus1 = min(vy)
	#print(sc)
	#print(c(minus,maxus,minus1,maxus1))
	#print(cbind(vx,vy))
	#print("ok")
	r=manifold_travel(m,lap,autov=c(v1,v2),con=0,dec,scale=sc,psvd=c(0,0,0,0,1),peig=c(0,0,0,0,1),zoom=c(minus,maxus,minus1,maxus1),main.text)
	#print("ok")
	# view graph
	if (1){
		#vx = m$Lw$eig$v[,v1]
		#vy = m$Lw$eig$v[,v2]
		coord = manifold2graph(vx,vy,inv,sc=1,tz);print(coord)
		#ppf5z$tk$plot1$lines2 = coord
		#tkplot.setcoords(ppf5z$tk$plot1$id,ppf5z$tk$plot1$lines2)
		tkplot.setcoords(m$tk$plot1$id,coord)
		tkplot.center(m$tk$plot1$id)
	}
	return(coord)
}

inside_vector = function(m,idn,vcol,vn,breakm=126,dy=c(1.1,0.05),text.cex=0.60,main.text){

	#tamus = length(m$Lw$svd$us[,vn])
	tamus = length(m[,vn])

	maxus = max(m[,vn])
	minus = min(m[,vn])
	
	x1 = 1
	x2 = breakm

	lin = ceiling(tamus/breakm)
	#print(lin)
	#print(tamus)
	par(mfrow=c(lin,1))
	par(mar=c(1,4,2,2))

	my = m[,vn]

	while (x1 <= tamus){

		#y = my[x1:x2,vn]
		y = my[x1:x2]

		rot = idn[x1:x2]
		mtext(main.text,3,0,cex=1.2)
		#x=barplot(y,ylim=c(minus*dy[1],maxus*dy[1]),col=vcol[x1:x2])#,main=paste("Inside Vector ",vn))#,mar=c(0,4,4,2))
		x=barplot(y,ylim=c(minus-dy[1]*dy[2],maxus+dy[1]*dy[2]),col=vcol[x1:x2])#,main=paste("Inside Vector ",vn))#,mar=c(0,4,4,2))
		ddy = c()

		for (i in 1:length(y)){
			if (y[i]>0) ddy[i] = y[i]+dy[2]
			else ddy[i] = y[i]-dy[2]
		}
		#idres = paste(id,res,sep="")
		text(x,ddy,labels=rot,cex=text.cex,srt=90)
		x1 = x2+1
		aux = x2+breakm
		if ((aux)<=tamus) x2 = aux
		else x2 = tamus

	}
}

inside_vectors = function(m,d,tam,vfrow,ncol,main.text,main.aux,sub.text,plotype,scale,xlab.text="",ylab.text=xlab.text,r=2,pos=4){

	#tami = dim(m)[1]
	#tamj = dim(m)[2]
	maxij = scale[1]*abs(max(m[,tam]))
	tamm = dim(m)[2]

	par(mfrow=vfrow)
	j = 1
	pi = "("
	po = ")"


	for (i in tam){

		k = tamm-i+1

		if (plotype==1){
			plot(m[,i],type="l",ylim=c(-maxij,maxij),main=paste(main.text,main.aux,i),col=ncol,xlab=xlab.text,ylab=ylab.text)
			mtext(paste(sub.text,round(d[i],r)),cex=0.6)
		}
		if (plotype==2){
			plot(1,xlim=c(0,scale[1]*tamm),ylim=c(-maxij,maxij),type="n",main=paste(main.text,main.aux,i),xlab=xlab.text,ylab=ylab.text)
			for (j in 1:tamm){
				arrows(j,0,j,m[j,i],length=0,col=ncol[j])
			}
			mtext(paste(sub.text,round(d[i],r)),cex=0.6)
		}
		if (plotype==3){
			#scale=1.25
			#maxij = scale*abs(max(m))
			#dy = scale*maxij
			#dy = 1.1*maxij
			dy = scale[2]*maxij
			#k = tamm-i+1
			#pi = "("
			#po = ")"
			plot(1,xlim=c(0,scale[1]*tamm),ylim=c(-maxij-dy,maxij),type="n",main=paste(main.text,main.aux,pi,i,k,po),xlab=xlab.text,ylab=ylab.text)			
			for (j in 1:tamm){
				arrows(j,0,j,m[j,i],length=0,col=ncol[j])
			}
			text(tamm+1,0,paste(main.aux,i),pos=pos,cex=0.8)
			for (j in 1:tamm){
				arrows(j,0-dy,j,m[j,k]-dy,length=0,col=ncol[j])
			}
			text(tamm+1,0-dy,paste(main.aux,k),pos=pos,cex=0.8)
			mtext(paste(sub.text,round(d[i],r),round(d[k],r)),cex=0.6)
		}
		if (plotype==4){
			#k = tamm-i+1
			plot(1,xlim=c(-maxij,maxij),ylim=c(-maxij,maxij),type="n",main=paste(main.text,main.aux,pi,i,k,po),xlab=paste(main.aux,i),ylab=paste(main.aux,k))
			for (j in 1:tamm){
				points(m[j,i],m[j,k],col=ncol[j],pch=1)
			}
			mtext(paste(sub.text,round(d[i],r),round(d[k],r)),cex=0.6)

		}

	}

}


#reduz a matrix
reduce_matrix = function(m,k,dec){

	if (dec=="EIG"){
		tamr = dim(m$v)[1]
		tamc = dim(m$v)[2]

		vin = solve(m$v)
		aux = m$v[1:tamr,1:k] %*% as.matrix(m$s[1:k,1:k]) %*% vin[1:k,1:tamr] 
	}

	if (dec=="SVD"){
		tamr = dim(m$u)[1]
		tamc = dim(m$u)[2]

		vin = t(m$v)
		aux = m$u[1:tamr,1:k] %*% as.matrix(m$s[1:k,1:k]) %*% vin[1:k,1:tamr] 
	}

	return(aux)
	
}

compute_error = function(m,dec,mdiv,tam,p,main.text,legxy,limy,labx="1:k",laby="Standard Deviation"){
	
	#m=ppf5z$Lw$eig
	#tam = dim(m)[1]
	auxsd = c()
	auxm = c()
	x1 = 1:mdiv
	x2 = (mdiv+1):tam
	for (i in 1:tam){
		if (dec=="EIG") x = round(reduce_matrix(m$eig,i,dec),p)
		if (dec=="SVD") x = round(reduce_matrix(m$svd,i,dec),p)
		mr = x - m$a
		mr1 = sd(mr[x1,x1])
		mr2 = sd(mr[x2,x1])
		mr3 = sd(mr[x1,x2])
		mr4 = sd(mr[x2,x2])
		#xr = round((x-m),5)
		auxsd = c(mr1,mr2,mr3,mr4)
		auxm = rbind(auxm,auxsd)
	}
	#return(auxm)
	if (limy[2]==0) limy=c(0,max(auxm))
	plot(auxm[,1],ylim=limy,type="l",main=paste(main.text,dec),xlab=labx,ylab=laby)
	legend(legxy[1],legxy[2],legend=1:4,lty=1:4)
	lines(auxm[,2],type="l",lty=2)
	lines(auxm[,3],type="l",lty=3)
	lines(auxm[,4],type="l",lty=4)
	return(auxm)
}

# GENERAL: Visualiza matriz m como um degrade de cores correspondente aos valores aij
#view_matrix_general = function(m,main.text="",tcol="BLACK",x.lab="",y.lab="",pcex=2.6,tcex=0.8,p.ch=22,ydimlim=100){
view_matrix_general = function(m,rownames,pcex=2.6,tcex=0.7,p.ch=22,ydimlim=10000,colbyc=F){

	#print(dim(m)[2]-1);readline()
	aux=m[,1:(dim(m)[2]-1)]
	if(!is.matrix(aux)){
		aux=t(as.matrix(aux))
	}
	#print(aux);readline()
	xdim = dim(aux)[2]
	ydim = dim(aux)[1]
	#if(!is.matrix(aux)){
		#print(1)
	#	aux = t(is.matrix(aux))
		#print(aux)
	#}

	ncol=max(aux)-min(aux)+1
	#ccol=rainbow(ncol)
	ccol=terrain.colors(ncol)

	#print(dim(aux))
	#xdim = dim(aux)[2]
	#ydim = dim(aux)[1]
	#plot(1,1,type="n",xlim=c(1,xdim),ylim=c(1,ydim),xlab=x.lab,ylab=y.lab)
	#title(main.text,col.main=tcol)
	#print(ydim)
	if (ydim<ydimlim){
		for (i in 1:xdim){
			if (colbyc){
				ncol=max(aux[,i])-min(aux[,i])+1
				#ccol=rainbow(ncol)
				ccol=terrain.colors(ncol)
			}
			for(j in 1:ydim){
				#points(i,j,bg=ccol[max(aux)-aux[j,i]+1],pch=p.ch,cex=pcex)
				if (!colbyc){
					points(i,j,bg=ccol[max(aux)-aux[j,i]+1],col=ccol[max(aux)-aux[j,i]+1],pch=p.ch,cex=pcex)
				} else {
					points(i,j,bg=ccol[max(aux[,i])-aux[j,i]+1],col=ccol[max(aux)-aux[j,i]+1],pch=p.ch,cex=pcex)
				}
				#points(i,j,bg=ccol[max(aux)-aux[j,i]+1],pch=".",cex=pcex)
				#text(i,j,aux[j,i],cex=tcex)
				#print(i);print(j);print(aux[i,j])	
			}
		}
	}
	i=xdim+1
	wcol="red"
	#auxcol = as.factor(somaum)
	#print(auxcol)
	#taml = length(levels(auxcol))
	#levels(auxcol)=red2blue_v1(taml)
	#wcol = as.factor(m[,i])
	#tamwcol = length(levels(wcol))
	#levels(wcol)=red2blue_v1(tamwcol)
	#print(wcol);readline()
	
	for(j in 1:ydim){
		wcol = color_pam_width_general(m[j,i])
		#points(i,j,bg=wcol,pch=p.ch,cex=pcex)
		points(i,j,bg=wcol,col=wcol,pch=p.ch,cex=pcex)
	}
	i=i+1
	#print(rownames(m));readline()
	#rname = rownames(m)
	rname = rownames
	for(j in 1:ydim){
		#wcol = color_pam_width_general(m[j,i])
		#points(i,j,bg=wcol,pch=as.character(m[j,i]),cex=pcex)
		text(i,j,rname[j],cex=tcex)
	}

}

#GENERAL: retorna fingerprint dado uma matriz

finger_print_matrix = function(m,v=0){

	tamcol = dim(m)[2]
	aux=c()
	#print(dim(m));readline()

	for (j in 1:tamcol){
		indv = m[,j]>v
		aux[j] = sum(indv)
	}
	return(aux)
}



#GENERAL: produz paleta de cores para silhouette width
 
 color_pam_width_general=function(w,mini=0.1,maxi=1.1,byi=0.1,pal=1){
 
	if (w<(-mini)){
		return("yellow")
	}
	if (w<mini){
		return("white")
	}
	interv = seq(mini,maxi,byi)
	ncol = length(interv)
	if (pal==1) {
		ccol = red2blue_v1(ncol)
	}
	indv = interv <= w
	id=sum(indv)
	if (id>0){
		return(ccol[id])
	}else{
		print("WARNING: no color found for pam width")
	}

 }

# Visualiza matriz m como um degrade de cores correspondente aos valores aij
view_matrix = function(m,pcex,minv,maxv,sym,pal,main.text="",x.lab="",y.lab=""){

	tamr = dim(m)[1]
	tamc = dim(m)[2]

	
	if (sym){
		if (abs(maxv)>abs(minv)) minv = -maxv
		else maxv = -minv
	}

	#plot(1:tamr,1:tamc,type="n",main=main.text)
	plot(1,1,type="n",main=main.text,xlim=c(1,tamc),ylim=c(1,tamr),xlab=x.lab,ylab=y.lab)
	#print("ok")
	#ncol = tamr*tamc
	ncol = tamr
	if (pal==0) vcol=red2blue_v1(ncol)
	if (pal==1) vcol=cm.colors(ncol)
	if (pal==2) vcol=red2blue_v2(ncol)
	if ((maxv-minv)==0) print("WARNING: possible division by zero")
	for (i in 1:tamr){
		for (j in 1:tamc){
			#aux = (m[i,j]+abs(minv))/(maxv+abs(minv))
			#icol = ((ncol-1)*(m[i,j]+abs(minv)))/(maxv+abs(minv))+1
			icol = (ncol*(m[i,j]-minv)+(maxv-m[i,j]))/(maxv-minv)
			icol = round(icol,0)
			if (icol>ncol) icol=ncol
	
			#points(i,j,pch=22,bg=vcol[icol],col=vcol[icol],cex=pcex)
			points(j,tamr-i+1,pch=22,bg=vcol[icol],col=vcol[icol],cex=pcex)
			#readline()
		}

	}
	

}

segment_vectors = function(id,v1,v2){

	v1p = id[v1>0]
	v1z = id[v1==0]
	v1n = id[v1<0]

	v2p = id[v2>0]
	v2z = id[v2==0]
	v2n = id[v2<0]

	segs = list (c(v1p),c(v2n),c(v1n),c(v2p),c(v1z,v2z))

	return(segs)

}

view_matrix_diff = function(ma,mdr,err,pcex,mark,col="GRAY"){

	tamr = dim(ma)[1]
	tamc = dim(ma)[2]

	for (i in 1:tamr){
		for (j in 1:tamc){
			if (abs(mdr[i,j])<=err){
				#points(i,j,pch=22,fg=col,cex=pcex)
				points(i,j,pch=mark,fg=col,cex=pcex)
			}
		}

	}


}

#Visualização intergrada das matrizes
compare_matrix = function(m,vn,lap,dec,main.text,p=5,pcex=2,err=0.001,mark=""){


	if (lap=="L"){
		ma = m$L$a
		if (dec=="EIG"){
			mr = m$L$eig
		}
		if (dec=="SVD"){
			mr = m$L$svd
		}		
	}
	if (lap=="Lw"){
		ma = m$Lw$a
		if (dec=="EIG"){
			mr = m$Lw$eig
		}
		if (dec=="SVD"){
			mr = m$Lw$svd
		}
	}

		#p=5
		#dec = "EIG"
		#lap = "Lw"
		x = round(reduce_matrix(mr,vn,dec),p)
		a = round(ma,p)
		xr = round((x-a),p)
		#pcex = 2
		par(mfrow=c(1,3))
		#view_matrix(a,vn,pcex,min(a),max(a),sym=1,pal=2,paste(main.text,lap,dec,"all"))
		view_matrix(a,pcex,min(a),max(a),sym=1,pal=2,paste(main.text,lap,dec,"all"))
		view_matrix_diff(a,xr,err,pcex,mark)
		#view_matrix(x,vn,pcex,min(a),max(a),sym=1,pal=2,paste(main.text,lap,dec,"k=",vn))
		view_matrix(x,pcex,min(a),max(a),sym=1,pal=2,paste(main.text,lap,dec,"k=",vn))
		view_matrix_diff(a,xr,err,pcex,mark)
		#view_matrix(xr,vn,pcex,min(a),max(a),sym=1,pal=2,paste(main.text,lap,dec,"k=",vn,"diff"))
		view_matrix(xr,pcex,min(a),max(a),sym=1,pal=2,paste(main.text,lap,dec,"k=",vn,"diff"))
		view_matrix_diff(a,xr,err,pcex,mark)


}

#retorna a norma dos vetores por linhas da matriz m
norm_vector = function(m,ntype="F"){
	
	tam = dim(m)[1]
	aux = c()

	for (i in 1:tam){
		aux= c(aux,norm(as.matrix(m[i,]),type=ntype))
	}
	return(aux)

}

# retorna angulo entre v1 e v2
angle = function(v1,v2,unit="degree",p=10){

	b = 180
	norms = (norm(as.matrix(v1),type="F")*norm(as.matrix(v2),type="F"))
	if (norms) {
		a = sum((v1*v2))/norms
	}
	else {
		print("WARNING: one norm = 0, angle = INF; but doing angle = 0")
		a = 0
	}
	a = round(a,p)
	#print(a)
	#print(acos(a))
	if (unit=="degree"){
		return(acos(a)*b/pi)
	}
	if (unit=="rad"){
		return(acos(a))
	}
}

# retorna norm e angulo em relacao a um vetor referencia
anglenorm = function(m,vref,unit="rad",p=10){

	ad = c()
	x = dim(m)[1]
	for (i in 1:x){
		#print(m[i,]);print(vref)
		a = angle(m[i,],vref,unit)
		n = norm(as.matrix(m[i,]),type="F")
		n = round(n,p)
		ad = rbind(ad,c(a,n))
		#print(ad);readline()
	}
	return(ad)
}

# função mais eficiente para anglenorm
anglenorm2 = function(m,vref,rc=1,unit="rad",p=10){

	a=apply(m,rc,angle,v2=vref,unit=unit)
	n=apply(m,rc,Norm)
	return(cbind(a,n))
}

#Transforma coordenada polar para cartesiana
pol2car = function(m,p=10){

	tam = dim(m)[1]

	aux = c()
	for (i in 1:tam){
		x = m[i,2]*cos(m[i,1])
		x = round(x,p)
		y = m[i,2]*sin(m[i,1])
		y = round(y,p)
		aux=rbind(aux,c(x,y))
	}
	return(aux)
}
# retorna um vetor perpendicular a v
perpendicular = function(v){

	tam = length(v)
	aux = rep(0,tam)
	
	for (i in 1:tam){
		#print(v[i])
		if (v[i]==0){
			aux[i]=1
			return(aux)
		}
	}
	print("WARNING: no perpendicular vector found")
	return(NA)
}

#projeta km centros sobre um vetor referencia vri
kmcenter2xy = function(km,vri){

	aux=c()
	tam = length(km$size)
	for (i in 1:tam){
		vr = km$center[i,]
		vr.xy = pol2car(anglenorm(t(as.matrix(vr)),vref=vri))
		aux = rbind(aux,vr.xy)
		#points(vr.xy,pch=csymb,col=colp)
	}
	return(aux)

}
# Tentativa de fazer um plot simétrico para os centroides e seus pontos
# Não deu muito certo...
kmcenterplan2xy = function(km,vri,vri2,p=10){

	aux=c()
	tam = length(km$size)
	for (i in 1:tam){
		a1 = angle(km$center[i,],vri,"rad")
		a2 = angle(km$center[i,],vri2,"rad")		
		x = Norm(km$center[i,])*cos(a1)
		y = Norm(km$center[i,])*cos(a2)
		x = round(x,p)
		y = round(y,p)
		vr.xy=c(x,y)
		#vr = km$center[i,]
		#vr.xy = pol2car(anglenorm(t(as.matrix(vr)),vref=vri))
		aux = rbind(aux,vr.xy)
		#points(vr.xy,pch=csymb,col=colp)
	}
	return(aux)

}

# Tentativa ligada a função anterior...
plan2car = function(m,vref1,vref2,p=10){

	tam = dim(m)[1]
	aux = c()
	a90 = pi/2
	a180 = pi 
	for (i in 1:tam){
		a1 = angle(m[i,],vref1,"rad")
		a2 = angle(m[i,],vref2,"rad")
		#print(a1);print(a2);readline()
		#if ((a1<=a90)&(a2<=a90)){
			x = Norm(m[i,])*cos(a1)
			y = Norm(m[i,])*cos(a2)
			x = round(x,p)
			y = round(y,p)
			aux = rbind(aux,c(x,y))
		#}
		#if ((a1<=a180)&(a1>a90)&(a2<=a90)){
		#	x = Norm(m[i,])*cos(a1)
		#	y = Norm(m[i,])*cos(a2)
		#	print("ok");print(x);print(y);readline()
		#	x = round(x,p)
		#	y = round(y,p)
		#	aux = rbind(aux,c(x,y))
		#}
		#if ((a1<=a180)&(a1>a90)&(a2<=a180)&(a2>a90)){
		#	x = -Norm(m[i,])*cos(a1)
		#	y = -Norm(m[i,])*cos(a2)
		#	x = round(x,p)
		#	y = round(y,p)
		#	aux = rbind(aux,c(x,y))
		#}
		#if ((a1<a90)&(a2<=a180)&(a2>a90)){
		#	x = Norm(m[i,])*cos(a1)
		#	y = -Norm(m[i,])*cos(a2)
		#	x = round(x,p)
		#	y = round(y,p)
		#	aux = rbind(aux,c(x,y))
		#}
		#x = Norm(m[i,])*cos(m[i,1])
		#x = round(x,p)
		#y = m[i,2]*sin(m[i,1])
		#y = round(y,p)
		#aux=rbind(aux,c(x,y))
	}
	return(aux)


}

# Um plot para visualizar clusters
#plotkm=plotkmcluster(u,ppf5zc$L$svd$km2u2,0,0,substr(ppf5zc$name,1,1),c("Real","Kmeans",paste("Laplacian - U",k,sep="")),c("comp_x vcenter ","comp_y vcenter "),ppf5zc$tk$plot1$id)
plotkmcluster = function(m,km,vi,plan,pattern,main.text,xylab,graphid,colp="red",csymb="c",pcex=0.8,sca=1.1,scb=3,mtl=0.3,bgp="grey60",ltys=1,lenp=0.1,p=5){

	plotkm = list()
	vr.xy=c()
	co=c(0,0)
	tam = length(km$size)
	if (vi){
		vri = km$center[vi,]
		xy = pol2car(anglenorm(m,vref=vri))
	}
	else {
		vrn = names(sort(apply(km$center,1,Norm),dec=T))
		vri = km$center[vrn[1],]
		#print(vri)
		#print(anglenorm(m,vref=vri))
		xy = pol2car(anglenorm(m,vref=vri))
	}
	#print(xy)
	
	if (plan){
		vri2 = perpendicular(vri)
		#print(vri2)
		#vri2 = rep(1,29)
		xy=plan2car(m,vref1=vri,vref2=vri2)
		#print(xy);readline()

	}

	yl = 	c(sca*min(xy),sca*max(xy))
	#print(yl)
	#yl = c(sca*min(xy[,2]),sca*max(xy[,2]))
	#xl = c(sca*min(xy[,1]),sca*max(xy[,1]))

	xl = yl
	#xl = c(-0.5,0.5)
	plotkm$xy = round(xy,p)
	par(mfrow=c(1,2))
	plot(xy,type="p",pch=pattern,ylim=yl,xlim=yl,main=main.text[1],xlab=paste(xylab[1],vi,sep=":"),ylab=paste(xylab[2],vi,sep=":"),cex=pcex)
	mtext(main.text[3],line=mtl)
	points(0,0,pch=21,cex=0.8,bg=bgp)
	abline(v=0,h=0,lty=3,col=bgp)
	if (plan){
		vr.xy = kmcenterplan2xy(km,vri,vri2)	
	} else {
		vr.xy = kmcenter2xy(km,vri)
	}
	#print(vr.xy);readline()
	plotkm$vr = round(vr.xy,p)
	for (i in 1:tam){
		text(vr.xy[i,1],vr.xy[i,2]-scb,labels=paste(csymb,i,sep=""),col=colp)
		arrows(0,0,vr.xy[i,1],vr.xy[i,2],lty=ltys,col=bgp,code=2,length=lenp)
	}
	plot(xy,type="p",pch=paste(km$cluster),ylim=yl,xlim=xl,main=paste(main.text[2],tam,sep="-"),xlab=paste(xylab[1],vi,sep=":"),ylab=paste(xylab[2],vi,sep=":"),cex=pcex)
	mtext(main.text[3],line=mtl)
	points(0,0,pch=21,cex=0.8,bg=bgp)
	abline(v=0,h=0,lty=3,col=bgp)
	for (i in 1:tam){
		text(vr.xy[i,1],vr.xy[i,2]-scb,labels=paste(csymb,i,sep=""),col=colp)
		arrows(0,0,vr.xy[i,1],vr.xy[i,2],lty=ltys,col=bgp,code=2,length=lenp)
	}
	#return(xy)
	if (graphid>0){
		# view graph
		#print("graph")
		vx = xy[,1]
		vy = xy[,2]
		coord = manifold2graph(vx,vy,0,1,0)
		#print(coord)
		tkplot.setcoords(graphid,coord)
		tkplot.center(graphid)
	}
	return(plotkm)
}

# Plot grafico com curva para numero de cluster em kmeans
howmanycluster = function(m,cint,imax,nseeds,main.text,yl=c(0,3000),verbose=1,maxyl=1){

	tam = dim(m)[1]
	# Determine number of clusters
	wss = (nrow(m)-1)*sum(apply(m,1,var))
	if (is.na(wss)) wss=0
	for (i in (cint[1]+1):cint[2]) {
		if (verbose) print(i)
		wss[i] = sum(kmeans(m,centers=i,iter.max=imax,nstart=nseeds)$withinss)
	}
	if (maxyl) yl=c(0,max(wss))
	#if (verbose) print(wss)
	plot(cint[1]:cint[2], wss, type="b", ylim=yl,main=main.text,xlab="Number of Clusters",ylab="Within groups sum of squares")
	#if (verbose) print(wss)

}

# Plot grafico para sondar correlação entre Norma de um vetor-linha e grau da linha
normaxgrau = function(m,diagL,symb,kmcluster,main.text,byrc,pcex=0.8,mtl=0.3,xla="Norm",yla="Degree"){


	auxn = apply(m,byrc,Norm)
	yli = c(0,1.1*max(diagL))
	xli = yli
	xli = c(0,1.1*max(auxn))
	#par(mfrow=c(1,2))
	plot(auxn,diagL,pch=symb,main=main.text[1],cex=pcex,xlab=xla,ylab=yla,xlim=xli,ylim=yli)
	mtext(main.text[3],line=mtl)
	#plot(auxn,diagL,pch=kmcluster,main=main.text[2],cex=pcex,xlab=xla,ylab=yla,xlim=xli,ylim=yli)
	#mtext(main.text[3],line=mtl)
}

##### FOURIER #####



sampledd = function(dd,win=10){

	i = 1
	aux = c()
	while((i+win)<=length(dd)){

		aux = c(aux,sample(dd[i:(i+win)],1))
		i = i+win
	}

	return(aux)

}


coding = function(dgt,type,Nu){

   N = length(dgt)
   s = c()
   if (type == "NRZ"){
      for (i in 1:N){
         if (dgt[i] == 1){
            s = c(s,rep(1,Nu))
         }
         if (dgt[i] == 0){
            s = c(s,rep(-1,Nu))
         }
         #print(s)
      }
   }
   if (type == "AM"){
      s = c(0)
      k = 10
      for (i in 1:N){
         if (dgt[i] == 1){
            tmp = sin(seq(0,2*Nu*pi,by=2*Nu*pi/(Nu*k)))
            s = c(s,tmp[2:(length(tmp))])
            #print(s)
         }
         if (dgt[i] == 0){
            tmp = 0.5*sin(seq(0,2*Nu*pi,by=2*Nu*pi/(Nu*k)))
            s = c(s,tmp[2:(length(tmp))])
            #print(s)
         }
      }
   }
   #print(s)
   return (s)
}


ifft=function(X,d=c(-1)){
   
   N = length(X)
   d = c(d,N-d)
   i = complex(i=1)
   x = c()
   for (n in 0:(N-1)){
      xn = 0
      for (k in 0:(N-1)){
         if (!(k %in% d))
            xn = xn + X[k+1]*exp(2*pi*n*i*k/N)
            #xn = xn + X[k+1]*(cos(2*pi*n*k/N)-i*sin(2*pi*n*k/N))
      }
      x[n+1] = xn

   }
   return(round(1/N*x,2))
}

matrixfft = function(D){


   X = fft(D)
   A = Mod(X)
   N=length(X)
   NuCom = N/2

   m = matrix(Re(ifft(X,c(1:NuCom))),nrow=1)
   s = matrix(m[1,],nrow=1)
   st = s

   j = 1

   while (j<=NuCom){#print("Ok")

      m = rbind(m,Re(ifft(X,c(0:(j-1),(j+1):(NuCom)))))
      st = st + m[j+1,]
      s = rbind(s,st)
      j = j + 1
   }

   mT = list(D=D,X=X,M=m,S=s,A=A)
   return(mT)
}

plotmatrixfft = function(T,NuCom,NuCyc,OnlyNZ,OverF,kyl,Ta,Bits,unit1,unit2){

   #library("gmp")
   amp = 0
   #cpnz = 0.2
   par(mar=c(1, 1, 1, 1))
   N=length(T$X)
   cpnz = 0.2
   lcpnz = 10
 #  if (OnlyNZ == 1){
 #     cpnz = 0.2
   k = 0

	#unit = "bit/s"
	#unit = "seg"
   for (i in 1:NuCom){
      if (sum(abs(T$M[i,(1:lcpnz)]))>cpnz)
         k = k + 1
     }
#print(k)
   if (OnlyNZ == 1){
      cpnz = 0.2
      par(mfrow=c(k+2,2))
   }
   else {
      cpnz = -1
      par(mfrow=c(NuCom+2,2))
   }

   #par(bg="BLACK",fg="YELLOW")
	par(bg="WHITE",fg="BLACK")
   lwdi = 2
   #yl = c(-1.5,1.6)
   yl = c(min(T$D)-kyl,max(T$D)+kyl)

   c = rainbow(NuCom+3)
   j = 1

	prec = 4

   #plot(rep(T$D,NuCyc),type="l",col="WHITE",ylim=yl,xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main="Funcao Original D(t)",col.main="WHITE")
	col1="WHITE"
	col1="BLACK"
   #plot(rep(T$D,NuCyc),type="l",col=col1,ylim=yl,xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main=paste("Funcao Original D(t): ",round((Bits*NuCyc)/Ta,2),unit),col.main=col1)
	plot(rep(T$D,NuCyc),type="l",col=col1,ylim=yl,xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main=paste("Funcao Original D(t): ",round(Ta,prec),unit1),col.main=col1)
   plot(rep(T$D,NuCyc),type="l",col=col1,ylim=yl,xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main="Funcao Original D(t)",col.main=col1)

   if (sum(abs(T$M[1,(1:lcpnz)]))>cpnz)
      plot(rep(T$M[1,],NuCyc),type="l",col=c[1],ylim=yl,xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main=paste("Harmonico","f",j-1,"=",0,unit2),col.main=c[1])
   while (j<NuCom){#print("Ok")

      if (sum(abs(T$M[j,(1:lcpnz)]))>cpnz){
         plot(rep(T$S[j,],NuCyc),col=c[j],ylim=yl,type="l",xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main=paste("F",j-1,"=",0,":",j-1),col.main=c[j])
         if (OverF == 1){
         
               #lines(rep(T$D,NuCyc),lty=3,col="WHITE")
					#lines(rep(T$D,NuCyc),lty=3,col="GREY")
					lines(rep(T$D,NuCyc),lty=3,col="BLACK")
         }
      }

      if (sum(abs(T$M[(j+1),(1:lcpnz)]))>cpnz){
         #plot(rep(T$M[(j+1),],NuCyc),col=c[j+1],ylim=yl,type="l",xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main=paste("Harmonico ", "f",j,"=",j,"/",N," Hz"),col.main=c[j+1])
         plot(rep(T$M[(j+1),],NuCyc),col=c[j+1],ylim=yl,type="l",xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main=paste("Harmonico ", "f",j,"=",round((j*NuCyc)/Ta,prec),unit2),col.main=c[j+1])
         if (OverF == 1){
#               lines(rep(T$D,NuCyc),lty=3,col="WHITE")
				lines(rep(T$D,NuCyc),lty=3,col="BLACK")
         }
      }     
      j = j + 1
   }

   if (sum(abs(T$M[j,1:lcpnz]))>cpnz){
      plot(rep(T$S[j,],NuCyc),col=c[j],ylim=yl,type="l",xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main=paste("F",j-1,"=",0,":",j-1),col.main=c[j])
      if (OverF == 1){#print("ok")
            #lines(rep(T$D,NuCyc),lty=3,col="WHITE")
				#lines(rep(T$D,NuCyc),lty=3,col="GREY")
				lines(rep(T$D,NuCyc),lty=3,col="BLACK")
      }

   }
	col2="YELLOW"
	col2="ORANGE"
   if (amp == 1){
      plot(T$A[1:(NuCom+1)],col=c[1:(NuCom+1)],type="h",xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main="Amplitudes Complexas",col.main=col2)
	}
   else
      plot(T$A[1:(N/2)],col=c[1:(NuCom+1)],type="h",xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main="Amplitudes Complexas",col.main=col2)


   plot(rep(T$D,NuCyc)-rep(T$S[j,],NuCyc),col="ORANGE",type="h",ylim=yl,xaxt="n",yaxt="n",xlab="",ylab="",lwd=lwdi,main="Residuos: D(t)-F(t)",col.main="ORANGE")

}

make_score_graph_superimposition_old2 = function(dist,g1,g2,a1,a2,w=c(1,1,1,1),consol=3,sep="-",infy=0.49){

	auxsc = list()
	auxc = names(dist)

	#print(dist)#;print(auxc);readline()

	#grau dos nos
	d1 = degree(g1) 
	d2 = degree(g2)

	# soma dos pesos das arestas de no divido pelo grau
	# TO DO: testar possibilidade de ser "soma"
	f1 = strength(g1)/d1
	f2 = strength(g2)/d2
	f1[is.nan(f1)]=0
	f2[is.nan(f2)]=0

	#print(strength(g1));print(strength(g2))
	#print(f1);print(f2)
	#print(d1);print(d2);readline()
	#print(f1/d1);print(f2/d2);readline()

	

	tam1 = length(d1)
	tam2 = length(d2)

	auxrd = c()
	auxra = c()
	auxrf = c()
	auxs = c()
	#auxsa = c()

	if (tam1<=tam2){ # sempre g1 < g2
		auxd1 = d1 # graus do grafo g1
		auxd2 = d2 # graus do grafo g2
		auxa1 = a1 # labels dos nos g1 
		auxa2 = a2 # labels dos nos g2
		auxf1 = f1 # labels dos vertices normalizados por grau de g1
		auxf2 = f2 # labels dos vertices normalizados por grau de g2
		auxg1 = g1 # grafo g1
		auxg2 = g2 # grafo g2
		tam = tam1 # quantidade de nos
	} else {
		auxd1 = d2
		auxd2 = d1
		auxa1 = a2
		auxa2 = a1
		auxf1 = f2
		auxf2 = f1
		auxg1 = g2
		auxg2 = g1
		tam = tam2
	}
	auxt1 = ends(auxg1,E(auxg1))
	auxt2 = ends(auxg2,E(auxg2))
	print(auxt1);print(auxt2)
	auxs = verify_edge_superimpostion(auxc,auxt1,auxt2)
	#auxe = 0
	par(frow=c(1,2))
	plot(g1);plot(g2)
	print(auxc);print(auxs)#;readline()

	#auxmaxd = dist/max(dist)
	#auxmaxd = dist/mean(dist)
	#auxmaxd = dist/median(dist)

	#auxmaxd = 1/(dist*dist) #usar infy=0.24 para dist minimo de 2.041
	auxmaxd = 1/(dist) #usar infy=0.49 para dist minimo 2.041
	#auxmaxd = 1/sqrt(dist)
	#auxmaxd[auxmaxd>infy]=infy
	auxmaxd[auxmaxd>infy]=1 #tetodist=1
	#print(auxmaxd)#;readline()

	for (i in 1:tam){
		#print(auxd1);print(auxd2);print(auxc)#;readline()
		#print(auxa1);print(auxa2);readline()
		auxi = as.numeric(unlist(strsplit(auxc[i],sep)))
		#print(auxi);readline()
		#auxdifa = abs(auxa1[auxi[1],auxi[1]]-auxa2[auxi[2],auxi[2]])
		#auxdifa = 1-min(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]])/max(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]]) #diff nos
		auxdifa = min(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]])/max(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]]) #diff nos
		#print(auxi);print(auxa1);print(auxa2);print(auxdifa);readline()
		#auxdifd = abs(auxd1[auxi[1]]-auxd2[auxi[2]])
		#auxdifd = 1-min(auxd1[auxi[1]],auxd2[auxi[2]])/max(auxd1[auxi[1]],auxd2[auxi[2]]) #diff graus
		auxdifd = min(auxd1[auxi[1]],auxd2[auxi[2]])/max(auxd1[auxi[1]],auxd2[auxi[2]]) #diff graus
		#auxdiff = 1-min(auxf1[auxi[1]],auxf2[auxi[2]])/max(auxf1[auxi[1]],auxf2[auxi[2]]) #diff pesos normalizados graus
		auxdiff = min(auxf1[auxi[1]],auxf2[auxi[2]])/max(auxf1[auxi[1]],auxf2[auxi[2]]) #diff pesos normalizados graus
		auxrd = c(auxrd,auxmaxd[i]*auxdifd)
		auxra = c(auxra,auxmaxd[i]*auxdifa)
		auxrf = c(auxrf,auxmaxd[i]*auxdiff)
		#print(auxra);readline()
		#auxs = c(auxs,w[1]*dist[i]+w[2]*auxdifd+w[3]*auxdifa+w[4]*auxdiff)
		auxs = c(auxs,w[1]*auxmaxd[i]+w[2]*auxdifd*auxmaxd[i]+w[3]*auxdifa*auxmaxd[i]+w[4]*auxdiff*auxmaxd[i])
	}
	
	#print(auxs);readline()
	#print(auxc);readline()
	auxsc$dist = dist # distancia euclideana nos
	auxsc$distn = auxmaxd # distancia elevada a 2	 
	auxsc$degv = auxrd # diff graus
	auxsc$nodv = auxra # diff nos (areas)
	auxsc$weiv = auxrf # diff pesos normalizados graus
	auxsc$scov = auxs # score final
	if (consol==1)	auxsc$scor = median(auxs)#+auxe
	if (consol==2)	auxsc$scor = mean(auxs)#+auxe
	if (consol==3)	auxsc$scor = sum(auxs)#+auxe ### dependente No de nós
	#auxsc$degs = sum(auxr)
	#print(auxsc);readline()
	return(auxsc) 
	

}

make_score_graph_superimposition_old1 = function(dist,g1,g2,a1,a2,w=c(1,1,1,1),consol=2,sep="-",type=1){

	auxsc = list()
	auxc = names(dist)

	#grau dos nos
	d1 = degree(g1) 
	d2 = degree(g2)

	# soma dos pesos das arestas de no divido pelo grau
	f1 = strength(g1)/d1
	f2 = strength(g2)/d2
	f1[is.nan(f1)]=0
	f2[is.nan(f2)]=0

	#print(strength(g1));print(strength(g2))
	#print(f1);print(f2)
	#print(d1);print(d2);readline()
	#print(f1/d1);print(f2/d2);readline()

	

	tam1 = length(d1)
	tam2 = length(d2)

	auxrd = c()
	auxra = c()
	auxrf = c()
	auxs = c()
	#auxsa = c()

	if (tam1<=tam2){
		auxd1 = d1 # graus do grafo g1
		auxd2 = d2 # graus do grafo g2
		auxa1 = a1 # labels dos nos g1 
		auxa2 = a2 # labels dos nos g2
		auxf1 = f1 # labels dos vertices normalizados por grau de g1
		auxf2 = f2 # labels dos vertices normalizados por grau de g2
		auxg1 = g1 # grafo g1
		auxg2 = g2 # grafo g2
		tam = tam1 # quantidade de nos
	} else {
		auxd1 = d2
		auxd2 = d1
		auxa1 = a2
		auxa2 = a1
		auxf1 = f2
		auxf2 = f1
		auxg1 = g2
		auxg2 = g1
		tam = tam2
	}
	#auxs = verify_edge_superimpostion(auxc,ends(auxg1,E(auxg1)),ends(auxg2,E(auxg2)))
	#auxe = 0
	#print();readline()

	#auxmaxd = dist/max(dist)
	#auxmaxd = dist/mean(dist)
	#auxmaxd = dist/median(dist)

	auxmaxd = (dist*dist)

	for (i in 1:tam){
		#print(auxd1);print(auxd2);print(auxc)#;readline()
		#print(auxa1);print(auxa2);readline()
		auxi = as.numeric(unlist(strsplit(auxc[i],sep)))
		#print(auxi);readline()
		#auxdifa = abs(auxa1[auxi[1],auxi[1]]-auxa2[auxi[2],auxi[2]])
		auxdifa = 1-min(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]])/max(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]]) #diff nos
		#auxdifd = abs(auxd1[auxi[1]]-auxd2[auxi[2]])
		auxdifd = 1-min(auxd1[auxi[1]],auxd2[auxi[2]])/max(auxd1[auxi[1]],auxd2[auxi[2]]) #diff graus
		auxdiff = 1-min(auxf1[auxi[1]],auxf2[auxi[2]])/max(auxf1[auxi[1]],auxf2[auxi[2]]) #diff pesos normalizados graus
		auxrd = c(auxrd,auxdifd)
		auxra = c(auxra,auxdifa)
		auxrf = c(auxrf,auxdiff)
		#print(auxra);readline()
		#auxs = c(auxs,w[1]*dist[i]+w[2]*auxdifd+w[3]*auxdifa+w[4]*auxdiff)
		auxs = c(auxs,w[1]*auxmaxd[i]+w[2]*auxdifd+w[3]*auxdifa+w[4]*auxdiff)
	}
	
	#print(auxs);readline()
	#print(auxc);readline()
	auxsc$dist = dist # distancia euclideana nos
	auxsc$distn = auxmaxd # distancia elevada a 2	 
	auxsc$degv = auxrd # diff graus
	auxsc$nodv = auxra # diff nos (areas)
	auxsc$weiv = auxrf # diff pesos normalizados graus
	auxsc$scov = auxs # score final
	if (consol==1)	auxsc$scor = median(auxs)#+auxe
	if (consol==2)	auxsc$scor = mean(auxs)#+auxe
	if (consol==3)	auxsc$scor = sum(auxs)#+auxe ### dependente No de nós
	#auxsc$degs = sum(auxr)
	#print(auxsc);readline()
	return(auxsc) 
	

}

#align_graph_pdb = function(command,args,idline=c(35,37,52)){
#align_graph_pdb = function(command,args,idline=c(31,33,48)){
align_graph_pdb_old = function(command,flags,script,args,idline=c(32,34,49),dline=19,v1=".pml",v2="a.pml"){

	auxrmsd = c()
	auxwin = c()
	auxrot = c()
	auxlist = list()
	
	#print(command);readline()
	#print(args);readline()
	args = c(flags,paste(script,v1,sep=""),args)
	#print(args);readline()
	auxout = system2(command,args=args,stdout=T)
	auxline = unlist(strsplit(auxout[idline[1]]," "))
	print(auxout)#;readline()
	#print(auxline);readline()
	if (auxline[1]=="RMSD"){
		auxnum = as.numeric(auxline[2])
		auxrot = auxout[idline[2]:idline[3]]
		auxrot = str_c(auxrot,collapse="")
		#print(nchar(auxrot))
		auxrot = str_sub(auxrot,2,(nchar(auxrot)-1))
		auxrot = as.numeric(str_split(auxrot,",")[[1]])
		#print(auxrot);readline()
		#auxrotlist = c(auxrotlist,auxrot) 
		#auxrotlist[[k]] = auxrot;k=k+1
		auxlist$rmsd = auxnum
		auxlist$rot = auxrot
	}else{
		#print(paste("WARNING: RMSD line not found in cealign output for",args))
		auxlist$rmsd = 0
		auxlist$rot = 0	
	}
	return(auxlist)

}

make_score_graph_superimposition_old3 = function(dist,g1,g2,a1,a2,w=c(1,1,1,1),consol=3,sep="-"){

	auxsc = list()
	auxc = names(dist)

	#print(dist)#;print(auxc);readline()

	#grau dos nos
	d1 = degree(g1) 
	d2 = degree(g2)

	# soma dos pesos das arestas de no divido pelo grau
	# TO DO: testar possibilidade de ser "soma"
	f1 = strength(g1)/d1
	f2 = strength(g2)/d2
	f1[is.nan(f1)]=0
	f2[is.nan(f2)]=0

	#print(strength(g1));print(strength(g2))
	#print(f1);print(f2)
	#print(d1);print(d2);readline()
	#print(f1/d1);print(f2/d2);readline()

	

	tam1 = length(d1)
	tam2 = length(d2)

	auxrd = c()
	auxra = c()
	auxrf = c()
	auxre = c()
	auxs = c()
	auxa = c()
	auxe = c()
	auxd = c()
	#auxsa = c()

	if (tam1<=tam2){ # sempre g1 < g2
		auxd1 = d1 # graus do grafo g1
		auxd2 = d2 # graus do grafo g2
		auxa1 = a1 # labels dos nos g1 
		auxa2 = a2 # labels dos nos g2
		auxf1 = f1 # labels dos vertices normalizados por grau de g1
		auxf2 = f2 # labels dos vertices normalizados por grau de g2
		auxg1 = g1 # grafo g1
		auxg2 = g2 # grafo g2
		tam = tam1 # quantidade de nos
	} else {
		auxd1 = d2
		auxd2 = d1
		auxa1 = a2
		auxa2 = a1
		auxf1 = f2
		auxf2 = f1
		auxg1 = g2
		auxg2 = g1
		tam = tam2
	}
	auxt1 = ends(auxg1,E(auxg1))
	auxt2 = ends(auxg2,E(auxg2))
	#print(auxt1);print(auxt2)
	auxe = verify_edge_superimpostion(auxc,auxt1,auxt2)
	names(auxe)=auxc
	#auxe = 0
	#par(mfrow=c(1,2))
	#plot(g1);plot(g2)
	#print(auxc);print(auxe);readline()

	#auxmaxd = dist/max(dist)
	#auxmaxd = dist/mean(dist)
	#auxmaxd = dist/median(dist)

	#auxmaxd = 1/(dist*dist) #usar infy=0.24 para dist minimo de 2.041
	#auxmaxd = 1/(dist) #usar infy=0.49 para dist minimo 2.041
	#auxmaxd = 1/sqrt(dist)
	#auxmaxd[auxmaxd>infy]=infy
	#infy = 1/dprox
	#auxmaxd[auxmaxd>infy]=1 #tetodist=1
	#auxmaxd[auxmaxd>infy]=1 #tetodist=1
	#print(auxmaxd)#;readline()
	auxmaxd = klogistic(dist,k=2,xo=3) # decaimento logistico...
	#auxmaxd = klogistic(dist,k=2,xo=2) # decaimento logistico...
	auxes = sum(auxe);print(auxd1)
	if (auxes) {
		#auxre = auxe/auxes
		auxre = auxe/auxd1
	}else{
		auxre = auxe
	}
	names(auxre)=auxc
	for (i in 1:tam){
		#print(auxd1);print(auxd2);print(auxc)#;readline()
		#print(auxa1);print(auxa2);readline()
		auxi = as.numeric(unlist(strsplit(auxc[i],sep)))
		#print(auxi);readline()
		#auxdifa = abs(auxa1[auxi[1],auxi[1]]-auxa2[auxi[2],auxi[2]])
		#auxdifa = 1-min(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]])/max(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]]) #diff nos
		auxdifa = min(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]])/max(auxa1[auxi[1],auxi[1]],auxa2[auxi[2],auxi[2]]) #diff nos
		#print(auxi);print(auxa1);print(auxa2);print(auxdifa);readline()
		#auxdifd = abs(auxd1[auxi[1]]-auxd2[auxi[2]])
		#auxdifd = 1-min(auxd1[auxi[1]],auxd2[auxi[2]])/max(auxd1[auxi[1]],auxd2[auxi[2]]) #diff graus
		auxdifd = min(auxd1[auxi[1]],auxd2[auxi[2]])/max(auxd1[auxi[1]],auxd2[auxi[2]]) #diff graus
		#auxdiff = 1-min(auxf1[auxi[1]],auxf2[auxi[2]])/max(auxf1[auxi[1]],auxf2[auxi[2]]) #diff pesos normalizados graus
		#auxdiff = min(auxf1[auxi[1]],auxf2[auxi[2]])/max(auxf1[auxi[1]],auxf2[auxi[2]]) #diff pesos normalizados graus
		#auxrd = c(auxrd,auxmaxd[i]*auxdifd)
		#auxdife = auxe[i]/auxes
		auxa = c(auxa,auxdifa)
		auxd = c(auxd,auxdifd)
		#auxre = c(auxre,auxmaxd[i]*auxdife[i])
		#auxrf = c(auxrf,auxmaxd[i]*auxdiff)
		#print(auxra);readline()
		#auxs = c(auxs,w[1]*dist[i]+w[2]*auxdifd+w[3]*auxdifa+w[4]*auxdiff)
		#auxs = c(auxs,w[1]*auxmaxd[i]+w[2]*auxdifd*auxmaxd[i]+w[3]*auxdifa*auxmaxd[i]+w[4]*auxdiff*auxmaxd[i])
		#auxs = c(auxs,w[1]*auxmaxd[i]+w[2]*auxdife[i]*auxmaxd[i]+w[3]*auxdifa*auxmaxd[i])#+w[4]*auxdiff*auxmaxd[i])
	}
	names(auxa)=auxc
	names(auxd)=auxc
	auxre = auxre*auxmaxd
	auxra = auxa*auxmaxd
	auxrd = auxd*auxmaxd
	auxs = w[1]*auxmaxd+w[2]*auxre+w[3]*auxrd+w[4]*auxra
	
	#print(auxs);readline()
	#print(auxc);readline()
	auxsc$dist = dist # distancia euclideana nos
	auxsc$distn = auxmaxd # distancia elevada a 2	 
	auxsc$edsi = auxe # edge superimposition
	auxsc$edsv = auxre # edge superimposition
	auxsc$degi = auxd # diff graus
	auxsc$degv = auxrd # diff graus
	auxsc$nodi = auxa # diff nos (areas)
	auxsc$nodv = auxra # diff nos (areas)

	#auxsc$weiv = auxrf # diff pesos normalizados graus
	auxsc$scov = auxs # score final
	
	if (consol==1)	auxsc$scor = median(auxs)#+auxe
	if (consol==2)	auxsc$scor = mean(auxs)#+auxe
	if (consol==3)	auxsc$scor = sum(auxs)#+auxe ### dependente No de nós
	#auxsc$degs = sum(auxr)
	#print(auxsc);readline()
	return(auxsc) 
	

}

transform_by_rot_matrix_old = function(xyz,rot){
	
	#yy=cbind(yy,rep(1,dim(yy)[1]))
	#rr=matrix(x[[7]]$rot,ncol=4,byrow=T)

	#print(xyz)
	auxyz = cbind(xyz,rep(1,dim(xyz)[1]))
	auxyz = t(auxyz)

	#print(auxyz);readline()
	
	#translation pre-rotation
	auxr1 = diag(4)
	auxr1[,4]=rot[4,]
	auxr =  auxr1 %*% auxyz

	#rotation
	auxr2 = rot
	auxr2[4,] = c(0,0,0,1)
	auxr2[,4] = c(0,0,0,1)
	auxr =  auxr2 %*% auxr

	#translation pos-rotation
	auxr3 = diag(4)
	auxr3[,4]=rot[,4]
	auxr =  auxr3 %*% auxr

	auxr = t(auxr)[,1:3]
	return(auxr)

}

#make_general_similarity_assessment = function(db,pdbid,clusid,minclusid=3,workdir="Pymol"){
make_general_similarity_assessment_old = function(db,dball,ids,w=c(1,1,1,1,1),r=3,verbose=T){

	auxidname = rownames(dball$info)
	auxpdbname = as.character(dball$info[,6])
	tamrow = length(auxpdbname)
	#print(tamrow);readline()
	#tamrow = c(1:51,76:77) # apenas peptidases, ate 1OYV, mais 5m2j
	#print(tamrow);readline()
	tamcol = length(ids)
	#print(tamcol);readline()
	auxm = matrix(rep(0,tamrow*tamcol),ncol=tamcol,nrow=tamrow) #tudo
	#auxm = matrix(rep(0,length(tamrow)*tamcol),ncol=tamcol,nrow=length(tamrow)) # apenas peptidases, ate 1OYV, mais 5m2j
	#print(auxm);readline()
	rownames(auxm) = auxidname
	#rownames(auxm) = auxidname[c(1:51,76:77)]# # apenas peptidases, ate 1OYV, mais 5m2j
	colnames(auxm) = auxidname[ids]
	
	auxr = list()
	auxtab = data.frame()

	k = 1
	#print(tamcol);print(tamrow);readline()
	#print(auxm);readline()

	#args = c("Pymol/1R0R_E_I_BSR_6_5.pdb","Pymol/1TEC_E_I_BSR_7_4.pdb",i)
	#command = "python"
	#script = "Pymol/script-python-v1.py"
	#allargs = c(script,args)
	#output = system2(command, args=allargs, stdout=TRUE)

	if (verbose) print(paste("Gerating supergraph alignment scores..."))

	#for (i in 1:(tam-1)){
	for (j in 1:tamcol){
	#for (i in 5){
		if (auxpdbname[ids[j]]!="none"){
			pdbid = as.numeric(dball$info[ids[j],5])
			clusid = as.numeric(dball$info[ids[j],3])
			xyz1 = db[[pdbid]]$exp$superclus[[clusid]]$geomc
			g1 = db[[pdbid]]$exp$superclus[[clusid]]$ga
			a1 = db[[pdbid]]$exp$superclus[[clusid]]$a
			sa1 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
			#print(pdbid);print(clusid);print(xyz1);readline()
			for (i in 1:tamrow){
			#for (i in tamrow){
			#for (i in 4:tamrow){
			#for (i in c(5,11)){ #(1PPF_E_I_BSRAA_1_5, 1ACB_E_I_BSR_2_5)
			#for (i in c(61)){ #5me5_4
			#for (i in c(53)){ #2AAI_A_B_BSRaa_10_4
			#for (i in c(11)){#1acb_5
			#for (i in c(48)){#1oyv_6
			#for (i in c(11,48)){#1acb_5,1oyv_6
				if (auxpdbname[i]!="none"){
					pdbid = as.numeric(dball$info[i,5])
					clusid = as.numeric(dball$info[i,3])
					xyz2 = db[[pdbid]]$exp$superclus[[clusid]]$geomc
					g2 = db[[pdbid]]$exp$superclus[[clusid]]$ga
					a2 = db[[pdbid]]$exp$superclus[[clusid]]$a
					sa2 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
					#sa = min(sa1,sa2)/max(sa1,sa2)					
					sa = sa2/sa1
					#print(pdbid);print(clusid);print(xyz2);readline()
					#x = align_graph_pdb(auxpdbname[i],auxpdbname[j])
					#print(dim(auxm));print(x);
					#print(i);print(j);print(dim(auxm));print(auxm[i,j]);readline()
					#auxm[i,j] = align_graph_pdb(auxpdbname[ids[j]],auxpdbname[i])
					#print(i);print(auxpdbname[i])
					auxr[[k]] = best_align_graph_pdb(auxpdbname[ids[j]],auxpdbname[i],xyz1,xyz2,g1,g2,a1,a2,w)
					#print(ids[j]);print(i);
					auxr[[k]]$names = c(auxpdbname[ids[j]],auxpdbname[i])					
					#print(auxpdbname[ids[j]]);print(auxpdbname[i]);print(auxm[i,j])
					#print(x);readline()
					#auxm[j,i] = auxm[i,j]
					#print(i);print(j);
					#print(sprintf("(%.0f,%.0f): %.3f",i,ids[j],auxm[i,j]))
					#print(auxr[[k]]);readline()
					#stat = as.numeric(round(summary(auxr[[k]]$score$scov),r))
					#stat = c(stat,sum(auxr[[k]]$score$scov))
					stat = c(auxr[[k]]$score$scor,sum(auxr[[k]]$score$pain),auxr[[k]]$score$scot,sa)#;print(stat);readline()
					#stat = c(stat,sa)
					#print(sprintf("%s - %s: %.3f - %.3f",auxpdbname[ids[j]],auxpdbname[i],auxr[[k]]$rmsd,stat[4]))
					#print(sprintf("%s - %s: %.3f - %.3f",auxpdbname[ids[j]],auxpdbname[i],auxr[[k]]$rmsd,stat[7]))
					print(sprintf("%s - %s: %.3f - %.3f",auxpdbname[ids[j]],auxpdbname[i],auxr[[k]]$rmsd,stat[3]))
					#print(ids[j]);print(i);print(k);print(auxpdbname[ids[j]]);print(auxpdbname[i]);print(auxr[[k]]$rmsd)
					#print(stat)
					#auxtab0 = data.frame(ids[j],i,k,auxpdbname[ids[j]],auxpdbname[i],auxr[[k]]$rmsd,stat[1],stat[3],stat[4],stat[6],stat[7],stat[8])
					auxtab0 = data.frame(ids[j],i,k,auxpdbname[ids[j]],auxpdbname[i],auxr[[k]]$rmsd,stat[1],stat[2],stat[3],stat[4])
					#print(auxtab0);readline()
					#print("ok")
					auxtab = rbind(auxtab,auxtab0)
					#auxr[[k]]$infoids = c(ids[j],i)
					#print(auxpdbname[i]);print(auxpdbname[j]);print(as.character(auxm[i,j]));readline()
					k=k+1
				}	
			}
		}
	}
	dball$align = auxr
	colnames(auxtab)=c("ref","inf","alg","n1","n2","rmsd","scorep","pain","scoret","sarea")
	rownames(auxtab)=NULL
	auxo = auxtab
	#auxscore = auxo[order(auxo[,9]),c(4:5,9)] #by mean
	#auxscore = auxo[order(auxo[,9],decreasing=T),c(4:5,9)] #by mean
	#auxscore = auxo[order(auxo[,11],decreasing=T),c(4:5,11:12)] #by sum
	auxscore = auxo[order(auxo[,9],decreasing=T),c(4:5,7:10)] #by sum
	#print(auxscore);readline()
	dball$tab = auxtab
	dball$prescore = auxscore
	dball$score = resume_score(score=dball$prescore)
	return(dball)

}


visualize_align_super_cluster_old = function(db,dball,pdbid,ids,clusid,score=1,xlim=NULL,ylim=NULL,cols=NULL,mfrow=NULL,sep="_",dxy=c(4,2),sub="",subid=5,k=3){


	#graphnames = as.character(dball$score[,1])[1]
	#graphnames = as.character(dball$group$fakename)[1]
	
	#print(graphnames);readline()

	#auxlim = get_xylim_from_model(db=db,dball=dball,graphnames[1],dxy=dxy,sep="_")#;print("ok")
	#xlim = auxlim$xlim
	#ylim = auxlim$ylim

	if (clusid==0){
		graphnames = as.character(dball$res[[pdbid]]$scorend$pdb2)[1]
		#print(graphnames);readline()
		auxlim = get_xylim_from_model(db=db,dball,graphnames[1],dxy=dxy,sep="_")#;print("ok")
		xlim = auxlim$xlim
		ylim = auxlim$ylim
		#print(auxlim)
		graphnames = as.character(dball$res[[pdbid]]$scorend$pdb2)[ids]
		auxn1 = as.factor(as.character(dball$res[[pdbid]]$scorend$pdb1))
		#print(auxn1)
		#auxc1 = pdbsplit(as.character(auxn1))[,6]
		auxc1 = pdbsplit(as.character(auxn1))[,"nclus"]
		#auxc1 = auxc1[,dim(auxc1)[2]]
		#print(auxc1);readline()
		#print(auxn1);readline()
		levels(auxn1) = 1:length(levels(auxn1))
		auxn1 = as.numeric(auxn1)
		#print(auxn1);readline()
		sub=round(dball$res[[pdbid]]$scorend$sclocal,k)
		sub=paste(sub,auxc1,sep=sep)#;print(sub);readline()
	}else{
		graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$score$pdbid)[1]
		auxlim = get_xylim_from_model(db=db,dball=dball,graphnames[1],dxy=dxy,sep="_")#;print("ok")
		xlim = auxlim$xlim
		ylim = auxlim$ylim
		if (score==1) {
			graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$score$pdbid)[ids]
			tam = length(graphnames)
			sub = round(dball$res[[pdbid]]$alignment[[clusid]]$score$sclocal,k)[ids]
		}
		else {
			graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$prescore$n2)[ids]
			tam = length(graphnames)
			sub1 = round(dball$res[[pdbid]]$alignment[[clusid]]$prescore$scorep,k)[ids]
			sub2 = round(dball$res[[pdbid]]$alignment[[clusid]]$prescore$scoret,k)[ids]
			sub = paste(sub1,sub2,sep=sep)
		}
		#graphnames = as.character(dball$res[[pdbid]]$alignment[[clusid]]$score$pdbid)[ids]
		#tam = length(graphnames)
		#sub = round(dball$res[[pdbid]]$alignment[[clusid]]$score$sclocal,k)
	}
	#print(auxn1);readline()
	#sub = round(dball$score[ids,subid],k)
	#print(sub)

	#aux1 = dba2[[i]]$exp$superclus[[5]]$geomc

	tam = length(graphnames)#;print(tam);readline()
	

	#xx=dba2[[1]]$exp$superclus[[5]]$geomc
	#yy=dba2[[13]]$exp$superclus[[4]]$geomc
	#rr=matrix(dball[[1]]$align[[51]]$rot,ncol=4,byrow=T)
	#zz=transform_by_rot_matrix(yy,rr)


	#pars1 = unlist(strsplit(graphnames[1],sep))
	#pdbid1 = as.numeric(pars1[5])
	#clusid1 = as.numeric(pars1[6])
	#auxyz1=db[[pdbid1]]$exp$superclus[[clusid1]]$geomc
	#xlim = range(auxyz1[,1])
	#ylim = range(auxyz1[,2])
	#if (!is.null(dxy)){
	#	xlim[1] = xlim[1]-dxy[1]
	#	xlim[2] = xlim[2]+dxy[1]
	#	ylim[1] = ylim[1]-dxy[2]
	#	ylim[2] = ylim[2]+dxy[2]
	#}
	#auxlim = get_xylim_from_model(db,dball,graphnames[1],dxy=dxy,sep=sep)
	#xlim = auxlim$xlim
	#ylim = auxlim$ylim
	if (!is.null(mfrow)) par(mfrow=mfrow)
	#visualize_super_cluster(db[[pdbid1]],i=clusid1,typeg="ga",resize=1,coordi=auxyz1,cols=NULL,xlim=xlim,ylim=ylim)
	#if (!is.null(cols)) plot3d(auxyz1,box=F,type="s",col=cols[1],add=T)
	#i = 2
	#while (i<=tam){
	for (i in 1:tam){
		#print(graphnames[i])		
		pars2 = unlist(strsplit(graphnames[i],sep))
		pdbid2 = as.numeric(pars2[5])
		clusid2 = as.numeric(pars2[6])
		#auxyz2=db[[pdbid2]]$exp$superclus[[clusid2]]$geomc
		idgrp = which(graphnames[i]==dball$group$fakename)
		if (clusid==0) {
			idali = auxn1[i]
		}else{
			idali = clusid
		}
		#if (clusid) {
		#	id = which(graphnames[i]==dbgrp$fakename)
		#}
		#else{
		#	id = auxn1[i] 
		#}
		#print(graphnames[i]);print(id)
		#rotid = id
		#rotid = which(graphnames[i]==dball$tab[,"n2"])
		#print(rotid)
		#if (length(rotid)==1){
		#	rotid = dball$tab[rotid,"alg"]
		#}else{
		#	print(paste("WARNING: some problem in rotid graphnames"))
		#}
		#auxrr=matrix(dball$align[[rotid]]$rot,ncol=4,byrow=T)#;print(auxrr)
		#auxzz = transform_by_rot_matrix(auxyz2,auxrr);print(auxzz);print(dball$aligndata[[id]]$supxyz);readline()
		auxzz = dball$res[[pdbid]]$alignment[[idali]]$aligndata[[idgrp]]$supxyz
		#print(dball$aligndata[[id]]$names)
		#print(idgrp);print(auxzz);readline()
		if (!is.null(cols)) plot3d(auxzz,box=F,type="s",col=cols[i],add=T)#;print(pdbid2);print(clusid2)
		visualize_super_cluster(db[[pdbid2]],i=clusid2,typeg="ga",resize=1,coordi=auxzz,cols=NULL,xlim=xlim,ylim=ylim,sub=sub[i])
		#i=i+1
	}



}

### WARNING: NAO SEI SE HA GARANTIA DE ORDENACAO DO SCOREND... algo a verificar...
consolidate_pdb_similarity_assessment_old = function(dbali,sep="_",type=2){

	tamj = length(dbali)
	tami = dim(dbali[[1]]$score)[1]
	auxn1 = c()
	auxr = data.frame()
	auxa = data.frame()
	aux1 = data.frame()

	ind = make_index_tab_of_assessment(dbali,sep)
	#print(ind)#;readline()

	#tami = dim(ind)[1]
	#tamj = dim(ind)[2]
	
	for (i in 1:tami){
		auxsc= c()
		aux1 = c()
		for (j in 1:tamj){
			
			auxsc=c(auxsc,dbali[[j]]$score[ind[i,j],"sclocal"])
		}
		auxmax = which(auxsc==max(auxsc))
		clu = auxmax[1]
		pdb = ind[i,auxmax[1]]
		auxr = dbali[[clu]]$score[pdb,]
		aux1 = dbali[[clu]]$score[1,1]
		#print(aux1);print(auxr);readline()
		auxr = data.frame(aux1,auxr)
		auxa = rbind(auxa,auxr)
		#print(auxsc);print(auxmax);print(auxa);readline()
	}
	#print(auxa[,1])
	auxr1 = as.factor(as.character(auxa[,1]))
	auxr1 = table(auxr1)
	#return(aux1)
	if (type==1){
		#auxr1 = as.factor(as.character(auxa[,1]))
		#auxr1 = table(auxr1)
		#return(aux1)
		clu1 = which(max(auxr1)==auxr1)[1]
		auxr = dbali[[clu1]]$score[1,]
		aux1 = dbali[[clu1]]$score[1,1]
		auxr = data.frame(aux1,auxr)
		auxa[1,] = auxr
	}
	if (type==2){
		auxr2 = sort(auxr1,decreasing=T)
		auxn1 = names(auxr1)
		auxn2 = names(auxr2)
		#print(auxn1);readline()
		clu1 = which(auxn2[1]==auxn1)
		auxr = dbali[[clu1]]$score[1,]
		aux1 = dbali[[clu1]]$score[1,1]
		auxr1 = data.frame(aux1,auxr)
		#print(auxr1)#;readline()
		if (length(auxr1)>1){
			clu2 = which(auxn2[2]==auxn1)#;print(auxn1);print(auxn2)
			auxr = dbali[[clu2]]$score[1,]
			aux2 = dbali[[clu2]]$score[1,1]
			auxr2 = data.frame(aux2,auxr)
			#print(auxr2);readline()
			auxa[1,] = auxr2
			auxa = rbind(auxr1,auxa)
		}else{
			auxa[1,] = auxr1
		}
		#print(auxa);readline()
		#clu1 = which(max(auxr1)==auxr1)[1]
		#auxr = dbali[[clu1]]$score[1,]
		#aux1 = dbali[[clu1]]$score[1,1]
		#auxr = data.frame(aux1,auxr)
		#auxa[1,] = auxr

	}
	#print(aux1);print(id1);readline()
	colnames(auxa)[1:2] = c("pdb1","pdb2")
	#print(auxa);readline()
	auxo = order(auxa[,"sclocal"],decreasing=T)
	return(auxa[auxo,])
}

find_nearest_nodes_old2 = function(xyz1,xyz2,id1=NULL,id2=NULL){

	tam1 = dim(xyz1)[1]
	tam2 = dim(xyz2)[1]
	auxr = c()

	if (is.null(id1)) id1 = 1:tam1
	if (is.null(id2)) id2 = 1:tam2


	if (tam1<=tam2){
		auxyz1 = xyz1
		auxyz2 = xyz2
		tam = tam1
	} else {
		auxyz1 = xyz2
		auxyz2 = xyz1
		tam = tam2
		temp = id1
		id1 = id2
		id2 = temp
	}
	for (i in 1:tam){
		auxd = dist.xyz(auxyz1[i,],auxyz2)
		#print(auxd)
		auxmin = min(auxd)
		auxid = which(auxd==auxmin,arr.ind=T)
		auxs = paste(id1[i],id2[auxid[2]],sep="-")
		names(auxmin)=auxs
		#print(auxmin);readline()
		auxr = c(auxr,auxmin)
	}
	#print(auxr);readline()
	return(auxr)
}

verify_edge_superimpostion_old = function(dist,e1,e2,sep="-"){

	#print(dist);print(e1);print(e2)
	auxc=as.numeric(unlist(strsplit(dist,sep)))
	#print(auxc)
	x = matrix(auxc,ncol=2,byrow=T)
	#print(x)#;readline()
	#tam1 = dim(e1)[1]
	tam1 = length(dist)
	rr = rep(0,tam1)

	for (i in 1:tam1){
		xi = x[i,1]
		me1 = e1[e1[,1]==xi,]
		#print(me1)
		if (!is.matrix(me1)){
			me1 = t(as.matrix(me1))
		}
		#print(me1)#;readline()
		tame1 = dim(me1)[1]
		if (tame1>0){
			for (j in 1:tame1){
				c1 = x[me1[j,1],2]
				c2 = x[me1[j,2],2]
				c12 = c(c1,c2)
				c21 = c(c2,c1)
				if (!( is.na(row.match(c12,e2))&is.na(row.match(c21,e2)) )){ #de Morgan...
					#n = n + 1
					rr[me1[j,1]]=rr[me1[j,1]]+1
					rr[me1[j,2]]=rr[me1[j,2]]+1
					#print(me1[j,1]);print(me1[j,2]);print(rr);readline()
				}
				#print(c12);print(c21);print(rr);readline()
			}
		}

	}
	#auxr=list()
	#auxr$esup = rr
	#auxr$esun = n 
	return(rr)
}

#visualize_align_supercluster_prescore(db=dba2,dball=dball[[i]],pdbname="1ACB",cols=cols,dxy=dxy)
visualize_align_supercluster_prescore_old = function(db,dball,pdbnames,cols,dxy){


	id1 = as.numeric(as.vector(dball$group[which(dball$group$pdb==pdbnames[1]),"id"]))
	id2 = as.numeric(as.vector(dball$group[which(dball$group$pdb==pdbnames[2]),"id"]))
	tam1 = length(id1)
	tam2 = length(id2)
	parid=c(id1[1],1,0)
	#print(id1);print(id2);print(parid);readline()

	ids=grep(pdbnames[2],dball$res[[parid[1]]]$alignment[[parid[2]]]$prescore$n2)
	visualize_align_super_cluster(db=dba2,dball=dball,pdbid=parid[1],clusid=parid[2],ids=c(1,ids),score=parid[3],cols=cols,mfrow=c(tam1,tam2+1),dxy=dxy)
	if (tam1>1){
		for (j in 2:tam1){
			parid=c(id1[1],j,0)
			#print(parid);readline()
			ids=grep(pdbnames[2],dball$res[[parid[1]]]$alignment[[parid[2]]]$prescore$n2)
			visualize_align_super_cluster(db=dba2,dball,pdbid=parid[1],clusid=parid[2],ids=c(1,ids),score=parid[3],cols=cols,mfrow=NULL,dxy=dxy)
		}
	}

}

find_nearest_nodes_old = function(xyz1,xyz2,tam1=NULL,tam2=NULL){


	#tam1 = dim(xyz1)[1]
	#tam2 = dim(xyz2)[1]

	if (is.null(tam1)) tam1 = dim(xyz1)[1]
	if (is.null(tam2)) tam2 = dim(xyz2)[1]
	auxr = c()

	if (tam1<=tam2){
		auxyz1 = xyz1
		auxyz2 = xyz2
		#tam = tam1
		#tam = dim(auxyz1)[1]
	} else {
		auxyz1 = xyz2
		auxyz2 = xyz1
		#tam = tam2
		#tam = dim(auxyz1)[1]
	}
	tam = dim(auxyz1)[1]
	for (i in 1:tam){
		auxd = dist.xyz(auxyz1[i,],auxyz2)
		#print(auxd)
		auxmin = min(auxd)
		auxid = which(auxd==auxmin,arr.ind=T)
		auxs = paste(i,auxid[2],sep="-")
		names(auxmin)=auxs
		#print(auxmin);readline()
		auxr = c(auxr,auxmin)
	}
	#print(auxr);readline()
	return(auxr)
}

find_asymetric_nodes_old = function(mx,my,mz){


	mz.upp = upper.tri(mz)
	mz.low = lower.tri(mz)
	mza = mz
	mzb = mz
	mza[mz.upp]=0#;print(za)
	mzb[mz.low]=0#;print(zb)

	vm = which(mza!=t(mzb),arr.ind=T)#;print(v)#;readline()
	
	mx.d = diag(mx)
	my.d = diag(my)

	vd = which(mx.d!=my.d)#;print(v);readline() 

	vd = t(mapply(rep,vd,2))

	res = rbind(vm,vd)
	
	res = lapply(apply(res,1,list),unlist)

	return(res)

}

make_matrix_symmetrics_old = function(mx,my,an,zero=0.001){

	#tam = dim(an)[1]
	msg = "WARNING: something strange in make matrix symmetrics:"

	while(length(an)){
		i=an[[1]][1]
		j=an[[1]][2]
		if (mx[i,j]>zero){
			print("ma=my")
			res = change_cel(ma=my,mb=mx,an=an,i=i,j=j)#;print(res);readline()
			mx = res$mb
			my = res$ma
			an = res$an
		}else if (my[i,j]>zero){#print("o")
			print("ma=mx")
			res = change_cel(ma=mx,mb=my,an=an,i=i,j=j)
			mx = res$ma
			my = res$mb
			an = res$an
		}else{
			print(paste(msg,3))
		}
	}
	res=list()
	res$mx = mx
	res$my = my
	return(res)
	
}

change_cel_old = function(ma,mb,an,i,j,go,zero=0.001){

	msg = "WARNING: something strange in make matrix symmetrics:"

	k = an[[1]][2]
	an.id = which(sapply(an,retvi,2,k))#;print(an.id);print(i);print(j)
	if (length(an.id)==2){
		if(is_valid_change_cel(i=j,m1=ma,m2=mb)){#;print(go);readline()
			ma[j,i] = mb[i,j]
			i = an[[an.id[2]]][1]
			j = an[[an.id[2]]][2]
			ma[i,j] = 0
			an[[an.id[2]]] = NULL
			an[[an.id[1]]] = NULL
		}else{
			temp = an[[an.id[1]]]
			an[[an.id[1]]] = an[[an.id[2]]]
			an[[an.id[2]]] = temp
		}
	}else
	if (length(an.id)==1){
		#print(mb);readline()
		mb[i,j] = 0
		an[[an.id[1]]] = NULL
		go = F
	}else{
		print(paste(msg,1))
	}
	res = list()
	res$ma = ma
	res$mb = mb
	res$an = an
	res$go = go
	#print(res);readline()
	return(res)

}

make_matrix_symmetrics_old3 = function(mx,my,an,go,zero=0.001){

	#tam = dim(an)[1]
	msg = "WARNING: something strange in make matrix symmetrics:"

	while(length(an)){
		i=an[[1]][1]
		j=an[[1]][2]
		if (mx[i,j]>zero){
			print("1- ma=my")
			res = change_cel(ma=my,mb=mx,an=an,i=i,j=j,go=go)#;print(res);readline()
			mx = res$mb
			my = res$ma
			an = res$an
			go = res$go
		}else
		if (my[j,i]>zero){
			print("2- ma=mx,i=j,j=i")
			res = change_cel(ma=mx,mb=my,an=an,i=j,j=i,go=go)#;print(res);readline()
			mx = res$ma
			my = res$mb
			an = res$an
			go = res$go
		}else
		if(my[i,j]>zero){#print("o")
			print("3- ma=mx")
			res = change_cel(ma=mx,mb=my,an=an,i=i,j=j,go=go)
			mx = res$ma
			my = res$mb
			an = res$an
			go = res$go
		}else
		if (mx[j,i]>zero){
			print("4- ma=my,i=j,j=i")
			res = change_cel(ma=my,mb=mx,an=an,i=j,j=i,go=go)#;print(res);readline()
			mx = res$mb
			my = res$ma
			an = res$an
			go = res$go
		}else{
			print(paste(msg,3))
		}
	}
	fin=list()
	fin$mx = mx
	fin$my = my
	fin$go = go
	return(fin)
	
}


change_cel_old2 = function(ma,mb,an,i,j,zero=0.001){

	msg = "WARNING: something strange in make matrix symmetrics:"
	res=list()

	if (ma[j,i]>zero){
		print(paste(msg,1))
	}else{
		ma[j,i] = mb[i,j]
		an[[1]]=NULL
		am = matrix(unlist(an),ncol=2,byrow=T)
		v = which(am[,1]==j)
		if (length(v)!=1){
			print(paste(msg,2))
		}else{
			i = an[[v]][1]
			j = an[[v]][2]
			ma[i,j] = 0
			an[[v]]=NULL
		}
	}
	res$ma = ma
	res$mb = mb
	res$an = an
	return(res)

}

make_superclus_to_fake_pdb_old = function(db,dball,type=1,minclusid=3,minpdb=16,verbose=T,seps="_",workdir="Pymol/",sufix=".pdb"){

	if (is.null(dball$info$fakename)){
		auxinfo = dball$info
		#print(auxinfo)
		tami = dim(auxinfo)[1]
		#print(tami)
		auxfilepdb = c()

	

		if (verbose) print(paste("Generating fake pdbs in workdir ",workdir))

		for (i in 1:tami){
			#tamj = length(db[[i]]$exp$superclus)
			j = as.numeric(auxinfo[i,3])#;j=1
			k = as.numeric(auxinfo[i,5])#;print(j);print(i);print(k);print(minclusid);readline()
			#print(is.data.frame(as.numeric(j)))
			#print(j<minclusid);readline()
			if (j<minclusid){
				#print("0k1")
				if (verbose) print(paste("WARNING: PDB",db[[k]]$pdbname,"has superclus length",j,"less than minclusid",minclusid))
				auxfilepdb = c(auxfilepdb,"none")#;print("0k")
			}else{
			#for (j in tamj){
				#print("0k2")
				auxname = paste(db[[k]]$matrixname,seps,k,seps,j,sep="")
				auxfilefake = paste(workdir,auxname,sufix,sep="")
				#print(auxfilefake);readline()
				if (file.exists(auxfilefake)){
					if (verbose) print(paste("WARNING: fake PDB",auxfilefake,"detected. It will not be overwritten"))
				}else{ 			
					auxyz = db[[k]]$exp$superclus[[j]]$geomc
					#print(auxyz)
					if (type==1){
						auxyz = increase_matrix_to_min_random(auxyz,minxyz=minpdb)
					}
					if (type==2){
						auxyz = increase_matrix_to_cubic(auxyz,minxyz=minpdb)
					}
					#print(auxyz);readline()		
					#auxname = paste(db[[k]]$matrixname,seps,k,seps,j,sep="")
					auxpdb = matrix2pdb(auxyz,remark=toupper(auxname))
					#auxpdbname = paste(auxname,sufix,sep="")
					#auxpdbname = auxname
					write.table(file=auxfilefake,auxpdb,quote=F,row.names=F,col.names=F)
				}
				auxfilepdb = c(auxfilepdb,auxname)
				#print(auxpdb);readline()
			}
		}
		auxcol = colnames(auxinfo)
		auxinfo = data.frame(auxinfo,auxfilepdb)
		colnames(auxinfo)=c(auxcol,"fakename")
		dball$info = auxinfo
		}
	else{
		print("column fakename already exists; so doing nothing with fakepdbs")
	}
	return(dball)

}

new_consolidate_pdb_similarity_assessment_old = function(db,dbres,dbsel,pdblist,debug=F){

	
	tami = length(pdblist)
	auxr = data.frame()
	auxa = data.frame()
	auxl = list()

	i=1
	#for (i in 1){
	#for (i in 2){
#		pdbnamei = pdblist[i]
#		#print(paste("NEW Consolidating similarity between",toupper(dbres$pdbname),"-",toupper(pdbnamei)))
#		aux1 = join_prescore_by_pdb(dbres,pdbnamei)#;print(auxr);readline()
#		aux1.len = dim(aux1)[1]
#		id = seq(1,aux1.len,as.integer(sqrt(aux1.len)))#;print(id)#;readline()
#		aux1 = aux1[id,]
		#auxo = order(auxr$scoret,decreasing=T)[1:id]
		#auxa = rbind(auxa,auxr[auxo,]);print(auxa);readline()
	#}
	#print(aux1)#;readline()

	aux1 = dbsel
	#tamj = dim(aux1)[1]
	#print(aux1)

	for (i in 1:tami){
	#for (i in 2){
		pdbnamei = pdblist[i]
		print(paste("NEW Consolidating similarity between",toupper(dbres$pdbname),"-",toupper(pdbnamei)))
		auxr = join_prescore_by_pdb(dbres,pdbnamei)#;print(auxr);readline()
		tamj = dim(auxr)[1]
		scorea = c()
		for (j in 1:tamj){
			#vid1 = auxr$n1==aux1[j,"n1"];print(vid1);readline()
			vid1 = as.character(aux1$n1)==as.character(auxr[j,"n1"])#;print(vid1)#;readline()
			vid2 = as.character(aux1$n2)==as.character(auxr[j,"n2"])#;print(vid2)#;
			scorea1 = aux1[vid1,"scoret"]
			scorea2 = aux1[vid2,"scoret"]
			#scorea = c(scorea,auxr[j,"scoret"]/min(scorea1,scorea2))#;print(scorea);readline()
			scor = (2*auxr[j,"scoret"])/(scorea1+scorea2)
			scorea = c(scorea,scor)#;print(scorea);readline()
			alig.id = auxr[j,"id1"]
			pre.id = auxr[j,"auxid"]
			dbres$alignment[[alig.id]]$prescore[pre.id,"scorea"]=scor
			#auxr = data.frame(auxr,scorea=scorea);
		}
		#auxr = data.frame(auxr,scorea=scorea)#;print(auxr);readline()
		auxr$id = NULL
		auxr$auxid = NULL
		auxr$scorea = scorea
		#colnames(auxr) = c("id0","id1","id2","n1","n2","local","global","scorea","kernel","scoret")
		#auxo = order(auxr$scoret,decreasing=T)[1]
		auxo = order(auxr$scorea,decreasing=T)[1]
		auxa = rbind(auxa,auxr[auxo,])#;print(auxa);readline()
	}
	#colnames(auxa) = c("id0","id1","id2","n1","n2","local","global","scorea","kernel","scoret")
	#auxo = order(auxa$scoret,decreasing=T)
	auxo = order(auxa$scorea,decreasing=T)
	auxa = auxa[auxo,]#;print(auxa);readline()
	dbres$scorend = auxa
	#print(dbres$scorend);readline()
	return(dbres)
}

# VERSAO NOVA QUE FICOU VELHA: usando funcao propria!!! ;-) de alinhamento para construir similaridade in memory
new_make_similarity_assessment_old = function(db,dball,id0,pdblist,dbhot=NULL,rotref=diag(4),self=F,w=c(1,1,1,1,1),r=3,verbose=T,outdir="Datalign/", sufix="align.csv",sep="_", comp=F){


	auxres = list()
	pdbname = toupper(pdblist[id0])
	auxres$pdbname = pdbname
	auxgrp = which(dball$group[,"pdb"]==pdbname)
	auxres$groupid = auxgrp
	auxres$alignment = list()
	id = auxgrp

	auxall=list()
	#id = 5
	id.len = length(id)

	if (self){
#		tami = id[1]
#		tamj = id
		if (id.len>1){
			tami = id[1]:id[(id.len-1)]
		}
		else{
			print(paste("WARNING: id length for self similarty is insufficient in",pdbname))
		}
	}else{
		tami = id
		tamj = 1:dim(dball$group)[1]#;tam=3
	}
	k = 1
	#id = 1 # test 1ppf 5
	#id=1
	#tami = 5
	#tamj = 3
	for (i in tami){
		#print(paste("i:",i))
		pdbid = as.numeric(as.vector(dball$group[i,"id"]))
		clusid = as.numeric(as.vector(dball$group[i,"nclus"]))#;print(pdbid);print(clusid)
		xyz1 = db[[pdbid]]$exp$superclus[[clusid]]$geomc#;print(xyz1)
		g1 = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print("o");print(g1);print("o")
		a1 = db[[pdbid]]$exp$superclus[[clusid]]$a#;print(a1)
		sa1 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
		auxname1 = as.character(dball$group[i,"fakename"])
		pdbname1 = as.character(dball$group[i,"pdb"])
		if (!is.null(dbhot)) hot1 = dbhot[[pdbid]]
		else hot1 = NULL

		auxr = list()
		auxtab = data.frame()
		id0 = i
		id1 = k
		if (self) tamj = c(i,i+1)
		for (j in tamj){
		#for (j in 22){ #test 1oyv 7
		#for (j in 5){ #test 1acb 5
		#for (j in 2){ #test 1ppf 3
			#print(paste("j:",j))
			id2 = j
			pdbid = as.numeric(as.vector(dball$group[j,"id"]))
			clusid = as.numeric(as.vector(dball$group[j,"nclus"]))
			xyz2 = db[[pdbid]]$exp$superclus[[clusid]]$geomc
			g2 = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print("o");print(g2);print("o")
			a2 = db[[pdbid]]$exp$superclus[[clusid]]$a
			sa2 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
			sa = sa2/sa1#;print("o")
			auxname2 = as.character(dball$group[j,"fakename"])#;print(auxname1);print(auxname2)
			pdbname2 = as.character(dball$group[j,"pdb"])
			if (!is.null(dbhot)) hot2 = dbhot[[pdbid]]
			else hot2 = NULL
			#if (pre){
			#auxpdbname = dball$group$pdb[id]
			#auxv = (auxlist[,"id1"]==i)&(auxlist[,"id2"]==j)#;print("oi")
			#print(sum(auxv));readline()
			#return(auxlist)
			#auxr[[j]] = pre_best_align_graph_pdb(auxlist[auxv,],xyz1,xyz2,g1,g2,a1,a2,rotref,w)
			#if (verbose) print(sprintf("%s - %s ",auxname1,auxname2))
			#print("o")
			auxr[[j]] = new_best_align_graph_pdb(xyz1,xyz2,g1,g2,a1,a2,pdbname1,pdbname2,hot1,hot2,rotref,w)#;print("o")	
			##auxr[[k]] = new_best_align_graph_pdb(xyz1,xyz2,g1,g2,a1,a2,pdbname1,pdbname2,rotref,w)#;print("o")			
			#print(auxr);readline()
			#}else{
			#	auxr[[i]] = best_align_graph_pdb(auxname1,auxname2,xyz1,xyz2,g1,g2,a1,a2,w)
			#}
			auxr[[j]]$names = c(auxname1,auxname2)
			##auxr[[k]]$names = c(auxname1,auxname2)
			stat=c()
			#print(auxr[[j]]$score$global$scor);readline()
			##stat = auxr[[k]]$score$local$scor[1]
			##stat = c(stat,auxr[[k]]$score$global$scor[1])
			##stat = c(stat,auxr[[k]]$score$scot)
			##stat = c(stat,auxr[[k]]$score$ker)
			stat = auxr[[j]]$score$local$scor[1]
			stat = c(stat,auxr[[j]]$score$global$scor[1])
			stat = c(stat,auxr[[j]]$score$scot)
			stat = c(stat,auxr[[j]]$score$ker)
			if (verbose) print(sprintf("%s - %s - %.3f - %.3f - %.3f - %.3f",auxname1,auxname2,stat[1],stat[2],stat[3],stat[4]))#;readline()
			auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,stat[1],stat[2],stat[3],stat[4])
			if (stat[3]==0) print(paste("WARNING: it was not possible to align",auxname1,"-",auxname2))
			auxtab = rbind(auxtab,auxtab0)
			#stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$scot,sa)#;print(stat);readline()
			#if (comp) {
			#	stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$gain,auxr[[j]]$score$scot,sa)
			#	if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[j]]$rmsd,stat[4]))
			#	auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,auxr[[j]]$rmsd,stat[1],stat[2],stat[3],stat[4],stat[5])
			#	auxtab = rbind(auxtab,auxtab0)
			#}else{
			#	stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$scot,sa)
			#	if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[j]]$rmsd,stat[3]))
			#	auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,auxr[[j]]$rmsd,stat[1],stat[2],stat[3],stat[4])
			#	auxtab = rbind(auxtab,auxtab0)
			#}
			#k = k+1
		}
	#}
		#dball$aligndata = auxr
		colnames(auxtab)=c("id0","id1","id2","n1","n2","local","global","scoret","kernel")
		rownames(auxtab)=NULL
		if (!self){
			auxo = order(auxtab[,8],decreasing=T)
			auxscore = auxtab[auxo,]#;prbyint(auxscore);readline()
		}else{
			auxscore = auxtab
		}
		#auxall[[1]] = list()
		#auxall[[1]]$aligndata = auxr
		#auxall[[1]]$prescore = auxscore
		#auxall[[1]]$prescore$scorea = 0
		auxall[[k]] = list()
		auxall[[k]]$aligndata = auxr
		auxall[[k]]$prescore = auxscore
		auxall[[k]]$prescore$scorea = 0
		
		#if (comp){
		#	colnames(auxtab)=c("id0","id1","id2","n1","n2","rmsd","scorep","pain","gain","scoret","sarea")
		#}else{
		#	colnames(auxtab)=c("id0","id1","id2","n1","n2","rmsd","scorep","pain","scoret","sarea")
		#}
		#rownames(auxtab)=NULL
		#auxo = auxtab
		#auxscore = auxo[order(auxo[,9]),c(4:5,9)] #by mean
		#auxscore = auxo[order(auxo[,9],decreasing=T),c(4:5,9)] #by mean
		#auxscore = auxo[order(auxo[,11],decreasing=T),c(4:5,11:12)] #by sum
		#auxscore = auxo[order(auxo[,"scoret"],decreasing=T),] #by sum
		#print(auxscore);readline()
		#dball$tab = auxtab
		#dball$tab = auxtab
		#dball$prescore = auxscore
		#dball$score = resume_score(score=dball$prescore)
		#auxall[[k]]$prescore = auxscore
		#auxall[[k]]$score = resume_score(score=auxall[[k]]$prescore)
		k=k+1
	}
	auxres$alignment = auxall
	#return(auxall)
	return(auxres)

}

make_group_table_old = function(dball,filename=NULL,outdir="Datalign/",sufix=".csv"){

	auxv = dball$info[,6]!="none"
	auxinf = dball$info[auxv,]
	auxr = pdbsplit(auxinf[,6])
	auxr=data.frame(auxinf[,c(1:2)],auxr,auxinf[,c(6)])
	colnames(auxr) = c("gacc","lacc","pdb","chain1","chain2","method","id","nclus","fakename")
	#print(auxr)
	dball$group = auxr
	if (!is.null(filename)){
		filename = paste(outdir,toupper(filename),sep,sufix,sep="")
		write.csv(auxr,filename)
	}
	return(dball)

}

resume_for_web_output_old2 = function (db,dball,id,k=0){

#> connetions_group_bsr[3]
#[[1]]
#  source source_volume target target_volume edge_sum      s_x      s_y      s_z      t_x      t_y      t_z
#1      1      5553.502      2     19821.187 1646.727 30.46871 69.08758 11.42233 27.80587 54.77508 12.79056
#2      1      5553.502      3      2085.385  147.141 30.46871 69.08758 11.42233 38.81296 61.25319 19.91191
#3      2     19821.187      3      2085.385  164.266 27.80587 54.77508 12.79056 38.81296 61.25319 19.91191

#auxp2 = pdbsplit(auxn1[j])#;print(auxp1);print(auxp2)
#pdbid2 = as.numeric(as.vector(auxp2[5]))
#clusid2 = as.numeric(as.vector(auxp2[6]))
#auxt1 = ends(auxg1,E(auxg1))
#dball[[35]]$res[[1]]$alignment[[3]]$aligndata[[3]]$supxyz
#matrix(dball[[1]]$align[[51]]$rot,ncol=4,byrow=T)
#find_nearest_nodes
#c(5,x[!(x %in% 5)])

	auxl = list()
	if (!is.null(dball$res[[id]]$scorend)){
		tami = dim(dball$res[[id]]$scorend)[1];print(tami)#;readline()
		#auxr = data.frame()
		#auxa = data.frame()
	
		print(paste("Resuming for web output..."))
		pdbname0 = dball$res[[id]]$pdbname
		i0 = grep(pdbname0,dball$res[[id]]$scorend$n2);print(i0);  
		al01 = dball$res[[id]]$scorend$id1[i0]
		al02 = dball$res[[id]]$scorend$id2[i0]
		xyz0 = dball$res[[id]]$alignment[[al01]]$aligndata[[al02]]$supxyz
		#xyz0 = transform_by_rot_matrix(xyz0,rotref)
		#i0 = 3
		ilist = 1:tami
		ilist = c(i0,ilist[!(ilist %in% i0)])
		#print(i0);print(al01);print(al02);print(xyz0);print(ilist);readline()
		#for (i in 1:tami){#print(dball$res[[id]]$scorend$n2[1])
		for (i in ilist){
			auxp = pdbsplit(dball$res[[id]]$scorend$n2[i])		
			#pdbid = as.numeric(as.vector(auxp[5]))
			#clusid = as.numeric(as.vector(auxp[6]))
			pdbid = as.numeric(as.vector(auxp["id"]))
			clusid = as.numeric(as.vector(auxp["nclus"]))
			a = round(db[[pdbid]]$exp$superclus[[clusid]]$a,k)
			g = db[[pdbid]]$exp$superclus[[clusid]]$ga
			e = ends(g,E(g))
			e = cbind(e,round(E(g)$weight,k))
			al1 = dball$res[[id]]$scorend$id1[i]
			al2 = dball$res[[id]]$scorend$id2[i]#;print(al1);print(al2)
			xyz = dball$res[[id]]$alignment[[al1]]$aligndata[[al2]]$supxyz
			ngroup = dim(xyz)[1]
			if (!is.null(db[[pdbid]]$exp$preclus)){
				ngroup = ngroup + length(db[[pdbid]]$exp$preclus$ids)
			}#;print(i)
			#xyz = transform_by_rot_matrix(xyz,rotref)
			rot = dball$res[[id]]$alignment[[al1]]$aligndata[[al2]]$rot
			#plot(g,layout=xyz);print(a);print(e);print(xyz);readline()
			tamj = dim(e)[1]
			auxa = data.frame()
			for (j in 1:tamj){
				auxr = data.frame(source=e[j,1],source_volume=a[e[j,1],e[j,1]],target=e[j,2],target_volume=a[e[j,2],e[j,2]],edge_sum=e[j,3])
				xyzj = c(xyz[e[j,1],],xyz[e[j,2],])
				auxr = data.frame(auxr,s_x=xyzj[1],s_y=xyzj[2],s_z=xyzj[3],t_x=xyzj[4],t_y=xyzj[5],t_z=xyzj[6])
				auxa = rbind(auxa,auxr)
				#print(auxa);readline()
			}
			#x = round(find_nearest_nodes(xyz0,xyz),k)
			#print(as.character(dball$res[[id]]$scorend[i,5]));print(x);readline()

			auxl[[i]]=list()
			auxl[[i]]$matrix = auxa
			auxl[[i]]$scorend = dball$res[[id]]$scorend[i,]
			auxl[[i]]$rot = matrix(rot,ncol=4,byrow=F) #ATENCAO!!! mudou byrow de T para F
			auxl[[i]]$ngroup = ngroup
			#auxl[[i]]$rotref = dball$rotref#;print(id);readline() 
			#auxl[[i]]$corr = names(x)

		}
		#print(auxl);readline()
	}
	return(auxl)
}

get_bsr_info_interactions_old = function(bsr,atom_names,atom_codes,k=1){#TERRIVELMENTE lenta...

	tami = dim(bsr)[1]
	tamj = dim(bsr)[2]

	#;print(atom_names);readline()
	#atom_codes = (atom_names,meso,aacode)#;print(pol);readline()

	bsr_tab = data.frame()

	for (i in 1:tami){
	
		for (j in i:tamj){
			bsr_value = round(bsr[i,j],k)
			if (bsr_value>0){
				source = atom_names[i]
				target = atom_names[j]
				if (atom_codes[i]=="a"){
					source_polarity = "APOLAR"
				}else{
					source_polarity = "POLAR"
				}
				if (atom_codes[j]=="a"){
					target_polarity = "APOLAR"
				}else{
					target_polarity = "POLAR"
				}
			
				bsr_line = data.frame(source=source,source_polarity=source_polarity,target=target,target_polarity=target_polarity,bsr_value=bsr_value)
				bsr_tab = rbind(bsr_tab,bsr_line)#;print(bsr_tab);readline()
			}
		}
	}
	return(bsr_tab)

}

make_hot_assessment_bad = function(db,dball,id,pdbname,rotref,w=c(1,1,1,1,1),r=3,verbose=T,outdir="Datalign/",sufix="align.csv",sep="_", comp=F){


	#filename = paste(outdir,toupper(pdbname),sep,sufix,sep="")
	#print(paste("Reading pre-computation file:",filename))
	#auxlist = read.table(filename,header=T,sep=",")
	#auxcoln = colnames(auxlist)
	#auxlist = as.matrix(auxlist)
	#colnames(auxlist) = auxcoln
	#print(auxlist[1,]);readline()
	#print(colnames(auxlist));readline()
	#print(auxlist[1,])
	#print(dim(auxlist));readline()

	#print(pdbname)#;readline()
	auxall=list()
#	tam = dim(dball$group)[1]#;tam=3
	tam = which(dball$group$pdb==pdbname)#;print(tam);readline()	
	k = 1
	#id = 1 # test 1ppf 5
	for (i in id){

		pdbid = as.numeric(as.vector(dball$group[i,"id"]))
		clusid = as.numeric(as.vector(dball$group[i,"nclus"]))#;print(pdbid);print(clusid)
		xyz1 = db[[pdbid]]$exp$superclus[[clusid]]$geomc#;print(xyz1)
		g1 = db[[pdbid]]$exp$superclus[[clusid]]$ga#;print(g1)
		a1 = db[[pdbid]]$exp$superclus[[clusid]]$a#;print(a1)
		sa1 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
		auxname1 = as.character(dball$group[i,"fakename"])

		auxr = list()
		auxtab = data.frame()
		id0 = i
		id1 = k
		for (j in tam){
		#for (j in 22){ #test 1oyv 7
		#for (j in 5){ #test 1acb 5
		#for (j in 3){ #test 1ppf 3
			id2 = j
			pdbid = as.numeric(as.vector(dball$group[j,"id"]))
			clusid = as.numeric(as.vector(dball$group[j,"nclus"]))
			xyz2 = db[[pdbid]]$exp$superclus[[clusid]]$geomc
			g2 = db[[pdbid]]$exp$superclus[[clusid]]$ga
			a2 = db[[pdbid]]$exp$superclus[[clusid]]$a
			sa2 = sum(db[[pdbid]]$exp$superclus[[1]]$a)
			sa = sa2/sa1
			auxname2 = as.character(dball$group[j,"fakename"])#;print(auxname1);print(auxname2)
			#if (pre){
			#auxpdbname = dball$group$pdb[id]
			#auxv = (auxlist[,"id1"]==i)&(auxlist[,"id2"]==j)#;print("oi")
			#print(sum(auxv));readline()
			#return(auxlist)
			#auxr[[j]] = pre_best_align_graph_pdb(auxlist[auxv,],xyz1,xyz2,g1,g2,a1,a2,rotref,w)
			#if (verbose) print(sprintf("%s - %s ",auxname1,auxname2))
			auxr[[j]] = new_best_align_graph_pdb(xyz1,xyz2,g1,g2,a1,a2,rotref,w)			
			#print(auxr);readline()
			#}else{
			#	auxr[[i]] = best_align_graph_pdb(auxname1,auxname2,xyz1,xyz2,g1,g2,a1,a2,w)
			#}
			auxr[[j]]$names = c(auxname1,auxname2)
			stat=c()
			#print(auxr[[j]]$score$global$scor);readline()
			stat = auxr[[j]]$score$local$scor[1]
			stat = c(stat,auxr[[j]]$score$global$scor[1])
			stat = c(stat,auxr[[j]]$score$scot)
			if (verbose) print(sprintf("%s - %s - %.3f - %.3f - %.3f",auxname1,auxname2,stat[1],stat[2],stat[3]))#;readline()
			auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,stat[1],stat[2],stat[3])
			auxtab = rbind(auxtab,auxtab0)
			#stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$scot,sa)#;print(stat);readline()
			#if (comp) {
			#	stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$gain,auxr[[j]]$score$scot,sa)
			#	if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[j]]$rmsd,stat[4]))
			#	auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,auxr[[j]]$rmsd,stat[1],stat[2],stat[3],stat[4],stat[5])
			#	auxtab = rbind(auxtab,auxtab0)
			#}else{
			#	stat = c(auxr[[j]]$score$scor,sum(auxr[[j]]$score$pain),auxr[[j]]$score$scot,sa)
			#	if (verbose) print(sprintf("%s - %s: %.3f - %.3f",auxname1,auxname2,auxr[[j]]$rmsd,stat[3]))
			#	auxtab0 = data.frame(id0,id1,id2,auxname1,auxname2,auxr[[j]]$rmsd,stat[1],stat[2],stat[3],stat[4])
			#	auxtab = rbind(auxtab,auxtab0)
			#}
		}
	#}
		#dball$aligndata = auxr
		colnames(auxtab)=c("id0","id1","id2","n1","n2","local","global","scoret")
		rownames(auxtab)=NULL
		auxo = order(auxtab[,8],decreasing=T)
		auxscore = auxtab[auxo,]#;prbyint(auxscore);readline()
		auxall[[k]] = list()
		auxall[[k]]$aligndata = auxr
		auxall[[k]]$prescore = auxscore
		
		#if (comp){
		#	colnames(auxtab)=c("id0","id1","id2","n1","n2","rmsd","scorep","pain","gain","scoret","sarea")
		#}else{
		#	colnames(auxtab)=c("id0","id1","id2","n1","n2","rmsd","scorep","pain","scoret","sarea")
		#}
		#rownames(auxtab)=NULL
		#auxo = auxtab
		#auxscore = auxo[order(auxo[,9]),c(4:5,9)] #by mean
		#auxscore = auxo[order(auxo[,9],decreasing=T),c(4:5,9)] #by mean
		#auxscore = auxo[order(auxo[,11],decreasing=T),c(4:5,11:12)] #by sum
		#auxscore = auxo[order(auxo[,"scoret"],decreasing=T),] #by sum
		#print(auxscore);readline()
		#dball$tab = auxtab
		#dball$tab = auxtab
		#dball$prescore = auxscore
		#dball$score = resume_score(score=dball$prescore)
		#auxall[[k]]$prescore = auxscore
		#auxall[[k]]$score = resume_score(score=auxall[[k]]$prescore)
		k=k+1
	}
	return(auxall)

}


verify_hot_nodes_old = function(db,dist,pdbname1,pdbname2,limit=c(1,0.9),sep="-",zero=0.1,k=3){

	
	dist.m=matrix(as.numeric(unlist(strsplit(names(dist),sep))),ncol=2,byrow=T)
	dist.m = as.data.frame(cbind(dist.m,dist))
	dist.m$a1 = 0 #;print(dist.m);readline()
	dist.m$a2 = 0
	dist.m$a = 0
	dist.id = which(dist.m[,3]<limit[1])
	tam = length(dist.id)
	if (tam>0){
		#dist.hot = tomatrix(dist.m[v,]);print(dist.hot)
		#tam = dim(dist.hot)[1]
		#print(pdbsplit(pdbname1));readline()
		#id1 = as.numeric(pdbsplit(pdbname1)[,5:6])#;print(id1)
		#id2 = as.numeric(pdbsplit(pdbname2)[,5:6])#;print(id2)
		id1 = as.numeric(pdbsplit(pdbname1)[,c("id","nclus")])#;print(id1)
		id2 = as.numeric(pdbsplit(pdbname2)[,c("id","nclus")])#;print(id2)
		d1 = diag(db[[id1[1]]]$exp$superclus[[id1[2]]]$a)#;print(d1)
		d2 = diag(db[[id2[1]]]$exp$superclus[[id2[2]]]$a)#;print(d2)
		dist.m$a1 = d1[dist.m[,1]]
		dist.m$a2 = d2[dist.m[,2]]
		for (i in dist.id){
			i1 = dist.m[i,1]
			i2 = dist.m[i,2]
			#d12 = abs(d1[i1]-d2[i2]);print(d12)
			d12 = round(min(d1[i1],d2[i2])/max(d1[i1],d2[i2]),k)
			if (d12>limit[2]){
				#dist.m[i,"a1"]=d1[i1]
				#dist.m[i,"a2"]=d2[i2]
				dist.m[i,"a"]=d12
			}else{
				print(paste("WARNING? there is doubt if these parameters (",i,dist.m[i,3],d1[i1],d2[i2],d12,") indicate a hot node in",pdbname1))
			}
		}
		#print(dist.m);readline()
		return(dist.m)
	}else{
		print(paste("WARNING: it was not found a hot nodes for",pdbname1))
		print(dist.m)#;readline()
		return(list(data.frame(V1=0,V2=0,dist=0,a1=0,a2=0,a=0,score=0)))
	}

}


toMongo_old = function(db, collection){
  
  library(mongolite)
  db_collection = mongo(url= "mongodb://localhost:27017", collection = collection, db = "ppi") # create connection, database and collection
  
  if(collection == "proteins"){
    # para cada pdb
    for(i in 1: length(db)){
      # chain chain
      if(length(db[[i]][[2]]$bsr_groups_nonpolar_nonpolar) == 0){
        # any lig
        if(length(db[[i]][[3]]$bsr_groups_nonpolar_nonpolar) == 0){
          db_collection$insert(db[[i]][[3]])    
          db_collection$insert(db[[i]][[2]])
          db_collection$insert(db[[i]][[1]])
        }else{
          db_collection$insert(db[[i]][[3]])    
          db_collection$insert(db[[i]][[1]])
          db_collection$insert(db[[i]][[2]])
        }
      }else{
        db_collection$insert(db[[i]][[2]])
        db_collection$insert(db[[i]][[1]])
        db_collection$insert(db[[i]][[3]])
      }
      
    }  
  }else{
    for(i in 1: length(db)){
      for(j in 1: length(db[[i]])){
        db_collection$insert(db[[i]][[j]])
      }
      
    }
  }
}


ToMongoDB_old2 = function(dbweb,collection){
	library(mongolite)
	db_collection = mongo(url= "mongodb://localhost:27017", collection = collection, db = "ppi") # create connection, database and collection
	
	for(i in 1 :length(dbweb)){
  
			item = dbweb[[i]]
			n1 = unlist(strsplit(as.character(item$scorend$n1), "_"))
			n2 = unlist(strsplit(as.character(item$scorend$n2), "_"))
			algorithm = substr(n1[4],1,3)
			source_polarity <- if(substr(n1[4],4,nchar(n1[4]))=="") "ALL" else substr(n1[4],4,nchar(n1[4]))
			target_polarity <- if(substr(n2[4],4,nchar(n2[4]))=="") "ALL" else substr(n2[4],4,nchar(n2[4]))

			#Monta a estrutura esperada no MongoDB
			s = structure(list(
  			source_pdb              = toupper(n1[1]), # nome do PDB que é a base de comparação
  			target_pdb              = toupper(n2[1]),
  			algorithm               = algorithm,
  			source_number_of_groups = n1[6],
  			target_number_of_groups = n2[6],
  			source_chains           = n1[2:3],
  			target_chains           = n2[2:3],
  			source_polarity         = source_polarity,
  			target_polarity         = target_polarity,
  			scoret                  = item$scorend$scoret,
  			kernel                  = item$scorend$kernel,
  			scorea                  = item$scorend$scorea,
  			matrix                  = item$matrix,
  			rot                     = item$rot,
  			hot						= item$hot,
			hotcols					= item$hotcols
			))
			
			if(collection == "hotspots"){
			  query = paste0('{"source_pdb.0" : "', toupper(n1[1]), '", "target_pdb.0" : "', toupper(n2[1]), '", "source_number_of_groups.0" : "', n1[6], '", "target_number_of_groups.0" : "', n2[6], '" } ')
			  hotspot_collection <- db_collection$find(query)
			  
			  #pode gravar pois não existe nada no banco e não vai gravar elementos duplicados
			  if(length(hotspot_collection)==0){
			    db_collection$insert(s)
			  }  
			}else{
			  #Persiste no banco de dados MongoDB
			  db_collection$insert(s)  
			}
	}
}

#Monta lista superclus onde cada cluster vira um nó. O rotulo do nó eh o volume das arestas internas
#e o rotulo das arestas do superclus eh o volume das arestas externas. Com isso, monta-se uma 
#matriz do supergrafo, que tambem pode ser vista como uma matriz de confusao. 
make_superclus_old = function(db,only_last=F,only_a=F,ma=NULL,k=3,verbose=T,geo=T,density=T){

	#aux = db
	#;print("o") 
	if (!is.null(db$A)){
		auxwclus = list()
		#;print("o") 
		#tam = 2:length(aux$exp$pam)

		auxlen = length(db$exp$pam)#;print(auxlen)
		if (only_last){
			tam = auxlen
			#tam = only_last
		}else{
			tam = 1:auxlen	
		}
		#auxa = aux$A$a
		;print(tam)#;readline() 

		for (i in tam){
		#for (i in 8){
			if (verbose) print(paste("Making superclus from pam with",i,"clusters"))
			auxids = db$exp$pam[[i]]$clustering
			#tamids = levels(auxids)
			#print(tamids);readline()
			#auxm = matrix(rep(0,i*i),i,i)
			auxwclus[[i]] = list()
			#auxwclus[[i]]$gi = 
			#for (j in 1:i){
			#auxj = which(auxids == j)
			#if (!is.null(ma)){
			#	auxsupa = ma
			#}else{		
			auxsupa = make_matrix_superclus(auxids,db$A$a)
			#}
			#print(auxsupa);readline()
			auxwclus[[i]]$a = auxsupa
			if (!only_a){
				#diag(auxsupa) = 0 E(aux$g)$weight
				auxg = graph.adjacency(auxsupa,weighted=TRUE,mode="undirected",diag=FALSE)
				#vertex_attr(auxg) = list(label=1:dim(auxsupa)[1])
				vertex_attr(auxg) = list(label=rep(1,dim(auxsupa)[1]))

				#print(i);print(E(auxg)$weight);readline()
				if (is.null(E(auxg)$weight)){
					diag(auxsupa)=1
					#print(auxsupa)
					auxg = graph.adjacency(auxsupa,weighted=NULL,mode="undirected")
					#vertex_attr(auxg) = list(label=1:dim(auxsupa)[1])
					vertex_attr(auxg) = list(label=rep(1,dim(auxsupa)[1]))		
				}# else {
				#print("ok")
				auxwclus[[i]]$ga = auxg
				#}
				if (length(levels(factor(auxsupa)))>1){
					auxcm = confusionMatrix(as.table(auxwclus[[i]]$a))
					auxwclus[[i]]$accall = round(t(as.matrix(auxcm$overal[1:2])),k)
					#print(as.matrix(auxcm$byClass));readline()
					if (!is.matrix(auxcm$byClass)){
						auxcm$byClass = round(t(as.matrix(auxcm$byClass[1:2])),k)
						#print(auxcm$byClass);readline()
						auxcm$byClass = rbind(auxcm$byClass,auxcm$byClass[1,]) 
						auxwclus[[i]]$accbyclass = auxcm$byClass
						#print(auxwclus[[i]]$accbyclass);readline()
					}else{
					#print(auxcm$byClass);readline()
						auxwclus[[i]]$accbyclass = round(auxcm$byClass[,1:2],k)
					}
				}else{
					auxwclus[[i]]$accall = t(as.matrix(c(Accuracy=1,Kappa=1)))
					auxwclus[[i]]$accbyclass = t(as.matrix(c(Sensitivity=1,Specificity=1)))
				}
				#auxwclus[[i]]$hot = rep(0,dim(auxsupa)[1])
				if (density){
					auxwclus[[i]]$density = calculate_superclus_density(db$A$a,db$element_name,auxids)
				}
			}
		}
		#print(auxwclus);readline()
		if (only_last){
			db$exp$superclus = c(db$exp$superclus,auxwclus)
		}else{
			db$exp$superclus = auxwclus
		}

		if (geo){
			db = add_geometric_center_per_super_clustering(db)
		}
		#print(db$exp$superclus);readline()
	}
	return(db)

}

#Monta lista superclus onde cada cluster vira um nó. O rotulo do nó eh o volume das arestas internas
#e o rotulo das arestas do superclus eh o volume das arestas externas. Com isso, monta-se uma 
#matriz do supergrafo, que tambem pode ser vista como uma matriz de confusao. 
make_superclus_old2 = function(db,only_last=F,only_a=F,ma=NULL,k=3,verbose=T,geo=T,density=T){

	#aux = db
	#;print("o") 
	if (!is.null(db$A)){
		auxwclus = list()
		auxlen = length(db$exp$pam)#;print(auxlen)
		if (only_last){
			auxids = db$exp$pam[[auxlen]]$clustering
			if (verbose) print(paste("Making superclus from pam for last",auxlen,"clusters"))
			auxwclus[[1]] = make_superclus_for(auxids,db$A$a,db$element_name,only_a,density,k)
			db$exp$superclus = c(db$exp$superclus,auxwclus)
		}else{
			tam = 1:auxlen
			for (i in tam){
			#for (i in 8){
				if (verbose) print(paste("Making superclus from pam with",i,"clusters"))
				auxids = db$exp$pam[[i]]$clustering
				auxwclus[[i]] = make_superclus_for(auxids,db$A$a,db$element_name,only_a,density,k)
			}
			db$exp$superclus = auxwclus
		}
		if (geo){
			db = add_geometric_center_per_super_clustering(db)
		}
		#print(db$exp$superclus);readline()
	}
	return(db)

}


calculate_superclus_density_old = function(m,e_name,ids,sep="[.]",ncol=3){

	res = list()#;print(ids)
	id = unique(ids)#;print(id)
	rex = c()
	rey = c()
	
	#for (i in id){
	for (i in 8){
		f = ids == i
		mi = m[f,f]
		gi = graph.adjacency(mi,weighted=TRUE,mode="undirected",diag=FALSE)#;print(sep_atom_name_id(e_name[f],sep=sep,ncol=ncol));readline()
		e_frame = sep_atom_name_id(e_name[f],sep=sep,ncol=ncol)
		e_table = table(e_frame[,1])
		if (length(e_table)>1){
			#parts = #;print(parts);readline()
			x = density_bipartite_graph(gi,e_table)
		}else{
			e_factor = factor(e_frame[,2])
			levels(e_factor) = 1:length(levels(e_factor))
			e_v = as.character(e_factor)
			e_table = table(e_v);print(e_table)#;readline()
			x = density_bipartite_graph(gi,e_table);print(x)#;readline()
		}
		y = edge_density(gi)
		#print(x);print(y);readline()
		rex = c(rex,x)
		rey = c(rey,y)
	}
	names(rex)=id
	names(rey)=id
	res$total = rey
	res$bipar = rex
	return(res)
}

make_all_metric_super_clustering_old = function(db,id,metricpar,lowcut,minclusid,maxclusid,join=F,only_last=0,minpar=4,pre=F,sep="_", continuum=F,dcut=c(0.2,0.01)){

	if (!is.null(db$A)){
		auxmet = make_all_super_clustering(db=db,type="met",par=c(metricpar,lowcut),id=id,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
		f = ((auxmet[,3]>=minclusid)&(auxmet[,3]<=maxclusid))
		auxmet = auxmet[f,]
		print(auxmet)#;readline()
		#print(lowcut);readline()
		lastline = dim(auxmet)[1]
		lastgo = F
		if (auxmet[lastline,3]>maxclusid){
			print(paste("WARNING: found the limit of maxclusid=",maxclusid,"for",db$pdbname,"Trimming it... "))
			if (lastline<maxclusid){
				#if ((auxmet[lastline,1]==1)&(auxmet[lastline,2]==1)){
					##auxmet = lastline_it(auxmet,lastline,minclusid,maxclusid)
				#	auxmet[lastline,3]=maxclusid
				#	rown = rownames(auxmet)[lastline]
				#	rown = unlist(strsplit(rown,sep))
				#	rown[length(rown)]=as.character(maxclusid)
				#	rown = paste(rown,collapse=sep)
				#	#print(rown)#;readline()
				#	rownames(auxmet)[2]=rown
				#	#print(auxmet);readline()
				#}else{
					auxmet = auxmet[1:lastline,]
				#}
			}else{
				auxmet = auxmet[1:maxclusid,]
			}#;print(auxmet);readline()
		}else{
			#print(auxmet[lastline,3]);minclusid;readline()
			if (auxmet[lastline,3]<minclusid) {
				print(paste("WARNING: it was not possible to reach the limit of minclusid",minclusid,"for",db$pdbname))
			}
			if (auxmet[lastline,3]==minclusid){
				print(paste("WARNING: found the limit of minclusid=",minclusid,"for",db$pdbname,"Trying overtake it... "))
				lowcuti = lowcut
				lastgo = T
				repeat{
					lowcuti[1]=lowcuti[1]-dcut[1]
					if (lowcuti[1]<0){
						lowcuti[2] = lowcuti[2]-dcut[2]
						lowcuti[1] = lowcut[1]
					}
					if (lowcuti[2]<0){
							print(paste("WARNING: could not ovetake the minclusid for",db$pdbname))
							lastgo = F
							break;
					}
					auxmet = make_all_super_clustering(db=db,type="met",par=c(metricpar,lowcuti),id=id,minpar=minpar, pre=pre,join=join,only_last=only_last,continuum=continuum)
					#print(auxmet);print(lowcuti);readline()
					lastline = dim(auxmet)[1]
					if (auxmet[lastline,3]>minclusid) {
						auxmet = auxmet[1:(minclusid+1),]
						break;
					}
				}
			}
		}
		#print(auxmet);readline()
	}else{
		auxmet = t(data.frame(c(rep(0,4),id)))
		rownames(auxmet) = paste(id,db$pdbname,"0",sep=sep)
		#print(auxmet);readline()
	}
	return(auxmet)
}

make_graph_alignment_old = function(db1,db2,t1,t2,vcores,minclus=4,maxclus=4,r=5,comp=F,frame=16,sep="-",zero=0.001){

	#print(t1);print(t2)#;readline()

	#print(paste())
	pdbname1 = as.character(t1[1,"pdbname"])#;print(pdbname1)
	xyz1 = db1[[t1[1,"pdbid"]]][[t1[1,"inter"]]]$clus[[t1[1,"pol"]]][[t1[1,"nclus"]]]$xyz
	a1 = db1[[t1[1,"pdbid"]]][[t1[1,"inter"]]]$clus[[t1[1,"pol"]]][[t1[1,"nclus"]]]$a
	g1 = graph.adjacency(a1,weighted=TRUE,mode="undirected",diag=FALSE)
	pdbname2 = as.character(t2[1,"pdbname"])#;print(pdbname2)
	xyz2 = db2[[t2[1,"pdbid"]]][[t2[1,"inter"]]]$clus[[t2[1,"pol"]]][[t2[1,"nclus"]]]$xyz
	a2 = db2[[t2[1,"pdbid"]]][[t2[1,"inter"]]]$clus[[t2[1,"pol"]]][[t2[1,"nclus"]]]$a
	g2 = graph.adjacency(a2,weighted=TRUE,mode="undirected",diag=FALSE)

	print(paste("Initiating align between",pdbname1," - ",pdbname2))

	p1 = list()
	p1$xyz = xyz1
	p1$a = a1
	p1$g = g1

	p2 = list()
	p2$xyz = xyz2
	p2$a = a2
	p2$g = g2

	#print(xyz1);print(xyz2);readline()

	rotlist = generate_all_new_alignments_par(target=xyz1,source=xyz2,minclus=minclus,maxclus=maxclus,frame=frame,pdb1=pdbname1,pdb2=pdbname2, ncores=vcores[1])

	tamlist = length(rotlist)
	##### ker.global = score_kernel(g1,g2)#;print("o") ####
	score.local = NULL
	base=0
	for (i in 1:tamlist){### tamlist sempre 1 ###
	#for (i in 2){
		#print("o")
		minclus = rotlist[[i]]$minclus
		framer = 1:frame
		framet = (frame+1):(frame+minclus)
		frames = (frame+minclus+1):(frame+2*minclus)
		#frame = frame+2*minclus
		rotmatrix = rotlist[[i]]$rotmatrix#;print(rotmatrix)
		tamatrix = dim(rotmatrix)[1]#;print(tamatrix);readline()
		### ESTUDAR PARALELISMO SOBRE tamatrix ###
		for (j in 1:tamatrix){
			rot0 = rotmatrix[j,framer]#;print(rot0)
			tid0 = rotmatrix[j,framet]#;print(tid0)
			sid0 = rotmatrix[j,frames]#;print(sid0)
			xyz1i = xyz1[tid0,]
			xyz2i = xyz2[sid0,]
			#xyz1i = scale(xyz1i,s=F)
			#xyz2i = scale(xyz2i,s=F)
			xyz2t0i = transform_by_rot_matrix(xyz2i,rot0)#;print("o")
			xyz2t0 = transform_by_rot_matrix(xyz2,rot0)
			#p2$xyz = xyz2t0
			#par(mfrow=c(1,2));limxyz=get_lim_xyz(p1$xyz,p2$xyz);plot_g(p1,limxyz);plot_g(p2,limxyz)
			dist0i = find_nearest_nodes(xyz1i,xyz2t0i)
			#print(dist0i)#;readline()
			if (!is.null(dist0i)){
				labels = fix_labels(labels=names(dist0i),id1=tid0,id2=sid0)
				names(dist0i) = labels;
				#print(dist0i);readline()
				dist0 = expand_to_nearest_global_nodes(dist0i,xyz1,xyz2t0,tid0,sid0)
				#print(dist0)#;readline()
				if(!is.null(dist0)){	
					score0 = make_score_graph_superimposition(dist0,g1,g2,a1,a2,xyz1,xyz2t0)
					#print(score0);readline()
					base0 = score0$scor[1]
					if (base0 > base){
						ijmax = c(i,j)
						rot = rot0
						#dist = dist0
						score = score0
						base = base0
						xyz2t = xyz2t0
					}
				}else{
					print(paste("WARNING: it was impossible to do a symmetric align with expand dist"))
				}
			}else{
				#print(paste("WARNING: it was impossible to do a symmetric align with dist"))
			}
		}
	}
	#print(score$dist);
	#print(score$scor);
	p2$xyz = xyz2t
	#par(mfrow=c(1,2));limxyz=get_lim_xyz(p1$xyz,p2$xyz);plot_g(p1,limxyz);plot_g(p2,limxyz)
	#readline()
	res = list()
	res$ijmax = ijmax
	res$score = score
	res$rot = rot
	#res$xyz2 = xyz2t
	res$t1 = t1
	res$t2 = t2
	res$p1 = p1
	res$p2 = p2
	print(paste("Finish align between",pdbname1," - ",pdbname2))
	return(res)
}


scan_rotations_old = function(j,vbin,rotmatrix,framer,framet,frames,p1,p2){

	rangej = which(vbin==j)

	xyz1 = p1$xyz
	a1 = p1$a
	g1 = p1$g

	xyz2 = p2$xyz
	a2 = p2$a
	g2 = p2$g

	base = 0

	for (i in rangej){#;print(i)
	#for (i in 1300){
		rot0 = rotmatrix[i,framer]#;print(rot0)
		tid0 = rotmatrix[i,framet]#;print(tid0)
		sid0 = rotmatrix[i,frames]#;print(sid0)
		xyz1i = xyz1[tid0,]
		xyz2i = xyz2[sid0,]
		#xyz1i = scale(xyz1i,s=F)
		#xyz2i = scale(xyz2i,s=F)
		xyz2t0i = transform_by_rot_matrix(xyz2i,rot0)#;print("o")
		xyz2t0 = transform_by_rot_matrix(xyz2,rot0)
		#p2$xyz = xyz2t0;par(mfrow=c(1,2));limxyz=get_lim_xyz(p1$xyz,p2$xyz);plot_g(p1,limxyz);plot_g(p2,limxyz)
		dist0i = find_nearest_nodes(xyz1i,xyz2t0i)#;print(dist0i);readline()
		print(dist0i)#;readline()
		if (!is.null(dist0i)){
			labels = fix_labels(labels=names(dist0i),id1=tid0,id2=sid0)#;print(labels)
			names(dist0i) = labels;
			#print(dist0i);readline()
			dist0 = expand_to_nearest_global_nodes(dist0i,xyz1,xyz2t0,tid0,sid0)
			print(dist0)#;readline()
			if(!is.null(dist0)){	
				score0 = make_score_graph_superimposition(dist0,g1,g2,a1,a2,xyz1,xyz2t0)
				#print(score0);readline()
				base0 = score0$scor[1]
				if (base0 > base){
					jmax = i
					rot = rot0
					score = score0
					base = base0
					xyz2t = xyz2t0
				}
			}else{
				print(paste("WARNING: it was impossible to do a symmetric align with expand dist"))
				#score0 = make_score_graph_superimposition(dist0i,g1,g2,a1,a2,xyz1,xyz2t0)
				#print(score0);readline()
			}
		}else{
			#print(paste("WARNING: it was impossible to do a symmetric align with dist"))
		}
	}
	#p2$xyz = xyz2t
	#par(mfrow=c(1,2));limxyz=get_lim_xyz(p1$xyz,p2$xyz);plot_g(p1,limxyz);plot_g(p2,limxyz)
	#readline()
	res = list()
	res$jmax = jmax
	res$score = score
	res$rot = rot
	res$xyz2 = xyz2t
	#res$t1 = t1
	#res$t2 = t2
	#res$p1 = p1
	#res$p2 = p2
	return(res)
}


scan_rotations_old2 = function(j,vbin,rotmatrix,framer,framet,frames,p1,p2){

	rangej = which(vbin==j)

	xyz1 = p1$xyz
	a1 = p1$a
	g1 = p1$g

	xyz2 = p2$xyz
	a2 = p2$a
	g2 = p2$g

	base = 0

	for (i in rangej){#;print(i)
	#for (i in 5){
		rot0 = rotmatrix[i,framer]#;print(rot0)
		tid0 = rotmatrix[i,framet]#;print(tid0)
		sid0 = rotmatrix[i,frames]#;print(sid0)
		xyz1i = xyz1[tid0,]
		xyz2i = xyz2[sid0,]
		#xyz1i = scale(xyz1i,s=F)
		#xyz2i = scale(xyz2i,s=F)
		xyz2t0i = transform_by_rot_matrix(xyz2i,rot0)#;print("o")
		xyz2t0 = transform_by_rot_matrix(xyz2,rot0)
		#p2$xyz = xyz2t0;par(mfrow=c(1,2));limxyz=get_lim_xyz(p1$xyz,p2$xyz);plot_g(p1,limxyz);plot_g(p2,limxyz)
		dist0i = find_nearest_nodes(xyz1i,xyz2t0i)#;print(dist0i);readline()
		#print(dist0i)#;readline()
		if (!is.null(dist0i)){
			labels = fix_labels(labels=names(dist0i),id1=tid0,id2=sid0)#;print(labels)
			dist1i = dist0i
			names(dist1i) = labels;
			#print(dist1i)#;readline()
			dist = expand_to_nearest_global_nodes(dist1i,xyz1,xyz2t0,tid0,sid0)
			#print(dist)#;readline()
			if(!is.null(dist)){
				local = F	
				score0 = make_score_graph_superimposition(dist,g1,g2,a1,a2,xyz1,xyz2t0)
				#print(score0);readline()
			}else{
				local = T
				print(paste("WARNING: it was impossible expand dist. Making only local dist..."))
				#xyz1i = xyz1[tid0,]
				#xyz2t0i = xyz2t0[sid0,]
				a1i = a1[tid0,tid0]
				a2i = a2[sid0,sid0]
				g1i = graph.adjacency(a1i,weighted=TRUE,mode="undirected",diag=FALSE)
				g2i = graph.adjacency(a2i,weighted=TRUE,mode="undirected",diag=FALSE)
				score0 = make_score_graph_superimposition(dist0i,g1i,g2i,a1i,a2i,xyz1i,xyz2t0i)
				#print(score0);readline()
			}
			base0 = score0$scor[1]
			if (base0 > base){
				jmax = i
				rot = rot0
				score = score0
				if (local) score = correct_score_labels(score,tid0,sid0)
				#print(score);readline()
				base = base0
				xyz2t = xyz2t0
			}
		}else{
			#print(paste("WARNING: it was impossible to do a symmetric align with dist"))
		}
	}
	#p2$xyz = xyz2t
	#par(mfrow=c(1,2));limxyz=get_lim_xyz(p1$xyz,p2$xyz);plot_g(p1,limxyz);plot_g(p2,limxyz)
	#readline()
	res = list()
	res$jmax = jmax
	res$score = score
	res$rot = rot
	res$xyz2 = xyz2t
	#res$t1 = t1
	#res$t2 = t2
	#res$p1 = p1
	#res$p2 = p2
	return(res)
}


get_info_connections_group_old = function(conn,local,global){

	tam = length(conn)
	res=list()

	for (i in 1:tam){
		
		res[[i]] = list()#;print(conn[[i]])
		res[[i]]$a = conn_to_adj_matrix(conn[[i]])
		res[[i]]$xyz = conn_to_xyz_matrix(conn[[i]])#;print(res);readline()
		res[[i]]$sensitivity = local[[i]]
		res[[i]]$accuracy = global[[i]]
		#res[[i]]$ga = graph.adjacency(res[[i]]$a,weighted=TRUE,mode="undirected",diag=FALSE)
		#vertex_attr(res[[i]]$ga) = list(label=rep(1,i))
		#print(res);readline()
	}

	return(res)

}

#' Evenly distributed n points on a sphere
#'
#' @description Distributes n points on a sphere in a relatively even fashion following the generalised Fibonacci algorithm, described at http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/
#'
#' @param n number of points to be placed on the sphere
#' @param r radius
#' @param out.xyz logical flag specifying whether to return Cartesian coordinates (default is TRUE)
#' @param out.sph logical flag specifying whether to return spherical coordinates (default is FALSE); \code{theta}=colatitute (=polar angle=angle from z-axis), \code{phi}=longitude (=azimuthal angle),
#'
#' @examples
#'
#' ## plot standard projections of a 1000-point Fibonacci sphere
#' xyz = fibonaccisphere()
#' plot(xyz, asp=1, pch=16, cex=0.5)
#'
#' @author Danail Obreschkow
#'
#' @seealso \code{\link{runif3}}
#'
#' @export

fibonaccisphere = function(n=1000, r=1, out.xyz=TRUE, out.sph=FALSE) {
  
  if (n<1 | round(n)!=n) stop('n must be a positive integer')
  if (!out.xyz & !out.sph) stop('either out.xyz and/or out.sph must be TRUE')
  
  goldenratio = (1+sqrt(5))/2
  i = seq(n)-0.5
  z = 1-2*i/n # z-coordinate for unit sphere
  theta = acos(pmax(-1,pmin(1,z))) # polar angle
  phi = (2*pi*i/goldenratio)%%(2*pi)
  
  if (out.xyz) {
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*z
    if (out.sph) {
      out = cbind(x=x,y=y,z=z,theta=theta,phi=phi)
    } else {
      out = cbind(x=x,y=y,z=z)
    }
  } else {
    out = cbind(theta=theta,phi=phi)
  }
  
  return(out)
  
}

	#auxlist = 


	#if (!is.matrix(auxlist)){
	#	auxlist = t(as.matrix(auxlist))
	#}
	#print(auxlist[1,])
	
	##script = paste(workdir,script,sep="")
	#command = paste(workdir,command,sep="")
	##name1 = paste(workdir,name1,sufix,sep="")
	##name2 = paste(workdir,name2,sufix,sep="")
#	auxrmsd = c()
#	auxwin = c()
#	auxrot = c()
#	auxsym = c()
#	auxr = c()
	#colrot = 5:20
	#auxrotlist = list()
	##auxlist = list()

	#i=minwin		
	#i = 13
	##flags = c("-cq")
	#args = c("-cp",script,name1,name2,i,guide) #ultimo se guide=0 (todos os atomos), guide=1 so CA
	##args = c(name1,name2,i,guide)
	#auxlist = align_graph_pdb(command,args)
	##print(command);print(flags);print(script);print(args);readline()
	##auxlist = align_graph_pdb(command,flags,script,args)
	#print(length(auxlist));print(auxlist);readline()
	#auxrmsd = auxlist[[1]]$rmsd
	#print(auxlist[1,])
	#print(colnames(auxlist))
	#auxrmsd = auxlist[1,"rmsd"]
	#print(auxrmsd);readline()
	#auxwin = auxlist[1,"win"]
	#auxrot = auxlist[1,colrot]
	#auxsym = auxlist[1,"sym"]
#	auxdist = 0
	#print(name1);print(name2)
#	if (crit=="rmsd") auxmin = auxrmsd
#	if (crit=="dist"){
		#print(xyz2);print(auxrot);print(matrix(auxrot,ncol=4,byrow=T))#;readline()
		#x = matrix(auxrot,ncol=4,byrow=T)
		#x = data.matrix(x)
		#print(x[4,]);readline()
#		xyz2a = fit_xyz_alignment(xyz1,xyz2)
#		xyz2t = transform_by_rot_matrix(xyz2,matrix(auxrot,ncol=4,byrow=T))#;print(i)
		#print(xyz2t);readline()
#		auxdist = find_nearest_nodes(xyz1,xyz2t)
		#auxscore = sum(auxdist)
		#print(xyz1);print(xyz2)
		#if(rgl.cur()!=0) rgl.close()
		#plot3d(xyz1,box=F,type="s",col="blue",add=T)
		#plot3d(xyz2t,box=F,type="s",col="red",add=T)#;readline()
#		auxscore = make_score_graph_superimposition(auxdist,g1,g2,a1,a2,w)
		#print(i);print(auxrmsd);print(auxscore);readline()
		#auxmin = auxscore$scor
#		auxbase = auxscore$scor
#	}
#	j=1#;print("oi")	
#	while(auxlist[j,"rmsd"] !=0 ){
		#i=i+1
		#args = c("-cq",script,name1,name2,i)
		#args[5] = i
		#print(args);
		##args[3] = i
		#print(args);readline()
		#auxlist = align_graph_pdb(command,args)
		#print(i)
		##auxlist = c(auxlist,align_graph_pdb(command,flags,script,args))
		#print(length(auxlist));readline()#;print(auxlist);readline()
#		j=j+1
		#print(auxlist[[j]]);readline()

#		if (auxlist[j,"rmsd"] == 0) break;
#		if (crit=="rmsd"){
#			if (auxlist[j,"rmsd"] < auxmin){
#				auxrmsd = auxlist[j,"rmsd"]
#				auxwin = auxlist[j,"win"]
#				auxrot = auxlist[j,colrot]
#				auxmin = auxlist[j,"rmsdj"]
#				auxsym = auxlist[j,"sym"]
#			}
#		}
#		if (crit=="dist"){
#			auxrmsd0 = auxlist[j,"rmsd"]
#			auxrot0 = auxlist[j,colrot]
#			auxsym0 = auxlist[j,"sym"]
#			xyz2t0 = transform_by_rot_matrix(xyz2,matrix(auxrot0,ncol=4,byrow=T))#;print(i)
#			auxdist0 = find_nearest_nodes(xyz1,xyz2t0)#;print("oi")
			#auxscore0 = sum(auxdist0)
			#auxscore0 = median(auxdist0)
#			auxscore0 = make_score_graph_superimposition(auxdist0,g1,g2,a1,a2,w)#;print(auxscore0);readline()
			#if(rgl.cur()!=0) rgl.close()
			#plot3d(xyz1,box=F,type="s",col="blue",add=T)
			#plot3d(xyz2t,box=F,type="s",col="red",add=T)#;readline()
			#print(j);print(auxrmsd0);print(auxscore0);print(auxmin);readline()
			#if (auxscore0$scor < auxmin){
#			if (auxscore0$scor > auxbase){
			#print(xyz2);print(auxrot);readline()
#				auxrmsd = auxrmsd0
#				auxwin = auxlist[j,"win"]
#				auxrot = auxrot0
#				auxdist = auxdist0
#				auxscore = auxscore0
				#auxmin = auxscore0$scor
#				auxbase = auxscore0$scor
#				auxsym = auxsym0
#				xyz2t = xyz2t0
				#xyz2t = transform_by_rot_matrix(xyz2,matrix(auxrot0,ncol=4,byrow=T))
				#plot3d(xyz1,box=F,type="s",col="blue",add=F)
				#plot3d(xyz2t,box=F,type="s",col="red",add=T);readline()
#			}
#		}
#	}
	#auxlist$dist = round(auxdist,r)
	#auxlist$score = round(auxscore,r)
	#auxlist$score = mean(auxlist$dist)
#	auxr$score = auxscore
#	auxr$rmsd = round(auxrmsd,r)
#	auxr$win = auxwin
#	auxr$rot = round(auxrot,r)
	#auxr$rotref = round(rotref,r)
#	auxr$sym = auxsym
#	auxr$score$pain = verify_penalty(auxscore$distn,a1,a2)
#	if (comp) auxr$score$gain = sum(auxscore$nocv) - sum(auxscore$nodv)
	#auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)
#	if (comp) auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)+auxr$score$gain
#	else auxr$score$scot = auxr$score$scor-sum(auxr$score$pain)
	#auxr$score$sco1 = auxr$score$scot/(length(auxr$score$scov)*4)
	#auxzz = transform_by_rot_matrix(auxyz2,auxrr)
#	if (!is.null(rotref)){
#		auxr$supxyz = transform_by_rot_matrix(xyz2t,rotref)
#	}else{
#		auxr$supxyz = xyz2t
#	} 
	#print(x);readline()
	#print(auxlist);readline()
#	return(auxr)
#}

#PRE#
#recover_transformation = function(out,id,r=5){

#	auxrot = out[id[2]:id[3]]
#	auxrot = str_c(auxrot,collapse="")
	#print(nchar(auxrot))
#	auxrot = str_sub(auxrot,2,(nchar(auxrot)-1))
#	auxrot = round(as.numeric(str_split(auxrot,",")[[1]]),r)
#	return(auxrot)
#}

#PRE#
#transrotation = function(rotall, rot, rot0){

#	auxrotall = rotall
#	auxrot = matrix(rot,ncol=4,byrow=T)
#	auxrotall[1:3,1:3] = auxrot[1:3,1:3] %*% auxrotall[1:3,1:3] 
	#auxrotall[,4] = auxrot[,4]
	#auxrotall[4,] = auxrot[4,]
#	auxrotall[,4] = rot0[,4]
#	auxrotall[4,] = rot0[4,]
	#print(auxrot)
#	return(auxrotall)

#}

#PRE#
#align_graph_pdb = function(command,args,idline=c(35,37,52)){
#align_graph_pdb = function(command,args,idline=c(31,33,48)){
#align_graph_pdb = function(command,flags,script,args,idline=c(32,34,49),dline=19,tam=10,v1=".pml",v2="a.pml",r=1,dr=0.001){
#align_graph_pdb = function(command,flags,script,args,idline=c(36,38,53),dline=19,tam=10,v1=".pml",v2="a.pml",r=1,dr=0.001){

#	auxrmsd = c()
#	auxwin = c()
#	auxrot = c()
#	auxlist = list()
	
	#print(command);readline()
	#print(args);readline()
#	args1 = c(flags,paste(script,v1,sep=""),args)
	#print(args);readline()
#	auxout = system2(command,args=args1,stdout=T)
#	auxline = unlist(strsplit(auxout[idline[1]]," "))
	#print(auxout);readline()
	#print(auxline);readline()
#	auxlist[[1]] = list()
#	if (auxline[1]=="RMSD"){
#		auxnum = as.numeric(auxline[2])
#		auxrot = recover_transformation(auxout,idline)
		#print(auxrot)#;readline()
#		idline2 = idline+dline
#		auxline = unlist(strsplit(auxout[idline2[1]]," "))
#		if (auxline[1]=="RMSD"){
#			auxnum2 = as.numeric(auxline[2])
#			auxrot2 = recover_transformation(auxout,idline2)
#			auxrot2a = matrix(auxrot2,ncol=4,byrow=T)[1:3,1:3]
			#print(auxrot2a)#;readline()
			#if (auxrot2a==diag(3)) {
#			auxtest = sum(abs(auxrot2a-diag(3)))
			#print(auxtest);readline()
#			if (auxtest<dr){
				#print("one");
#				auxlist[[1]]$rmsd = auxnum
#				auxlist[[1]]$rot = auxrot
#				auxlist[[1]]$sym = 0
#			}else{
				#print("two");
#				i=1
				#print(auxrot);
#				args2 = c(flags,paste(script,v2,sep=""),args)
#				auxout = system2(command,args=args2,stdout=T)
				#print(auxout);readline()
				#idline = idline+8
#				auxrotall = diag(4)
#				auxline = unlist(strsplit(auxout[idline[1]]," "))
#				auxnum = as.numeric(auxline[2])
#				auxrot = recover_transformation(auxout,idline)
#				auxrotall = matrix(auxrot,ncol=4,byrow=T)
				#auxrotall = transrotation(auxrotall,auxrot)
				#transform_by_rot_matrix(y,matrix(r,ncol=4,byrow=T))
				#print(auxrotall);readline()
#				auxlist[[i]]=list()
#				auxlist[[i]]$rmsd = auxnum
#				auxlist[[i]]$rot = as.vector(t(auxrotall))
#				auxlist[[i]]$sym = i
				#auxrotall2 = diag(4)
#				idline = idline+dline
#				auxrot0 = auxrotall
#				for (i in 2:tam){
					#print(idline);readline()
#					auxline = unlist(strsplit(auxout[idline[1]]," "))
					#print(auxline)
#					auxnum = as.numeric(auxline[2])
#					auxrot = recover_transformation(auxout,idline)
#					auxrotall = transrotation(auxrotall,auxrot,auxrot0)
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
#					auxlist[[i]]=list()
#					auxlist[[i]]$rmsd = auxnum
#					auxlist[[i]]$rot = as.vector(t(auxrotall))
#					auxlist[[i]]$sym = i
#					idline = idline+dline
					
#				}
#			}
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
#		}else{
#			print(paste("WARNING: something wrong in cealign"))
#		}
#	}else{
		#print(paste("WARNING: RMSD line not found in cealign output for",args))
#		auxlist[[1]]$rmsd = 0
#		auxlist[[1]]$rot = as.vector(diag(4))
#		auxlist[[1]]$sym = 0	
#	}
#	return(auxlist)

#}

#PRE#
#pre_align_graph_pdb = function(name1,name2,id1,id2,command="pymol",script="script-align-v2",workdir="Pymol/",minwin=3,sufix=".pdb",r=5,guide=0){

#	script = paste(workdir,script,sep="")
	#command = paste(workdir,command,sep="")
#	name1 = paste(workdir,name1,sufix,sep="")
#	name2 = paste(workdir,name2,sufix,sep="")
#	auxrmsd = c()
#	auxwin = c()
#	auxrot = c()
#	auxsym = c()
	#auxrotlist = list()
#	auxlist = list()

#	auxa=c()

#	i=minwin		
	#i = 13
#	flags = c("-cq")
	#args = c("-cp",script,name1,name2,i,guide) #ultimo se guide=0 (todos os atomos), guide=1 so CA
#	args = c(name1,name2,i,guide)
	#auxlist = align_graph_pdb(command,args)
	#print(command);print(flags);print(script);print(args);readline()


#	while(1){

#		auxlist = align_graph_pdb(command,flags,script,args)
		#if (auxlist[[1]]$rmsd == 0) break;
		#print(auxlist);print(unlist(auxlist));readline()
#		auxlen = length(auxlist)
#		if (auxlen>1){
#			for (j in 1:auxlen){
#				auxr = c(id1=id1,id2=id2,win=i,unlist(auxlist[j]))
				#print(auxr);readline()
#				auxa = rbind(auxa,auxr)
#			}

#		}else{
			#if (auxlist[[1]]$rmsd == 0) break;
#			auxr = c(id1=id1,id2=id2,win=i,unlist(auxlist))
#			auxa = rbind(auxa,auxr)
#			if (auxlist[[1]]$rmsd == 0) break;
#		}
		#auxr=unlist(auxlist)
		#print(length(auxlist));readline()
		#auxa = rbind(auxa,auxr)
		#print(auxa);readline()
#		args[3] = i
#		i=i+1

#	}
	#print(auxa);readline()
#	return(auxa)
#}

#PRE#
#pre_computation_of_alignments = function(db,dball,generalpath,version){

#	auxall=list()
	#tamk = dim(dball$group)[1]#;tam=3
	#tam = c(1)
	#pdblist = levels(dball$group$pdb)
#	pdblist = unique(dball$group$pdb)
#	tami = length(pdblist)
	#print(pdblist);readline()
	#auxa = c()

#	for (i in 1:tami){
#		pdbi = pdblist[i]
		#pdbname = as.character(pdbi)
		#tamj = which(dball$group$pdb==pdbi)
		#auxa = c()
#		print(paste("Pre-computation alignment for",pdbname))
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

#pre_computation_of_alignments = function(db,dball,outdir="Datalign/",sufix="align.csv",sep="_"){
#pre_computation_of_alignments_old = function(db,dball,outdir="Datatest/",sufix="align.csv",sep="_"){

#	auxall=list()
#	tamk = dim(dball$group)[1]#;tam=3
	#tam = c(1)
	#pdblist = levels(dball$group$pdb)
#	pdblist = unique(dball$group$pdb)
#	tami = length(pdblist)
	#print(pdblist);readline()
	#auxa = c()

#	for (i in 1:tami){
#		pdbi = pdblist[i]
#		pdbname = as.character(pdbi)
#		tamj = which(dball$group$pdb==pdbi)
#		auxa = c()
#		print(paste("Pre-computation alignment for",pdbname))
#		for (j in tamj){
#			auxnamei = as.character(dball$group[j,"fakename"])			
#			for (k in 1:tamk){
				#j=5
#				auxnamek = as.character(dball$group[k,"fakename"])
#				#print(auxnamei);print(auxnamek)
#				#print(sprintf("%s - %s: %.3f - %.3f",auxnamei,auxnamek,j,k));readline()
#				auxr=pre_align_graph_pdb(auxnamei,auxnamek,j,k)
#				#print(auxr[dim(auxr)[1],])
#				#print(dim(auxr)[1]);readline()
#				auxa = rbind(auxa,auxr)
#				print(sprintf("%s - %s: %i - %i - %i",auxnamei,auxnamek,j,k,dim(auxa)[1]))#;readline()
#			}
#		}
#		filename = paste(outdir,toupper(pdbname),sep,sufix,sep="")
#		#print(filename);readline()
#		write.table(auxa,filename,row.names=F,col.names=T,sep=",")
#		#print(filename);readline()
#	
#	}
#
#}

#make_all_self_similarity_assessment = function(db,dball,i,pdblist,dbhot=NULL,rotref=diag(4),self=T){
#
#	auxr = list()
#	pdbname = toupper(pdblist[i])
#	auxr$pdbname = pdbname
#	auxgrp = which(dball$group[,"pdb"]==pdbname)
#	auxr$groupid = auxgrp
#	auxr$alignment = list()
#	auxr$alignment = new_make_similarity_assessment(db=db,dball=dball,dbhot=dbhot,id=auxgrp,pdbname=pdbname,rotref=rotref,w=c(1,1),self=self)
#	return(auxr)
#
#}

