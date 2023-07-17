
proper=function(x) {
	paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

}

no_word = function(p=NULL,ext=NULL,enc=NULL){

	if (!is.null(ext)){
		p = c(p,ext)
		return (p)
	}
	prep = c("no","na","de","do","da","em","nos","nas","dos","das","uma","um","umas","uns","lhe","lhes")
	art = c("a","o","e","à","às")
	pron = c("esse","este","essa","esta","isso","isto","aquele","aquela","aquilo","esses","estes","essas","estas","aqueles","aquelas")
	pron2 = c("sua","seu","suas","seus","nesse","nessa","nisso","nisto","nesses","nessas","teu","tua")
	conj = c("que","se","mas","por","pelo","pela","com","para","como")
	other = c("ao","aos","figura","the")
	p = c(prep,art,pron,pron2,conj,other)
	p2 = sapply(p,proper)
	names(p2)=NULL
	p = c(p,p2)
	if (!is.null(enc)){
		Encoding(p)=enc
	}
	#print(p);readline()
	return(p)

}


filter_word = function(v,p=NULL){

	if (is.null(p)){
		p = no_word()
	} else {
		p = no_word()
	}
	auxv = v %in% p
	#auxr = v[auxv]
	return(auxv)

}

filter_cloud = function(m,inter,p=NULL){

	auxv = filter_word(m[,2],p)
	
	auxm = m[!auxv,]
	#auxv = filter_word(aux[,2],p)
	#print(auxv);readline()
	#return(aux[!auxv,])
	return(auxm[inter,])

}

norm_word = function(m,nor){

	auxv1 = m[,2] %in% nor[,1]
	auxv2 = nor[,1] %in% m[,2]
	#print(auxv)
	m[auxv1,1] = m[auxv1,1] - nor[auxv2,2]
	return(m)

}

norm_cloud_lower = function(m){

	auxs = tolower(m[,2])
	auxm = data.frame(m[,1],auxs)
	tam = length(auxs)
	auxr = data.frame(a=0,b="")

	for (i in 1:tam){
		#print(i)
		#print(auxs[i])
		#auxv = auxs[i] %in% auxm[,2]
		auxv = auxs %in% auxs[i]
		if (auxm[i,1]!=0){
			auxx = sum(auxm[auxv,1])
			#print(auxv)
			auxids = which(auxv==T)
			auxy = auxm[auxids[1],2]
			#print(auxids)
			auxm[auxids[1:length(auxids)],1]=0
			#print(auxx)
			auxdf = data.frame(a=auxx,b=auxy)
			#print(auxdf)
			auxr = rbind(auxr,auxdf)
			#print(auxr);readline()
		}
	}
	auxr = auxr[2:dim(auxr)[1],]
	auxids = order(auxr[,1],decreasing=T)
	return(auxr[auxids,])
}


remove_words = function(m,w,qte=NULL,ord=T){
	
	if (!is.null(qte)){
		w.df = data.frame(w,qte)
		aux = norm_word(m,w.df)
	}else{
		auxv1 = m[,2] %in% w
		m[auxv1,1] = m[auxv1,1] - m[auxv1,1]
		aux = m
	}
	#aux = norm_word(m,w.df)
	if (ord){
		auxo = order(aux[,1],decreasing=T)
		auxr = aux[auxo,]
		return(auxr)
	}else{
		return(aux)
	}

}






