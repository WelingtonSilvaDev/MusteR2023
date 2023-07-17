
#####################################################################
#Biblioteca funções para geração em massa de scripts Pymol
#####################################################################

# OBSOLETE: Nova versão da pymol_save_enz e pymol_save_inb
pymol_save = function(name,version="-v1",append="",file=T){
	#salvando imagem
	dir = "/home/profcaveira/Documents/UFMG/Bioinfo/Nilma/"
	img = "Imagens/ALL/"
	str01 = "set ray_shadows, off\n"
	str02 = "ray 1024,768\n"
	str03 = "png "
	pym01 = paste(str03,dir,img,name,append,version,".png","\n",sep="")
	script = paste(str01,str02,pym01,"\n",sep="")
	str04 ="Scripts-Pymol/"
	if (file)
		write(script,file=paste(dir,str04,"script-img-all-int-",name,append,".pml",sep=""))
	else
		return(script)
}

# OBSOLETE: Nova versão pymol_select_enz e pymol_select_inb unificada
pymol_select = function(map,mpo,name,type,color_ap,color_po,color_al,model,append="",file=T){

	#select = c()
	#script = c()
	res = "resi "
	n = "n. "
	op = "("
	cp = ")"
	plus = "+"
	and = " and "
	not = " not "
	hif = "-"
	bar ="/"
	dir = "/home/profcaveira/Documents/UFMG/Bioinfo/Nilma/"
	sub01 = "PDB/"
	sub02 = paste(type,bar,sep="")
	sub03 = "Scripts-Pymol/"
	chainX = substr(type,1,1)
	int = "INT"
	spc = " "
	comma = ", "
	nil = "NIL"
	hide = "hide "
	lines = "lines "
	load = "load "
	bg_color = "bg_color "	
	create = "create "
	chain = " chain "
	solvent = "solvent"
	delete = "delete "
	color = "color "
	rebuild = "rebuild "
	set = "set "
	rayob = "ray_opaque_background, 0"
	pdbsuf = ".pdb"
	select = "select "
	apolar = "apolar "
	polar = "polar "
	all = "all "
	no = "no"
	stickrad = "stick_radius, 0.25 "
	sphere = "sphere "
	stick = "stick "
	both = "both "
	show = "show "
	run = "run "
	py = ".py"
	spectrumany = "spectrumany"
	param03 = "b "
	pml = ".pml"
	scriptpre = "script-"
	spheretrans = "sphere_transparency=0.5 "
	#, white red, " 	
	
	white = "white "
	cyan = "cyan "
	green = "green "
	#str00 = paste(name,"-",type,"-NIL",sep="")
	#str00 = paste(name,"-INB-NIL",sep="")
	

	pdbrad = paste(name,hif,type,sep="")
	pdball = paste(pdbrad,pdbsuf,sep="")
	#pdb = paste(str00,".pdb",sep="")
	#str01 = "load "

	#pym00 = paste(str01,dir,"PDB/INB/",pdb,"\n",sep="")
	pym00 = paste(load,dir,sub01,sub02,pdball,"\n",sep="")
	
	pym01 = paste(bg_color,white,"\n",sep="")

	param01 = paste(and,chain,chainX,and,not,solvent,sep="")
	#ename = paste(chainX,hif,int,hif,name,sep="")
	ename = paste(pdbrad,hif,int,sep="")
	pym02 = paste(create,ename,comma,pdbrad,and,chain,chainX,and,not,solvent,"\n",sep="")

	pym03 = paste(delete,pdbrad,"\n",sep="")

	pym04 = paste(hide,lines,comma,ename,"\n",sep="")
	
	pym05 = paste(set,rayob,"\n",sep="")

	#str05 = "color cyan, "
	pym06 = paste(color,color_al,comma,ename,"\n",sep="")

	#str06 = "show spheres, "
	#str06 = "show dots, "
	#pym06 = paste(str06,ename,"\n",sep="")
	
	pym07 = paste(rebuild,"\n",sep="")	

	scriptin = paste(pym00,pym01,pym02,pym03,pym04,pym05,pym06,pym07,"\n",sep="")

	selectap = paste(op,sep="")
	selectal = paste(op,sep="")
	lenap = dim(map)[1]
	for (i in 1:lenap){
		if (i<lenap){
			selectap = paste(selectap,op,res,map[i,2],and,n,map[i,3],cp,plus,sep="")
			selectal = paste(selectal,op,res,map[i,2],cp,plus,sep="")
		}
		else{ 
			selectap = paste(selectap,op,res,map[i,2],and,n,map[i,3],cp,cp,"\n",sep="")
			selectal = paste(selectal,op,res,map[i,2],cp,plus,sep="")
		}
	}
	
	#inameap = paste("ifrInb",name,"ap",sep="")
	nameap = paste(ename,hif,apolar,sep="")
	pym08 = paste(select,nameap,comma,ename,and,selectap,"\n",sep="")
	
	selectpo = paste(op,sep="")
	lenpo = dim(mpo)[1]
	for (i in 1:lenpo){
		if (i<lenpo){
			selectpo = paste(selectpo,op,res,mpo[i,2],and,n,mpo[i,3],cp,plus,sep="")
			selectal = paste(selectal,op,res,mpo[i,2],cp,plus,sep="")
		}
		else{ 
			selectpo = paste(selectpo,op,res,mpo[i,2],and,n,mpo[i,3],cp,cp,"\n",sep="")
			selectal = paste(selectal,op,res,mpo[i,2],cp,cp,"\n",sep="")
		}
	}
	namepo = paste(ename,hif,polar,sep="")
	pym09 = paste(select,namepo,comma,ename,and,selectpo,"\n",sep="")

	nameall = paste(ename,hif,all,sep="")
	pym10 = paste(select,nameall,comma,ename,and,selectal,"\n",sep="")

		
	noname = paste(ename,hif,no,sep="")
	#str07b = " and not "
	#pym09a = paste(str07,noname,", ",ename,str07b,inameap,str07b,inamepo,"\n",sep="")
	pym11 = paste(select,noname,comma,ename,and,not,nameap,and,not,namepo,"\n",sep="")

	pym12 = paste(set,stickrad,comma,ename,"\n",sep="")

	if (model == strwrap(stick)) pym13 = paste(show, stick,comma, nameall,"\n",sep="")
	if (model == strwrap(sphere)) pym13 = paste(show, sphere,comma, nameall,"\n",sep="")
	if (model == strwrap(both)){
		pym13 = paste(show, stick,comma, nameall,"\n",sep="")
		pym13 = paste(pym13, show, sphere,comma, nameall,"\n",sep="")
		pym13 = paste(pym13, set,spheretrans,comma,nameall,"\n",sep="")

	}
	#pym09c = paste(str07c,ename,str07b,noname,"\n",sep="")
	#pym09c = paste(str07c,ename,str07b,inameal,"\n",sep="")
	#str07d = "show spheres, "
	#pym09d = paste(str07d,ename,str07b,noname,"\n",sep="")
	#pym09e = paste("set sphere_transparency=0.7, ", ename) 

	scriptsel = paste(pym08,pym09,pym10,pym11,pym12,pym13,"\n",sep="")
	

	pym10 = paste(run,dir,sub03,spectrumany,py,"\n",sep="")
	#str10 = "spectrumany b, white red, "

	pym11 = paste(spectrumany,spc,param03,comma,color_ap[1],spc,color_ap[2],comma,nameap,"\n",sep="")
	#str11 = "color yellow, "

	pym12 = paste(color,color_po,comma,namepo,"\n",sep="")

	scriptfin = paste(pym10,pym11,pym12)
	script = paste(scriptin,scriptsel,scriptfin,"\n",sep="")
	#str12 ="Scripts-Pymol/INB/"
	if (file)
		write(script,file=paste(dir,sub03,type,bar,scriptpre,ename,append,pml,sep=""))
	else
		#return(cat(script))
		return(script)
}

# Inicinaliza script Pymol conforme dirpath e pdbname
pymol_load_init = function(dirpath,pdbname,sufix){

	line = list()
	delete = "delete "
	all = "all"
	end = "\n"
	i = 1

	line[[i]] = paste(delete,all,sep="");i=i+1

	load = "load "
	line[[i]] = paste(load,dirpath,pdbname,sufix,sep="");i=i+1

	remove = "remove "
	hydrogens = "hydrogens " 
	solvent = "solvent " 
	line[[i]] = paste(remove, hydrogens,sep="");i=i+1
	line[[i]] = paste(remove, solvent,sep="");i=i+1
	
	bg = "bg_color "
	white = "white "
	line[[i]] = paste(bg,white,sep="")
	lineall = paste(line[1:i],end,sep="")
	return(lineall)

}


#Transforma [atom]name no padrão eR146.C para E/146/C
pymol_transname = function(v){

	tam = length(v)
	aux = c()
	auxbar = "/"
	auxsp = " "
	#auxplus = " + "

	for (i in 1:tam){
	
		#auxcnia1 = splitatom(v[i-1])
		auxcnia = splitatom(v[i])
		auxpaste = paste(auxsp,auxcnia[1],auxbar,auxcnia[3],auxbar,auxcnia[4],auxsp,sep="")
		#auxpaste2 = paste(auxcnia2[1],auxcnia2[3],auxcnia2[4])
		#auxpaste = paste(auxpaste1,auxpaste2)
		#auxcat = strcat(auxpaste,collapse="/")
		#print(auxpaste)
		#readline()
		aux = c(aux,auxpaste)
		
	}
	return(aux)

}

#select 1R0R.fe.grafo, (id 250+252+444+455+456+458+459+673+674+675+676+697+709+719+725+726+727+728+729+730+731+748+750+881+887+889+1064+1078+1079+1322+1323+1324+1325+1326+1327+1328+1524+1525+1526+1527+1532+1544+1551+1559+1560+2007+2008+2009+2010+2023+2024+2025+2026+2032+2033+2034+2039+2050+2052+2057+2058+2059+2060+2065+2066+2067+2074+2075+2076+2077+2078+2079+2080+2090+2206 and 1R0R);

# Monta Pymol select conforme selecname e vetor de nomes v
pymol_select = function(v,selectname,pdbname){

	select = "select "
	comma = " , "
	op = "("
	cl = ")"
	and = " and "
	end = "\n"
	#tam = length(v)
	auxv = strcat(v,collapse="+")
	#auxv = ""
	line = list()
	i = 1

	#for (i in 1:(tam-1)){
	#	aux = paste(aux,v[i],
#
	#}	
	#pymol01 = paste(pdbname,selectname,sep="")
	line[[i]] = paste(select,selectname,comma,op,auxv,and,pdbname,cl,sep="")
	lineall = paste(line[1:i],end,sep="")
	return(lineall)

}

# Configura spheres para selectname
pymol_spheres = function(selectname,pdbname,colors,spherescale){

	
	line = list()
	
	end = "\n"
	comma = " , "
	i = 1 

	#line[[i]] = paste("color ",colors[1],comma,pdbname,end,sep="");i=i+1
	#line[[i]] = paste("hide all",end,sep="");i=i+1
	line[[i]] = paste("color ",colors,comma,selectname,sep="");i=i+1
	line[[i]] = paste("show spheres, ",selectname,sep="");i=i+1
	line[[i]] = paste("set sphere_scale, ",spherescale,comma,selectname,sep="") 
	lineall = paste(line[1:i],end,sep="")

	return(lineall)

}

#distance dist, (e/35/CB), (e/35/CG)
#hide labels, dist
#set dash_color, gray, dist
#set dash_gap, 0, dist
#set dash_width, 1.0, dist

# Monta parte do script Pymol responsável pelas arestas conforme matriz m
pymol_set_graph_edges = function(m,transname,edgename,edgecolor,edgegap,edgewidth){

	tami = dim(m)[1]
	tamj = dim(m)[2]

	comma = " , "
	op = " ( "
	cl = " ) "
	end = "\n"
 
	k = 1
	line = list()

	for (i in 1:tami){
		for (j in i:tamj){
			#print(m[1,1])
			#readline()
			if (m[i,j]){
				line[[k]] = paste("distance ",edgename,comma,op,transname[i],cl,comma,op,transname[j],cl,sep="");k=k+1
			}
		}
	}
	line[[k]] = paste("hide labels, ", edgename);k=k+1
	line[[k]] = paste("set dash_color, ", edgecolor,comma,edgename);k=k+1
	line[[k]] = paste("set dash_gap, ", edgegap,comma,edgename);k=k+1
	line[[k]] = paste("set dash_width, ", edgewidth,comma,edgename)#;k=k+1
	#print(line[[k]])
	#readline()
	lineall = paste(line[1:k],end,sep="")
	return(lineall)
	
}

# Permite uma visão da alça inibitório no contexto geral do complexo enzima-inibidor
# DO TODO: precisa ser adaptado para gerar script Pymol correto quando tem mais de uma cadeia.
pymol_global = function(resname,pdbname,vertexsufix,chainnames,eivector,loopname="loop",basename="base",stick_trans=0.7,colors=c("yelloworange ","white "),sursufix="sur"){

	select = "select "
	comma = " , "
	op = "("
	cl = ")"
	and = " and "
	no = "no"
	not = " not "
	end = "\n"
	resi = " resi "
	chain = " chain "
	show = " show "
	sticks = " sticks "
	surface = " surface "
	set_trans = " set_bond stick_transparency "
	color = " color " 
	create = " create " 
	 
	#tam = length(v)
	#auxv = strcat(v,collapse="+")
	#auxv = ""
	line = list()
	i = 1

	auxname = find_inibitory_loop(resname,eivector)
	#print(auxname$loop)
	auxv = strcat(auxname$loop,collapse="+")
	auxloopname = paste(pdbname,loopname,sep="")
	#print(chainnames)
	#readline()
	line[[i]] = paste(select,auxloopname,comma,resi,auxv,and,chain,chainnames[2],and,pdbname,sep="");i=i+1

	if (auxname$noloop[1]!=0){
		auxv = strcat(auxname$noloop,collapse="+")
		noloopname = paste(pdbname,no,loopname,sep="")
		line[[i]] = paste(select,noloopname,comma,resi,auxv,and,chain,chainnames[2],and,pdbname,sep="");i=i+1
	}
	auxv = strcat(auxname$base,collapse="+")
	auxbasename = paste(pdbname,basename,sep="")
	line[[i]] = paste(select,auxbasename,comma,resi,auxv,and,chain,chainnames[1],and,pdbname,sep="");i=i+1

	line[[i]] = paste(show,sticks,comma,auxloopname,sep="");i=i+1
	if (auxname$noloop[1]!=0) line[[i]] = paste(show,sticks,comma,noloopname,sep="");i=i+1
	line[[i]] = paste(show,sticks,comma,auxbasename,sep="");i=i+1

	line[[i]] = paste(set_trans,comma,stick_trans,comma,auxloopname,sep="");i=i+1
	if (auxname$noloop[1]!=0) line[[i]] = paste(set_trans,comma,stick_trans,comma,noloopname,sep="");i=i+1
	line[[i]] = paste(set_trans,comma,stick_trans,comma,auxbasename,sep="");i=i+1

	tam = length(vertexsufix)
	auxodd = seq(1,tam,by=2)
	auxeven = seq(2,tam,by=2)
	auxsufix.i = vertexsufix[auxeven]
	auxsufix.e = vertexsufix[auxodd]
	auxvertex.i = paste(pdbname,auxsufix.i,sep="-")
	auxvertex.i = paste(auxvertex.i,comma)
	auxvertex.i = strcat(auxvertex.i)
	auxvertex.e = paste(pdbname,auxsufix.e,sep="-")
	auxvertex.e = paste(auxvertex.e,comma)
	auxvertex.e = strcat(auxvertex.e)
	#print(auxvertex.e)

	line[[i]] = paste(color,colors[1],comma,auxloopname,and,not,op,auxvertex.i,cl,sep="");i=i+1
	line[[i]] = paste(color,colors[2],comma,pdbname,and,chain,chainnames[1],and,not,op,auxvertex.e,cl,sep="");i=i+1

	auxsur = paste(pdbname,sursufix,sep="")
	#print(chainnames)
	#readline()
	line[[i]] = paste(create,auxsur,comma,pdbname,and,chain,chainnames[1],sep="");i=i+1
	line[[i]] = paste(show,surface,comma,auxsur,sep="")

	lineall = paste(line[1:i],end,sep="")
	return(lineall)


}



# Monta script com load, select por indicator vectors, configura spheres
pymol_view_by_vector_indicators = function(m,vlist,dirpath,resname,pdbname,filename,transname,edgename,colorlist,vpart,vertexsufix,chainnames,eivector,end="\n",spherescale=0.25,pymolglobal=TRUE){

	tam = length(vlist)

	hif="-"	

	line = list()
	j = 1

	line[[j]] = strcat(pymol_load_init(dirpath,pdbname,sufix=".pdb"));j=j+1
	line[[j]] = strcat("hide all",end);j=j+1

	for (i in 1:tam){

		selectname = paste(pdbname,hif,vertexsufix[i],sep="")
		line[[j]] = strcat(pymol_select(transname[vlist[[i]]],selectname,pdbname));j=j+1
		line[[j]] = strcat(pymol_spheres(selectname,pdbname,colorlist$color[colorlist$id[i]],spherescale=spherescale));j=j+1

	}
	
	if (vpart[1]==1) edv = make_edge_partition(vpart[2:length(vpart)])
	else edv = vpart[2:length(vpart)]
	#print(edv)
	tamedv = length(edv)
	if (tamedv < 2) print("WARNING: partition unitary in pymol edge constructor")
	for (i in 2:tamedv){
		if (i != tamedv){
			mref = (m >= edv[i]) & (m < edv[i-1])
		} else {
			mref = (m >= edv[i]) & (m <= edv[i-1])
		}

		line[[j]] = strcat(pymol_set_graph_edges(mref,transname,paste(pdbname,edgename,edv[i],hif,edv[i-1],sep=""),"gray",0,round(edv[i-1]/edv[1],2)));j=j+1
	}
	#j = j - 1
	# TO DO: compatibilizar pymol_global (encontro da alça inibitória) com qualquer tipo de cadeia-cadeia
	#print(resname)
	if (pymolglobal) line[[j]] = strcat(pymol_global(resname,pdbname,vertexsufix,chainnames,eivector))
	else j = j - 1
	lineall = paste(line[1:j],end,sep="")
	return(lineall)
	#write(lineall,file=filename)
}

#pymol_set_graph_edges = function(m,transname,edgename,edgecolor,edgegap,edgewidth)







