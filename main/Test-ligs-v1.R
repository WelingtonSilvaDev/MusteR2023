# 
# ############## TESTs #################
# 
# source("../bib/fun-base-v45.R")
# 
# if(1){
#   
#   i = 4
#   x = list()
#   x$a = as.matrix(dba2[[1]][[i]]$A$a)
#   diag(x$a) = 15000
#   x$n = lig.table$obs[ceiling(i/2)]
#   x$g = dba2[[1]][[i]]$g
#   x$xyz = 12*dba2[[1]][[i]]$pdb$xyz
#   x$e = E(x$g)$weight
#   
#   elen = dba2[[1]][[i]]$element_name
#   
#   x.tab = sep_atom_name_id(elen, sep="[.]", ncol = 3) %>% as_tibble()
#   x.tab = x.tab %>% extract(V2, c("resn","resi"),"(^.{0,2}\\D)(\\d*$)",
#                             remove=F)
#   x.tab = x.tab %>% mutate(resc=aa321(resn)) 
#   x.tab = x.tab %>% mutate(resC = case_when(
#     resc == "X" ~ resn,
#     TRUE ~ resc)
#   )
#   x.tab = x.tab %>% mutate(resCi = paste0(resC,resi))
#   x.tab = x.tab %>% relocate(any_of(c("resc","resC","resCi")), .after = resn)
#   x.tab = x.tab %>% mutate(resatom = paste0(resc,resi))
#   
#   col = rep("cornflowerblue", dim(x$a)[1])
#   f = grepl("B.", elen)
#   col[f] = "tomato"
#   f = grepl(".HOH", elen)
#   col[f] = "aquamarine"
#   
#   #crp = colorRampPalette(c("black", "white"))
#   crp = colorRamp(c("white", "black"))
#   values = x$e
#   rr <- range(values)
#   svals <- (values-rr[1])/diff(rr)
#   x$ecol = rgb(crp(svals)/255, alpha = 0.4)
#   
#   limxyz = get_lim_xyz(x$xyz,x$xyz)
#   #limxyz$xlim[1] = limxyz$xlim[1] + 5
#   limxyz$ylim[1] = limxyz$ylim[1] - 50
#   limxyz$ylim[2] = limxyz$ylim[2] + 50
#   
#   plot_g(x,limxyz,col = col, vlabel = x.tab$resatom, ecolor = x$ecol, 
#          title = x$n)
#   
#   #o = order(values)
#   #plot(values[o], col = x$ecol[o], pch = 19)
#   
#   
#   
# }