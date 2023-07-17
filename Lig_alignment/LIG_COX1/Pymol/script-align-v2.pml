python
#import sys 
#graph1="test1e"
#graph2="test2e"
#output="out"

#protin=sys.argv[3]
#prot3="energy_54.pdb"
#prot3name=prot3.split(".pdb")[0]
#print prot3, prot3name
#staticStructurePath.split('/')[-1].split('.')[0]
graph1p = sys.argv[3]
graph1n = graph1p.split('/')[-1].split('.')[0]
#graph1n = graph1.split(".pdb")[0]
graph2p = sys.argv[4]
#graph2n = graph2.split(".pdb")[0]
graph2n = graph2p.split('/')[-1].split('.')[0]

#parw = int(float(sys.argv[5]))
parw = sys.argv[5]
parg = sys.argv[6]

print sys.argv
print graph1p
print graph2p
print graph1n
print graph2n
print parw
print parg
python end

#cmd.load("%s.pdb"%graph1n)
#cmd.load("%s.pdb"%graph2n)

#cmd.load("%s"%graph1p,"%s"%graph1n)
#cmd.load("%s"%graph2p,"%s"%graph2n)

#print "%s"%graph1n,"%s"%graph2n,"window=%s"%parw

cmd.cealign("%s"%graph1n,"%s"%graph2n,window=parw,quiet=-1,guide=parg)
#align=cmd.cealign("%s"%graph1n,"%s"%graph2n,window=parw,quiet=-1,guide=parg)
#print(align)
#aligns = str(align)
#rmsd = aligns.split("'RMSD':")[1].split(',')[0]
#leng = aligns.split("'alignment_length':")[1].split(',')[0]
#rot = aligns.split("'rotation_matrix':")[1]
#ini = '['
#end = ']'
#rot = rot[rot.find(ini)+len(ini):rot.rfind(end)]

#print(rmsd)
#print(leng)
#print(rot)

#alignp = aligns.split(":")
#print(alignp[0])
#print(alignp[1])

#python
#writefile=open('%s-%s.out'%(graph1n,graph2n),"a")
#writefile.write('%s %s\n'%(graph1n, graph2n))
#writefile.write(' '.join('%s' % x for x in align))
#writefile.write('%s'%(rmsd))
#writefile.write('\n')
#writefile.write('%s'%(leng))
#writefile.write('\n')
#writefile.write('%s'%(rot))
#writefile.write('\n')
#python end


