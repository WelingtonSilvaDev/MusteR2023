delete all
load C:/temp/1PPF.pdb
remove hydrogens 
remove solvent 
bg_color white 


hide all


select 1PPF-e1 , ( E/42/SG + E/58/SG  and 1PPF)


color red , 1PPF-e1
show spheres, 1PPF-e1
set sphere_scale, 0.25 , 1PPF-e1


select 1PPF-i1 , ( I/19/CB  and 1PPF)


color pink , 1PPF-i1
show spheres, 1PPF-i1
set sphere_scale, 0.25 , 1PPF-i1


select 1PPF-e2 , ( E/57/CD2 + E/99/CD2  and 1PPF)


color yellow , 1PPF-e2
show spheres, 1PPF-e2
set sphere_scale, 0.25 , 1PPF-e2


select 1PPF-i2 , ( I/17/CG2  and 1PPF)


color paleyellow , 1PPF-i2
show spheres, 1PPF-i2
set sphere_scale, 0.25 , 1PPF-i2


select 1PPF-e3 , ( E/147/CG  and 1PPF)


color blue , 1PPF-e3
show spheres, 1PPF-e3
set sphere_scale, 0.25 , 1PPF-e3


select 1PPF-i3 , ( I/29/CE  and 1PPF)


color lightblue , 1PPF-i3
show spheres, 1PPF-i3
set sphere_scale, 0.25 , 1PPF-i3


hide labels,  1PPFedA4-5
set dash_color,  gray  ,  1PPFedA4-5
set dash_gap,  0  ,  1PPFedA4-5
set dash_width,  1  ,  1PPFedA4-5


hide labels,  1PPFedA3-4
set dash_color,  gray  ,  1PPFedA3-4
set dash_gap,  0  ,  1PPFedA3-4
set dash_width,  0.8  ,  1PPFedA3-4


hide labels,  1PPFedA2-3
set dash_color,  gray  ,  1PPFedA2-3
set dash_gap,  0  ,  1PPFedA2-3
set dash_width,  0.6  ,  1PPFedA2-3


select 1PPFloop ,  resi 17+18+19 and  chain I and 1PPF
select 1PPFnoloop ,  resi 29 and  chain I and 1PPF
select 1PPFbase ,  resi 42+57+58+99+147 and  chain E and 1PPF
 show  sticks  , 1PPFloop
 show  sticks  , 1PPFnoloop
 show  sticks  , 1PPFbase
 set_bond stick_transparency  , 0.7 , 1PPFloop
 set_bond stick_transparency  , 0.7 , 1PPFnoloop
 set_bond stick_transparency  , 0.7 , 1PPFbase
 color yelloworange  , 1PPFloop and  not (1PPF-i1  , 1PPF-i2  , 1PPF-i3  , )
 color white  , 1PPF and  chain E and  not (1PPF-e1  , 1PPF-e2  , 1PPF-e3  , )
 create 1PPFsur , 1PPF and  chain E
 show  surface  , 1PPFsur


