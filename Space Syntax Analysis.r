
# Space Syntax Analysis plugin for R Project - R script version
#
# Copyright (C) Burak Beyhan
# 
# This plugin constructs adjusted graphs, calculate geodesic distance between the
# features and basic space syntax parameters. For a user guide and further explanations
# please refer to the following papers and book if you use this plugin in your studies;
#
# Beyhan, B. (2011) "Developing Space Syntax Tools for Free and Open Source Software for GIS", 
# in Proceedings of the 19th International Conference on Geoinformatics (Geoinformatics 2011), 
# Shanghai, China.
#
# Beyhan, B. (2012) "Developing Graph Theoretic Analysis Tools in FOSS4GIS: An Experiment in 
# OpenJUMP with a Specific Focus on Space Syntax", FOSS4G-CEE and Geoinformatics, Prague 2012.
#
# Beyhan, B. (2012) "A simple installation and user's guide for the plugins and scripts 
# developed to conduct space syntax analysis (SSA) in FOSS4GIS: OpenJUMP, gvSIG, OrbisGIS,
# Quantum GIS, OpenEv, Thuban, MapWindow GIS, SAGA, and R Project", 
# http://mekandizim.mersin.edu.tr/. 
#
# Beyhan, B. (2012) "Plugins and Scripts Developed to Conduct Space Syntax Analysis in FOSS4GIS: 
# OpenJUMP, gvSIG, OrbisGIS, Quantum GIS, OpenEv, Thuban, MapWindow GIS, SAGA, and R Project", 
# http://mekandizim.mersin.edu.tr/, forthcoming.
#
# This program is free software; you can redistribute it and/or modify it under the terms of 
# the GNU General Public License as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see http://www.gnu.org/licenses or write to the Free 
# Software Foundation,Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA


library(shapefiles)
library(tcltk)
library(sna)
library(maptools)

shpfilec = tkgetOpenFile(title = "Choose a shp file", filetypes = "{{GIS Files} {.shp}} {{All files} *}")
shpfiled = tclvalue(shpfilec)

if (shpfiled != "")
{

shpfiler = substr(shpfiled,1,nchar(shpfiled)-4)

dosya = read.shapefile(shpfiler)
shps = dosya$shp$shp
dbfs = dosya$dbf$dbf
objs = nrow(dbfs)

sna = array(0, c(objs,objs))
atr = array(0, c(objs))
var1 = array(0, c(objs))
var2 = array(0, c(objs))

colsim = colnames(dbfs)


# main script constructing adjusted graph, calculating geodesic and space syntax parameters

 # construction of adjusted graph

 for(i in 1:objs) {
   shp1 = shps[[i]]
   box1 = shp1$box
   xmin = box1[1]
   ymin = box1[2]
   xmax = box1[3]
   ymax = box1[4]
   row1 = shp1$num.points
   cntr = 0

   for(j in 1:objs) {
     shp2 = shps[[j]]
     box2 = shp2$box
     umin = box2[1]
     vmin = box2[2]
     umax = box2[3]
     vmax = box2[4]
     row2 = shp2$num.points

     if (xmax < umin || ymin > vmax) bonuc = 0
     else if (xmin > umax || ymin > vmax) bonuc = 0
     else if (xmax < umin || ymax < vmin) bonuc = 0
     else if (xmin > umax || ymax < vmin) bonuc = 0
     else if (i == j) bonuc = 0
     else bonuc = 1

     if (bonuc == 1) {

       for (p1 in 1:(row1-1)) {
         x1 = shp1$points[p1,1]
         y1 = shp1$points[p1,2]
         x2 = shp1$points[p1+1,1]
         y2 = shp1$points[p1+1,2]
         k = 0
         
        if(is.na(x1) || length(x1)==0 || is.nan(x1)) x1 = 0
        if(is.na(x2) || length(x2)==0 || is.nan(x2)) x2 = 0
        if(is.na(y1) || length(y1)==0 || is.nan(y1)) y1 = 0
        if(is.na(y2) || length(y2)==0 || is.nan(y2)) y2 = 0
         
         #if (identical(x1,x2)){ 
         if (x1 == x2){ 
           k = 1
         } else {
           b1 = (y2 - y1)/(x2 - x1)
           a1 = y1 - b1*x1 
         }

         for (p2 in 1:(row2-1)) {
           u1 = shp2$points[p2,1]
           v1 = shp2$points[p2,2]
           u2 = shp2$points[p2+1,1]
           v2 = shp2$points[p2+1,2]
           
           if(is.na(u1) || length(u1)==0 || is.nan(u1)) u1 = 0
           if(is.na(u2) || length(u2)==0 || is.nan(u2)) u2 = 0
           if(is.na(v1) || length(v1)==0 || is.nan(v1)) v1 = 0
           if(is.na(v2) || length(v2)==0 || is.nan(v2)) v2 = 0           
           
           l = 0
           sonuc = 0
           if (u2 == u1) l = 1
           else {
             b2 = (v2 - v1)/(u2 - u1)
             a2 = v1 - b2*u1
           }
           if ((k == 0) && (l == 0)) {
             #if (!identical(b1,b2)) {
             if (b1 != b2) {
               xi = 0 - (a1-a2)/(b1-b2) 
               yi = a1 + b1*xi
               if ((x1-xi)*(xi-x2)>=0 && (u1-xi)*(xi-u2)>=0 && (y1-yi)*(yi-y2)>=0 && (v1-yi)*(yi-v2)>=0) sonuc = 1
             }
           }
           if ((k == 1) && (l == 0)) {
             t1 = a2 + b2*x1
             if ((((u1 >= x1) && (u2 <= x1)) || ((u2 >= x1) && (u1 <= x1))) && (((y1 >= t1) && (y2 <= t1)) || ((y2 >= t1) && (y1 <= t1)))) sonuc = 1
           }
           if ((k == 0) && (l == 1)) { 
              t1 = a1 + b1*u1
              if ((((x1 >= u1) && (x2 <= u1)) || ((x2 >= u1) && (x1 <= u1))) && (((v1 >= t1) && (v2 <= t1)) || ((v2 >= t1) && (v1 <= t1)))) sonuc = 1
           }
           if (((x1 == u1) && (y1 == v1)) || ((x2 == u1) && (y2 == v1)) || ((x1 == u2) && (y1 == v2)) || ((x2 == u2) && (y2 == v2))) sonuc = 1
		
           if (sonuc == 1) break
         }
         if (sonuc == 1) {
           sna[i,j] = 1
           cntr = cntr + 1
           break
         }
       }
     }

   }
 atr[i] = cntr
 }


# calculation of geodesic by employing sna package (library)

gd = geodist(sna)
dst = gd$gdist


OnOK=function(){

 shpfilecs = tkgetSaveFile(title = "File to save results")
 shpfileds = tclvalue(shpfilecs)

# creation of new fields for the parameters

dbs = array(0, c(objs,8))
names = colnames(dbs)

names[1] = "Lineno"
names[2] = "Connectivity"
names[3] = "TotalDepth"
names[4] = "MeanDepth"
names[5] = "GlobalInteg"
names[6] = "LocalDepth"
names[7] = "LocalInteg"
names[8] = "Control"

colnames(dbs) = names


# calculation of space syntax parameters

 for(s1 in 1:objs) {
   dbs[s1,1] = s1 

   td = 0
   locd = as.integer(tclvalue(LRVal))
   deg3 = 0
   ld3 = 0
   cntrl = 0
   cnt = 0

   dvl = (2*(objs*(log((objs+2)/3,2)-1)+1))/((objs-1)*(objs-2))

   for(vas in 1:objs) {
     td = td + dst[s1,vas]
     if (dst[s1,vas] == 1) {
       cnt = cnt + 1
       cntrl = cntrl + 1/atr[vas]
     }
     if (dst[s1,vas] < locd) {
       ld3 = ld3 + dst[s1,vas]
       deg3 = deg3 + 1
     }
   }

   md = td/(objs-1)
   ra = 2*(md-1)/(objs-2)
   rr = ra/dvl
   gint = 1/rr

   lmd = ld3/(deg3-1)
   lra = 2*(lmd-1)/(deg3-2)
   ldvl = (2*(deg3*(log((deg3+2)/3,2)-1)+1))/((deg3-1)*(deg3-2))

   lrr = lra/ldvl	
   lint = 1/lrr

   dbs[s1,2] = cnt 
   dbs[s1,3] = td
   dbs[s1,4] = md
   dbs[s1,5] = gint
   dbs[s1,6] = ld3
   dbs[s1,7] = lint
   dbs[s1,8] = cntrl

   var1[s1] = cnt
   var2[s1] = gint

 }

dosya$dbf$dbf = dbs


# creation of new shp file storing the information for space syntax parameters

write.shapefile(dosya, shpfileds, T)


# production of the thematic map for global integration

dosya = readShapeLines(shpfileds)
dbfs = dosya@data
bb = as.matrix(dbfs[5])
bb = 1-bb/max(bb)
if(is.na(bb) || length(bb)==0 || is.nan(bb)) bb = 0
plot(dosya, col=rgb(bb, 0.2, 0.4), lwd=3)


tkdestroy(tt)


# calculation of intellegibility value

if(tclvalue(IVVal)=="1")
{
tkmessageBox(title = "Result", message = paste("Intelligibility value: ", cor(var1,var2)), icon = "info", type = "ok")
}

}


# script for GUI

tclServiceMode(FALSE)

tt <- tktoplevel()
kk <- tkframe(tt)
mm <- tkframe(tt)

tkwm.state(tt,"withdrawn") 

tktitle(tt) = "Space Syntax Analysis"

RF <- tklabel(kk, text = "Radius for local:")
LRVal <- tclVar("3")
LR <- tkentry(kk, textvariable=LRVal, width=3)

IV <- tkcheckbutton(tt, text="Calculate intellegibility value.")
IVVal=tclVar("0")
tkconfigure(IV,variable=IVVal)

BOS <- tklabel(mm, text = "    ")
HOS <- tklabel(tt, text = "    ")

tamam <- tkbutton(mm,text='OK', command=OnOK)
iptal <- tkbutton(mm,text='Cancel', command=function() tkdestroy(tt))

tkpack(RF, side='left')
tkpack(LR, side='right')
tkpack(kk, anchor="w")
tkpack(IV, anchor="w")
tkpack(HOS, anchor="w")
tkpack(tamam, side='left')
tkpack(BOS, side='left')
tkpack(iptal)
tkpack(mm)

tclServiceMode(TRUE)

tkwm.state(tt,"normal")
}
