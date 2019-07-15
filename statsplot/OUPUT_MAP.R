

for(i in 1:7){
  if(i==1){
    integrated.maplist_P1[[i]] <- integrated.maplist_P1[[i]][-1046,]
  } else if (i==4){
    integrated.maplist_P1[[i]] <- integrated.maplist_P1[[i]][-659,]
  }else{
    integrated.maplist_P1[[i]] <- integrated.maplist_P1[[i]]
  }

  }

integrated.maplist_P2[[2]]

for(i in 1:7){
  if(i==2){
    integrated.maplist_P2[[i]] <- integrated.maplist_P2[[i]][-1,]
  } else if (i==4){
    integrated.maplist_P2[[i]] <- integrated.maplist_P2[[i]][-1137,]
  }else{
    integrated.maplist_P2[[i]] <- integrated.maplist_P2[[i]]
  }

}



for(i in 1:7){
  if(i==1){
    integrated.maplist[[i]] <- integrated.maplist[[i]][-1909,]
  } else if (i==2){
    integrated.maplist[[i]] <- integrated.maplist[[i]][-1,]
  }else if (i==4){
    integrated.maplist[[i]] <- integrated.maplist[[i]][-c(1,2,1399),]
  }else {
    integrated.maplist[[i]] <- integrated.maplist[[i]]
  }

  }



alldt<-getwd()

for(i in 1:7){
  if (i==1){
    #dir.create(paste("integrated.maplist_P1"))
    #setwd(paste(alldt,"/integrated.maplist_P1",sep = ""))
    write.csv(integrated.maplist_P1[[i]],
                file=paste("integrated.maplist_P1_LG",i,".cvs",sep = ""),quote = F,row.names = F)
    #dir.create(paste("integrated.maplist_P2"))
    #setwd(paste(alldt,"/integrated.maplist_P2",sep = ""))
    write.csv(integrated.maplist_P2[[i]],
                file=paste("integrated.maplist_P2_LG",i,".csv",sep = ""),quote = F,row.names = F)
    #dir.create(paste("integrated.maplist_all"))
   # setwd(paste(alldt,"/integrated.maplist_all",sep = ""))
    write.csv(integrated.maplist_P1[[i]],
                file=paste("integrated.maplist_LG",i,".csv",sep = ""),quote = F,row.names = F)
  }else {
    #dir.create(paste("integrated.maplist_P1"))
    #setwd(paste(alldt,"/integrated.maplist_P1",sep = ""))
    write.csv(integrated.maplist_P1[[i]],
                file=paste("integrated.maplist_P1_LG",i,".csv",sep = ""),quote = F,row.names = F)
    #dir.create(paste("integrated.maplist_P2"))
    #setwd(paste(alldt,"/integrated.maplist_P2",sep = ""))
    write.csv(integrated.maplist_P2[[i]],
                file=paste("integrated.maplist_P2_LG",i,".csv",sep = ""),quote = F,row.names = F)
    #dir.create(paste("integrated.maplist_all"))
    #setwd(paste(alldt,"/integrated.maplist_all",sep = ""))
    write.csv(integrated.maplist_P1[[i]],
                file=paste("integrated.maplist_LG",i,".csv",sep = ""),quote = F,row.names = F)

  }

}


#dir.create(paste("integrated.maplist_P1"))
#setwd(paste(getwd(),"/integrated.maplist_P1",sep = ""))

all_P1<-do.call(rbind,integrated.maplist_P1)
write.csv(all_P1,
            file= "integrated.maplist_P1_all.csv",quote = F,row.names = F)

#dir.create(paste("integrated.maplist_P2"))
#setwd(paste(getwd(),"/integrated.maplist_P2",sep = ""))
all_P2<-do.call(rbind,integrated.maplist_P2)
write.csv(all_P1,
            file= "integrated.maplist_P2_all.csv",quote = F,row.names = F)
#dir.create(paste("integrated.maplist_all"))
#setwd(paste(getwd(),"/integrated.maplist",sep = ""))
all<-do.call(rbind,integrated.maplist_P2)
write.csv(all,
            file= "integrated.maplist_all.csv",quote = F,row.names = F)



####output genetype
map<-read.csv("polymap.map.tetra.csv",header=T)

ALL_dos<-read.csv("polymap.20k.tetra.csv",header=T)


###FATHER MAP
MarN<-NULL;MapD<-NULL;AvD<-NULL;MaxG<-NULL;#Gap5cM<-NULL
Position_P1<-list();Genetype_P1<-list()
for(i in 1:7){
  MarN[i]<-dim(integrated.maplist_P1[[i]])[1]
  MapD[i]<-max(integrated.maplist_P1[[i]][,2])
  AvD[i]<-max(integrated.maplist[[i]][,2])/dim(integrated.maplist[[i]])[1]
  MaxG[i]<-max(integrated.maplist_P1[[i]][,2])-min(integrated.maplist_P1[[i]][,2])
  #Gap5cM
  Position_P1[[i]]<-map[which(ALL_dos$SNP %in% integrated.maplist_P1[[i]][,1]),]
  Genetype_P1[[i]]<-ALL_dos[which(ALL_dos$SNP %in% integrated.maplist_P1[[i]][,1]),]
  }

phymap<-do.call(rbind,Position_P1)
write.table(phymap,file="phymap_P1.txt",
           quote = F,row.names = F,col.names = T)

Genetype<-do.call(rbind,Genetype_P1)

write.table(Genetype,file="Genetype_P1.txt",
            quote = F,row.names = F,col.names = T)


all_stats<-as.matrix(cbind(MarN,MapD,AvD,MaxG))

colnames(all_stats)<-c("Mar Num","Map Distance","Aver Distance","Max Gap")

for(i in 1:7){
  rownames(all_stats)[i]<-paste("LG",i,sep = "")
}
write.csv(all_stats,file = "stats_P1.csv",quote=F,row.names=T)


####MOTHER MAP
MarN<-NULL;MapD<-NULL;AvD<-NULL;MaxG<-NULL;#Gap5cM<-NULL
Position_P2<-list();Genetype_P2<-list()
for(i in 1:7){
  MarN[i]<-dim(integrated.maplist_P2[[i]])[1]
  MapD[i]<-max(integrated.maplist_P2[[i]][,2])
  AvD[i]<-max(integrated.maplist[[i]][,2])/dim(integrated.maplist[[i]])[1]
  MaxG[i]<-max(integrated.maplist_P2[[i]][,2])-min(integrated.maplist_P2[[i]][,2])
  #Gap5cM
  Position_P2[[i]]<-map[which(ALL_dos$SNP %in% integrated.maplist_P2[[i]][,1]),]
  Genetype_P2[[i]]<-ALL_dos[which(ALL_dos$SNP %in% integrated.maplist_P2[[i]][,1]),]
}

phymap<-do.call(rbind,Position_P2)
write.table(phymap,file="phymap_P2.txt",
            quote = F,row.names = F,col.names = T)

Genetype<-do.call(rbind,Genetype_P2)

write.table(Genetype,file="Genetype_P2.txt",
            quote = F,row.names = F,col.names = T)


all_stats<-as.matrix(cbind(MarN,MapD,AvD,MaxG))

colnames(all_stats)<-c("Mar Num","Map Distance","Aver Distance","Max Gap")

for(i in 1:7){
  rownames(all_stats)[i]<-paste("LG",i,sep = "")
}
write.csv(all_stats,file = "stats_P2.csv",quote=F,row.names=T)


###INTEGRATE MAP
MarN<-NULL;MapD<-NULL;AvD<-NULL;MaxG<-NULL;#Gap5cM<-NULL
Position_PP<-list();Genetype_PP<-list()
for(i in 1:7){
  MarN[i]<-dim(integrated.maplist[[i]])[1]
  MapD[i]<-max(integrated.maplist[[i]][,2])
  AvD[i]<-max(integrated.maplist[[i]][,2])/dim(integrated.maplist[[i]])[1]
  MaxG[i]<-max(integrated.maplist[[i]][,2])-min(integrated.maplist[[i]][,2])
  #Gap5cM
  Position_PP[[i]]<-map[which(ALL_dos$SNP %in% integrated.maplist[[i]][,1]),]
  Genetype_PP[[i]]<-ALL_dos[which(ALL_dos$SNP %in% integrated.maplist[[i]][,1]),]
}

phymap<-do.call(rbind,Position_PP)
write.table(phymap,file="phymap_PP.txt",
            quote = F,row.names = F,col.names = T)

Genetype<-do.call(rbind,Genetype_PP)

write.table(Genetype,file="Genetype_PP.txt",
            quote = F,row.names = F,col.names = T)


all_stats<-as.matrix(cbind(MarN,MapD,AvD,MaxG))

colnames(all_stats)<-c("Mar Num","Map Distance","Aver Distance","Max Gap")

for(i in 1:7){
  rownames(all_stats)[i]<-paste("LG",i,sep = "")
}
write.csv(all_stats,file = "stats_PP.csv",quote=F,row.names=T)
