library(openxlsx)
library(igraph)
library(ggplot2)

#change the path according to where these files are.
#genes names and node ids
nodes_genes<-read.csv('../data/node_gene.csv',header = T,stringsAsFactors = F)

#GO term associations with genes - ground truth
tf2<-read.xlsx(sheet=1,'../data/gt_all_solution.xlsx',colNames = F)

gene<-tf2$X1 #go term genes only, some gene do not have go term associations.
tf2$X1<-NULL

n_c<-ncol(tf2)

nm<-data.frame()
ra<-data.frame()

ns = 100 # number of samples

nm[c(1:ns),1]<-1:ns
ra[c(1:ns),1]<-1:ns

all_partitions <- read.xlsx(sheet = 1,'../data/all_partitions.xlsx')

for (cluster_column in colnames(all_partitions)[c(3:14)]){
  #each partition
  
  tf <- subset(all_partitions,select = c('node','gene',cluster_column))
  tf<-setNames(tf,c('node','gene','V2'))
  
  
  ra1<-1:ns
  nm1<-1:ns
  
  for(j in c(1:ns))
  {
    tf2$XR<-0
    for (i in c(1:nrow(tf2))){
      
      tf2[i,n_c+1]<-as.numeric(sample(tf2[i,which(!is.na(tf2[i,1:n_c])),],1))
      
    }
    
    tf1<-cbind(gene,tf2)
    tf1<-subset(tf1,select=c('gene','XR'))
    
    tf1<-merge(tf,tf1,by='gene',all.x=T)
    tf1<-tf1[which(!is.na(tf1$XR)),]
    tf1<-tf1[which(tf1$V2>0),]
    
    smallcluster = TRUE #include clusters with sizes < 3 ?
    
    if (smallcluster){
      
      nm1[j]=compare(comm1 = tf1$V2,comm2 = tf1$XR,method='adjusted.rand')
      
    } else {
      
      tf1_2<-tf1
      V2freq<-data.frame(table(tf1_2$V2))
      V2freq<-setNames(V2freq,c("V2","V2freq"))
      XRfreq<-data.frame(table(tf1_2$XR))
      XRfreq<-setNames(XRfreq,c("XR","XRfreq"))
      tf1_2<-merge(tf1_2,V2freq,by="V2",all.x=T)
      tf1_2<-merge(tf1_2,XRfreq,by="XR",all.x=T)
      
      tf1_2<-tf1_2[which(tf1_2$V2freq>2 & tf1_2$XRfreq>2),] #cluster sizes larger than 2
      nm1[j]=compare(comm1 = tf1_2$V2,comm2 = tf1_2$XR,method='adjusted.rand')
    }
    
    
  }
  
  for (k in c(1:length(nm1))){
    ra1[k] <- sum(nm1[1:k])/k
  }
  
  
  #plot(ra1,type='l',xlab="samples",ylab="running avg. ARI")
  
  nm[cluster_column]<-nm1
  ra[cluster_column]<-ra1
  
}


#Plot distribution with CLR network
ggplot(nm)+geom_histogram(aes(cluster_Q_f),alpha=0.6,bins=500)+geom_histogram(aes(cluster_Qx_f),alpha=0.6,fill='blue',bins=500)+
  geom_histogram(aes(cluster_MCL_f),alpha=0.6,fill='red',bins=500)+
  geom_histogram(aes(cluster_Qg_f),alpha=0.6,fill='green',bins=500)+
  xlab('Sample ARI')+theme_bw()
ggsave("out/ARI_distribution.pdf", width = 4, height = 4)

#write.csv(nm,'../out/ensemble_ARI.csv',row.names = F,quote = F)

eARI<-data.frame(matrix(NA, nrow = 12, ncol = 3))

eARI$X1<-as.vector(colnames(nm[2:13]))

for(i in c(2:ncol(nm))){
  
  eARI[i-1,2] <- mean(nm[,i])
  eARI[i-1,3] <- sd(nm[,i])
  
}

eARI<-eARI[order(eARI$X1),]

xlabels = eARI$X1
ggplot(eARI, aes(x=X1, y=X2,fill=X1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width = 0.5) +
  geom_errorbar(aes(ymin=X2-X3, ymax=X2+X3), width=.1,
                position=position_dodge(2.))+coord_flip()+theme_classic()+ylab("Average ARI score")+xlab("Clustering method")+
  theme(legend.position = "none")+scale_x_discrete(breaks=eARI$X1,labels=xlabels)
ggsave("../out/avg_ARI.pdf", width = 4, height = 4)
#Specify the directory and file name.

#ggplot(ra)+geom_line(aes(V1,cluster_Q_f))+
#  geom_line(aes(V1,cluster_MCL_f),color='red')+
#  geom_line(aes(V1,cluster_Qx_f),color='green')+
#  geom_line(aes(V1,cluster_Qg_f),color='purple')+xlab('#Samples')+ylab('Running avg. ARI')+theme_bw()

# One partition only
#ggplot(nm)+
#  geom_density(aes(cluster_Qx_f),alpha=0.6,fill='green',bins=100)+xlab('Sample ARI')+theme_bw()

#ggplot(ra)+
#  geom_line(aes(V1,cluster_Qx_f),color='green')+xlab('#Samples')+ylab('Running avg. ARI')+theme_bw()
