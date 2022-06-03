library(igraph)
#as a function for a igraph object with attributes "name" and "state".  State should be a numeric with 0 indicating uninfected nodes.  Infected nodes should have a state of 1.

ego.k <- function(graph, k, type=c("both","all","out","in"),node.type=NULL,par=T){
   if (!is.igraph(graph)) {stop("Not a igraph object")}
   if (is.null(V(graph)$name)) {stop("No 'name' attribute in igraph object")}
   if (!is.character(V(graph)$name)) {stop("Name attribute in iGraph object must be a character.  Use 'as.character()' to modify")}
   if (is.null(V(graph)$state)) {stop("No V(graph)$state attribute. infecteds must have non-zero value")}
   if(missing(type)) {type <- "both"}
   if(missing(k)) {stop("Number of steps (k) not specified")}
   if (is.directed(graph)==T & type=="all") {stop("Using all with a directed graph will ignore edge direction.  Use 'both' for type")}
  #node.type refers to if you want to assess ego.k for a specific type of node.  It should be given as the attribute name.  There should a vertex attribute of "node type"
   
  if(is.null(node.type)==T){#if no node type specified, run k-statistic for all nodes
   TB<- data.frame(name = as.character(V(graph)$name),state=as.numeric(V(graph)$state))
   TB.id <- subset(TB, state>0 & is.na(state)==F)
   TB.names <- as.character(TB.id$name)
  }
  
  if(is.null(node.type)==F){# run k-statistic for all nodes on specified node type
    TB<- data.frame(name = as.character(V(graph)$name),state=as.numeric(V(graph)$state),v.type=V(graph)$node.type)
    TB.id <- subset(TB, state>0 & is.na(state)==F & v.type== node.type)
    TB.names <- as.character(TB.id$name)
  }
  
   #output.k <- data.frame(Infected.node=as.character(rep("n",length(TB.names))),k.steps=rep(k,length(TB.names)),mode=rep(type,length(TB.names)),infecteds.k = rep(0,length(TB.names)),stringsAsFactors=F)
   
   library(snowfall);library(parallel)
   #for(i in 1:length(TB.names)){
  #procedure for k = 1
  if(k==1){
    type2 <- ifelse(type=="both","total",type)
    sub <- subgraph(graph,v=V(graph)$name[V(graph)$state>0 & is.na(V(graph)$state)==F])
    deg <- degree(simplify(sub),v=TB.names,mode=type2)
    output.k <- data.frame(Infected.node=names(deg),
                           k.steps=k,
                           mode=type,
                           infecteds.k=deg)
    }
  
  #procedure for k>1 (where we can't just do degree) 
  if(k>1){
    f <- function(x){  
     out.id <- V(graph)$name[V(graph)$name==TB.names[x]]
      
      if(type != "both"){
         ego.i <- graph.neighborhood(graph,k,out.id,mode=type)
         if(length(ego.i)==1){
            ego.i <- ego.i[[1]]
            
            k.stat <- length((V(ego.i)$state)[(V(ego.i)$state != "0") & is.na(V(ego.i)$state)==F])#Total infecteds within k steps.
            is.out <- V(graph)$state[match(TB.names,V(graph)$name)]
            k.stat <- ifelse(is.out[x]>0 & is.na(is.out[x])==F, k.stat-1,k.stat) #-1 if ego is outbreak
            
         }
         
      }
      if(type=="both"){
         ego.i.in <- graph.neighborhood(graph,k,out.id,mode="in")
         ego.i.out <- graph.neighborhood(graph,k,out.id,mode="out")
         if(length(ego.i.in)==1){
            e.u <- graph.union(ego.i.in[[1]],ego.i.out[[1]])
            #need to merge attibutes    
            V(e.u)$state <- as.numeric(V(graph)$state[match(V(e.u)$name,V(graph)$name)])
            
            k.stat <- length((V(e.u)$state)[(V(e.u)$state != "0") & is.na(V(e.u)$state)==F]) #Total infecteds within k steps.
            is.out <- V(graph)$state[match(TB.names,V(graph)$name)]
            k.stat <- ifelse(is.out[x]>0 & is.na(is.out[x])==F, k.stat-1,k.stat) #-1 if ego is outbreak
            
         }
      }
   
   
      row <- data.frame(Infected.node=TB.names[x],k.steps=k,mode=type,infecteds.k=k.stat)
      
      #output.k$Infected.node[i] <- TB.names[i]
      #print(i)
   }


  if(par==T){
    cores.n <- detectCores()-1
    sfInit(parallel=par,cpus=cores.n)
    sfLibrary(igraph)
    sfExport("TB.names","graph","k","type","node.type")
  }
  
  output.k <- sfLapply(1:length(TB.names),f)
  sfStop()
  output.k <- do.call("rbind",output.k)
  }
  
  output.k[is.na(output.k)==T]<-0
    
     
   return(output.k)
}




#if node type is specified, than randomization occurs within node types

k.test <- function(graph, k, iterations=1000, bin=c("degree","hood.size","full"),type=c("all","out","in","both"),node.type=NULL,par=TRUE){
   if (!is.igraph(graph)) {stop("Not a graph object")}
   if (is.null(V(graph)$name)) {stop("No V(graph)$name attribute")}
   if(missing(type)) {type <- "both"}
   if (is.null(V(graph)$state)) {stop("Not V(graph)$state attribute. infecteds must have non-zero value")}
   if(missing(k)) {stop("Number of steps (k) not specified")}
   if(missing(bin)) {stop("Type of randomization (bin=) not specified. 
                          Choose 'degree' for randomization within degree bins of width 25. 
                          Choose 'hood.size' for randomization among neighborhood size bins with width 200. 
                          Or choose 'full' for complete randomization. 
                          Full is the default")}
  #node.type refers to if you want to assess ego.k for a specific type of node.  It should be given as the attribute name.  There should a vertex attribute of "node type"
   library(reshape2)
   TB<- data.frame(name = as.character(V(graph)$name),outbreak=as.numeric(V(graph)$state))
   TB.id <- subset(TB, outbreak>0 & is.na(outbreak)==F)
   TB.names <- as.character(TB.id$name)
   
   o <- ego.k(graph,k,type,node.type=node.type)
   
   random.infecteds <- length(TB.names) #how big of a sample should be taken
   
   output <- data.frame(iter=rep(0,iterations),k.steps=rep(k,iterations),mode=rep(type,iterations),infecteds.median = rep(0,iterations),infecteds.mean = rep(0,iterations),stringsAsFactors=F)
   
   #for creating degree and hood.size bins to randomize within
   V(graph)$degree <- degree(graph)
   V(graph)$hood.size <- neighborhood.size(graph,2)
   t <- data.frame(V(graph)$name,V(graph)$state,V(graph)$degree,V(graph)$hood.size)
   cutpoints.d <-  unique(quantile(t[,3]+1))
   t$bin.d <- cut(t[,3],c(0,cutpoints.d))
   cutpoints <-   unique(quantile(t[,4]+1))
   t$bin.hs <- cut(t[,4],c(0,cutpoints))
   colnames(t)<- c("name","state","degree","hood.size","bin.d","bin.hs")
   t$bin.full <- 1
   if(is.null(node.type)==F){
     t$v.type <- V(graph)$node.type[match(t$name,V(graph)$name)]
   } else{ 
     t$v.type <- "any"
   }
   t.tb<- subset(t,state>0)  
   
   library(snowfall);library(parallel)
   func.iter <- function(x){
     i <- x
      graph.i <- graph
      #select which type of randomization will e performed (full, within degree bins(width=25), within neighborhood bins(width=200))
      if(bin=="degree"){
        #make a vecteor TB.names for use in graph.neighborhood.... length 58.  By drawing n from each bin
        bin.tab <- melt(table(t.tb$v.type,t.tb$bin.d))#get the number of observed infecteds per bin
        names.v <- c()
        for(j in 1:nrow(bin.tab)){
          v <- as.character(sample(t$name[t$bin.d==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
          names.v <- c(names.v,v)
        }
        V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
      }else if(bin=="hood.size"){
        #make a vecteor TB.names for use in graph.neighborhood.... length 58.  By drawing n from each bin
        bin.tab <- melt(table(t.tb$v.type,t.tb$bin.hs))#get the number of observed infecteds per bin
        names.v <- c()
        for(j in 1:nrow(bin.tab)){
          v <- as.character(sample(t$name[t$bin.hs==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
          names.v <- c(names.v,v)
        }
        V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
      }else{
        bin.tab <- melt(table(t.tb$v.type,t.tb$bin.full))#get the number of observed infecteds per bin
        names.v <- c()
        for(j in 1:nrow(bin.tab)){
          v <- as.character(sample(t$name[t$bin.full==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
          names.v <- c(names.v,v)
        }
        V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
      }
      
      r <- ego.k(graph.i,k,type,node.type=node.type,par=F)
      #output$iter[i] <- i
      #output$infecteds.mean[i] <- mean(r$infecteds.k)
      #output$infecteds.median[i] <- median(r$infecteds.k)
      row <- data.frame(iter=i,k.steps=k,mode=type,infecteds.median=median(r$infecteds.k),infecteds.mean=mean(r$infecteds.k))
      return(row)
      
   }
   
   cores.n <- detectCores()-1 
   sfInit(parallel=par,cpus=cores.n)
   sfLibrary(igraph);sfLibrary(reshape2)
   sfExport("ego.k","graph","bin","t.tb","t","k","type","node.type")
   
   output <- sfLapply(1:iterations,func.iter)
   sfStop()
   
   output <- do.call("rbind",output)
   
   
   
   m <- mean(o$infecteds.k)
   med.out <- median(o$infecteds.k)
   outbreak.med.p <- nrow(subset(output,infecteds.median>=med.out))/iterations
   outbreak.mean.p <- nrow(subset(output,infecteds.mean>=m))/iterations
   
   
   pvals <- data.frame(test=c("P-value (median infecteds within k steps)","P-value (mean infecteds within k steps)"),obs.vals=c(med.out,m),pval=c(outbreak.med.p,outbreak.mean.p))
   
   library(ggplot2)
   g <- ggplot(output,aes(x=infecteds.mean))+geom_density(fill="lightblue",adjust=1.5)+ geom_vline(xintercept=m[1],colour="red")+scale_x_continuous("Mean number of infecteds within k steps of a case",limits=c(0,max(c(m[1],output$infecteds.mean))+1))+theme_bw()+
      theme(legend.position="none",
            axis.title.y=element_text(size=16,vjust=1),
            axis.title.x=element_text(size=16,vjust=0),
            axis.text=element_text(size=12,colour="black"),
            panel.grid=element_blank())#median is somewhat worthless because median is usually 0 (meaning your pvalue can be 0 when all medians are 0(both observed and random))
   return(list(output,pvals,g))
   }






#significance of clustering by PATHS. If we randomly distribute at equal number of infecteds, what's the probability that the observed pathlenghts are shorter than if randomly distributed
#if you assign path.test to object, will get output of iterations as a dataframe.
#if network is directed, then you should use mode="in".  Otherwise, all arrows will be ignored (and paths will be calculated between nodes even if they aren't actually reachable in a directed network)
#Thus, when "in" is used, pathlengths are nearest paths that are incoming into a node (all arrows in path must flow towards the node)
#Both not an optiono here because else all shortest paths will be counted twice
#node type not implemented here accept for randomization purposes.
path.test <- function(graph, weight,  iterations=1000, bin=c("degree","hood.size","full"),type=c("all","out","in"),node.type=NULL,par=T){#types refers shortest paths function.  Weight should refer to an edge attribute
   if (!is.igraph(graph)) {stop("Not a graph object")}
   if (is.null(V(graph)$name)) {stop("No V(graph)$name attribute")}
   if (is.null(V(graph)$state)) {stop("Not V(graph)$state attribute. Infecteds must have non-zero value")}
   if(missing(type)) {type <- "in"}
   if (is.directed(graph)==T & type %in% c("all","out","both")) {stop("Type should be 'in'  for directed graphs")}
   if(missing(weight)) {weight <- NA}
   if(missing(bin)) {stop("Type of randomization (bin=) not specified. 
                          Choose 'degree' for randomization within degree bins of width 25. 
                          Choose 'hood.size' for randomization among neighborhood size bins with width 200. 
                          Or choose 'full' for complete randomization. 
                          Full is the default")}
   library(reshape2)
   TB<- data.frame(name = as.character(V(graph)$name),outbreak=as.numeric(V(graph)$state))
   TB.id <- subset(TB, outbreak>0 & is.na(outbreak)==F)
   TB.names <- as.character(TB.id$name)
   
   diam<-diameter(graph)
   
   #get nearest Infecteds in observed net (both weighted and unweighted)
   p.w <- shortest.paths(graph,v=TB.names,to=TB.names,weights=1/weight,mode=type)
   p.w <- ifelse(is.finite(p.w)==F,diam+1,p.w)#if unreachable, given the pathlength one longer of the longest geodesic in the graph
   nearest.obs.w <- c()
   for(i in 1:nrow(p.w)){
      p.i <- p.w[i,]
      i.min <- p.i[which(rank(p.i,ties.method="first")==2)]#find 2nd because first is diagonanol of matrix
      nearest.obs.w <- c(nearest.obs.w,i.min)
   }
   
   p <- shortest.paths(graph,v=TB.names,to=TB.names,mode=type)
   p <- ifelse(is.finite(p)==F,diam+1,p)
   nearest.obs <- c()
   for(i in 1:nrow(p)){
      p.i <- p[i,]
      i.min <- p.i[which(rank(p.i,ties.method="first")==2)]#find 2nd because first is diagonanol of matrix
      nearest.obs <- c(nearest.obs,i.min)
   }

   
   random.outbreaks <- length(TB.names) #how big of a sample should be taken
   
   output <- data.frame(iter=rep(0,iterations),mode=rep(type,iterations),pathlength.mean = rep(0,iterations),pathlength.w.mean = rep(0,iterations),stringsAsFactors=F)
   
   #for creating degree and hood.size bins to randomize within
   V(graph)$degree <- degree(graph)
   V(graph)$hood.size <- neighborhood.size(graph,2)
   t <- data.frame(as.character(V(graph)$name),V(graph)$state,V(graph)$degree,V(graph)$hood.size)
   cutpoints.d <-  unique(quantile(t[,3]+1))
   t$bin.d <- cut(t[,3],c(0,cutpoints.d))
   cutpoints <-   unique(quantile(t[,4]+1))
   t$bin.hs <- cut(t[,4],c(0,cutpoints))
   colnames(t)<- c("name","state","degree","hood.size","bin.d","bin.hs")
   t$bin.full <- 1
   if(is.null(node.type)==F){
     t$v.type <- V(graph)$node.type[match(t$name,V(graph)$name)]
   } else{ 
     t$v.type <- "any"
   }
   t.tb<- subset(t,state>0)   
   
   library(snowfall);library(parallel)
   func.iter <- function(x){
     i <- x
     graph.i <- graph
      #select which type of randomization will e performed (full, within degree bins(width=25), within neighborhood bins(width=200))
     if(bin=="degree"){
       #make a vecteor TB.names for use in graph.neighborhood.... length 58.  By drawing n from each bin
       bin.tab <- melt(table(t.tb$v.type,t.tb$bin.d))#get the number of observed infecteds per bin
       names.v <- c()
       for(j in 1:nrow(bin.tab)){
         v <- as.character(sample(t$name[t$bin.d==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
         names.v <- c(names.v,v)
       }
       V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
     }else if(bin=="hood.size"){
       #make a vecteor TB.names for use in graph.neighborhood.... length 58.  By drawing n from each bin
       bin.tab <- melt(table(t.tb$v.type,t.tb$bin.hs))#get the number of observed infecteds per bin
       names.v <- c()
       for(j in 1:nrow(bin.tab)){
         v <- as.character(sample(t$name[t$bin.hs==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
         names.v <- c(names.v,v)
       }
       V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
     }else{
       bin.tab <- melt(table(t.tb$v.type,t.tb$bin.full))#get the number of observed infecteds per bin
       names.v <- c()
       for(j in 1:nrow(bin.tab)){
         v <- as.character(sample(t$name[t$bin.full==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
         names.v <- c(names.v,v)
       }
       V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
     }
      
      p.w <- shortest.paths(graph,v=names.v,to=names.v,weights=1/weight,mode=type)
      p.w <- ifelse(is.finite(p.w)==F,diam+1,p.w)#if unreachable, given the pathlength one longer of the longest geodesic in the graph
      nearest.w <- c()
      for(g in 1:nrow(p.w)){
         p.i <- p.w[g,]
         i.min <- p.i[which(rank(p.i,ties.method="first")==2)]#find 2nd because first is diagonanol of matrix
         nearest.w <- c(nearest.w,i.min)
      }
      p <- shortest.paths(graph,v=names.v,to=names.v,mode=type)
      p <- ifelse(is.finite(p)==F,diam+1,p)
      nearest <- c()
      for(h in 1:nrow(p)){
         p.i <- p[h,]
         i.min <- p.i[which(rank(p.i,ties.method="first")==2)]#find 2nd because first is diagonanol of matrix
         nearest <- c(nearest,i.min)
      }

      #output$iter[i] <- i
      #output$pathlength.mean[i] <- mean(nearest)
      #output$pathlength.w.mean[i] <- mean(nearest.w)
      row <- data.frame(iter=i , mode=type, pathlength.mean=mean(nearest),pathlength.w.mean=mean(nearest.w) )
      return(row)
   }
   
   cores.n <- detectCores()-1 
   sfInit(parallel=T,cpus=cores.n)
   sfLibrary(igraph);sfLibrary(reshape2)
   sfExport("graph","bin","t.tb","t","type","node.type","weight","diam")
   
   output <- sfLapply(1:iterations,func.iter)
   sfStop()
   
   output <- do.call("rbind",output)
   
   nearest.obs.mean <- mean(nearest.obs)
   nearest.w.obs.mean <- mean(nearest.obs.w)
   
   nearest.w.p <- nrow(subset(output,pathlength.w.mean<=nearest.w.obs.mean))/iterations
   nearest.p <- nrow(subset(output,pathlength.mean<=nearest.obs.mean))/iterations
   
   pvals <- data.frame(test=c("P-value (Mean pathlength to nearest outbreak)","P-value (Mean weighted pathlength to nearest outbreak)"),obs.vals=c(nearest.obs.mean,nearest.w.obs.mean),pval=c(nearest.p,nearest.w.p))
   
   library(ggplot2)
   g <- ggplot(output,aes(x=pathlength.mean))+geom_density(fill="lightblue",adjust=1.5)+ geom_vline(xintercept=nearest.obs.mean,colour="red")+scale_x_continuous("Mean pathlength to nearest case")+
      theme(legend.position="none",
            axis.title.y=element_text(size=16,vjust=1),
            axis.title.x=element_text(size=16,vjust=0),
            axis.text=element_text(size=12,colour="black"),
            panel.grid=element_blank())#median is somewhat worthless because median is usually 0 (meaning your pvalue can be 0 when all medians are 0(both observed and random))
   g1 <- ggplot(output,aes(x=pathlength.w.mean))+geom_density(fill="lightblue",adjust=1.5)+ geom_vline(xintercept=nearest.w.obs.mean,colour="red")+scale_x_continuous("Mean weighted pathlength to nearest case")+
      theme(legend.position="none",
            axis.title.y=element_text(size=16,vjust=1),
            axis.title.x=element_text(size=16,vjust=0),
            axis.text=element_text(size=12,colour="black"),
            panel.grid=element_blank())#median is somewhat worthless because median is usually 0 (meaning your pvalue can be 0 when all medians are 0(both observed and random))
   return(list(output,pvals,g,g1))
   }













#significance of clustering. If we randomly distribute at equal number of outbreaks, what's the probability that the observed outbreaks have a greater value than if randomly distributed.  Caculates JOINT DISTRIBUTION (spatial within a threshold distance of outbreak)
#if you assign k.test to object, will get output of iterations as a dataframe
#Should do "both" if directed.  All will ignore directions of arrows
#if node.type is specifried, than randomizaiton occurs within node types.
ks.test <- function(graph, k, iterations=1000, bin=c("degree","hood.size","full"),type=c("all","out","in","both"),threshold=5,node.type=NULL,par=T,plot.ks=F){
   if (!is.igraph(graph)) {stop("Not a graph object")}
   if(missing(type)) {type <- "both"}
   if (is.directed(graph)==T & type=="all") {stop("Using all with a directed graph will ignore edge direction.  Use 'both' for type")}
   if (is.null(V(graph)$name)) {stop("No V(graph)$name attribute")}
   if (is.null(V(graph)$state)) {stop("No V(graph)$state attribute. Infecteds must have non-zero value")}
   if (is.null(V(graph)$lat)) {stop("No V(graph)$lat and long attributes")}
   if(missing(k)) {stop("Number of steps (k) not specified")}
   if(missing(threshold)) {threshold <- 10} #default threshold 
   if(missing(bin)) {stop("Type of randomization (bin=) not specified. 
                          Choose 'degree' for randomization within degree bins of width 25. 
                          Choose 'hood.size' for randomization among neighborhood size bins with width 200. 
                          Or choose 'full' for complete randomization. 
                          Full is the default")}
   library(reshape2)
   #get number outbreaks within k steps of observed outbreaks
   o <- ego.k(graph,k,type,node.type=node.type)
   #get numer outbreaks within treshold distance of observed outbreaks
   TB<- data.frame(name = as.character(V(graph)$name),outbreak=as.numeric(V(graph)$state))
   TB.id <- subset(TB, outbreak>0 & is.na(outbreak)==F)
   TB.names <- as.character(TB.id$name)
   
   TB.pairs <- merge(TB.names,TB.names)
   TB.pairs <- subset(TB.pairs,x != y)
   TB.pairs$long1 <-  V(graph)$long[match(TB.pairs$x,V(graph)$name)]
   TB.pairs$long2 <- V(graph)$long[match(TB.pairs$y,V(graph)$name)]
   TB.pairs$lat1 <- V(graph)$lat[match(TB.pairs$x,V(graph)$name)]
   TB.pairs$lat2 <- V(graph)$lat[match(TB.pairs$y,V(graph)$name)]
   TB.pairs$dist <- sqrt((TB.pairs$long1 - TB.pairs$long2)^2 + (TB.pairs$lat1 - TB.pairs$lat2)^2)/1000
   TB.pairs.threshD <- subset(TB.pairs,dist<threshold)
   
   if(is.null(node.type)==F){# run  for only nodes on specified node type.  Otherwise, run for all (default)
     TB.pairs.threshD$v.type <- V(graph)$node.type[match(TB.pairs.threshD$x,V(graph)$name)]
     TB.pairs.threshD <- subset(TB.pairs.threshD,  v.type== node.type)
   }
   
   infecteds.d <- tapply(TB.pairs.threshD$dist,as.character(TB.pairs.threshD$x),length)
   mean.out.d <- sum(infecteds.d,na.rm=T)/length(infecteds.d)
   
   random.infecteds <- length(TB.names) #how big of a sample should be taken
   
   output <- data.frame(iter=rep(0,iterations),k.steps=rep(k,iterations),mode=rep(type,iterations),infecteds.median = rep(0,iterations),infecteds.mean = rep(0,iterations),stringsAsFactors=F)
   
   #for creating degree and hood.size bins to randomize within
   V(graph)$degree <- degree(graph)
   V(graph)$hood.size <- neighborhood.size(graph,2)
   t <- data.frame(V(graph)$name,V(graph)$state,V(graph)$degree,V(graph)$hood.size,V(graph)$long)
   cutpoints.d <-  unique(quantile(t[,3]+1))
   t$bin.d <- cut(t[,3],c(-1,cutpoints.d))
   cutpoints <-   unique(quantile(t[,4]+1))
   t$bin.hs <- cut(t[,4],c(-1,cutpoints))
   colnames(t)<- c("name","state","degree","hood.size","long","bin.d","bin.hs")
   t <- subset(t,is.na(long)==F)
   t$bin.full <- 1
   if(is.null(node.type)==F){
     t$v.type <- V(graph)$node.type[match(t$name,V(graph)$name)]
   } else{ 
     t$v.type <- "any"
   }
    
   t.tb<- subset(t,state>0)  
   
   
   library(snowfall);library(parallel)
   func.iter <- function(x){
     i <- x
      graph.i <- graph
      #select which type of randomization will e performed (full, within degree bins(width=25), within neighborhood bins(width=200))
      if(bin=="degree"){
         #make a vecteor TB.names for use in graph.neighborhood.... length 58.  By drawing n from each bin
         bin.tab <- melt(table(t.tb$v.type,t.tb$bin.d))#get the number of observed infecteds per bin
         names.v <- c()
         for(j in 1:nrow(bin.tab)){
            v <- as.character(sample(t$name[t$bin.d==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
            names.v <- c(names.v,v)
         }
         V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
      }else if(bin=="hood.size"){
         #make a vecteor TB.names for use in graph.neighborhood.... length 58.  By drawing n from each bin
        bin.tab <- melt(table(t.tb$v.type,t.tb$bin.hs))#get the number of observed infecteds per bin
        names.v <- c()
        for(j in 1:nrow(bin.tab)){
          v <- as.character(sample(t$name[t$bin.hs==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
          names.v <- c(names.v,v)
         }
         V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
      }else{
        bin.tab <- melt(table(t.tb$v.type,t.tb$bin.full))#get the number of observed infecteds per bin
        names.v <- c()
        for(j in 1:nrow(bin.tab)){
          v <- as.character(sample(t$name[t$bin.full==bin.tab$Var2[j] & t$v.type==bin.tab$Var1[j]],bin.tab$value[j]))
          names.v <- c(names.v,v)
        }
         V(graph.i)$state <- ifelse(V(graph.i)$name %in% names.v,2011,0)
      }
      
      r <- ego.k(graph.i,k,type,node.type=node.type,par=F)
      #output$iter[i] <- i
      #output$infecteds.mean[i] <- mean(r$infecteds.k)
      #output$infecteds.median[i] <- median(r$infecteds.k)
      
      ####RUN STATS HERE FOR SPATIAL CULSTERING
      names.v <- as.character(names.v)
      r.pairs <- merge(names.v,names.v)
      r.pairs <- subset(r.pairs,x != y)
      r.pairs$long1 <- V(graph)$long[match(r.pairs$x,V(graph)$name)]
      r.pairs$long2 <- V(graph)$long[match(r.pairs$y,V(graph)$name)]
      r.pairs$lat1 <- V(graph)$lat[match(r.pairs$x,V(graph)$name)]
      r.pairs$lat2 <- V(graph)$lat[match(r.pairs$y,V(graph)$name)]
      r.pairs$dist <- sqrt((r.pairs$long1 - r.pairs$long2)^2 + (r.pairs$lat1 - r.pairs$lat2)^2)/1000
      r.pairs.threshD <- subset(r.pairs,dist<threshold)
      
      if(is.null(node.type)==F){# run  for only nodes on specified node type.  Otherwise, run for all (default)
        r.pairs.threshD$v.type <- V(graph.i)$node.type[match(r.pairs.threshD$x,V(graph.i)$name)]
        r.pairs.threshD <- subset(r.pairs.threshD,  v.type== node.type)
      }
      
      r.infecteds.d <- tapply(r.pairs.threshD$dist,as.character(r.pairs.threshD$x),length)
      #output$infecteds.mean.dist[i] <- sum(r.infecteds.d,na.rm=T)/length(names.v)
      
      row <- data.frame(iter=i,k.steps=k,mode=type,infecteds.median=median(r$infecteds.k),infecteds.mean=mean(r$infecteds.k),infecteds.mean.dist=sum(r.infecteds.d,na.rm=T)/length(names.v))
      return(row)
      
   }
   
   cores.n <- detectCores()-1 
   sfInit(parallel=T,cpus=cores.n)
   sfLibrary(igraph);sfLibrary(reshape2)
   sfExport("ego.k","graph","bin","t.tb","t","k","type","node.type","threshold")
   
   output <- sfLapply(1:iterations,func.iter)
   sfStop()
   
   output <- do.call("rbind",output)
   
   m <- mean(o$infecteds.k)
   med.out <- median(o$infecteds.k)
   outbreak.med.p <- nrow(subset(output,infecteds.median>=med.out))/iterations
   outbreak.mean.p <- nrow(subset(output,infecteds.mean>=m))/iterations
   
   output$infecteds.mean.dist <- ifelse(is.nan(output$infecteds.mean.dist)==T,0,output$infecteds.mean.dist)
   outbreak.dist.p <- nrow(subset(output,infecteds.mean.dist>=mean.out.d))/iterations
   pvals <- data.frame(test=c("P-value (median infecteds within k steps)","P-value (mean infecteds within k steps)","P-value (mean infecteds within threshold distance)"),obs=c(med.out,m,mean.out.d),pval=c(outbreak.med.p,outbreak.mean.p,outbreak.dist.p))
   
   
   library(ggplot2)
   g <- ggplot(output,aes(x=infecteds.mean))+geom_density(fill="lightblue",adjust=1.5)+ geom_vline(xintercept=m[1],colour="red")+scale_x_continuous("Mean number of infecteds within k steps of a case",limits=c(0,max(c(m[1],output$infecteds.mean))+1))+theme_bw()+
      theme(legend.position="none",
            axis.title.y=element_text(size=16,vjust=1),
            axis.title.x=element_text(size=16,vjust=0),
            axis.text=element_text(size=12,colour="black"),
            panel.grid=element_blank()) #median is somewhat worthless because median is usually 0 (meaning your pvalue can be 0 when all medians are 0(both observed and random))
   gs <- ggplot(output,aes(x=infecteds.mean.dist))+geom_density(fill="lightblue",adjust=1.5)+ geom_vline(xintercept=mean.out.d,colour="red")+scale_x_continuous("Mean number of infecteds within d distance of a case",limits=c(0,max(c(mean.out.d,output$infecteds.mean.dist))+1)) +theme_bw()+
      theme(legend.position="none",
            axis.title.y=element_text(size=16,vjust=1),
            axis.title.x=element_text(size=16,vjust=0),
            axis.text=element_text(size=12,colour="black"),
            panel.grid=element_blank()) #median is somewhat worthless because median is usually 0 (meaning your pvalue can be 0 when all medians are 0(both observed and random))
   
   #joint plot for ks test generic with ggplot
   if(plot.ks==T){   
     j <- data.frame(x=output$infecteds.mean,y=output$infecteds.mean.dist)
         #Nans are from cases where there are 0 farms within d
         j$y <- ifelse(is.nan(j$y)==T,0,j$y)
         j$x <- ifelse(is.nan(j$x)==T,0,j$x)
         g.ks <- NULL
         library(MASS)
         bw <- function(x){ #modified from bandwidth.nrd
            r <- quantile(x, c(0.25, 1))
            h <- (r[2] - r[1])/1.34
            4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
         }
         
         
            r <- c(range(j$x),range(j$y))
            r[c(1,3)] <- r[c(1,3)]-1
            r[c(2,4)] <- r[c(2,4)]+.1
            r <- ifelse(r<0,0,r)
            bw.x <- bw(j$x)
            bw.x <- ifelse(bw.x==0,.5,bw.x)
            bw.y <- bw(j$y)
            bw.y <- ifelse(bw.y==0,.5,bw.y)
            kk <- kde2d(j$x,j$y,h=c(bw.x,bw.y),n=100,lims=r)
            library(reshape2)
            dimnames(kk$z)<- list(kk$x,kk$y)
            dc <- melt(kk$z)
            g.ks <- ggplot(dc,aes(x=Var1,y=Var2))+
               geom_tile(aes(fill=value))+ 
               #xlim(min(c(j$x,pvals$obs[2])),max(c(j$x,pvals$obs[2])))+
               #ylim(min(c(j$y,pvals$obs[3])),max(c(j$y,pvals$obs[3])))+
               scale_fill_gradient(low="white",high="gray10")+
               theme_bw()+geom_contour(aes(z=value),bins=20,color="gray40",alpha=.4) +
               geom_contour(aes(z=value),breaks=quantile(dc$value,probs=c(.95)),col="black",size=1.5) + xlab("Mean number of infecteds within k steps of a case")+
               ylab("Mean number of infecteds within d distance")+ 
               theme(legend.position="none",
                     axis.title.y=element_text(size=16,vjust=1),
                     axis.title.x=element_text(size=16,vjust=0),
                     axis.text=element_text(size=12),
                     panel.grid=element_blank()) + 
               geom_point(data=pvals,aes(x=obs[2],y=obs[3]),pch="*",size=10)
         } else {g.ks=NULL}
   
   return(list(output,pvals,g,gs,g.ks))
   }








