library(circlize)
library(viridis)

layout(matrix(1:2,1,2))

df1<-read.table("SNPden_raw.snpden",sep="\t",header=T)
df2<-read.table("SNPden_fit.snpden",sep="\t",header=T)

circle_plot <- function(df)
{
  df<-df[,c(1,2,4)]
  colnames(df)<-c("Chr","X","Y")
  df$Chr<-as.factor(df$Chr)
  levels(df$Chr)<-c(as.character(c(1:13)),rev(as.character(c(14:26))))
  df$X<-df$X/1000000
  options(scipen=999)
  
  col<-c(viridis(13,0.2,0.4,alpha=0.8),rev(viridis(13,0.8,1,alpha=0.8)))
  circos.initialize(factors=df$Chr,x=df$X)
  circos.trackPlotRegion(factors=df$Chr,y=df$Y,track.height = 0.05)
  for(i in 1:26){
    highlight.sector(sector.index = as.character(i),col=col[i])
    if(i<=13){
      circos.text(mean(get.cell.meta.data("xlim",sector.index = as.character(i))), 
                  mean(get.cell.meta.data("ylim")),
                  labels = paste0("A",i),
                  sector.index = as.character(i),cex=0.8)
    }
    else{
      circos.text(mean(get.cell.meta.data("xlim",sector.index = as.character(i))), 
                  mean(get.cell.meta.data("ylim")),
                  labels = paste0("D",i-13),
                  sector.index = as.character(i),cex=0.8)
    }
  }
  col<-c(viridis(13,0.2,0.4,alpha=0.5),rev(viridis(13,0.8,1,alpha=0.5)))
  circos.trackPlotRegion(factors=df$Chr,y=df$Y,track.height = 0.05)
  for(i in 1:26){
    highlight.sector(sector.index = as.character(i),track.index = 2,col=col[i])
    if(i<=13){
      circos.text(mean(get.cell.meta.data("xlim",sector.index = as.character(i))), 
                  mean(get.cell.meta.data("ylim")),
                  labels = i,
                  sector.index = as.character(i),cex=0.8)
    }
    else{
      circos.text(mean(get.cell.meta.data("xlim",sector.index = as.character(i))), 
                  mean(get.cell.meta.data("ylim")),
                  labels = i,
                  sector.index = as.character(i),cex=0.8)
    }
  }
  col<-c(viridis(13,0.2,0.4,alpha=0.8),viridis(13,0.8,1,alpha=0.8))
  circos.trackPlotRegion(factors=df$Chr,y=df$Y,track.height = 0.3, ylim = c(0,40))
  circos.trackLines(df$Chr,df$X,df$Y,col=col)
  
  draw.sector(get.cell.meta.data("cell.start.degree", sector.index = "1"), 
              get.cell.meta.data("cell.end.degree", sector.index = "13"), 
              rou1 = get.cell.meta.data("cell.bottom.radius", track.index = 3), col = col[6])
  draw.sector(get.cell.meta.data("cell.end.degree", sector.index = "13"), 
              get.cell.meta.data("cell.start.degree", sector.index = "1"), 
              rou1 = get.cell.meta.data("cell.bottom.radius", track.index = 3), col = col[19])
  circos.clear()
}
circle_plot(df1)
circle_plot(df2)
