library(tnet)
library(lubridate)
library(igraph)
library(lpSolve)
library(tidyverse)
library(fable)

source(file="Functions.R")

args=(commandArgs(TRUE))
print(args)
index <- as.numeric(args)


data(tnet)
dat <- tnet::OnlineSocialNetwork.n1899.lnet
head(dat)
dat$t <-  as_datetime(dat$t)
dat$date <- as_date(dat$t)
length(which(dat$i == dat$j))
inds <- which(dat$i == dat$j)
dat2 <- dat[-inds, ]

dat3 <- dat2[ ,c("i", "j", "date")]
len <-  length(unique(dat3$date))
unique_dates <- unique(dat3$date)
num_networks <- length(unique_dates)

graphlist <- list()
for(i in 1:(index+11)){ #
  nn <- unique_dates[i]
  inds <- which( dat3$date == nn )
  datwin <- dat3[inds, 1:2]

  # get the vertices present up to that point
  indsvert <- which(dat3$date <= nn)
  datvert <- dat3[indsvert, 1:2]
  vertix <- unique(sort(c(datvert[ ,1], datvert[ ,2])))

  # make a graph using those vertices
  gr <- graph_from_data_frame(datwin, directed = FALSE, vertices = vertix)
  gr <- delete_vertex_attr(gr, 'name')
  gradj <- as_adjacency_matrix(gr)

  # As we're still doing an unweighted graph put the weights back to 1
  gradj[gradj > 1] <- 1

  # make a graph from that
  gr2 <- graph_from_adjacency_matrix(gradj)
  graphlist[[i]] <- gr2
}


for(h in 1:10){
  time_step <- h
  forgr <- predict_graph(graphlist[1:index],
                         formulation = 2,
                         weights_opt = 4,
                         conf_level1 = 80,
                         conf_level2 = 90,
                         weights_param = 0.001,
                         h = time_step)

  # Evaluation metrics
  eval1 <- eval_metrics(graphlist[[(index+time_step)]], forgr$graph_mean)

  eval2 <- eval_metrics(graphlist[[(index+time_step)]], forgr$graph_lower)

  eval3 <- eval_metrics(graphlist[[(index+time_step)]], forgr$graph_upper)

  evaldf <- rbind(unlist(eval1), unlist(eval2), unlist(eval3))
  evaldf <- evaldf %>%
    as.data.frame() %>%
    mutate(pred_nodes = c(vcount(forgr$graph_mean), vcount(forgr$graph_lower), vcount(forgr$graph_upper)),
           actual_nodes = vcount(graphlist[[(index+time_step)]]),
           pred_edges = c(ecount(forgr$graph_mean), ecount(forgr$graph_lower), ecount(forgr$graph_upper)),
           actual_edges = ecount(graphlist[[(index+time_step)]]),
           exp = c('mean', 'lower', 'upper'),
           h = time_step) %>%
    relocate(exp, h, actual_nodes, pred_nodes, actual_edges, pred_edges)

  if(h == 1){
    evaldf_all <- evaldf
  }else{
    evaldf_all <- bind_rows(evaldf_all, evaldf)
  }

}


filename <- paste("Data_Output/UCI/Evaluation_Metrics_UCI_index_", index,  ".csv", sep="" )
write.csv(evaldf_all, filename, row.names = FALSE)








