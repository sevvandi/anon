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


dat <- read.csv("Data_Input/email-Eu-core-temporal.txt", header = FALSE, sep = " ")
head(dat)
colnames(dat) <- c("From", "To", "Time")
sort(unique(c(dat$From, dat$To)))
range(dat$Time)
range(dat$Time)[2]/(60*60*24)


# First rename the nodes
# Change vertex names
all_nodes <- unique(c(t(as.matrix(dat[ ,1:2]))))
df_rename <- data.frame(actual_name = all_nodes, new_name = 1:length(all_nodes))

# renaming the From column to new_From
dff <- df_rename %>%
  mutate(From = actual_name) %>%
  right_join(dat) %>%
  select(new_name, To, Time) %>%
  rename(From = new_name)

# renaming the To column to new To
dff2 <- df_rename %>%
  mutate(To = actual_name) %>%
  right_join(dff) %>%
  select( From, new_name, Time) %>%
  rename(To = new_name) %>%
  arrange(Time)
head(dff2)

df3 <- dff2 %>%
  mutate(
    day = ceiling(Time/(60*60*24)) 
  )

df3[1, 4] <- 1 # update day of time 0 to be 1
num_networks <- length(unique(df3$day))
days <- sort(unique(df3$day))
days



graphlist <- list()
for(win in 1:(index+11)){
  inds <- which(df3$day %in% days[win])
  datwin <- df3[inds, 1:2]

  # get the vertices present up to that point
  indsvert <- which(df3$day <= days[win])
  datvert <- df3[indsvert, 1:2]
  vertix <- unique(sort(c(datvert[ ,1], datvert[ ,2])))


  gr <- graph_from_data_frame(datwin, directed = FALSE, vertices = vertix )
  gradj <- as_adjacency_matrix(gr)

  # As we're still doing an unweighted graph put the weights back to 1
  gradj[gradj > 1] <- 1

  # make a graph from that
  gr2 <- graph_from_adjacency_matrix(gradj, mode = 'undirected')
  gr2 <- delete_vertex_attr(gr2, 'name')
  graphlist[[win]] <- gr2

}


for(h in 1:10){
  time_step <- h
  forgr <- predict_graph(graphlist[1:index],
                         formulation = 2,
                         weights_opt = 2,
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


filename <- paste("Data_Output/EU_Email/Evaluation_Metrics_EU_Email_index_", index,  ".csv", sep="" )
write.csv(evaldf_all, filename, row.names = FALSE)














