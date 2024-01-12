library(tnet)
library(lubridate)
library(igraph)
library(lpSolve)
library(tidyverse)
library(fable)

source(file="Functions.R")
# http://www.sociopatterns.org/datasets/hypertext-2009-dynamic-contact-network/
# https://www.sciencedirect.com/science/article/pii/S0022519310006284

args=(commandArgs(TRUE))
print(args)
array_index <- as.numeric(args)
print(array_index)

index <- rep(20:30, each = 5)
time_step <- rep(1:5, 11)
para_df <- data.frame(index = index, time_step = time_step)
index <- para_df[array_index, 1]
time_step <- para_df[array_index, 2]
print(index)
print(time_step)

dat <- read.delim('Data_Input/ht09_contact_list.dat', header = FALSE)
head(dat)
colnames(dat) <- c("Seconds", "From", "To")
diff(range(dat$Seconds))/(60*60)

dat2 <- dat %>%
  mutate(Time = ceiling((Seconds)/(60*60)) ) %>%
  relocate(From, To, Time, Seconds) %>%
  select(-Seconds)
head(dat2)


# First rename the nodes
# Change vertex names
all_nodes <- unique(c(t(as.matrix(dat2[ ,1:2]))))
df_rename <- data.frame(actual_name = all_nodes, new_name = 1:length(all_nodes))

# renaming the From column to new_From
dff <- df_rename %>%
  mutate(From = actual_name) %>%
  right_join(dat2) %>%
  select(new_name, To, Time) %>%
  rename(From = new_name)

# renaming the To column to new To
dff2 <- df_rename %>%
  mutate(To = actual_name) %>%
  right_join(dff) %>%
  select( From, new_name, Time) %>%
  rename(To = new_name) %>%
  arrange(Time, From)
head(dff2)


# Rename time
unique(dff2$Time)
time_rename <- data.frame(oriTime = unique(dff2$Time), Time = 1:length(unique(dff2$Time)))
dff3 <- dff2 %>%
  rename(oriTime = Time) %>%
  full_join(time_rename)

length(unique(c(dff3$From, dff3$To)))
dim(dff3)



graphlist <- list()
for(win in 1:(index+11)){
  inds <- which(dff3$Time == win)
  datwin <- dff3[inds, 1:2]

  # get the vertices present up to that point
  indsvert <- which(dff3$Time <= win)
  datvert <- dff3[indsvert, 1:2]
  vertix <- unique(sort(c(datvert[ ,1], datvert[ ,2])))


  gr <- igraph::graph_from_data_frame(datwin, directed = FALSE, vertices = vertix )
  gradj <- igraph::as_adjacency_matrix(gr)

  # As we're still doing an unweighted graph put the weights back to 1
  gradj[gradj > 1] <- 1

  # make a graph from that
  gr2 <- igraph::graph_from_adjacency_matrix(gradj, mode = 'undirected')
  gr2 <- igraph::delete_vertex_attr(gr2, 'name')
  graphlist[[win]] <- gr2

}


forgr <- predict_graph(graphlist[1:index],
                       formulation = 2,
                       weights_opt = 4,
                       conf_level1 = NULL,
                       conf_level2 = 90,
                       weights_param = 0.001,
                       h = time_step)

# Evaluation metrics
eval1 <- eval_metrics(graphlist[[(index+time_step)]], forgr$graph_mean)

evaldf <- as.data.frame(eval1)
evaldf <- evaldf %>%
  as.data.frame() %>%
  mutate(pred_nodes = vcount(forgr$graph_mean),
         actual_nodes = vcount(graphlist[[(index+time_step)]]),
         pred_edges = ecount(forgr$graph_mean),
         actual_edges = ecount(graphlist[[(index+time_step)]]),
         exp = 'mean',
         h = time_step,
         graphs_in_training = index) %>%
  relocate(exp, graphs_in_training, h, actual_nodes, pred_nodes, actual_edges, pred_edges)


filename <- paste("Data_Output/HT09/Evaluation_Metrics_HT09_array_index_", array_index,  ".csv", sep="" )
write.csv(evaldf, filename, row.names = FALSE)






