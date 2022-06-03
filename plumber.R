# plumber.R
source("Cluster.test.source.parallel.R")

#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  plumber::forward()
}

#* Echo back the input
#* @param msg The message to echo
#* @get /echo
function(msg="") {
  list(msg = paste0("The message is: '", msg, "'"))
}

#* Plot a histograma
#* @serializer png
#* @get /plot
function() {
  rand <- rnorm(100)
  hist(rand)
}

#* Return the sum of two numbers
#* @param a The first number to add
#* @param b The second number to add
#* @post /sum
function(a, b) {
  as.numeric(a) + as.numeric(b)
}

#* Return the k-test result on sample data
#* @get /sample
function() {
  e.A <- read.csv("data/edges.A_full.csv", stringsAsFactors = F)
  attr1 <- read.csv("data/attr.csv", stringsAsFactors = F)
  netA <- graph.data.frame(e.A, directed = F)
  V(netA)$state <- as.integer(attr1$state[match(V(netA)$name,attr1$name)])
  test2 <- k.test(netA,k=1,type="both",bin="full",iterations=20)
  test2[c(1,2)]
}

#* Return a k-test statistics within a window
#* @param window The week to snapshot infected cases
#* @get /k-test
function(type, window) {
  if(type == "fs") {
    edges <- read.csv("data/flight_salient_edges.csv", stringsAsFactors = F) 
    cases <- read.csv("data/flight_cases.csv", stringsAsFactors = F)   
  } else if(type == "ps") {
    edges <- read.csv("data/phone_salient_edges.csv", stringsAsFactors = F) 
    cases <- read.csv("data/phone_cases.csv", stringsAsFactors = F)   
  }else {
    edges <- read.csv("data/flight_salient_edges.csv", stringsAsFactors = F) 
    cases <- read.csv("data/flight_cases.csv", stringsAsFactors = F) 
  }

  cases$state <- ifelse(cases$arrival < window, 1, 0)
  cases$name <- cases$id
  net <- graph.data.frame(edges, directed = F)
  V(net)$state <- as.integer(cases$state[match(V(net)$name,cases$name)])
  test <- k.test(net,k=1,type="both",bin="full",iterations=20)
  test[c(1,2)]
}
