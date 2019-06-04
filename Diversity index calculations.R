#===================================================================================
#      alpha diversity calculations ("shannon diversity index")
#===================================================================================

# shannon diversity calculation -> H = - sum pi log(pi)
# pi = ni/N
# ni = number of individuals of the ith species
# N  = total number of individuals 

# therefore: 
# shannon diversity -> H = -sum ((ni/N)*ln(ni/N))

library(vegan)

#====main R code 
# diversity(x, index = "shannon", MARGIN = 1, base = exp(1))
#
# x  =  vector or matrix or data to be used 
# index = diversity index to be used  ( "shannon" or "simpson" or "invsimpson" etc )
# Margin = 	Margin for which the index is computed = column (2) row (1)
# base = The logarithm base used in shannon.

#data 
# need at least two coloumns to determine alpha diversity for each sample 
# format 
# col 1 = species -> " unique clonotype"
# col 2 = species abundance -> abundance of each clontype in a sample 
# data must be numeric "remove first column"

# eventually will use 
setwd(here())
getwd()

#need to decide on location and folder 
#targetDir<-'../Results/InterrogateRunReport/'

#files<-list.files(targetDir,pattern = "^Run.*xlsx")


diversity(x, index = "shannon", MARGIN = 1, base = exp(1))


#======================
# test 
#======================

# test 1-> made a vector 
communityI <- c(10, 1, 1, 1, 1)
communityII <- c(5, 5, 5, 5, 5)
communities <- rbind(communityI,communityII)
#default diveristy index
# shannon index does it by row 
diversity(communities)


# upload excel file 
TEST <-read.csv(file.choose())

# test 2 -> calculation the alpha diversity by rows 
#-1 removes the first column 
diversity (TEST)

#test 3 -> calculation by colunm , Margin = 2 
diversity(TEST[,-1], index = "shannon", MARGIN = 2, base = exp(1))

#test 4 -> calculating a specifc column 
diversity(TEST$dog4, index = "shannon", MARGIN = 2, base = exp(1))

