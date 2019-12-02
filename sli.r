

#=========== SLI PAPER DATA-ANALYSIS ==============#

# By Component 5 team members: Reza & rudy

# Access Reza's suites of new programs:
source("https://raw.githubusercontent.com/rnorouzian/m/master/m.r")

# Access the SLI paper Data:
D <- read.csv("https://raw.githubusercontent.com/rnorouzian/i/master/SLI.csv", h = T)


# total participants:
total <- nrow(D)         

# non-responses:
none <- nrow(na.omit(D)) 


#Percetange of data missing:
mis <- noquote(paste0(round(1 - none/total, 4)*1e2, "%")) # 12.68% of responses missing


# Use Reza's newly developed program to recover the missing responses:

D <- impute(D)   # impute the missing using median responses

pre <- D[,1:8]
pos <- D[,9:16]

# Use Reza's newly developed program to measure latent constructs (pre-post):

a <- (efa(pre, factors = 1, scores = "reg")$score)$Fa  # measure pre latent construct (educational leadership)
b <- (efa(pos, factors = 1, scores = "reg")$score)$Fa  # measure pos latent construct (educational leadership)

test <- t.test(a, b, paired = TRUE)  # run a paired t-test on the latent contstructs before & after SLI rather then items

t.val <- unname(test$statistic)

N <- total

# Use Reza's newly developed program to measure effect size (pre-post):

cohen.d  <- t2d(t.val, N) # measure effect size of change after SLI
# a d of "0.5896104"

# Use Reza's newly developed program to interpret effect size in percentages (pre-post):

dint.norm(cohen.d) # 'd' shows SLI has made 22.23% improvement.


# Show the conceptual model:


m2 <- " 
LBpr = ~Q1_a+Q6_a
BAVpr= ~Q2_a+Q3_a
OSpr = ~Q4_a+Q5_a
EFpr = ~Q7_a+Q8_a

LBps = ~Q1_b+Q6_b
BAVps= ~Q2_b+Q3_b
OSps = ~Q4_b+Q5_b
EFps = ~Q7_b+Q8_b

SLIpr =~ LBpr + BAVpr + OSpr+ EFpr

SLIps =~ LBps + BAVps + OSps+ EFps

"


library(semPlot)
library(lavaan)
library(ReporteRs)


# fit the conceptual model:

fit2 <- cfa(m2, data = D)

# Create a figure of the conceptual model
semPaths(fit2, layout = "spring", residuals = F)


# print the conceptual model for publication:

G2 <- function() semPaths(fit2, layout = "spring", residuals = F)


doc2 <- addPlot(docx(), fun = G2, vector.graphic = TRUE, width = 3.7, height = 4.5,  
                par.properties = parCenter(), editable = TRUE)


writeDoc(doc2, file = "CFA2.docx")
