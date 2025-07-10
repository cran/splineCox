## ----setup, echo = TRUE-------------------------------------------------------
library(splineCox)
library(joint.Cox)  # Required for example data

## ----example-data-------------------------------------------------------------
# Load the dataset
data(dataOvarian)

# Display the first few rows
head(dataOvarian)

## ----fit-predefined-model-----------------------------------------------------
# Define variables
t.event <- dataOvarian$t.event
event <- dataOvarian$event
Z <- dataOvarian$CXCL12
M <- c("constant", "increase", "decrease")

# Fit the model
reg2 <- splineCox.reg2(t.event, event, Z, model = M, plot = TRUE)

# Display the results
print(reg2)

## ----fit-custom-model---------------------------------------------------------
# Define custom numeric vectors for baseline hazard shapes
custom_models <- list(c(0.1, 0.2, 0.3, 0.2, 0.2), c(0.2, 0.3, 0.3, 0.1, 0.1))

# Fit the model
reg2_custom <- splineCox.reg2(t.event, event, Z, model = custom_models, plot = TRUE)

# Display the results
print(reg2_custom)

## ----display-predefined-results-----------------------------------------------
# Print a summary of the results
print(reg2)

## ----display-custom-results---------------------------------------------------
# Print a summary of the results
print(reg2_custom)

## ----copula-density-plot, message=FALSE, warning=FALSE------------------------
library(ggplot2)

N <- 50
u <- v <- seq(from = 0, to = 1, length.out = N)
U <- rep(u, N)
V <- rep(v, each = N)

# Positive Exchangeable
c.data <- data.frame(
  U = U, V = V,
  C = spline.copula(U, V, R = "PE1", density = TRUE, mat = FALSE)
)

ggplot(aes(x=U, y=V), data=c.data) +
  geom_contour(aes(x=U,y=V,z=C,colour=after_stat(level)),
               data=c.data,bins=25)+xlab("u")+ylab("v")

# Negative Exchangeable
c.data <- data.frame(
  U = U, V = V,
  C = spline.copula(U, V, R = "NE3", density = TRUE, mat = FALSE)
)

ggplot(aes(x=U, y=V), data=c.data) +
  geom_contour(aes(x=U,y=V,z=C,colour=after_stat(level)),
               data=c.data,bins=25)+xlab("u")+ylab("v")

## ----copula-dependence-measures-----------------------------------------------
# Compute Kendall's tau and Spearman's rho for preset R matrix "PE1"
out <- spline.copula(U, V, R = "PE1", Kendall = TRUE, Spearman = TRUE)

# Display the results
out$Kendall_tau
out$Spearman_rho

