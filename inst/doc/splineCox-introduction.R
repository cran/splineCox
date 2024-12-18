## ----setup, echo = TRUE-------------------------------------------------------
library(splineCox)
library(joint.Cox)  # Required for example data

## ----example-data-------------------------------------------------------------
# Load the dataset
data(dataOvarian)

# Display the first few rows
head(dataOvarian)

## ----fit-model----------------------------------------------------------------
# Define variables
t.event <- dataOvarian$t.event
event <- dataOvarian$event
Z <- dataOvarian$CXCL12
M <- c("constant", "increase", "decrease")

# Fit the model
reg2 <- splineCox.reg2(t.event, event, Z, model = M)

# Display the results
print(reg2)

## ----display-results----------------------------------------------------------
# Print a summary of the results
print(reg2)

