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
reg2 <- splineCox.reg2(t.event, event, Z, model = M)

# Display the results
print(reg2)

## ----fit-custom-model---------------------------------------------------------
# Define custom numeric vectors for baseline hazard shapes
custom_models <- list(c(0.1, 0.2, 0.3, 0.2, 0.2), c(0.2, 0.3, 0.3, 0.1, 0.1))

# Fit the model
reg2_custom <- splineCox.reg2(t.event, event, Z, model = custom_models)

# Display the results
print(reg2_custom)

## ----display-predefined-results-----------------------------------------------
# Print a summary of the results
print(reg2)

## ----display-custom-results---------------------------------------------------
# Print a summary of the results
print(reg2_custom)

