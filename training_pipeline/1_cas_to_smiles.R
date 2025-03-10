# Loading libraries
library(webchem)
library(readxl)

# Reading training data and setting up values by extracting specific columns and storing them in separate variables
my_data <- read_excel("hc50traningdata.xlsx", sheet = 5)
cas_numbers <- my_data$CAS
hc50_values <- my_data$HC50

# Setting up an empty character vector to store generated SMILES notations for each CAS number
smiles_numbers <- character(length(cas_numbers))

# Generating SMILES notations and storing them in a previously generated character vector
results <- cir_query(cas_numbers[1:length(cas_numbers)], 
                     from = "cas", 
                     to = "smiles")
smiles_numbers <- results$smiles

# Generating a dataframe to store CAS numbers, SMILES notations and HC50
results_dataframe <- data.frame(CAS = cas_numbers, 
                                SMILES = smiles_numbers, 
                                HC50 = hc50_values)

# Saving previously generated dataframe as a CSV (comma-separated values) file
write.csv(results_dataframe, 
          file = "smiles_hc50.csv", 
          row.names = TRUE, 
          quote = FALSE)