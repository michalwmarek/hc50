library(webchem)
library(readxl)

my_data <- read_excel("hc50traningdata.xlsx", sheet = 5)
cas_numbers <- my_data$CAS
hc50_values <- my_data$HC50

smiles_numbers <- character(length(cas_numbers))

results <- cir_query(cas_numbers[1:length(cas_numbers)], 
                     from = "cas", 
                     to = "smiles")
smiles_numbers <- results$smiles

results_dataframe <- data.frame(CAS = cas_numbers, 
                                SMILES = smiles_numbers, 
                                HC50 = hc50_values)

write.csv(results_dataframe, 
          file = "smiles_hc50.csv", 
          row.names = TRUE, 
          quote = FALSE)