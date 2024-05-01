
# NB: Download original/desired/latest file from http://www.healthdata.org/covid/data-downloads or historic predictions from https://www.healthdata.org/node/8787 and save in data folder. 

## Select all filenames from folder where information on hospitalizations is saved
file_list <- list.files("./data/reference_hospitalization_all_locs_timeline")
# Store as vector
file_list <- as.vector(file_list)
# Select only those that include ihme
file_list <- file_list[grepl("ihm", file_list) & grepl("csv", file_list)]

#file_name <- file_list[15:17] #use for debugging to find out which datasets give any errors

# Loop over the file list and process each file
for (file_name in file_list) {
  # Load the CSV file
  data <- read.csv(paste0("./data/reference_hospitalization_all_locs_timeline/", file_name), header=TRUE, stringsAsFactors=FALSE)
  
  # Subset the data to only include rows where location_name is "United States of America"
  data_subset <- data[data$location_name == "United States of America", ]
  
  # Rename 'date_reported' column to 'date' if it exists
  # This is needed as there are some variations of this column across the datasets
  if ("date_reported" %in% colnames(data_subset)) {
    colnames(data_subset)[which(colnames(data_subset) == "date_reported")] <- "date"
  }
  
  # Keep only the columns we want
  data_subset <- data_subset[, c("date", "admis_mean", "admis_upper", "admis_lower")]
  
  # Save specific data subset
  
  
  # Save the subsetted data as an RData file
  save(data_subset, file = paste0("./data/reference_hospitalization_all_locs_timeline/", gsub(".csv", ".RData", file_name)))
}




# Create an empty data frame to store the results
ihme_range <- data.frame(file_name = character(),
                      first_date = character(),
                      last_date = character(),
                      stringsAsFactors = FALSE)

for (file_name in file_list) {
  # Load the RData file
  load(paste0("./data/reference_hospitalization_all_locs_timeline/", gsub(".csv", ".RData", file_name)))
  
  # Assign the loaded data to the file name without the extension
  assign(gsub(".csv", "", file_name), data_subset)
  
  # Get the first and last dates
  first_date <- as.character(range(data_subset$date)[1])
  last_date <- as.character(range(data_subset$date)[2])
  
  # Add the results to the data frame
  ihme_range <- rbind(ihme_range, data.frame(file_name = file_name,
                                       first_date = first_date,
                                       last_date = last_date,
                                       stringsAsFactors = FALSE))
  
}
ihme_range

#Note: 2021_12_22 contained the wrong data, and 2021_04_30 contained no data
