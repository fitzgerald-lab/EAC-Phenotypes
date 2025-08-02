#=====================================================================# 
######  There are multiple .xxxxx.csv files in the list
######  The codes will keep the most recent one for each patient   
#=====================================================================#
# Extract patient ID and date from each path

# file_list <- file_list_BO
get_latest_files <- function(file_list) {
  file_df <- data.frame(
    file_path = file_list,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      file_name = basename(file_path),
      patient_id = str_extract(file_name, "^[A-Z]{2}_[0-9]{3}"),
      file_date = str_extract(file_name, "\\d{4}-\\d{2}-\\d{2}"),
      file_date = as.Date(file_date)
    )
  
  # Keep the latest file per patient
  latest_files <- file_df %>%
    group_by(patient_id) %>%
    slice_max(order_by = file_date, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Result: only the most recent .csv per patient
  cleaned_file_list <- latest_files$file_path
  
  return(cleaned_file_list)
}