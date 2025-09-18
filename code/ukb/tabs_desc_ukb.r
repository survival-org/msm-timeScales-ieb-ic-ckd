library(dplyr)
library(data.table)
library(lubridate)

dir_out <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb/"
file_pheno_out <- "pheno_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical.txt"
pheno <- fread(file.path(dir_out, file_pheno_out)) %>%
  mutate(date_at_exam = as.Date(date_at_exam))

# avg_gap_years <- function(df) {
#   df %>%
#     arrange(id, date_at_exam) %>%
#     group_by(id) %>%
#     mutate(gap_days = as.numeric(difftime(lead(date_at_exam), date_at_exam, units = "days"))) %>%
#     summarise(id_avg_years = mean(gap_days / 365.25, na.rm = TRUE), .groups = "drop") %>%
#     filter(!is.nan(id_avg_years)) %>%
#     summarise(mean_years = mean(id_avg_years),
#               n_ids = n())
# }

# res_overall <- pheno %>%
#   avg_gap_years() %>%
#   mutate(scenario = "overall")

# res_00 <- pheno %>%
#   filter(flag1 == 0, flag2 == 0) %>%
#   avg_gap_years() %>%
#   mutate(scenario = "flag1=0 & flag2=0")

# res_10 <- pheno %>%
#   filter(flag1 == 1, flag2 == 0) %>%
#   avg_gap_years() %>%
#   mutate(scenario = "flag1=1 & flag2=0")

# res_01plus <- pheno %>%
#   filter(flag2 == 1) %>%
#   avg_gap_years() %>%
#   mutate(scenario = "flag2=1")

# summary_gaps <- bind_rows(res_overall, res_00, res_10, res_01plus) %>%
#   select(scenario, mean_years, n_ids)

# summary_gaps

median_gap_years <- function(df) {
  df %>%
    dplyr::arrange(id, date_at_exam) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(gap_days = as.numeric(difftime(lead(date_at_exam), date_at_exam, units = "days"))) %>%
    dplyr::summarise(id_avg_years = mean(gap_days / 365.25, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(!is.nan(id_avg_years)) %>%
    dplyr::summarise(median_years = median(id_avg_years),
                     min_years = min(id_avg_years),
                     max_years = max(id_avg_years),
                     n_ids = n())
}

res_overall <- pheno %>%
  median_gap_years() %>%
  mutate(scenario = "overall")

res_00 <- pheno %>%
  filter(flag1 == 0, flag2 == 0) %>%
  median_gap_years() %>%
  mutate(scenario = "flag1=0 & flag2=0")

res_10 <- pheno %>%
  filter(flag1 == 1, flag2 == 0) %>%
  median_gap_years() %>%
  mutate(scenario = "flag1=1 & flag2=0")

res_01plus <- pheno %>%
  filter(flag2 == 1) %>%
  median_gap_years() %>%
  mutate(scenario = "flag2=1")

summary_gaps <- bind_rows(res_overall, res_00, res_10, res_01plus) %>%
  select(scenario, median_years, min_years, max_years, n_ids)

summary_gaps
