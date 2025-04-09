source("wis37138/kidney_survival/code/helpers.r")

library(ggplot2)

dir_data <- "/wis37138/kidney_survival/data/ukb/"
dir_fig <- "/wis37138/kidney_survival/results/figures/"
file_pheno <- "pheno_ukb.txt"
file_pheno_noaki6months <- "pheno_ukb_noAKI6months.txt"
file_events <- "events_ukb.txt"
file_events_noaki6months <- "events_ukb_noAKI6months.txt"
file_events_withBackTransitions <- "events_withBackTransitions_ukb.txt"
file_events_withBackTransitions_noaki6months <- "events_withBackTransitions_ukb_noAKI6months.txt"

pheno <- read_file(file.path(dir_data, file_pheno))
events <- read_file(file.path(dir_data, file_events))
events_withBackTransitions <- read_file(file.path(dir_data, file_events_withBackTransitions))

nrow(pheno)
length(unique(pheno$id))

## distribution of number of eGFR assessments per individual
egfr_counts <- pheno %>%
    group_by(id) %>%
    summarise(count = n())

summary(egfr_counts$count)

egfr_counts2 <- egfr_counts %>%
    mutate(count = ifelse(count > 50, 50, count))

visit_counts <- ggplot(egfr_counts2, aes(x = count)) +
    geom_histogram(binwidth = 1, fill = "blue", color = "black") +
    labs(title = "Distribution of Number of eGFR Assessments per Individual",
             x = "Number of eGFR Assessments",
             y = "Frequency") +
    theme_minimal()

ggsave(file.path(dir_fig, "visit_counts.png"), visit_counts, width = 6, height = 4)

## distribution of time between eGFR assessments
pheno <- pheno %>%
    arrange(id, event_dt)

pheno <- pheno %>%
    group_by(id) %>%
    mutate(time_diff = as.numeric(difftime(event_dt, lag(event_dt), units = "days"))/365.25) %>%
    ungroup()

time_diff <- pheno %>%
    filter(!is.na(time_diff))

time_diff_plot <- ggplot(time_diff, aes(x = time_diff)) +
    geom_histogram(binwidth = 1, fill = "blue", color = "black") +
    labs(title = "Distribution of Time Between eGFR Assessments",
         x = "Time Between Assessments (years)",
         y = "Frequency") +
    scale_x_continuous(limits = c(0, 20)) +
    theme_minimal()

ggsave(file.path(dir_fig, "time_diff_plot.png"), time_diff_plot, width = 6, height = 4)

# events
table(events$transition)

a <- pheno %>% filter(eGFRcrea<30 &is.na(event_dt_eskd)) %>% pull(id) %>% unique

# 2 -> 0 back transitions
idx <- events_withBackTransitions %>% filter(transition == "2 -> 0") %>% pull(id) %>% unique
idx

## spaghetti plots for the individuals in idx
p_14 <- ggplot(pheno %>% filter(id %in% idx), aes(x = age, y = eGFRcrea, color = factor(id))) +
    geom_line(aes(group = id)) +
    labs(x = "Age", y = "eGFRcrea", color = "Individual ID") +
    theme_minimal()
p_14

aki_age <- pheno %>%
    filter(id %in% idx & !is.na(event_dt_aki)) %>%
    select(dateOfBirth, event_dt_aki) %>%
    slice(1) %>%
    mutate(age_at_aki = as.numeric(difftime(event_dt_aki, dateOfBirth, units = "days"))/365.25) %>%
    pull(age_at_aki)
p_14_aki <- ggplot(pheno %>% filter(id %in% idx & !is.na(event_dt_aki)), aes(x = age, y = eGFRcrea, color = factor(id))) +
    geom_line(aes(group = id)) +
    labs(x = "Age", y = "eGFRcrea", color = "Individual ID") +
    theme_minimal() +
    geom_vline(xintercept = aki_age, linetype = "dashed", color = "black")
p_14_aki

# (back) transitions of aki individuals
idx_aki <- pheno %>% filter(!is.na(event_dt_aki)) %>% pull(id) %>% unique

p_aki <- ggplot(pheno %>% filter(id %in% idx_aki), aes(x = age, y = eGFRcrea, color = factor(id))) +
    geom_line(aes(group = id)) +
    labs(x = "Age", y = "eGFRcrea") +
    theme_minimal() +
    theme(legend.position = "none")
p_aki

events_withBackTransitions %>% filter(id %in% idx_aki) %>% pull(transition) %>% table

# transitions of aki individuals if excluding +- 6 months around AKI
events_withBackTransitions_noaki6months <- read_file(file.path(dir_data, file_events_withBackTransitions_noaki6months))
events_withBackTransitions_noaki6months %>% filter(id %in% idx_aki) %>% pull(transition) %>% table
