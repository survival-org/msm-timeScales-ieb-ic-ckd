library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
# library(lubridate)
library(zoo)
library(stringr)
library(mvtnorm)
library(doParallel)
library(foreach)
library(openxlsx)
library(pammtools)
library(ggplot2)
library(viridis)
library(kableExtra)

# data prep ----

read_file <- function(a_infile, header="auto",sep="auto",fill=FALSE, ...){
	tblInRaw <- fread(a_infile, header=header, sep=sep)
	tblIn <- as.data.frame(tblInRaw,fill=fill)
	return(tblIn)
}

extract_date <- function(pheno, code){
    if(code != "emr_death"){
        out <- pheno %>%
            filter(record_type == "emr_diagnostic_code" & diagnostic_code == code) %>%
            select(id, date_at_exam)
    } else if(code == "Diabetes") {
        out <- pheno %>%
            filter(record_type %in% c("emr_diagnostic_code", "sc_diagnostic_code") & diagnostic_code == code) %>%
            select(id, date_at_exam) %>%
            distinct()
    } else {
        out <- pheno %>%
            filter(record_type == code) %>%
            select(id, date_at_exam)
    }

    return(out)
}

extract_riskFactor <- function(pheno, code){
    out <- pheno %>%
        filter(record_type %in% code) %>%
        select(id, date_at_exam, record_value)

    return(out)
}

create_tvc <- function(pheno, risk_factors, diagnostic_codes) {

  dt <- as.data.table(pheno)

  # --- 1. Process risk_factors from 'record_type' ---
  dt_risk <- dt[
    record_type %in% risk_factors,
    .(id, date_at_exam, date_of_birth, record_type, record_value)
  ][,
    age := as.numeric(as.yearmon(date_at_exam) - as.yearmon(date_of_birth)),
    by = id
  ]

  tvc_risk <- dcast(
    dt_risk,
    id + date_at_exam + date_of_birth + age ~ record_type,
    value.var = "record_value"
  )

  # --- 2. Process diagnostic_codes from 'diagnostic_code' ---
  dt_diag <- dt[
    diagnostic_code %in% diagnostic_codes,
    .(id, date_at_exam, date_of_birth, diagnostic_code)
  ][,
    age := as.numeric(as.yearmon(date_at_exam) - as.yearmon(date_of_birth)),
    by = id
  ][,
    record_value := 1  # set binary 1 to indicate presence
  ] %>%
    unique()

  tvc_diag <- dcast(
    dt_diag,
    id + date_at_exam + date_of_birth + age ~ diagnostic_code,
    value.var = "record_value"
  )

  # --- 3. Merge the two wide tables ---
  tvc <- merge(tvc_risk, tvc_diag,
               by = c("id", "date_at_exam", "date_of_birth", "age"),
               all = TRUE)

  # --- 4. Drop date columns ---
  tvc[, c("date_at_exam", "date_of_birth") := NULL]

  return(tvc)
}

merge_tvc <- function(ped, tvc, vars, fill_before_first = TRUE, fill_subsequent = TRUE) {

  ped_dt <- as.data.table(ped)
  tvc_dt <- as.data.table(tvc)

  # Compute ped range for each ID
  ped_range <- ped_dt[, .(min_start = min(tstart), max_start = max(tstart)), by = id]
  tvc_dt <- ped_range[tvc_dt, on = "id", nomatch = 0L]

  # 4a) collapse early observations to first time point
  early <- tvc_dt[age < min_start]
  early2 <- if (nrow(early) > 0) {
    early[order(id, age), .SD[.N], by = id][, age := min_start]
  } else {
    tvc_dt[0]
  }

  # 4b) keep observations within window
  mid <- tvc_dt[age >= min_start & age <= max_start]
  tvc_clean <- rbind(early2, mid, use.names = TRUE)[, c("min_start", "max_start") := NULL]

  # Set key for rolling joins
  setkey(tvc_clean, id, age)

  # 6) Rolling join on tstart
  merged <- tvc_clean[
    ped_dt,
    on = .(id, age <= tstart),
    mult = "last",
    nomatch = NA,
    j = {
      ped_vals <- mget(names(ped_dt))
      tvc_vals <- mget(vars)
      c(ped_vals, tvc_vals, list(tvc_age = age))
    }
  ]

  # 7) Fill before first observed value (per‐var approach)
  if (fill_before_first) {
    # Ensure tvc_clean is ordered by id, then age
    setorder(tvc_clean, id, age)

    # For each var separately, find its first non-NA age & value
    for (v in vars) {
      # 1) Extract only the rows where v is not NA
      first_info_v <- tvc_clean[
        !is.na(get(v)),
        .(first_age = age[1L], first_val = get(v)[1L]),  # take top row for each id
        keyby = id
      ]

      # 2) If there are any first observations for this var, merge them
      if (nrow(first_info_v) > 0L) {
        # Rename columns so they won't clash
        setnames(first_info_v, c("first_age", "first_val"),
                            c(paste0("first_age_", v), paste0("first_val_", v)))

        # 3) Left-join into merged by id
        merged <- merge(
          merged,
          first_info_v,
          by    = "id",
          all.x = TRUE
        )

        # 4) Fill v wherever tstart < first_age_v AND v is still NA
        merged[
          get(paste0("first_age_", v)) > tstart & is.na(get(v)),
          (v) := get(paste0("first_val_", v))
        ]

        # 5) Remove the helper columns for this var
        merged[, c(paste0("first_age_", v), paste0("first_val_", v)) := NULL]
      }
    }
  }

  # 8) Forward-fill subsequent values
  if (fill_subsequent) {
    setorder(merged, id, tstart)
    merged[, (vars) := lapply(.SD, function(x) zoo::na.locf(x, na.rm = FALSE)),
           by = id, .SDcols = vars]
  }

  merged[, tvc_age := NULL]

  return(merged[])

}

CKDEpi.creat.rf <- function (creatinine, sex, age) {

  crea_mgdl_to_umoll <- 88.4

  if (!is.null(creatinine) & !is.null(sex) & !is.null(age)) {
    creatinine <- as.numeric(creatinine) / crea_mgdl_to_umoll
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    n <- length(creatinine)

    if (length(sex) == n & length(age) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age)]

      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10

      # CKD-Epi equation
      k <- a <- numeric(n)
      k[sex == 0] <- 0.7
      k[sex == 1] <- 0.9
      a[sex == 0] <- -0.241
      a[sex == 1] <- -0.302
      one <- rep(1, n)
      eGFR <- apply(cbind(creatinine/k, one), 1, min, na.rm = T)^a * apply(cbind(creatinine/k, one), 1, max, na.rm = T)^-1.200 * 0.9938^age
      eGFR[sex == 0] <- eGFR[sex == 0] * 1.012

      # Restore missing data at the indexed positions
      eGFR[idx] <- NA

      # Output
      142 * eGFR
    } else
      stop("Different number of observations between variables")
  } else
    stop("Some variables are not defined")
}

create_eGFR_final_cutoff <- function(pheno, cutoff){ # requires columns id, eGFRcrea, date_at_exam
    out <- pheno %>%
        group_by(id) %>%
        filter(any(eGFRcrea < cutoff)) %>%  # Keep only individuals with any eGFRcrea < 15
        filter(eGFRcrea < cutoff) %>%        # Filter rows with eGFRcrea < 15
        slice_min(order_by = date_at_exam, n = 1, with_ties = FALSE) %>%  # Get the first date when eGFRcrea < 15
        select(id, date_at_exam) %>%    # Keep only id and date_at_exam columns
        ungroup()

    return(out)
}

add_status <- function(pheno, cutoffs){ # requires columns id, eGFRcrea, date_at_exam, date_at_exam_eskd

	pheno <- pheno[order(id, date_at_exam)]
	pheno[, flag1 := as.integer(eGFRcrea < cutoffs[1])]
	pheno[, flag2 := as.integer(eGFRcrea < cutoffs[2])]
	pheno[, status := 0]
	pheno[date_at_exam == date_at_exam_eskd, status := 3]

	# pheno[, status_withBackTransitions := {
	# s <- status  # Initialize status vector for the group

	# # Handle eGFRcrea < cutoffs[2] for at least 2 consecutive measurements (status = 2)
	# for (i in seq_len(.N)) {
	# 	if (flag2[i] == 1 && i > 1 && flag2[i - 1] == 1) {
	# 	# Backfill status to the start of the consecutive sequence
	# 	start_i <- i - 1
	# 	while (start_i >= 1 && flag2[start_i] == 1) {
	# 		if (s[start_i] < 2) s[start_i] <- 2
	# 		start_i <- start_i - 1
	# 	}
	# 	if (s[i] < 2) s[i] <- 2
	# 	}
	# }

	# # Handle eGFRcrea < cutoffs[1] for at least 2 consecutive measurements (status = 1)
	# for (i in seq_len(.N)) {
	# 	if (flag1[i] == 1 && s[i] == 0 && i > 1 && flag1[i - 1] == 1) {
	# 	# Backfill status to the start of the consecutive sequence
	# 	start_i <- i - 1
	# 	while (start_i >= 1 && flag1[start_i] == 1) {
	# 		if (s[start_i] < 1) s[start_i] <- 1
	# 		start_i <- start_i - 1
	# 	}
	# 	if (s[i] < 1) s[i] <- 1
	# 	}
	# }

	# s  # Return the updated status vector
	# }, by = id]

	pheno[, status := {
		s <- status  # Initialize status vector with current statuses (0 or 3)

		# Assign status 2 where eGFRcrea < cutoffs[2] for at least 2 consecutive measurements
		idx2 <- which(flag2 == 1 & shift(flag2, type = "lag") == 1)
		s[idx2] <- pmax(s[idx2], 2)
		s[idx2 - 1] <- pmax(s[idx2 - 1], 2)

		# Assign status 1 where eGFRcrea < cutoffs[1] for at least 2 consecutive measurements
		idx1 <- which(flag1 == 1 & shift(flag1, type = "lag") == 1)
		s[idx1] <- pmax(s[idx1], 1)
		s[idx1 - 1] <- pmax(s[idx1 - 1], 1)

		# Ensure that status does not decrease over time (excluding transitions to status 3)
		for (i in 2:.N) {
			if (s[i - 1] == 3) {
				s[i] <- 3  # Once in status 3, remain in status 3
			} else if (s[i] != 3) {
				s[i] <- max(s[i - 1], s[i])  # Enforce non-decreasing status
			}
		}
		s  # Return the updated status vector
	}, by = id]

	return(pheno)
}

add_vars <- function(pheno) { # requires columns id, date_at_exam, date_at_exam_eskd, status
    pheno[, date_at_exam_last := max(date_at_exam), by = id]
    pheno[, need_adjustment := !is.na(date_at_exam_eskd) & (date_at_exam_eskd != date_at_exam_last), by = id]
    pheno[, age0_original := age_at_exam[which.min(date_at_exam)], by = id]
    pheno[, age0 := as.numeric(as.yearmon(min(date_at_exam)) - as.yearmon(date_of_birth)), by = id]
    pheno[, age := as.numeric(as.yearmon(date_at_exam) - as.yearmon(date_of_birth)), by = id]
    pheno[, age_eskd := as.numeric(as.yearmon(date_at_exam_eskd) - as.yearmon(date_of_birth)), by = id]
    # pheno[, t_years := (as.numeric(date_at_exam) - min(as.numeric(date_at_exam))) / as.numeric(365.25), by = id]
    # pheno[, t := (as.numeric(date_at_exam) - min(as.numeric(date_at_exam))) / as.numeric(30.44), by = id]
    # pheno[, t_eskd := as.numeric(date_at_exam_eskd - min(date_at_exam))/as.numeric(365.25), by = id]
    # pheno[, t := as.numeric(as.yearmon(date_at_exam) - as.yearmon(min(date_at_exam))), by = id]
    # pheno[, t_eskd := as.numeric(as.yearmon(date_at_exam_eskd) - as.yearmon(min(date_at_exam))), by = id]
    # pheno[, t_last := max(t), by = id]
    pheno[, run_id := rleid(status), by = id]
    # pheno[, run_id_withBackTransitions := rleid(status_withBackTransitions), by = id]

    return(pheno)
}

create_events <- function(pheno, snp_cols, adjust_ids, age_eskd_dt, selected_cols = NULL, sim = FALSE, death_info = TRUE) {
    # requires columns: status, t, age0, eGFRcrea, run_id; and external vectors: snp_cols, adjust_ids, age_eskd_dt
    # additionally, pheno is assumed to include date_at_exam_death and date_at_exam_last (the latter being the date of the last observed event)

    events <- pheno[, c(
        .(
            from = as.character(first(status)),
            agestart = first(age),
            agestop_raw = last(age),
            age0 = first(age0),       # time-constant
            eGFRcrea = first(eGFRcrea), # value at tstart of interval
            date_at_exam = last(date_at_exam) # event date at the transition time point
        ),
        if (!sim) .(
            sex = as.factor(first(sex)),          # time-constant
            diabetes = as.factor(first(diabetes)) # value at tstart of interval
        ),
        # Add SNP columns
        lapply(.SD, first)
        ), by = .(id, run_id), .SDcols = snp_cols]

    events[, to := as.character(shift(from, type = "lead")), by = id]
    events[, is_last_event := seq_len(.N) == .N, by = id]
    events[, need_adjustment := id %in% adjust_ids]

    # For last events, assign terminal state or censoring.
    events[is_last_event == TRUE, to := ifelse(
        need_adjustment,
        3,
        ifelse(from %in% c("0", "1", "2"), "cens", NA)
    )]

    # Create transition column
    events[, transition := paste0(from, "->", to)]

    # Remove transitions from state 3 without a defined "to"
    events <- events[!(from == "3" & is.na(to))]

    # Merge in age_eskd_dt info (as in your original function)
    events <- merge(events, age_eskd_dt, by = 'id', all.x = TRUE)

    # Adjust stop times for individuals needing adjustment
    events[need_adjustment == TRUE & is_last_event == TRUE, agestop := age_eskd]

    # Assign episode numbers
    events[, episode := seq_len(.N), by = .(id, transition)]

    # Set default status: if the event is censoring, status = 0; otherwise 1.
    events[, status := ifelse(to == "cens", 0, 1)]

    # Determine agestop: if state is "3" then use age_eskd, else use the next agestart
    events[, agestop := ifelse(
        to == "3",
        age_eskd,
        shift(agestart, type = "lead")
    ), by = id]
    events[is.na(agestop), agestop := agestop_raw]

    if(death_info){
        death_df <- unique(pheno[, .(id, date_at_exam_death, date_at_exam_last)])
        events <- merge(events, death_df, by = "id", all.x = TRUE)

        events[from %in% c("0", "1", "2") & to == "cens" & # Update transitions originally set as censoring from states 0, 1, or 2.
               !is.na(date_at_exam_death) & date_at_exam_death > date_at_exam_last, # For these rows, if a death date is observed and it occurs after the last event date,
            `:=`(to = "death", status = 1)] # then update: set to = "death", status = 1
        events[from %in% c("0", "1", "2") & to == "death",
            transition := paste0(from, "->", to)] # and update the transition string
    }

    setorder(events, id, agestart)
    if (is.null(selected_cols)) selected_cols <- colnames(events)
    events <- events[, ..selected_cols]

    # Remove rows where from equals to (invalid transitions)
    events <- events[from != to]

    return(events)
}

adjust_events <- function(events, pheno) {

  # Ensure date columns are of type Date
  events <- events %>%
    mutate(date_at_exam = as.Date(date_at_exam))

  pheno <- pheno %>%
    mutate(date_at_exam = as.Date(date_at_exam))

  # Merge in the age value from pheno (renamed as agestop_lag)
  events <- events %>%
    left_join(pheno %>% select(id, date_at_exam, agestop_lag = age),
              by = c("id", "date_at_exam"))

  # Arrange events by id and agestart and then update sequentially
  events_updated <- events %>%
    arrange(id, agestart) %>%
    group_by(id) %>%
    group_modify(~ {
      df <- .x
      n <- nrow(df)
      # Create columns to hold the new (adjusted) values
      df$new_agestart <- rep(NA_real_, n)
      df$new_agestop  <- rep(NA_real_, n)

      for (i in seq_len(n)) {
        if (i == 1) {
          # First transition: leave agestart unchanged
          df$new_agestart[i] <- df$agestart[i]
          # For agestop, if this single event is final with to=="3", leave unchanged; otherwise update
          if (n == 1 && df$to[i] == "3") {
            df$new_agestop[i] <- df$agestop[i]
          } else {
            df$new_agestop[i] <- (df$agestop[i] + df$agestop_lag[i]) / 2
          }
        } else {
          # For later transitions: set new agestart to the updated agestop of the previous transition
          df$new_agestart[i] <- df$new_agestop[i - 1]
          # And update agestop as the average, except if this is the last transition and to=="3"
          if (i == n && df$to[i] == "3") {
            df$new_agestop[i] <- df$agestop[i]
          } else {
            df$new_agestop[i] <- (df$agestop[i] + df$agestop_lag[i]) / 2
          }
        }
      }

      return(df)
    }) %>% ungroup()

  # Replace the original agestart and agestop with the updated ones
  events_updated <- events_updated %>%
    mutate(
      agestart_unadjusted = agestart,
      agestop_unadjusted = agestop,
      agestart = new_agestart,
      agestop = new_agestop
    ) %>%
    select(-new_agestart, -new_agestop)

  return(events_updated)
}

impute_events <- function(events, pheno, cutoffs){

    formula <- "eGFRcrea ~ age"

    # 1. Compute idx_exclude:
    #    a) Identify individuals with a 1->3 or 0->3 transition in events.
    idx_candidates <- events %>% filter(transition %in% c("1->3", "0->3")) %>% pull(id) %>% unique()

    #    b) Get each individual’s last row in pheno (by maximum age)
    pheno_last <- pheno %>% group_by(id) %>% slice_max(age, with_ties = FALSE) %>% ungroup()

    #    c) For those individuals in idx_candidates, exclude if their last pheno row has eGFRcrea > cutoffs[3], e.g. 15
    idx_exclude <- pheno_last %>% filter(id %in% idx_candidates & eGFRcrea >= cutoffs[3]) %>% pull(id)

    # 2. Select individuals with jumps that are not in idx_exclude
    idx <- events %>% filter(transition %in% c("0->2", "1->3", "0->3")) %>% pull(id) %>% unique()
    idx <- setdiff(idx, idx_exclude)

    events_idx <- events %>% filter(id %in% idx)
    pheno_idx <- pheno %>% filter(id %in% idx)

    new_events_list <- list()

    for(i in idx){

        transitions_i <- events_idx %>% filter(id == i) %>% pull(transition)
        if (sum(transitions_i %in% c("0->2", "1->3", "0->3")) > 1) {
            stop(paste("Individual", i, "has multiple jumps!"))
        }

        if(!"0->3" %in% transitions_i){

            if("0->2" %in% transitions_i){
                from_state <- as.character(0)
                to_state <- as.character(2)
                imputed_state <- as.character(1)
                cutoff_egfr <- cutoffs[1] - 1 # e.g. 59
            } else {
                from_state <- as.character(1)
                to_state <- as.character(3)
                imputed_state <- as.character(2)
                cutoff_egfr <- cutoffs[2] - 1 # e.g. 29
            }

            events_i <- events_idx %>% filter(id == i)
            if("0->2" %in% transitions_i){
                cutoff_age_i <- events_i %>% filter(from == to_state) %>% pull(agestart)
                pheno_i <- pheno_idx %>% filter(id == i & age <= cutoff_age_i) %>% arrange(date_at_exam)
            } else {
                pheno_i <- pheno_idx %>% filter(id == i) %>% arrange(date_at_exam)
            }

            from_age_i <- events_i %>% filter(from == from_state) %>% pull(agestart)
            to_age_i <- events_i %>% filter(from == from_state) %>% pull(agestop)

            mod <- lm(formula, data=pheno_i)

            if(nrow(pheno_i) >= 2){
                coefs <- coef(mod)
                beta0 <- coefs[1]
                beta1 <- coefs[2]
                transition_age_i_calculated <- (cutoff_egfr - beta0) / beta1
            }

            if(!is.na(mod$coefficients["age"]) & mod$coefficients["age"] < 0 & transition_age_i_calculated > from_age_i & transition_age_i_calculated < to_age_i){
                transition_age_i <- transition_age_i_calculated
            } else {
                transition_age_i <- from_age_i + (to_age_i - from_age_i)/2
            }

            row_to_split <- events_i %>% filter(from == from_state, to == to_state)
            row_original <- row_to_split
            row_original$to <- imputed_state
            row_original$transition <- paste0(from_state, "->", imputed_state)
            row_original$agestop <- transition_age_i
            new_row <- row_to_split
            new_row$from <- imputed_state
            new_row$transition <- paste0(imputed_state, "->", to_state)
            new_row$agestart <- transition_age_i
            new_row$eGFRcrea <- cutoff_egfr

            events_i <- events_i %>%
                filter(!(from == from_state & to == to_state)) %>%
                bind_rows(row_original, new_row)

            new_events_list[[i]] <- events_i

        } else{

            from_state <- as.character(0)
            to_state <- as.character(3)
            imputed_states <- c(as.character(1), as.character(2))
            cutoff_egfrs <- cutoffs[1:2] - 1 # e.g. 59, 29

            events_i <- events_idx %>% filter(id == i)
            pheno_i <- pheno_idx %>% filter(id == i) %>% arrange(date_at_exam) # no specification of cutoff age necessary because ESKD is last event

            from_age_i <- events_i %>% filter(from == from_state) %>% pull(agestart)
            to_age_i <- events_i %>% filter(from == from_state) %>% pull(agestop)

            mod <- lm(formula, data=pheno_i)

            # if(nrow(pheno_i)==1) stop(paste("Individual", i, "has only one observation!"))
            if(nrow(pheno_i) >= 2){
                coefs <- coef(mod)
                beta0 <- coefs[1]
                beta1 <- coefs[2]
                transition_ages_i_calculated <- c((cutoff_egfrs[1] - beta0) / beta1, (cutoff_egfrs[2] - beta0) / beta1)
            }

            if(!is.na(mod$coefficients["age"]) & mod$coefficients["age"] < 0 & transition_ages_i_calculated[1] < transition_ages_i_calculated[2] & transition_ages_i_calculated[1] > from_age_i & transition_ages_i_calculated[2] < to_age_i){
                transition_ages_i <- transition_ages_i_calculated
            } else {
                transition_ages_i <- c(from_age_i + (to_age_i - from_age_i)/3, from_age_i + (to_age_i - from_age_i)/3*2)
            }

            row_to_split <- events_i %>% filter(from == from_state, to == to_state)
            row_original <- row_to_split
            row_original$to <- imputed_states[1]
            row_original$transition <- paste0(from_state, "->", imputed_states[1])
            row_original$agestop <- transition_ages_i[1]
            new_row_1 <- row_to_split
            new_row_1$from <- imputed_states[1]
            new_row_1$to <- imputed_states[2]
            new_row_1$transition <- paste0(imputed_states[1], "->", imputed_states[2])
            new_row_1$agestart <- transition_ages_i[1]
            new_row_1$agestop <- transition_ages_i[2]
            new_row_1$eGFRcrea <- cutoff_egfrs[1]
            new_row_2 <- row_to_split
            new_row_2$from <- imputed_states[2]
            new_row_2$transition <- paste0(imputed_states[2], "->", to_state)
            new_row_2$agestart <- transition_ages_i[2]
            new_row_2$eGFRcrea <- cutoff_egfrs[2]

            events_i <- events_i %>%
                filter(!(from == from_state & to == to_state)) %>%
                bind_rows(row_original, new_row_1, new_row_2)

            new_events_list[[i]] <- events_i

        }
    }

    # Bind the imputed events for individuals not in idx_exclude
    events_modified <- bind_rows(new_events_list)

    # For individuals in idx_exclude, remove their last event row.
    events_exclude <- events %>% filter(id %in% idx_exclude)
    events_exclude_modified <- events_exclude %>%
      group_by(id) %>%
      arrange(agestart) %>%
      filter(row_number() != n()) %>%
      ungroup()

    # For individuals who were not processed (i.e. not in idx or idx_exclude)
    events_rest <- events %>% filter(!id %in% c(idx, idx_exclude))

    # Combine all events and order by id and agestart.
    events_new <- bind_rows(events_rest, events_modified, events_exclude_modified) %>%
      arrange(id, agestart)

    return(events_new)

}

# because r() uses round half to even, which is not what we want here
round_up_half <- function(x, digits = 0) {
  floor(x * 10^digits + 0.5) / 10^digits
}

# merge_events <- function(events, age_scale=FALSE) {
#     # for all rows where to != "cens" & the row with tstart==tstop is the individual's last row, we need to merge the 2 last rows:
#     # from 2nd last row: id, tstart, from, status, episode, age0, sex, eGFRcrea, diabetes, snp_cols,
#     # from last row: tstop_raw, tstop, to
#     # recompute transition for these rows
#     last_rows <- events[, .I[.N], by = id]$V1
#     second_last_rows <- events[, .I[.N - 1], by = id]$V1

#     if(!age_scale) {
#         rows_to_merge <- events[, .I[tstart == tstop & .I != .I[1]], by = id]$V1

#         while (length(rows_to_merge) > 0) {
#             last_row <- rows_to_merge[1]
#             second_last_row <- last_row - 1

#             if (second_last_row %in% events[, .I]) {
#                 # Save updated "to" value
#                 new_to <- events[last_row]$to

#                 # Update the second-last row
#                 events[second_last_row, `:=`(
#                     tstop_raw = events[last_row]$tstop_raw,
#                     tstop = events[last_row]$tstop,
#                     to = new_to,
#                     transition = paste0(from, "->", new_to)
#                 )]

#                 # Delete the last row
#                 events <- events[-last_row]
#             }

#             # Recalculate rows after each merge
#             last_rows <- events[, .I[.N], by = id]$V1
#             second_last_rows <- events[, .I[.N - 1], by = id]$V1
#             rows_to_merge <- events[, .I[tstart == tstop & .I != .I[1]], by = id]$V1
#         }
#     } else {
#         rows_to_merge <- events[, .I[agestart == agestop & .I != .I[1]], by = id]$V1

#         while (length(rows_to_merge) > 0) {
#             last_row <- rows_to_merge[1]
#             second_last_row <- last_row - 1

#             if (second_last_row %in% events[, .I]) {
#                 # Save updated "to" value
#                 new_to <- events[last_row]$to

#                 # Update the second-last row
#                 events[second_last_row, `:=`(
#                     # tstop_raw = events[last_row]$agestop_raw,
#                     # tstop = events[last_row]$agestop,
#                     to = new_to,
#                     transition = paste0(from, "->", new_to)
#                 )]

#                 # Delete the last row
#                 events <- events[-last_row]
#             }

#             # Recalculate rows after each merge
#             last_rows <- events[, .I[.N], by = id]$V1
#             second_last_rows <- events[, .I[.N - 1], by = id]$V1
#             rows_to_merge <- events[, .I[agestart == agestop & .I != .I[1]], by = id]$V1
#         }
#     }

#     return(events)
# }

round_events <- function(events, r=2, age_scale=TRUE) {

    setDT(events)
    # remove all rows of events where tstart == tstop and to == "cens"
    # remove all rows of events where tstart == tstop and it's the first row per id
    if(!age_scale) {

        # round tstart, tstop_raw, tstop
        events[, tstart := round_up_half(tstart, r)]
        events[, tstop_raw := round_up_half(tstop_raw, r)]
        events[, tstop := round_up_half(tstop, r)]

        events[, age0 := round_up_half(age0, r)]
        events[, agestart := round_up_half(agestart, r)]
        events[, agestop_raw := round_up_half(agestop_raw, r)]
        events[, agestop := round_up_half(agestop, r)]

        events <- events[!(tstart == tstop)]

    } else {

        # round tstart, tstop_raw, tstop
        events[, age0 := round_up_half(age0, r)]
        events[, agestart := round_up_half(agestart, r)]
        events[, agestop_raw := round_up_half(agestop_raw, r)]
        events[, agestop := round_up_half(agestop, r)]

        events <- events[!(agestart == agestop)]
    }

    # # merge all rows of events where tstart == tstop and it's not the first row per id
    # events <- merge_events(events, age_scale=age_scale)

    return(events)
}

censor_events <- function(events, c) {
    events <- events %>%
        filter(tstart <= c) %>%
        mutate(
            tstop = ifelse(tstop > c, c, tstop),
            to = ifelse(tstop > c, "cens", to),
            transition = ifelse(tstop > c, paste0(from, "->", to), transition),
            status = ifelse(tstop > c, 0, status)
        )

    return(events)
}

# reset_clock <- function(events) {
#     events <- events %>%
#         mutate(
#             agestart_old = agestart,
#             agestop_old = agestop,
#             agestart_bl = ifelse(from == 0, 0, agestart),
#             agestop = ifelse(from == 0, agestop, agestop - agestart),
#             agestart = ifelse(from == 0, agestart, 0)
#         )

#     # events <- events %>%
#     #     group_by(id) %>%
#     #     arrange(agestart) %>%
#     #     mutate(duration_ckd = if_else(from == 2 & lag(transition) == "1->2", lag(agestop) - lag(agestart), 0)) %>%
#     #     ungroup()

#     setDT(events)                     # convert to data.table if not already
#     setorder(events, id, agestart)    # sort by id and agestart

#     events[, duration_ckd := ifelse(
#     from == 2 & shift(transition, type = "lag") == "1->2",
#     shift(agestop, type = "lag") - shift(agestart, type = "lag"),
#     0
#     ), by = id]

#     return(events)
# }

add_timeScales <- function(ped_events) {

    saved_attrs <- attributes(ped_events)

    out <- ped_events %>%
        # For each individual, determine the onset and progression ages:
        group_by(id) %>%
        mutate(
            age_onset =
                if(any(from == 1)) {
                    first(tstart[from == 1])
                } else if(any(to == 1)){
                    last(tend[to == 1])
                } else if(any(from== 2)){ # for those whose first transition is 2->
                    first(tstart[from == 2])
                } else {
                    0
                },
            age_progression =
                if(any(from == 2)) {
                    first(tstart[from == 2])
                } else if(any(to == 2)){
                    last(tend[to == 2])
                } else {
                    0
                },
        ) %>%
        ungroup() %>%
        # Now, compute the onset columns,
        # forcing them to 0 for rows with 0-> transitions.
        mutate(
            age_onset = case_when(
                from == 0 ~ 0,
                TRUE ~ age_onset
            ),
            tstart_onset = case_when(
                age_onset == 0 ~ 0,
                from == 0 ~ 0,
                TRUE ~ tstart - age_onset
            ),
            tend_onset = case_when(
                age_onset == 0 ~ 0,
                from == 0 ~ 0,
                TRUE ~ tend - age_onset
            ),
            # And compute the progression columns:
            # For individuals with no 2-> transition (i.e. age_progression == 0)
            # or for rows that occur before the progression event,
            # we set progression values to 0.
            age_progression = case_when(
                from %in% c(0, 1) ~ 0,
                TRUE ~ age_progression
            ),
            tstart_progression = case_when(
                age_progression == 0 ~ 0,
                from %in% c(0, 1) ~ 0,
                TRUE ~ tstart - age_progression
            ),
            tend_progression = case_when(
                age_progression == 0 ~ 0,
                from %in% c(0, 1) ~ 0,
                TRUE ~ tend - age_progression
            ),
            time_until_progression = case_when(
                from == 2 & age_progression > 0 & age_onset > 0 ~ age_progression - age_onset,
                TRUE ~ 0
            )
        )

    # Reapply all original attributes (excluding names, row.names, and class which mutate preserves)
    attr_names <- setdiff(names(saved_attrs), c("names", "row.names", "class"))
    for (a in attr_names) {
    attr(out, a) <- saved_attrs[[a]]
    }

    return(out)
}

split_ped <- function(events_ped) {

    out <- list()

    out$transition_01 <- events_ped %>% filter(from == 0)
    out$transition_12 <- events_ped %>% filter(from == 1)
    out$transition_23 <- events_ped %>% filter(from == 2)

    return(out)
}

# simulations ----

# createAgeTrajectories <- function(lengths, mu1, sd1, t2, t3=0, t4=0, t5=0) {
#   n <- sum(lengths)
#   result_list <- vector("list", n)
#   lengths_vec <- rep(1:5, lengths)
#   num_timepoints <- sum(lengths != 0) # assuming no "gaps" in the lengths vector

#   start_index <- 1
#   for (i in 1:num_timepoints) {
#     end_index <- start_index + lengths[i] - 1
#     for (j in start_index:end_index) {
#       result_list[[j]] <- numeric(i)
#       result_list[[j]][1] <- rnorm(1, mu1, sd1)
#       if (i > 1) result_list[[j]][2] <- result_list[[j]][1] + t2
#       if (i > 2) result_list[[j]][3] <- result_list[[j]][2] + t3
#       if (i > 3) result_list[[j]][4] <- result_list[[j]][3] + t4
#       if (i > 4) result_list[[j]][5] <- result_list[[j]][4] + t5
#     }
#     start_index <- end_index + 1
#   }

#   result_list <- setNames(result_list, seq_along(result_list))
#   return(result_list)
# }

create_phenoData <- function(trajectories, sigmaRE, MAF, intercept, ageEffect, trueInterceptEffect, trueSlopeEffect, sigma=1, n=NULL){

    if(is.null(n)) n <- length(trajectories)

    # sample trajectories
    subjects <- sample(names(trajectories), n)
    sampleTrajectories <- trajectories[subjects]
    m <- length(unlist(sampleTrajectories))

    # random effects
    REs <- mvtnorm::rmvnorm(n, mean = c(0,0), sigma = sigmaRE)

    G <- rbinom(n, 2, MAF)

    dataList <- lapply(seq_along(sampleTrajectories), function(i) {

        subjectTrajectory <- sampleTrajectories[[i]]
        len <- nrow(subjectTrajectory)
        age_c50_10 <- (subjectTrajectory$age - 50)/10

        data.frame(
            id = rep(subjects[i], len),
            date_at_exam = subjectTrajectory$date_at_exam,
            eGFRcrea = (intercept + REs[i, 1] + G[i] * trueInterceptEffect) + age_c50_10 * (ageEffect + REs[i, 2] + G[i] * trueSlopeEffect) + rnorm(len, 0, sigma),
            age = subjectTrajectory$age, # original-scale age
            date_of_birth = subjectTrajectory$date_of_birth,
            G = rep(G[i], len)
        )
    })
    data <- do.call(rbind, dataList)

    # # Set up parallel backend: Use available cores minus one (adjust as needed)
    # cores <- parallel::detectCores() - 30
    # cl <- makeCluster(cores)
    # registerDoParallel(cl)

    # # Parallel loop using foreach: iterate over subjects (sampleTrajectories)
    # dataList <- foreach(i = seq_along(sampleTrajectories), .combine = rbind, .packages = "mvtnorm") %dopar% {

    #     subjectTrajectory <- sampleTrajectories[[i]]
    #     len <- nrow(subjectTrajectory)
    #     # Center and scale age (example: subtract 50 and divide by 10)
    #     age_c50_10 <- (subjectTrajectory$age - 50) / 10

    #     # Create a data frame for this subject
    #     data.frame(
    #     id = rep(subjects[i], len),
    #     date_at_exam = subjectTrajectory$date_at_exam,
    #     eGFRcrea = (intercept + REs[i, 1] + G[i] * trueInterceptEffect) +
    #                 age_c50_10 * (ageEffect + REs[i, 2] + G[i] * trueSlopeEffect) +
    #                 rnorm(len, 0, sigma),
    #     age = subjectTrajectory$age,        # original-scale age
    #     date_of_birth = subjectTrajectory$date_of_birth,
    #     G = rep(G[i], len),
    #     stringsAsFactors = FALSE
    #     )
    # }

    # # Stop the parallel cluster
    # stopCluster(cl)

    rownames(data) <- NULL

    return(data)
}

create_eventData <- function(data, rounding = NULL, age_scale = TRUE, censoring = NULL) {

    eGFR15 <- create_eGFR15(data)
    data <- data %>%
        left_join(
            eGFR15 %>% group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>% rename(date_at_exam_eskd = date_at_exam),
            by = "id")

    setDT(data)
    data[, date_at_exam := as.IDate(date_at_exam, format="%Y-%m-%d")]
    data <- data[, .SD[date_at_exam <= date_at_exam_eskd | is.na(date_at_exam_eskd)], by = id]

    data <- add_status(data)
    data <- add_vars(data)

    adjust_ids <- data[need_adjustment == TRUE, unique(id)]
    age_eskd_dt <- data[!is.na(age_eskd), .(age_eskd = unique(age_eskd)), by = id]
    snp_cols <- c("G")
    selected_cols <- c(
    "id", "tstart", "tstop_raw", "tstop", "from", "to", "transition",
    "status", "episode", "age0", "agestart", "agestop_raw", "agestop", "eGFRcrea",
    snp_cols
    )
    events <- create_events(data, snp_cols=snp_cols, selected_cols=selected_cols, adjust_ids=adjust_ids, age_eskd_dt=age_eskd_dt, sim=TRUE)

    ### rounding (due to too many event time points otherwise)
    if(!is.null(rounding)) events <- round_events(events, r=rounding, age_scale=age_scale)

    ### administrative censoring
    if(!is.null(censoring) & !age_scale) events <- censor_events(events, c=censoring)

    return(events)
}

simulate <- function(
    nSims, ncores,
    datatype = "cohort", # "cohort" or "empirical"
    lengths=NULL, mu1=NULL, sd1=NULL, t2=NULL, t3=NULL, t4=NULL, t5=NULL,
    empiricalTrajectories=NULL, n=NULL,
    sigmaRE, MAF, intercept, ageEffect, trueSlopeEffect, trueInterceptEffect, sigma=1,
    dir_out="", return_results=FALSE){

    # results loop
    registerDoParallel(cores = ncores)

    results <- foreach(i = 1:nSims, .combine = 'rbind', .packages = c("dplyr", "lme4", "lmerTest", "mvtnorm")) %dopar% {
        if(datatype == "cohort") {
            ageTrajectories <- createAgeTrajectories(lengths, mu1, sd1, t2, t3, t4, t5)
        } else if(datatype == "empirical") {
            ageTrajectories <- sample(empiricalTrajectories, n, replace = FALSE)
        }
        df <- createData(
            trajectories = ageTrajectories, sigmaRE = sigmaRE, MAF = MAF, intercept = intercept, ageEffect = ageEffect,
            trueSlopeEffect = trueSlopeEffect, trueInterceptEffect = trueInterceptEffect, sigma = sigma, n = NULL)
        res_list <- lapply(names(models), function(model) {
            runModel(df = df, model = model, f = models[[model]], trueSlopeEffect = trueSlopeEffect)})
        res <- do.call(rbind, res_list)
        return(res)}

    stopImplicitCluster()

    # summarize results

    summary_df <- results %>%
        group_by(model) %>%
        summarise(
            avg_coef = mean(coef, na.rm = TRUE),
            se_avg_coef = sd(coef, na.rm = TRUE) / sqrt(n()),
            ci_avg_coef = paste0(round(avg_coef - 1.96 * se_avg_coef, 3), "; ", round(avg_coef + 1.96 * se_avg_coef, 3)),
            power = sum(p < alpha, na.rm = TRUE) / n(),
            power_ci_lower = binom.test(sum(p < alpha), n(), 0.5)$conf.int[1],
            power_ci_upper = binom.test(sum(p < alpha), n(), 0.5)$conf.int[2],
            power_ci = paste0(power_ci_lower, "; ", power_ci_upper)
        ) %>%
        arrange(match(model, names(models)))

    path_out_summary <- file.path(
        dir_out,
        paste0(
            "powerSims_nSims", nSims,
            # "_mu1", mu1, "_sd1", sd1, "_t2", t2, "_t3", t3, "_t4", t4, "_t5", t5,
            datatype,
            "_sigmaRE", sigmaRE[1,1], "+", sigmaRE[2,2], "+", sigmaRE[2,1], "_MAF", MAF, "_intercept", intercept, "_ageEffect", ageEffect,
            "_trueSlopeEffect", trueSlopeEffect, "_trueInterceptEffect", trueInterceptEffect, "_sigma", sigma,
            "_alpha", alpha,
              ".txt"))

    path_out_raw <- file.path(
        dir_out,
        paste0(
            "powerSims_nSims", nSims,
            # "_mu1", mu1, "_sd1", sd1, "_t2", t2, "_t3", t3, "_t4", t4, "_t5", t5,
            datatype,
            "_sigmaRE", sigmaRE[1,1], "+", sigmaRE[2,2], "+", sigmaRE[2,1], "_MAF", MAF, "_intercept", intercept, "_ageEffect", ageEffect,
            "_trueSlopeEffect", trueSlopeEffect, "_trueInterceptEffect", trueInterceptEffect, "_sigma", sigma,
            "_alpha", alpha,
            "_rawData",
              ".txt"))

    fwrite(summary_df, path_out_summary, sep="\t")
    fwrite(results, path_out_raw, sep="\t")

    if(return_results) return(summary_df)

}

add_transVars <- function(ped) {
    out <- ped %>%
        mutate(
            transition_to_death = factor(ifelse(to == 4, 1, 0), levels = c(0, 1)),
            transition_after_onset = factor(
                case_when(
                    from == 0 ~ "none",
                    transition %in% c("1->2", "2->3") ~ "progression",
                    transition %in% c("1->4", "2->4") ~ "death"
                ),
                levels = c("none", "progression", "death"),
                ordered = TRUE
            ),
            transition_after_progression = factor(
                case_when(
                    from %in% c(0, 1) ~ "none",
                    transition == "2->3" ~ "eskd",
                    transition == "2->4" ~ "death"
                ),
                levels = c("none", "eskd", "death"),
                ordered = TRUE
            ),
            transition_after_onset_strat = factor(
                case_when(
                    from == 0 ~ "none",
                    TRUE ~ transition
                ),
                levels = c("none", "1->2", "2->3", "1->4", "2->4"),
                ordered = TRUE
            )
        )

    return(out)
}

make_newped <- function(ped, mod, ci = TRUE) {

  from_states <- sort(unique(ped$from))

  obligatory_cols <- c("tstart", "tend", "intlen", "from", "to", "transition", "age_onset", "tend_onset")

  mod_vars <- attr(terms(mod), "term.labels")

  if ("age_progression" %in% mod_vars) {
    obligatory_cols <- c(
      obligatory_cols,
      "age_progression",
      "tend_progression",
      mod_vars[!mod_vars %in% obligatory_cols & !grepl(":", mod_vars)]
    )
  }

  if("age_progression" %in% mod_vars) {
    ped_new_in <- ped %>%
      make_newdata(
        transition              = unique(transition),
        tend                    = sort(unique(tend)),
        age_onset               = sort(unique(tend)) - 1, # hard coded!!!
        age_progression         = sort(unique(tend)) - 1 # hard coded!!!
      )
  } else {
    ped_new_in <- ped %>%
      make_newdata(
        transition              = unique(transition),
        tend                    = sort(unique(tend)),
        age_onset               = sort(unique(tend)) - 1 # hard coded!!!
      )
  }

  ped_new_in <- ped_new_in %>%
    mutate(
      from = as.integer(sub("([0-9])->.*","\\1", transition)),
      to   = as.integer(sub(".*->([0-9])","\\1", transition)),
      tstart = tend - 1, # hard coded!!!
      intlen = tend - tstart,
      age_onset = ifelse(from == 0, 0, age_onset),
      tend_onset = ifelse(from == 0, 0, tend - age_onset),
      age_progression = ifelse(from %in% c(0,1), 0, age_progression),
      tend_progression = ifelse(from %in% c(0,1), 0, tend - age_progression)
    ) %>%
    add_transVars()

  ped_new <- bind_rows(
    lapply(from_states, function(from) {
      ped_new_trans_in <- ped_new_in %>% filter(from == !!from)

      if (from == "0") {
        ped_new_trans <- ped_new_trans_in %>%
          distinct() # for transition 0->1, for a given tend, all rows are identical (because age_onset = 0 always); same for 0->4

      } else if (from == "1") {
        ped_new_trans <- ped_new_trans_in %>%
          distinct() %>%
          filter(tend_onset >= 0 & tend_onset <= 23) %>%
          filter(tstart >= age_onset) %>%
          filter(age_onset < 79)

      } else if (from == "2") {
        ped_new_trans <- ped_new_trans_in %>%
          distinct() %>%
          filter(tend_onset >= 0 & tend_onset <= 23) %>%
          filter(tstart >= age_onset) %>%
          filter(tend_progression >= 0 & tend_progression <= 16) %>%
          filter(age_progression == 0 | age_progression > age_onset) %>%
          filter(tstart >= age_progression) %>%
          filter(age_progression < 79) %>%
          filter(!(tend_onset == 23 & tend_progression == 1))
      }

      return(ped_new_trans)
    })
  )

  ped_new <- ped_new %>%
    select(one_of(obligatory_cols)) %>%
    group_by(age_progression, age_onset, transition) %>%
    arrange(age_progression, age_onset, transition, tend) %>%
    add_hazard(mod, type = "link", ci = ci) %>%
    rename(loghazard = hazard, loghazard_lower = ci_lower, loghazard_upper = ci_upper) %>%
    add_hazard(mod, ci = ci) %>%
    rename(hazard_lower = ci_lower, hazard_upper = ci_upper) %>%
    add_cumu_hazard(mod, ci = ci) %>%
    add_trans_prob(mod, ci = ci) %>%
    ungroup() %>%
    mutate(
      transition = factor(transition, levels = c("0->1", "1->2", "2->3", "0->4", "1->4", "2->4")),
      from = as.character(sub("([0-9])->.*","\\1", transition)),
      to   = as.character(sub(".*->([0-9])","\\1", transition))
    )

  return(ped_new)
}

extract_ranges <- function(ped_new_1, ped_new_2, transitions = c("0->1", "0->4", "1->2", "1->4", "2->3", "2->4"), scales = c("loghazard", "hazard", "cumu_hazard")) {

  ranges <- list()

  for (scale in scales) {
  for (transition in transitions) {
    # Create variable name
    transition_code <- gsub("->", "", transition)
    var_name <- paste0("fit_range_", scale, "_", transition_code)

    # Calculate range across both datasets
    values_1 <- ped_new_1[[scale]][ped_new_1$transition == transition]
    values_2 <- ped_new_2[[scale]][ped_new_2$transition == transition]

    min_val <- min(min(values_1, na.rm = TRUE), min(values_2, na.rm = TRUE))
    max_val <- max(max(values_1, na.rm = TRUE), max(values_2, na.rm = TRUE))

    ranges[[var_name]] <- c(min_val, max_val)
  }
  }

  return(ranges)
}


create_2d_plots <- function(ped_new, model, trans, ranges) {

  ped_new_trans <- ped_new %>% filter(transition == trans)
  digs       <- regmatches(trans, regexec("([0-9]+)->([0-9]+)", trans))[[1]]
  from_digit <- digs[2]; to_digit <- digs[3]
  fit_range_loghazard <- ranges[[paste0("fit_range_loghazard_", from_digit, to_digit)]]
  fit_range_hazard <- ranges[[paste0("fit_range_hazard_", from_digit, to_digit)]]
  fit_range_cumu_hazard <- ranges[[paste0("fit_range_cumu_hazard_", from_digit, to_digit)]]

  p_loghazard <- ggplot(ped_new_trans, aes(x = tend, y = loghazard)) +
    # geom_ribbon(aes(ymin = loghazard_lower, ymax = loghazard_upper), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = fit_range_loghazard) +
    labs(title = paste("Transition", trans), x = "Age", y = "Log-hazard")

  p_hazard <- ggplot(ped_new_trans, aes(x = tend, y = hazard)) +
    # geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = fit_range_hazard) +
    labs(title = paste("Transition", trans), x = "Age", y = "Hazard")

  p_cumu <- ggplot(ped_new_trans, aes(x = tend, y = cumu_hazard)) +
    # geom_ribbon(aes(ymin = cumu_lower, ymax = cumu_upper), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = fit_range_cumu_hazard) +
    labs(title = paste("Transition", trans), x = "Age", y = "Cumulative Hazard")

  p_tp <- ggplot(ped_new_trans, aes(x = tend, y = trans_prob)) +
    # geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by=0.2)) +
    labs(title = paste("Transition", trans), x = "Age", y = "Transition Probability")

  ggsave(file.path(dir_out, "figures", model, sprintf("line_loghazard_age_%s_%s.png", from_digit, to_digit)),
    plot = p_loghazard, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("line_hazard_age_%s_%s.png", from_digit, to_digit)),
    plot = p_hazard, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("line_cumu_age_%s_%s.png", from_digit, to_digit)),
    plot = p_cumu, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("line_tp_age_%s_%s.png", from_digit, to_digit)),
    plot = p_tp, width = 10, height = 10, units = "cm")

  out <- list(
    p_loghazard = p_loghazard,
    p_hazard = p_hazard,
    p_cumu = p_cumu,
    p_tp = p_tp
  )

  return(out)

}

create_3d_plots <- function(ped_new, model, trans, ranges, time_scale, age_onset_slices, prog_after_onset_slices = NULL, max_time_since_onset = 23, max_time_since_progression = 16) {

  ped_new_trans <- ped_new %>% filter(transition == trans)
  digs <- regmatches(trans, regexec("([0-9]+)->([0-9]+)", trans))[[1]]
  from_digit <- digs[2]; to_digit <- digs[3]
  fit_range_loghazard <- range(ped_new_trans$loghazard, na.rm = TRUE)
  fit_range_hazard <- range(ped_new_trans$hazard, na.rm = TRUE)
  fit_range_loghazard <- ranges[[paste0("fit_range_loghazard_", from_digit, to_digit)]]
  fit_range_hazard <- ranges[[paste0("fit_range_hazard_", from_digit, to_digit)]]
  fit_range_cumu_hazard <- ranges[[paste0("fit_range_cumu_hazard_", from_digit, to_digit)]]

  if(time_scale == "age") {
    time <- "tend"
    xlab <- "Age"
  } else if (time_scale == "time") {
    if(from_digit == 1) time <- "tend_onset" else if(from_digit == 2) time <- "tend_progression"
    if(from_digit == 1) xlab <- "Time since CKD onset" else if(from_digit == 2) xlab <- "Time since CKD progression"
  } else {
    stop("x_scale must be either 'age' or 'time'")
  }

  if(is.null(prog_after_onset_slices) | trans %in% c("1->2", "1->4")) {

    # contour plots
    p_contour_loghazard <- ggplot(ped_new_trans, aes(x = !!as.symbol(time), y = age_onset, z = loghazard)) +
      geom_tile(aes(fill = loghazard)) +
      stat_contour(aes(z = loghazard), colour = "grey30", bins = 6) +
      scale_fill_gradient2(name = "Log-hazard", low = "blue", mid = "grey", high = "red", midpoint = 0, limits = fit_range_loghazard) +
      labs(title = paste("Transition", trans), x = xlab, y = "Age at CKD onset")

    p_contour_hazard <- ggplot(ped_new_trans, aes(x = !!as.symbol(time), y = age_onset, z = hazard)) +
      geom_tile(aes(fill = hazard)) +
      stat_contour(aes(z = hazard), colour = "grey30", bins = 6) +
      scale_fill_gradient(name = "Hazard", low = "blue", high = "red", limits = fit_range_hazard) +
      labs(title = paste("Transition", trans), x = xlab, y = "Age at CKD onset")

    p_contour_cumu <- ggplot(ped_new_trans, aes(x = !!as.symbol(time), y = age_onset, z = cumu_hazard)) +
      geom_tile(aes(fill = cumu_hazard)) +
      stat_contour(aes(z = cumu_hazard), colour = "grey30", bins = 6) +
      scale_fill_gradient(name = "Cumulative Hazard", low = "blue", high = "red", limits = fit_range_cumu_hazard) +
      labs(title = paste("Transition", trans), x = xlab, y = "Age at CKD onset")

    p_contour_tp <- ggplot(ped_new_trans, aes(x = !!as.symbol(time), y = age_onset, z = trans_prob)) +
      geom_tile(aes(fill = trans_prob)) +
      stat_contour(aes(z = trans_prob), colour = "grey30", bins = 6) +
      scale_fill_gradient(name   = "Transition\nProbability", low = "blue", high = "red", limits = c(0,1)) +
      labs(title = paste("Transition", trans), x = xlab, y = "Age at CKD onset")

    # slice plots
    slice_df <- ped_new_trans %>% filter(age_onset %in% age_onset_slices)

    p_slice_loghazard <- ggplot(slice_df, aes(x = !!as.symbol(time), y = loghazard, group = factor(age_onset), colour = factor(age_onset), fill = factor(age_onset))) +
      geom_ribbon(aes(ymin = loghazard_lower, ymax = loghazard_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 1) +
      scale_colour_brewer(name = "Age at CKD onset", palette = "RdBu") +
      scale_fill_brewer(name = "Age at CKD onset", palette = "RdBu") +
      labs(title  = paste("Transition", trans), x  = xlab, y = "Log-hazard")

    p_slice_hazard <- ggplot(slice_df, aes(x = !!as.symbol(time), y = hazard, group = factor(age_onset), colour = factor(age_onset), fill = factor(age_onset))) +
      geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 1) +
      scale_colour_brewer(name = "Age at CKD onset", palette = "RdBu") +
      scale_fill_brewer(name = "Age at CKD onset", palette = "RdBu") +
      labs(title  = paste("Transition", trans), x  = xlab, y = "Hazard")

    p_slice_cumu <- ggplot(slice_df, aes(x = !!as.symbol(time), y = cumu_hazard, group = factor(age_onset), colour = factor(age_onset), fill = factor(age_onset))) +
      geom_ribbon(aes(ymin = cumu_lower, ymax = cumu_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 1) +
      scale_colour_brewer(name = "Age at CKD onset", palette = "RdBu") +
      scale_fill_brewer(name = "Age at CKD onset", palette = "RdBu") +
      labs(title  = paste("Transition", trans), x  = xlab, y = "Cumulative Hazard")

    p_slice_tp <- ggplot(slice_df, aes(x = !!as.symbol(time), y = trans_prob, group = factor(age_onset), colour = factor(age_onset), fill = factor(age_onset))) +
      geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 1) +
      scale_colour_brewer(name = "Age at CKD onset", palette = "RdBu") +
      scale_fill_brewer(name = "Age at CKD onset", palette = "RdBu") +
      scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by=0.2)) +
      labs(title  = paste("Transition", trans), x  = xlab, y = "Transition Probability")

  } else if(!is.null(prog_after_onset_slices) & trans %in% c("2->3", "2->4")) {

    # contour plots
    contour_df <- ped_new_trans %>% filter(age_onset %in% age_onset_slices)

    p_contour_loghazard <- ggplot(contour_df, aes(x = !!as.symbol(time), y = age_progression, z = loghazard)) +
      geom_tile(aes(fill = loghazard)) +
      stat_contour(aes(z = loghazard), colour = "grey30", bins = 6) +
      scale_fill_gradient2(name = "Log-hazard", low = "blue", mid = "grey", high = "red", midpoint = 0, limits = fit_range_loghazard) +
      facet_wrap(~ age_onset, ncol = 2, labeller = labeller(age_onset = function(x) paste0("Age at CKD onset: ", x))) +
      labs(title = paste("Transition", trans), x = xlab, y = "Age at CKD onset") +
      theme(panel.spacing = unit(1.5, "lines"))

    p_contour_hazard <- ggplot(contour_df, aes(x = !!as.symbol(time), y = age_progression, z = hazard)) +
      geom_tile(aes(fill = hazard)) +
      stat_contour(aes(z = hazard), colour = "grey30", bins = 6) +
      scale_fill_gradient(name = "Hazard", low = "blue", high = "red", limits = fit_range_hazard) +
      facet_wrap(~ age_onset, ncol = 2, labeller = labeller(age_onset = function(x) paste0("Age at CKD onset: ", x))) +
      labs(title = paste("Transition", trans), x = xlab, y = "Age at CKD onset") +
      theme(panel.spacing = unit(1.5, "lines"))

    p_contour_cumu <- ggplot(contour_df, aes(x = !!as.symbol(time), y = age_progression, z = cumu_hazard)) +
      geom_tile(aes(fill = cumu_hazard)) +
      stat_contour(aes(z = cumu_hazard), colour = "grey30", bins = 6) +
      scale_fill_gradient(name = "Cumulative Hazard", low = "blue", high = "red", limits = fit_range_cumu_hazard) +
      facet_wrap(~ age_onset, ncol = 2, labeller = labeller(age_onset = function(x) paste0("Age at CKD onset: ", x))) +
      labs(title = paste("Transition", trans), x = xlab, y = "Age at CKD onset") +
      theme(panel.spacing = unit(1.5, "lines"))

    p_contour_tp <- ggplot(contour_df, aes(x = !!as.symbol(time), y = age_progression, z = trans_prob)) +
      geom_tile(aes(fill = trans_prob)) +
      stat_contour(aes(z = trans_prob), colour = "grey30", bins = 6) +
      scale_fill_gradient(name   = "Transition\nProbability", low = "blue", high = "red", limits = c(0,1)) +
      facet_wrap(~ age_onset, ncol = 2, labeller = labeller(age_onset = function(x) paste0("Age at CKD onset: ", x))) +
      labs(title = paste("Transition", trans), x = xlab, y = "Age at CKD onset") +
      theme(panel.spacing = unit(1.5, "lines"))

    # slice plots
    slices <- expand.grid(
      age_onset           = age_onset_slices,
      prog_after_onset    = prog_after_onset_slices
    ) %>%
      mutate(age_progression = age_onset + prog_after_onset)

    slice_df <- ped_new_trans %>%
      inner_join(
        slices,
        by = c("age_onset", "age_progression")
      )

    # --- CHANGE 1: Create a dataframe to expand the x-axis limits ---
    expansion_df <- slice_df %>%
      distinct(age_onset) %>%
      mutate(!!time := age_onset)

    p_slice_loghazard <- ggplot(slice_df, aes(x = !!sym(time), y = loghazard, group = factor(prog_after_onset), colour = factor(prog_after_onset), fill = factor(prog_after_onset))) +
      geom_ribbon(aes(ymin = loghazard_lower, ymax = loghazard_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 1) +
      # --- CHANGE 2: Add an invisible geom_blank layer ---
      geom_blank(data = expansion_df, inherit.aes = FALSE, aes(x = !!sym(time))) +
      scale_colour_brewer(name = "Time until progression", palette = "RdBu") +
      scale_fill_brewer(name = "Time until progression", palette = "RdBu") +
      facet_wrap(~ age_onset, ncol = 2, labeller = labeller(age_onset = function(x) paste0("Age at CKD onset: ", x))) +
      theme(legend.position = "bottom", legend.box = "horizontal") +
      labs(title = paste("Transition", trans), x = xlab, y = "Log-hazard") +
      theme(panel.spacing = unit(1.5, "lines"))

    p_slice_hazard <- ggplot(slice_df, aes(x = !!sym(time), y = hazard, group = factor(prog_after_onset), colour = factor(prog_after_onset), fill = factor(prog_after_onset))) +
      geom_ribbon(aes(ymin = hazard_lower, ymax = hazard_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 1) +
      geom_blank(data = expansion_df, inherit.aes = FALSE, aes(x = !!sym(time))) +
      scale_colour_brewer(name = "Time until progression", palette = "RdBu") +
      scale_fill_brewer(name = "Time until progression", palette = "RdBu") +
      facet_wrap(~ age_onset, ncol = 2, labeller = labeller(age_onset = function(x) paste0("Age at CKD onset: ", x))) +
      theme(legend.position = "bottom", legend.box = "horizontal") +
      labs(title = paste("Transition", trans), x = xlab, y = "Hazard") +
      theme(panel.spacing = unit(1.5, "lines"))

    p_slice_cumu <- ggplot(slice_df, aes(x = !!sym(time), y = cumu_hazard, group = factor(prog_after_onset), colour = factor(prog_after_onset), fill = factor(prog_after_onset))) +
      geom_ribbon(aes(ymin = cumu_lower, ymax = cumu_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 1) +
      geom_blank(data = expansion_df, inherit.aes = FALSE, aes(x = !!sym(time))) +
      scale_colour_brewer(name = "Time until progression", palette = "RdBu") +
      scale_fill_brewer(name = "Time until progression", palette = "RdBu") +
      facet_wrap(~ age_onset, ncol = 2, labeller = labeller(age_onset = function(x) paste0("Age at CKD onset: ", x))) +
      theme(legend.position = "bottom", legend.box = "horizontal") +
      labs(title = paste("Transition", trans), x = xlab, y = "Cumulative Hazard") +
      theme(panel.spacing = unit(1.5, "lines"))

    p_slice_tp <- ggplot(slice_df, aes(x = !!sym(time), y = trans_prob * 100, group = factor(prog_after_onset), colour = factor(prog_after_onset), fill = factor(prog_after_onset))) +
      geom_ribbon(aes(ymin = trans_lower * 100, ymax = trans_upper * 100), alpha = 0.2, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 1) +
      geom_blank(data = expansion_df, inherit.aes = FALSE, aes(x = !!sym(time))) +
      scale_colour_brewer(name = "Time until progression", palette = "RdBu") +
      scale_fill_brewer(name = "Time until progression", palette = "RdBu") +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20), labels = function(x) paste0(x, "%")) +
      facet_wrap(~ age_onset, ncol = 2, labeller = labeller(age_onset = function(x) paste0("Age at CKD onset: ", x))) +
      theme(legend.position = "bottom", legend.box = "horizontal") +
      labs(title = paste("Transition", trans), x = xlab, y = "Transition probability (%)") +
      theme(panel.spacing = unit(1.5, "lines"))

  }

  ggsave(file.path(dir_out, "figures", model, sprintf("contour_loghazard_%s_%s_%s.png", time_scale, from_digit, to_digit)),
    plot = p_contour_loghazard, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("contour_hazard_%s_%s_%s.png", time_scale, from_digit, to_digit)),
    plot = p_contour_hazard, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("contour_cumu_%s_%s_%s.png", time_scale, from_digit, to_digit)),
    plot = p_contour_cumu, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("contour_tp_%s_%s_%s.png", time_scale, from_digit, to_digit)),
    plot = p_contour_tp, width  = 10, height = 10, units  = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("slice_loghazard_%s_%s_%s.png", time_scale, from_digit, to_digit)),
    plot = p_slice_loghazard, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("slice_hazard_%s_%s_%s.png", time_scale, from_digit, to_digit)),
    plot = p_slice_hazard, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("slice_cumu_%s_%s_%s.png", time_scale, from_digit, to_digit)),
    plot = p_slice_cumu, width = 10, height = 10, units = "cm")
  ggsave(file.path(dir_out, "figures", model, sprintf("slice_tp_%s_%s_%s.png", time_scale, from_digit, to_digit)),
    plot = p_slice_tp, width = 10, height = 10, units = "cm")

  out <- list(
    p_contour_loghazard = p_contour_loghazard,
    p_contour_hazard = p_contour_hazard,
    p_contour_cumu = p_contour_cumu,
    p_contour_tp = p_contour_tp,
    p_slice_loghazard = p_slice_loghazard,
    p_slice_hazard = p_slice_hazard,
    p_slice_cumu = p_slice_cumu,
    p_slice_tp = p_slice_tp
  )

  return(out)

}


extract_risk_factor_effects <- function(model, term_main) {

  su <- summary(model)
  vc <- vcov(model)
  all_coef_names <- names(su$p.coeff)

  if (!term_main %in% all_coef_names) {
    stop(paste("Main term '", term_main, "' not found in model coefficients.",
               "\nPlease check your model formula and the spelling of the term."))
  }

  coef_main <- su$p.coeff[term_main]
  se_main   <- su$se[term_main]
  p_main    <- su$p.pv[term_main]

  get_interaction_effect <- function(transition_level) {
    possible_names <- all_coef_names[
      grepl(term_main, all_coef_names, fixed = TRUE) &
      grepl(transition_level, all_coef_names, fixed = TRUE)
    ]
    term_int <- setdiff(possible_names, term_main)

    if (length(term_int) == 0) {
      warning(paste(
        "Interaction term for '", term_main, "' and transition '", transition_level,
        "' not found. Returning NA. Check your model formula."
      ))
      return(c(coef = NA, se = NA, p = NA))
    }
    if (length(term_int) > 1) {
       warning(paste(
         "Multiple possible interaction terms found for '", term_main,
         "' and transition '", transition_level, "'. Using the first one:", term_int[1]
       ))
    }
    term_int <- term_int[1]

    coef_combined <- coef_main + su$p.coeff[term_int]
    se_combined   <- sqrt(se_main^2 + su$se[term_int]^2 + 2 * vc[term_main, term_int])
    p_combined    <- 2 * pnorm(abs(coef_combined / se_combined), lower.tail = FALSE)

    return(c(coef = coef_combined, se = se_combined, p = p_combined))
  }

  effects_12 <- get_interaction_effect("1->2")
  effects_23 <- get_interaction_effect("2->3")

  tibble(
    term = c("0->1", "1->2", "2->3"),
    coef = c(coef_main, effects_12[1], effects_23[1]),
    se   = c(se_main,   effects_12[2], effects_23[2]),
    p    = c(p_main,    effects_12[3], effects_23[3])
  )
}


convert_to_latex_risk_factors <- function(final_table) {

  table_out <- final_table %>%
    dplyr::filter(Effect %in% c("G", "1→2", "2→3")) %>%
    dplyr::mutate(
      Coef = sprintf("%.3f", as.numeric(Coef)),
      SE = sprintf("%.3f", as.numeric(SE)),
      p_value = dplyr::if_else(round(as.numeric(p_value), 3) == 0, "$<$0.001", sprintf("%.3f", as.numeric(p_value))),
      Scenario = dplyr::case_when(
        grepl("BMI \\+ uACR \\+ Smoking \\+ eGFR", Scenario) ~ "\\shortstack[l]{G + PGS + Diabetes + BMI\\\\ + uACR + Smoking + eGFR}",
        grepl("BMI \\+ uACR \\+ Smoking", Scenario) ~ "\\shortstack[l]{G + PGS + Diabetes + BMI\\\\ + uACR + Smoking}",
        grepl("BMI \\+ uACR", Scenario) ~ "\\shortstack[l]{G + PGS + Diabetes + BMI\\\\ + uACR}",
        grepl("BMI", Scenario) ~ "\\shortstack[l]{G + PGS + Diabetes + BMI}",
        grepl("G \\+ PGS \\+ Diabetes", Scenario) ~ "\\shortstack[l]{G + PGS + Diabetes}",
        TRUE ~ Scenario
      ),
      Effect = dplyr::recode(Effect, "G" = "0→1")
    ) %>%
    dplyr::rename(
      `Risk Factors` = Scenario,
      Transition = Effect,
      Coefficient = Coef,
      `P-value` = p_value
    ) %>%
    dplyr::select(`Risk Factors`, Transition, Coefficient, SE, `P-value`)

  header <- c(
    "\\begin{table}[!ht]",
    "\\caption{\\captionukbmodelsriskfactors}",
    "\\label{tab:ukb-models-risk-factors}",
    "\\centering",
    "\\begin{tabular}[t]{lcrrr}",
    "\\toprule",
    "Risk Factors & Transition & Coefficient & SE & P-value\\\\"
  )

  body_lines <- c()
  for (i in seq(1, nrow(table_out), by = 3)) {
    is_shaded_block <- (ceiling(i / 3) %% 2 != 0)
    color_prefix <- if (is_shaded_block) "\\rowcolor{gray!10} " else ""
    risk_factor_text <- table_out$`Risk Factors`[i]

    other_cells_r1 <- paste(table_out[i, 2:ncol(table_out)], collapse = " & ")
    body_lines <- c(body_lines, paste0(color_prefix, " & ", other_cells_r1))

    other_cells_r2 <- paste(table_out[i + 1, 2:ncol(table_out)], collapse = " & ")
    body_lines <- c(body_lines, paste0(color_prefix, " & ", other_cells_r2))

    if (is_shaded_block) {
      multirow_cell <- paste0("\\multirow{-3}{*}{\\cellcolor{gray!10}", risk_factor_text, "}")
    } else {
      multirow_cell <- paste0("\\multirow{-3}{*}{", risk_factor_text, "}")
    }
    other_cells_r3 <- paste(table_out[i + 2, 2:ncol(table_out)], collapse = " & ")
    body_lines <- c(body_lines, paste0(color_prefix, multirow_cell, " & ", other_cells_r3))
  }

  body <- paste(body_lines, collapse = " \\\\\n")
  body <- paste0(body, " \\\\")

  footer <- c(
    "\\bottomrule",
    "\\end{tabular}",
    "\\parbox{\\linewidth}{\\footnotesize",
    "\\begin{itemize}[leftmargin=*, noitemsep, label=\\textcolor{white}{\\textbullet}]",
    " \\item G: Genetic variant rs77924615.",
    " \\item 0: Healthy; 1: CKD; 2: Severe CKD; 3: ESKD.",
    "\\end{itemize}",
    "}",
    "\\end{table}"
  )

  final_latex_code <- paste(c(header, "\\midrule", body, footer), collapse = "\n")

  return(final_latex_code)
}

create_summary_for_latex <- function(events_ps) {

  risk_factors <- c(
    "pgs_cross_594_umod", "diabetes", "sc_BMI",
    "sc_UACR", "smoking", "eGFRcrea"
  )
  risk_factor_labels <- c(
    "PGS", "Diabetes", "BMI",
    "uACR", "Smoking", "eGFR"
  )
  risk_factor_map <- setNames(risk_factor_labels, risk_factors)

  df_prep <- events_ps %>%
    dplyr::filter(if_all(all_of(risk_factors), ~ !is.na(.))) %>%
    dplyr::mutate(
      genotype = dplyr::case_when(
        genotype == 0 ~ "AA",
        genotype == 1 ~ "AG/GA",
        genotype == 2 ~ "GG"
      ),
      from = factor(
        dplyr::case_when(
          from == "State 0" ~ "Healthy",
          from == "State 1" ~ "CKD",
          from == "State 2" ~ "Severe CKD"
        ),
        levels = c("Healthy", "CKD", "Severe CKD")
      ),
      diabetes = as.numeric(as.character(diabetes)),
      smoking = as.numeric(as.character(smoking))
    ) %>%
    dplyr::select(id, genotype, from, all_of(risk_factors)) %>%
    dplyr::filter(!is.na(from)) %>%
    dplyr::distinct(id, from, .keep_all = TRUE)

  ns <- list(
    n_aa = prettyNum(sum(df_prep$genotype == "AA", na.rm = TRUE), big.mark = ","),
    n_ag = prettyNum(sum(df_prep$genotype == "AG/GA", na.rm = TRUE), big.mark = ","),
    n_gg = prettyNum(sum(df_prep$genotype == "GG", na.rm = TRUE), big.mark = ","),
    n_total = prettyNum(nrow(df_prep), big.mark = ",")
  )

  summarise_risk_factors <- function(df) {
    df %>%
      dplyr::summarise(
        pgs_cross_594_umod = sprintf("\\makecell{%.0f \\\\ (%.0f-%.0f)}", median(pgs_cross_594_umod, na.rm = TRUE), min(pgs_cross_594_umod, na.rm = TRUE), max(pgs_cross_594_umod, na.rm = TRUE)),
        sc_BMI = sprintf("\\makecell{%.1f \\\\ (%.1f-%.1f)}", median(sc_BMI, na.rm = TRUE), min(sc_BMI, na.rm = TRUE), max(sc_BMI, na.rm = TRUE)),
        sc_UACR = sprintf("\\makecell{%.1f \\\\ (%.1f-%.1f)}", median(sc_UACR, na.rm = TRUE), min(sc_UACR, na.rm = TRUE), max(sc_UACR, na.rm = TRUE)),
        eGFRcrea = sprintf("\\makecell{%.1f \\\\ (%.1f-%.1f)}", median(eGFRcrea, na.rm = TRUE), min(eGFRcrea, na.rm = TRUE), max(eGFRcrea, na.rm = TRUE)),
        diabetes = {
          n_positive <- sum(diabetes == 1, na.rm = TRUE)
          n_total <- sum(!is.na(diabetes))
          percentage <- ifelse(n_total > 0, (n_positive / n_total) * 100, 0)
          sprintf("%s (%s\\%%)", prettyNum(n_positive, big.mark = ","), sprintf("%.1f", percentage))
        },
        smoking = {
          n_positive <- sum(smoking == 1, na.rm = TRUE)
          n_total <- sum(!is.na(smoking))
          percentage <- ifelse(n_total > 0, (n_positive / n_total) * 100, 0)
          sprintf("%s (%s\\%%)", prettyNum(n_positive, big.mark = ","), sprintf("%.1f", percentage))
        },
        .groups = "drop"
      )
  }

  summary_by_genotype <- df_prep %>% dplyr::group_by(genotype, from) %>% summarise_risk_factors()
  summary_total <- df_prep %>% dplyr::group_by(from) %>% summarise_risk_factors() %>% dplyr::mutate(genotype = "Total")
  combined_summary <- dplyr::bind_rows(summary_by_genotype, summary_total)

  final_table <- combined_summary %>%
    tidyr::pivot_longer(cols = all_of(risk_factors), names_to = "Risk Factor", values_to = "Summary_Value") %>%
    tidyr::pivot_wider(names_from = "genotype", values_from = "Summary_Value") %>%
    dplyr::select(`Risk Factor`, State = from, AA, `AG/GA`, GG, Total) %>%
    dplyr::mutate(`Risk Factor` = dplyr::recode(`Risk Factor`, !!!risk_factor_map)) %>%
    dplyr::mutate(`Risk Factor` = factor(`Risk Factor`, levels = risk_factor_labels)) %>%
    dplyr::arrange(`Risk Factor`, State)

  n_table <- df_prep %>%
    dplyr::count(from, genotype, name = "count") %>%
    tidyr::pivot_wider(names_from = "genotype", values_from = "count") %>%
    dplyr::mutate(Total = AA + `AG/GA` + GG) %>%
    dplyr::mutate(dplyr::across(-from, ~prettyNum(.x, big.mark = ",")))

  return(list(ns = ns, table = final_table, n_table = n_table))
}

convert_to_latex_risk_factor_distributions <- function(distributions) {

  ns <- distributions$ns
  final_table <- distributions$table
  n_table <- distributions$n_table

  header_lines <- c(
    "\\begin{table}[!ht]",
    "\\footnotesize",
    "\\caption{\\captionukbriskfactordistributions}",
    "\\label{tab:ukb-risk-factor-distributions}",
    "\\centering",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}[t]{llcccc}",
    "\\toprule",
    "& & \\multicolumn{3}{c}{Genotype of rs77924615} & \\\\",
    "\\cmidrule(lr){3-5}",
    sprintf(
      "Risk Factor & State & AA & AG/GA & GG & Total \\\\ & & (n=%s) & (n=%s) & (n=%s) & (n=%s) \\\\",
      ns$n_aa, ns$n_ag, ns$n_gg, ns$n_total
    )
  )

  body_lines <- c()
  for (i in seq(1, nrow(final_table), by = 3)) {

    is_shaded_block <- (ceiling(i / 3) %% 2 != 0)
    risk_factor_text <- final_table$`Risk Factor`[i]

    # Helper function to apply cell coloring if the block is shaded
    apply_color <- function(cell_content) {
      if (is_shaded_block) {
        return(paste0("\\cellcolor{gray!10}", cell_content))
      }
      return(cell_content)
    }

    # Row 1
    cells_r1 <- c(apply_color(""), apply_color("0"), sapply(final_table[i, 3:ncol(final_table)], apply_color))
    body_lines <- c(body_lines, paste(cells_r1, collapse = " & "))

    # Row 2
    cells_r2 <- c(apply_color(""), apply_color("1"), sapply(final_table[i + 1, 3:ncol(final_table)], apply_color))
    body_lines <- c(body_lines, paste(cells_r2, collapse = " & "))

    # Row 3 (with multirow)
    if (is_shaded_block) {
      multirow_cell <- paste0("\\multirow{-3}{*}{\\cellcolor{gray!10}", risk_factor_text, "}")
    } else {
      multirow_cell <- paste0("\\multirow{-3}{*}{", risk_factor_text, "}")
    }
    cells_r3 <- c(multirow_cell, apply_color("2"), sapply(final_table[i + 2, 3:ncol(final_table)], apply_color))
    body_lines <- c(body_lines, paste(cells_r3, collapse = " & "))
  }

  body_lines <- c(body_lines, "\\midrule\\midrule")

  # 'n' block does not get shaded
  n_cells_r1 <- paste(n_table[1, 2:ncol(n_table)], collapse = " & ")
  body_lines <- c(body_lines, paste(" & 0 &", n_cells_r1))

  n_cells_r2 <- paste(n_table[2, 2:ncol(n_table)], collapse = " & ")
  body_lines <- c(body_lines, paste(" & 1 &", n_cells_r2))

  multirow_n_cell <- "\\multirow{-3}{*}{n}"
  n_cells_r3 <- paste(n_table[3, 2:ncol(n_table)], collapse = " & ")
  body_lines <- c(body_lines, paste(multirow_n_cell, " & 2 &", n_cells_r3))

  body <- paste(body_lines, collapse = " \\\\\n")
  body <- paste0(body, " \\\\")

  footer_lines <- c(
    "\\bottomrule",
    "\\end{tabular}",
    "\\normalsize",
    "\\end{table}"
  )

  final_latex_code <- paste(c(header_lines, "\\midrule", body, footer_lines), collapse = "\n")

  return(final_latex_code)
}
