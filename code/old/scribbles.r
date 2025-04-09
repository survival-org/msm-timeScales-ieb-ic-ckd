




## ESKD ----
eskd <- pheno %>% filter(!is.na(event_dt_eskd)) %>% select(id, event_dt=event_dt_eskd) %>% unique()

## CKD ----
setDT(pheno)
setorder(pheno, id, event_dt)
pheno[, event_dt := as.Date(event_dt)]
pheno[, eGFRcrea_less60 := eGFRcrea < 60]
ckd <- pheno[, {
  N <- .N
  spell_starts <- integer()
  spell_ends <- integer()
  i <- 2
  while (i <= N) {
    if (eGFRcrea_less60[i-1] == TRUE && eGFRcrea_less60[i] == TRUE &&
        (i == 2 || eGFRcrea_less60[i-2] == FALSE)) {
      # Start of CKD spell at position i-1
      start_pos <- i - 1
      # Spell continues until eGFRcrea >= 60
      end_pos <- i
      while (end_pos + 1 <= N && eGFRcrea_less60[end_pos + 1] == TRUE) {
        end_pos <- end_pos + 1
      }
      # Record spell starts and ends
      spell_starts <- c(spell_starts, start_pos)
      spell_ends <- c(spell_ends, end_pos)
      # Move to the next position after the spell ends
      i <- end_pos + 1
    } else {
      i <- i + 1
    }
  }
  # Output the CKD spells for this individual without the explicit 'id' column
  if (length(spell_starts) > 0) {
    data.table(
      event_dt_start = event_dt[spell_starts],
      event_dt_end = event_dt[spell_ends]
    )
  } else {
    NULL
  }
}, by = id]

## severe CKD ----
setDT(pheno)
setorder(pheno, id, event_dt)
pheno[, event_dt := as.Date(event_dt)]
pheno[, eGFRcrea_less30 := eGFRcrea < 30]
severe_ckd <- pheno[, {
  N <- .N
  spell_starts <- integer()
  spell_ends <- integer()
  i <- 1
  while (i <= N) {
    if (eGFRcrea_less30[i] == TRUE) {
      # Start of severe CKD spell at position i
      start_pos <- i
      # Spell continues until eGFRcrea >= 30 or end of data
      end_pos <- i
      while (end_pos + 1 <= N && eGFRcrea_less30[end_pos + 1] == TRUE) {
        end_pos <- end_pos + 1
      }
      # Record spell starts and ends
      spell_starts <- c(spell_starts, start_pos)
      spell_ends <- c(spell_ends, end_pos)
      # Move to the next position after the spell ends
      i <- end_pos + 1
    } else {
      i <- i + 1
    }
  }
  # Output the severe CKD spells for this individual
  if (length(spell_starts) > 0) {
    data.table(
      event_dt_start = event_dt[spell_starts],
      event_dt_end = event_dt[spell_ends]
    )
  } else {
    NULL
  }
}, by = id]

## create status ----
idx <- sort(unique(pheno$id))
for(i in idx) {
    ckd_times <- ckd[id == i, .(t = event_dt_start, status = 1)]
    severe_ckd_times <- severe_ckd[id == i, .(t = event_dt_start, status = 2)]
    eskd_time <- eskd[id == i, .(t = event_dt, status = 3)]
    max_eGFR_time <- pheno[id == i, .(t = max(.SD$t))]
    all_times <- rbind(ckd_times, severe_ckd_times, eskd_time)
}


setDT(pheno)
setDT(ckd)
setDT(severe_ckd)
setDT(eskd)

# Ensure event dates are of Date type
pheno[, event_dt := as.Date(event_dt)]
ckd[, `:=`(event_dt_start = as.Date(event_dt_start), event_dt_end = as.Date(event_dt_end))]
severe_ckd[, `:=`(event_dt_start = as.Date(event_dt_start), event_dt_end = as.Date(event_dt_end))]
eskd[, event_dt := as.Date(event_dt)]

# Compute baseline date for each individual (earliest event_dt in pheno)
baseline_dates <- pheno[, .(baseline_date = min(event_dt)), by = id]

# Merge baseline_dates with datasets and compute time since baseline (t)
ckd <- merge(ckd, baseline_dates, by = "id", all.x = TRUE)
ckd[, `:=`(
  start_t = as.numeric(event_dt_start - baseline_date) / 365.25,
  end_t = as.numeric(event_dt_end - baseline_date) / 365.25
)]

severe_ckd <- merge(severe_ckd, baseline_dates, by = "id", all.x = TRUE)
severe_ckd[, `:=`(
  start_t = as.numeric(event_dt_start - baseline_date) / 365.25,
  end_t = as.numeric(event_dt_end - baseline_date) / 365.25
)]

eskd <- merge(eskd, baseline_dates, by = "id", all.x = TRUE)
eskd[, t := as.numeric(event_dt - baseline_date) / 365.25]

# Compute time since baseline (t) in pheno
pheno <- merge(pheno, baseline_dates, by = "id", all.x = TRUE)
pheno[, t := as.numeric(event_dt - baseline_date) / 365.25]

# Compute max_t (end of observation) for each individual
max_t <- pheno[, .(t = max(t)), by = id]

# Collect all relevant time points for each individual
time_points_list <- list(
  # Baseline (t = 0)
  pheno[, .(id, t = 0)],
  # CKD spell start and end times
  ckd[, .(id, t = start_t)],
  ckd[, .(id, t = end_t)],
  # Severe CKD spell start and end times
  severe_ckd[, .(id, t = start_t)],
  severe_ckd[, .(id, t = end_t)],
  # ESKD event times
  eskd[, .(id, t)],
  # End of observation (max_t)
  max_t
)

# Combine and get unique time points
all_times <- rbindlist(time_points_list, use.names = TRUE)
all_times <- all_times[order(id, t)]
all_times <- all_times[, .(t = unique(t)), by = id]

# Create intervals by shifting time points
all_times[, `:=`(
  start = t,
  stop = shift(t, type = "lead")
), by = id]

# Remove the last interval where stop is NA
all_times <- all_times[!is.na(stop)]

assign_status <- function(dt) {
  id_val <- dt$id[1]  # Assuming this never returns an empty vector

  # Get individual's spells and events
  eskd_t <- eskd_t <- eskd[id == id_val, .(t)]
  severe_spells <- severe_spells <- severe_ckd[id == id_val, .(start_t, end_t)]
  ckd_spells <- ckd[id == id_val, .(start_t, end_t)]

  # Initialize status vector
  statuses <- integer(nrow(dt))

  for (i in seq_len(nrow(dt))) {
    s <- dt$start[i]
    e <- dt$stop[i]

    # Initialize status
    status <- 0

    # Status 3: ESKD event (highest priority)
    if (!is.null(eskd_t) && any(eskd_t == e)) {
        status <- 3
    # Status 2: Severe CKD spell
    } else if (any(severe_spells$start_t <= e & severe_spells$end_t >= s)) {
        status <- 2
    # Status 1: CKD spell
    } else if (any(ckd_spells$start_t <= e & ckd_spells$end_t >= s)) {
        status <- 1
    # Status 0: Default
    } else {
        status <- 0
    }

    statuses[i] <- status
  }

  return(statuses)
}



# Apply the function and assign statuses
all_times[, status := assign_status(.SD), by = id]


# Determine transitions
all_times[, next_status := shift(status, type = "lead"), by = id]
all_times[, transition := paste(status, "->", next_status)]
all_times[is.na(next_status), transition := "censoring"]

# Remove intervals where status does not change
all_times <- all_times[status != next_status | is.na(next_status)]

# Assign episodes
all_times[, episode := seq_len(.N), by = .(id, transition)]

# Assemble the eventData dataset
eventData <- all_times[, .(id, start, stop, status, transition, episode)]

# View the resulting dataset
print(eventData)





a <- pheno %>% filter(!is.na(event_dt_eskd)) %>% pull(id) %>% unique
length(a)
table(events$transition)
b <- events %>% filter(transition %in% c("0 -> 3", "1 -> 3", "2 -> 3")) %>% pull(id) %>% unique
length(b)
diff <- setdiff(a,b)
