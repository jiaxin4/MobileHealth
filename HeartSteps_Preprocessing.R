
rm(list = ls())
jbslot <- read.csv("https://raw.githubusercontent.com/klasnja/HeartStepsV1/main/data_files/jbsteps.csv")
gfslot <- read.csv("https://raw.githubusercontent.com/klasnja/HeartStepsV1/main/data_files/gfsteps.csv")
users <- read.csv("https://raw.githubusercontent.com/klasnja/HeartStepsV1/main/data_files/users.csv")
suggest <- read.csv("https://raw.githubusercontent.com/klasnja/HeartStepsV1/main/data_files/suggestions.csv")

library(tidyverse)
library(xtable)

# minutes of interest
total_T <- 60

CREATE_PLOTS <- FALSE # if TRUE, will make sanity check plots
PRINT_SANITY_CHECK <- TRUE # if TRUE, will print sanity check results

jbslot <- as_tibble(jbslot)

suggest <- as_tibble(suggest)

# remove the baseline covariates (decision.index = NA) -------------------------------------------
jbslot <- jbslot %>% filter(!is.na(decision.index))
# the publicly available dataset has baseline information recorded (i.e. at decision.index = NA)
# we will remove these to keep the form of the datasets to be consistent with our analysis.


# Remove travel dates in suggest
var_from_suggest <- c("user.index", "decision.index.nogap",  "decision.index", "avail", "send",
                          "send.active", "send.sedentary", "sugg.decision.utime", 
                          "dec.location.category", "jbsteps30pre", "jbsteps30pre.zero")
suggest_selectedvar <- suggest[, var_from_suggest]

suggest_selectedvar <- suggest_selectedvar %>% 
  mutate(sugg.decision.udate = date(ymd_hms(suggest_selectedvar$sugg.decision.utime)))
tmp_df <- data.frame()

for (i in 1:nrow(users)) {
  
  user_i_dta <- filter(suggest_selectedvar, user.index == users$user.index[i])
  if (users$travel.start[i] != "") {
    nrow_before <- nrow(user_i_dta)
    user_i_dta <- filter(user_i_dta, (sugg.decision.udate < users$travel.start[i]) | (sugg.decision.udate > users$travel.end[i]))
    nrow_after <- nrow(user_i_dta)
    cat(paste0("user ", users$user.index[i], ": travel.start ", users$travel.start[i], ", travel.end ", users$travel.end[i]), "\n")
    cat(paste0("        nrow_before ", nrow_before, ", nrow_after ", nrow_after, ", nrow_removed ",
               nrow_before - nrow_after, "\n"))
  }
  tmp_df <- rbind(tmp_df, user_i_dta)
}
suggest_tdr <- tmp_df

# user 1: travel.start 2015-08-12, travel.end 2015-08-31 
# nrow_before 278, nrow_after 178, nrow_removed 100
# user 3: travel.start 2015-08-13, travel.end 2015-08-20 
# nrow_before 255, nrow_after 215, nrow_removed 40
# user 6: travel.start 2015-08-10, travel.end 2015-08-15 
# nrow_before 212, nrow_after 182, nrow_removed 30
# user 13: travel.start 2015-09-01, travel.end 2015-09-05 
# nrow_before 215, nrow_after 192, nrow_removed 23
# user 14: travel.start 2015-10-12, travel.end 2015-10-22 
# nrow_before 265, nrow_after 210, nrow_removed 55
# user 16: travel.start 2015-09-20, travel.end 2015-09-22 
# nrow_before 214, nrow_after 199, nrow_removed 15
# user 31: travel.start 2015-12-15, travel.end 2016-01-08 
# nrow_before 309, nrow_after 184, nrow_removed 125

# Construct a new decision.index.nogap variable in suggest dataset, then paste it into jbslot -------------------------------------------
## Issue: Why do we need to do this? b/c somehow there are observations that have decision.index.nogap skipped
## e.g. user 13's orginal decision.index.nogap has labeled wrong at j = 23, 24 (they got skipped)
## Solution: construct a new decision.index.nogap variable

suggest_tdr <- suggest_tdr %>% 
  group_by(user.index) %>% 
  mutate(max.decision.index.nogap = n(),
         decision.index.nogap.new = 0:(unique(max.decision.index.nogap) - 1))
a <- which(suggest_tdr$decision.index.nogap != suggest_tdr$decision.index.nogap.new)
b <- unique(suggest_tdr[a, ]) 

# Paste the newly constructed decision.index.nogap.new into jbslot
jbslot <- suggest_tdr %>% 
  dplyr::select("user.index", "decision.index", "decision.index.nogap.new", 
         "sugg.decision.utime", "sugg.decision.udate") %>% 
  right_join(jbslot, by = c("user.index", "decision.index"))


# remove travel dates in jbslot-------------------------------------------

## Issue: At i  = 1, j = 203, the suggest dataset has notification sent when they are in travel dates. 
## But the steps are recorded. And the steps recorded time is not in travel dates (the morning they got back). 
## So there are NAs when joining these two datasets. We can remove these observations since their notifications are sent at travel dates.
## Solution: we use the sugg.decision.udate from suggest dataset (to filter the observations in travel dates) to avoid discrepancy between suggest and jbslot.

# travel day is stored in "users" data.frame, with variable names "travel.start" and "travel.end"
tmp_df <- data.frame()
for (i in 1:nrow(users)) {
  
  user_i_dta <- filter(jbslot, user.index == users$user.index[i])
  if (users$travel.start[i] != "") {
    nrow_before <- nrow(user_i_dta)
    user_i_dta <- filter(user_i_dta, (sugg.decision.udate < users$travel.start[i]) | (sugg.decision.udate > users$travel.end[i]))
    nrow_after <- nrow(user_i_dta)
    cat(paste0("user ", users$user.index[i], ": travel.start ", users$travel.start[i], ", travel.end ", users$travel.end[i]), "\n")
    cat(paste0("        nrow_before ", nrow_before, ", nrow_after ", nrow_after, ", nrow_removed ",
               nrow_before - nrow_after, "\n"))
  }
  tmp_df <- rbind(tmp_df, user_i_dta)
}
jbslot_tdr <- tmp_df

# user 1: travel.start 2015-08-12, travel.end 2015-08-31 
# nrow_before 7965, nrow_after 7542, nrow_removed 423
# user 3: travel.start 2015-08-13, travel.end 2015-08-20 
# nrow_before 4611, nrow_after 4441, nrow_removed 170
# user 6: travel.start 2015-08-10, travel.end 2015-08-15 
# nrow_before 7364, nrow_after 7094, nrow_removed 270
# user 13: travel.start 2015-09-01, travel.end 2015-09-05 
# nrow_before 3389, nrow_after 2866, nrow_removed 523
# user 14: travel.start 2015-10-12, travel.end 2015-10-22 
# nrow_before 9199, nrow_after 7828, nrow_removed 1371
# user 16: travel.start 2015-09-20, travel.end 2015-09-22 
# nrow_before 7041, nrow_after 6496, nrow_removed 545
# user 31: travel.start 2015-12-15, travel.end 2016-01-08 
# nrow_before 5357, nrow_after 5108, nrow_removed 249


#  construct min.after.decision variable in the dataset -------------------------------------------
jbslot_tdr <- jbslot_tdr %>% 
  mutate(steps.utime = ymd_hms(steps.utime),
         sugg.decision.utime = ymd_hms(sugg.decision.utime),
         min.after.decision = floor(as.numeric(difftime(steps.utime, sugg.decision.utime, units = 'mins'))) + 1
         )

jbslot_tdr <- jbslot_tdr[order(jbslot_tdr$user.index, 
                               jbslot_tdr$decision.index.nogap.new, 
                               jbslot_tdr$min.after.decision), ]


if (PRINT_SANITY_CHECK) {
  # make sure primary keys are unique
  jbslot_tdr %>% 
    count(user.index, decision.index.nogap.new, min.after.decision) %>% 
    filter(n > 1) # no duplicate!
}

if (PRINT_SANITY_CHECK) {
  # make sure all user has decision.index.nogap from 0 to max consecutively, no gap in between
  print(jbslot_tdr %>% group_by(user.index) %>%
          summarise(n_dp = length(decision.index.nogap.new),
                    dp_min = min(decision.index.nogap.new),
                    dp_max = max(decision.index.nogap.new),
                    no_gap_in_dp = ((dp_max - dp_min + 1) == n_dp)), n = Inf)
}

#Note: User.index = 4 has no observation recorded at decision.index 0, because his/her connection failed.
# # A tibble: 37 × 5
#   user.index  n_dp dp_min dp_max no_gap_in_dp
# <int> <int>  <int>  <int> <lgl>       
# 1          1  7542      0    169 FALSE       
# 2          2  5020      0    207 FALSE       
# 3          3  4441      0    213 FALSE       
# 4          4 11715      1    217 FALSE       
# 5          5  4202      0    215 FALSE       
# 6          6  7094      0    181 FALSE       
# 7          7  6841      0    214 FALSE       
# 8          8  7898      0    219 FALSE       
# 9          9  5296      0    205 FALSE       
# 10         10  6937      0    213 FALSE       
# 11         11  6488      0    216 FALSE       
# 12         12  7744      0    242 FALSE       
# 13         13  2866      0    184 FALSE       
# 14         14  7828      0    208 FALSE       
# 15         15  7426      0    254 FALSE       
# 16         16  6496      0    197 FALSE       
# 17         17  7618      0    208 FALSE       
# 18         18  5226      0    210 FALSE       
# 19         19  6168      0    208 FALSE       
# 20         20  4779      0    233 FALSE       
# 21         21  8558      0    233 FALSE       
# 22         22  4733      0    213 FALSE       
# 23         23  5198      0    223 FALSE       
# 24         24  6215      0    173 FALSE       
# 25         25 10346      0    222 FALSE       
# 26         26  6289      0    211 FALSE       
# 27         27  3897      0    210 FALSE       
# 28         28  8327      0    208 FALSE       
# 29         29  2575      0     83 FALSE       
# 30         30  7494      0    208 FALSE       
# 31         31  5108      0    182 FALSE       
# 32         32  3922      0    223 FALSE       
# 33         33  7813      0    226 FALSE       
# 34         34  3039      0    211 FALSE       
# 35         35  5797      0    208 FALSE       
# 36         36  9155      0    201 FALSE       
# 37         37  5451      0    228 FALSE       

# Construct new template data set with "user", "decision.index.nogap", "min.after.decision" -----------------------------------
## make a combined keyword, then separate out (because expand.grid only accepts vector arguments)
user.di_unique <- unite(suggest_tdr, "user.di", "user.index", "decision.index.nogap.new")$user.di

dta_template <- expand.grid(user.di = user.di_unique, min.after.decision = 0:(total_T - 1)) %>%
  separate(user.di, into = c("user", "decision.index.nogap")) %>%
  arrange(user, decision.index.nogap, min.after.decision)
dta_template$user <- as.numeric(dta_template$user)
dta_template$decision.index.nogap <- as.numeric(dta_template$decision.index.nogap)
dta_template <- as_tibble(dta_template)

## sanity checks on the number of decision points and the number of minutes
if (PRINT_SANITY_CHECK) {
  summary(dta_template)
  print(dta_template %>% group_by(user) %>% summarise(n.dp = length(decision.index.nogap) / total_T), n = Inf)
  summary(dta_template %>% group_by(user, decision.index.nogap) %>%
            summarise(n.minute = length(min.after.decision),
                      minute.max = max(min.after.decision),
                      minute.min = min(min.after.decision),
                      n.unique.minute = length(unique(min.after.decision))))
}

# # A tibble: 37 × 2
# user  n.dp
# <dbl> <dbl>
#   1     1   178
# 2     2   209
# 3     3   215
# 4     4   219
# 5     5   217
# 6     6   182
# 7     7   216
# 8     8   221
# 9     9   207
# 10    10   215
# 11    11   220
# 12    12   244
# 13    13   192
# 14    14   210
# 15    15   256
# 16    16   199
# 17    17   210
# 18    18   211
# 19    19   210
# 20    20   235
# 21    21   235
# 22    22   214
# 23    23   225
# 24    24   175
# 25    25   223
# 26    26   212
# 27    27   211
# 28    28   210
# 29    29   211
# 30    30   210
# 31    31   184
# 32    32   225
# 33    33   228
# 34    34   212
# 35    35   213
# 36    36   202
# 37    37   230

jbslot_tdr_selectedvar <- jbslot_tdr %>% 
  dplyr::select(c("user.index", "decision.index.nogap.new", 
                  "steps", "steps.utime", "study.day.nogap", "min.after.decision"))
suggest_tdr_selectedvar <- suggest_tdr %>% 
  dplyr::select(c("user.index", "avail",
                  "send", "send.active", "send.sedentary", "dec.location.category", 
                  "sugg.decision.utime", "sugg.decision.udate", 
                  "max.decision.index.nogap", "decision.index.nogap.new",
                  "jbsteps30pre", "jbsteps30pre.zero"))
  
jbslot60 <- dta_template %>%
  left_join(suggest_tdr_selectedvar, by = c("user" = "user.index",
                                "decision.index.nogap" = "decision.index.nogap.new")) %>% 
  left_join(jbslot_tdr_selectedvar, by = c("user" = "user.index",  
                                       "decision.index.nogap" = "decision.index.nogap.new",
                                       "min.after.decision" = "min.after.decision")) 

# Recoding some variables to be used in our analysis
jbslot60$steps[is.na(jbslot60$steps)] <- 0
jbslot60$steps.log <- log(jbslot60$steps + 0.5)
jbslot60$jbsteps30pre.log <- log(jbslot60$jbsteps30pre.zero + 0.5)

jbslot60$location.other <- !(jbslot60$dec.location.category %in% c("home", "work"))
jbslot60 <- jbslot60 %>% 
  mutate(send.active = case_when(send.active == "True" ~ TRUE,
                                 send.active == "False"|send.active == "" ~ FALSE),
         send.sedentary = case_when(send.sedentary == "True" |send.sedentary == "" ~ TRUE,
                                    send.sedentary == "False" ~ FALSE),
         avail = case_when(avail == "True" ~ TRUE,
                           avail == "False" ~ FALSE),
         decision.utime = ymd_hms(sugg.decision.utime),
         send = case_when(send.active == FALSE & send.sedentary == FALSE ~ FALSE,
                           send.active == TRUE | send.sedentary == TRUE ~ TRUE)
         )
saveRDS(jbslot60, file = "jbslot_public_60min.RDS")

