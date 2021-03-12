library(tidyverse)
library(assertthat)

# choose stochastically from a vector of probabilities
p_chooser <- function(probs) {
  p <- runif(1, 0, 1)
  upper_probs <- cumsum(probs)
  lower_probs <- lag(upper_probs, default = 0)
  which(map2_lgl(.x = lower_probs, .y = upper_probs, ~between(p, .x, .y)))
}

# transform estimates according to softmax
softmax <- function(estimates, temperature) {
  exps <- exp(estimates / temperature)
  exps / sum(exps)
}

# select using softmax
choose_softmax <- function(estimates, temperature) {
  probs <- softmax(estimates, temperature = temperature)
  p_chooser(probs)
}

# draw a single reward from binomial distribution (1 = rewarded, 0 = unrewarded) with a given probability
draw_reward_from_prob <- function(prob) {
  rbinom(1, 1, prob)
}

# initialize the memory of each agent and update the reward probabilities
# fill out the given data table with 'memories' column in initialized blank state
sim_initialize <- function(data_tbl) {
  
  if (data_tbl %has_name% "memories") print("Column 'memories' was overwritten")
  data_tbl %>%
  mutate(memories = pmap(list(prob, l_rate, temperature),
                         ~list(loc = 1:12, mask = mask, probs = mask * ..1, estimates = mask * 0,
                               reward = 0, loc_visited = 0, n_visited = mask * 0,
                               l_rate = ..2, temperature = ..3)))
}

# choose where to go, evaluate and update memory
update_memo <- function(memories, old) { 
  # memories is a list consisting of:
  # loc(ation), mask (rewarding status),
  # prob(abilitie)s, estimates (best to be initialized at 0), 
  # reward (status of last visit 1/0), loc_visited (last visited location),
  # n_visited (vector of number of visits per location)
  # l_rate (learning rate), temperature (exploration parameter in softmax)
  loc <- choose_softmax(memories$estimates, memories$temperature) # where to go
  prob_vec <- case_when( # determine the current state of rewarding or non-rewarding locations
    over_target_visits(memories$n_visited, rep(1, 12), target_visits) == 0 ~ memories$probs, 
    over_target_visits(memories$n_visited, rep(1, 12), target_visits) == 1 ~ memories$probs * 0,
    over_target_visits(memories$n_visited, rep(1, 12), target_visits) == 2 ~ mask2 * max(memories$probs)
  )
  reward <- draw_reward_from_prob(prob_vec[loc])
  memories$reward <- reward # update memo
  # update estimate using gamma (l_rate) updating
  memories$estimates[loc] <- (1 - memories$l_rate) * memories$estimates[loc] + memories$l_rate * reward
  # tally visits
  memories$n_visited[loc] <- memories$n_visited[loc] + 1
  # remember where last visit happened
  memories$loc_visited <- loc
  memories
}


# check if the sum of visits for a given subset of options (mask == 1)  has reached the target number
# the target is a vector of maximal visit numbers in an ascending order, e.g. c(50, 100, 200)
# returns the number of thresholds that has been surpassed (0 if none, 1 after the first threshold, etc.)
over_target_visits <- function(vis_vec, mask, target) {
  sum(sum(vis_vec[mask == 1]) >= target)
}
# over_target_visits(c(2, 2, 4), c(1, 1, 1), c(3, 5, 7))

# expand table with necessary number of animals and decisions per animal
exp_tbl <- function(data_tbl, n_pokes = 100, n_inds = 100) {
  obligatory_columns <- c("cond", "experiment")
  assert_that(data_tbl %has_name% obligatory_columns,
              msg = glue::glue("The data table must have columns called '{obligatory_columns[1]}' and '{obligatory_columns[2]}'!"))
  
  # helper function to make function vectorized
  exp_hlp <- function(data_tbl, n) {
    data_tbl %>%
      slice(rep(1:n(), each = n))
  }
  
  data_tbl %>%
    exp_hlp(n_pokes) %>%
    group_by(cond, experiment) %>%
    mutate(n_poke = 1:n()) %>%
    exp_hlp(n_inds) %>%
    group_by(cond, experiment, n_poke) %>%
    mutate(ind = 1:n()) %>%
    arrange(ind, experiment, cond, n_poke) %>%
    ungroup()
}

# label figures more easily and, clearly
learning_label <- function(string) {
  paste0("learning rate: ", string)
}
temperature_label <- function(string) {
  paste0("exploration: ", string)
}
prob_label <- function(string) {
  paste0("probability: ", string)
}

global_labeller <- labeller(
  l_rate = learning_label,
  temperature = temperature_label,
  prob = prob_label,
  .default = label_both
)

######## test of plasticity of sampling

# rewarding flowers
mask <- c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0) # mask showing the rewarding and non-rewarding options
mask2 <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0) # the new mask after the reversal
rewarding <- tibble(loc = 1:12, mask = mask, mask2 = mask2) # table of rewarding status information
target_visits <- c(10000) # a high number ensures this experiment has only one condition and no (mid-session) reversal

plast_test <- cross_df(list(prob = c(0.83, 0.5, 0.3), experiment = 1,
                          l_rate = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
                          temperature = c(0.0125, 0.025, 0.05, 0.1))) %>%
  mutate(cond = paste0(prob, l_rate, temperature)) %>% # created cond column
  # slice(1) %>%
  exp_tbl(n_pokes = 501) %>% # expand table with 100 animals and 500 choices per animal 
  # (501 because the first one is the initial state of the memory)
  # for each probability condition, learning rate and exploration (temperature) initialize memories
  # initiate column with blank memories, first row is before any visit happened; it is later removed
  # mutate(memories = pmap(list(prob, l_rate, temperature),
  #                        ~list(loc = 1:12, mask = mask, probs = mask * ..1, estimates = mask * 0,
  #                              reward = 0, loc_visited = 0, n_visited = mask * 0,
  #                              l_rate = ..2, temperature = ..3))) %>%
  sim_initialize() %>%
  group_by(ind, prob, experiment, l_rate, temperature) %>%
  mutate(memories = accumulate(memories, .f = update_memo)) # run the simulations

plast_res <- plast_test %>%
  filter(n_poke > 1) %>%
  mutate(loc = map_dbl(memories, ~pluck(., "loc_visited")),
         rewarded = map_dbl(memories, ~pluck(., "reward")),
         estimate = map2_dbl(memories, loc, ~pluck(.x, "estimates") %>% pluck(.y))) %>%
  left_join(rewarding)

folder <- "localoutput/"
write.table(plast_res %>% select(-memories), file = paste0(folder, "PlasticitySimulation.csv"), sep = ";", row.names = FALSE)

# plast_res <- read.csv2(paste0(folder, "PlasticitySimulation.csv", dec = "."))

# all.equal(plast_res %>% select(-memories), plast_res2)
samp_sim <- plast_res %>%
  filter(!is.na(estimate), n_poke > 100) %>%
  group_by(prob, ind, l_rate, temperature) %>%
  summarise(correct = mean(mask, na.rm = TRUE),
            sampling = 1 - correct,
            n_flowers = n_distinct(loc))

write.table(samp_sim, file = paste0(folder, "PlasticitySimSummary.csv"), sep = ";", row.names = FALSE)


samp_sim <- read.csv2(paste0(folder, "PlasticitySimSummary.csv"), dec = ".") %>% 
  mutate(ind = factor(ind))

samp_sim %>%
  ggplot(aes(prob, sampling, color = as.factor(temperature))) +
  geom_jitter(alpha = 0.3) +
  stat_summary(size = 0.7) +
  theme_bw() +
  facet_grid(. ~ l_rate, labeller = global_labeller) +
  labs(title = "Effect of learning rate and exploration on sampling at different probabilities",
       x = "probability", color = "exploration")

slopes_sim <- samp_sim %>%
  ungroup() %>%
  nest_by(ind, l_rate, temperature) %>%
  mutate(slopes = glm(sampling ~ prob, family = binomial, data = data) %>%
           coef() %>%  pluck(2),
         lin_slopes = lm(sampling ~ prob, data = data) %>%
           coef() %>%  pluck(2)) %>%
  unnest()

slopes_sim %>%
  ggplot(aes(lin_slopes, slopes)) +
  geom_point() +
  stat_smooth()


# take differences rather than slopes
samp_sim %>%
  ggplot(aes(prob, sampling, color = as.factor(temperature), group = temperature)) +
  geom_point() +
  stat_smooth(method = lm) +
  facet_grid(. ~ l_rate)

diffs <- samp_sim %>%
  ungroup() %>% 
  filter(prob != 0.5) %>%
  select(-correct) %>%
  pivot_wider(id = c(ind, l_rate, temperature), names_from = prob,
              values_from = sampling, names_prefix = "prob") %>%
  mutate(diff = prob0.83 - prob0.3) %>%
  select(ind, l_rate, temperature, diff)

# diff is basically the same as lin_slope
slopes_sim %>%
  left_join(diffs) %>%
  ggplot(aes(diff, lin_slopes, color = as.factor(temperature))) +
  geom_point() +
  facet_grid(. ~ temperature)

slopes_sim %>%
  ungroup() %>%
  nest_by(ind, l_rate, temperature, lin_slopes) %>%
  select(-data) %>%
  ggplot(aes(l_rate, -lin_slopes, color = as.factor(temperature))) +
  geom_point() +
  facet_grid(. ~ temperature)

# flexibility test as with bats in real experiment
# rewarding flowers

target_visits <- c(100, 101) # if a middle phase needed, use c(100, 150)

flex_test <- cross_df(list(prob = c(0.83, 0.5, 0.3), experiment = 1,
                           l_rate = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
                           temperature = c(0.0125, 0.025, 0.05, 0.1))) %>%
  mutate(cond = paste0(prob, l_rate, temperature)) %>% # created cond column
  exp_tbl(n_pokes = 301) %>% # expand table with 100 animals and 300 choices per animal
  # initiate column with blank memories, first row is before any visit happened; it is later removed
  sim_initialize() %>%
  group_by(ind, prob, experiment, l_rate, temperature) %>%
  mutate(memories = accumulate(memories, .f = update_memo))

flex_res <- flex_test %>%
  filter(n_poke > 1) %>%
  mutate(loc = map_dbl(memories, ~pluck(., "loc_visited")),
         rewarded = map_dbl(memories, ~pluck(., "reward")),
         estimate = map2_dbl(memories, loc, ~pluck(.x, "estimates") %>% pluck(.y)),
         phase = case_when(
           n_poke < target_visits[1] + 1 ~ 1,
           between(n_poke, target_visits[1] + 1, target_visits[2] + 1) ~ 2,
           n_poke > target_visits[2] + 1 ~ 3
         )) %>%
  left_join(rewarding)

flex_sim <- flex_res %>%
  filter(!is.na(mask), prob == 0.5, phase > 2) %>%
  group_by(cond, prob, ind, l_rate, temperature) %>%
  mutate(correct_cumsum = cumsum(mask2),
         reward_cumsum = cumsum(rewarded)) %>%
  filter(reward_cumsum > 0) %>%
  mutate(n_phase2 = 1:n()) %>%
  filter(n_phase2 < 61) %>%
  summarise(correct = sum(mask2, na.rm = TRUE),
            perseverations = sum(mask, na.rm = TRUE),
            flexibility = 60 - perseverations)

write.table(flex_sim, file = paste0(folder, "FlexibilitySimSummary.csv"), sep = ";", row.names = FALSE)

flex_sim <- read.csv2(paste0(folder, "FlexibilitySimSummary.csv"))

bin_size <- 20
binned_flex <- flex_res %>%
  filter(!is.na(mask)) %>%
  group_by(prob, ind, l_rate, temperature, cond) %>%
  mutate(n_poke = n_poke - 1,
         vis_bins = cut(n_poke, breaks = seq(0, 500, bin_size))) %>%
  group_by(prob, ind, l_rate, temperature, vis_bins, cond) %>%
  summarise(n_rewarding = sum(mask),
            n_new_rewarding = sum(mask2),
            n_total = n(),
            phase = mean(phase)) %>%
  group_by(prob, ind, l_rate, temperature, cond) %>%
  mutate(mean_rewarding = mean(n_rewarding, na.rm = TRUE))

annotation_table_lcf <- tibble(label = c("phase 1", "reversal", "chance"),
                               x = c(5, 20, 25),
                               y = c(5, 5, 2.5),
                               l_rate = 0.1,
                               temperature = 0.025,
                               prob = 0.83)

write.table(binned_flex, file = paste0(folder, "FlexibilitySimBinned.csv"), sep = ";", row.names = FALSE)

binned_flex <- read.csv2(paste0(folder, "FlexibilitySimBinned.csv"))

binned_flex %>%
  # filter(ind == 38) %>%
  group_by(prob, l_rate, temperature, cond, vis_bins) %>%
  summarise(n_rewarding = mean(n_rewarding),
            n_new_rewarding = mean(n_new_rewarding)) %>%
  ggplot() +
  geom_line(aes(as.numeric(vis_bins), n_rewarding, color = as.factor(l_rate), group = cond)) +
  geom_line(aes(as.numeric(vis_bins), n_new_rewarding, color = as.factor(l_rate), group = cond), linetype = 2) +
  facet_grid(prob ~ temperature, labeller = global_labeller) +
  labs(title = "Mean flexibility learning curves as a function of exploration, learning rate, and reward probability",
       x = paste("Block of", bin_size, "visits"), y = "visits to rewarding flowers",
       color = "learning rate") +
  theme_bw() +
  geom_vline(xintercept = 100/bin_size, linetype = 2) +
  geom_text(data = annotation_table_lcf, aes(x = x, y = y, label = label)) +
  geom_hline(yintercept = 2/12*bin_size, linetype = 3)
