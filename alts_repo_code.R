library(tidyverse)
library(proxy)
library(reshape2)
library(car)
library(broom)
library(lme4)
library(pbapply)
library(mipfp)
library(caret)


#load data
data_raw <- read_csv('alts_data.csv')

#load csv with list of participants who did study at both T1 and T2
valid_parts <- read_csv("valid_participants.csv")



#DATA EXCLUSION ONLY: creates df including RT and with median RT column, 
#then plots boxplot showing distribution of median rt for each part_time

median_rt_data <- data_raw %>% 
  mutate(expName = str_replace(expName, 'idiolect_t1', 't1'),
         expName = str_replace(expName, 'idiolect_t2', 't2'),
         selection = if_else(response.keys == 'right', right_ident, left_ident),
         part_time = paste(participant, expName, sep = '_')) %>% 
  select(participant, expName, part_time, sent_ident, selection, response.rt) %>% 
  filter(!is.na(sent_ident)) %>% 
  pivot_wider(id_cols = c(participant, expName, part_time), names_from = sent_ident, 
              values_from = response.rt) %>% 
  subset(duplicated(participant) | duplicated(participant, fromLast=TRUE)) %>% 
  select(-participant, -expName) %>% 
  rowwise() %>%
  mutate(median_rt = median(c_across(where(is.numeric))))

ggplot(median_rt_data) +
  aes(x = median_rt) +
  geom_histogram(fill = '#00BFC4', binwidth = 0.4) +
  theme_minimal() +
  labs(x = 'Median RT per response (seconds)', y = 'Count')

#Low RTs of outlier participant at T2
data_raw %>% 
  filter(participant == '5fefbee',
         expName == 'idiolect_t2',
         !is.na(sent_ident)) %>% 
  select(response.rt) %>% 
  print(n = 40) -> outlier_RTs

mean(outlier_RTs$response.rt)




#--------ANALYSIS 1A--------

#Logistic regression to measure within-speaker consistency 

#reorganise into analysable form with binary responses for logistic regression
data_raw %>% 
  filter(!is.na(sent_ident),
         participant %in% valid_parts$x,
         participant != "5fefbee") %>%
  mutate(left_ident_binary = if_else(expName == 'idiolect_t1', 0, 1), 
         right_ident_binary = if_else(expName == 'idiolect_t1', 1, 0),
         selection_binary = if_else(response.keys == 'right', 
                                    right_ident_binary, left_ident_binary)) %>% 
  select(participant, expName, sent_ident, selection_binary) %>%
  pivot_wider(names_from = expName, values_from = selection_binary) %>% 
  rename('item' = 'sent_ident',
         'resp_T1' = 'idiolect_t1',
         'resp_T2' = 'idiolect_t2') -> results



#create matrix showing raw numbers of response combinations
conf_matrix <- confusionMatrix(as.factor(results$resp_T2), as.factor(results$resp_T1))
conf_matrix

#turn into data frame
tbl_resp_combs <- data.frame(T2_0 = c(1183, 533),
                             T2_1 = c(733, 991),
                             row.names = c('T1_0', 'T1_1'))



#mixed effects logistic regression model with by-participant and by-item varying slopes
model <- glmer(resp_T2 ~ resp_T1 + (1 + resp_T1|participant) + (1 + resp_T1|item),
               data = results, family = 'binomial')
summary(model)


#get fixed effects coefficients table
model_fixed_effects_coefs <- as.data.frame(coef(summary(model)))

model_fixed_effects_coefs %>% 
  mutate(Estimate = round(Estimate, 2),
         `Std. Error` = round(`Std. Error`, 2),
         `z value` = round(`z value`, 2),
         `Pr(>|z|)` = round(`Pr(>|z|)`, 2)) -> tbl_fixed_effects

tbl_fixed_effects


#get probabilities table
intercept <- fixef(model)[1]
slope <- fixef(model)[2]

topright <- plogis(intercept + slope * 0)
bottomright <- plogis(intercept + slope * 1)
topleft <- 1 - plogis(intercept + slope * 0)
bottomleft <- 1 - plogis(intercept + slope * 1)

mainprobs <- data.frame(T2_0 = c(topleft, bottomleft),
                        T2_1 = c(topright, bottomright),
                        row.names = c('T1_0', 'T1_1'))

mainprobs %>% 
  mutate(T2_0 = round(T2_0, 2),
         T2_1 = round(T2_1, 2)) -> tbl_main_probs

tbl_main_probs



#get random effects coefficients table
model_random_effects_coefs <- as.data.frame(VarCorr(model))

model_random_effects_coefs %>% 
  unite(Model, c(var1, var2), sep = '+', na.rm = TRUE) %>% 
  rename('Groups' = 'grp') %>% 
  mutate(vcov = round(vcov, 2),
         sdcor = round(sdcor, 2)) -> tbl_random_effects

tbl_random_effects

#Exploring random effects
coef(model)
ranef(model)

#showing SP6 results as example of unbalanced variable
results %>% 
  filter(item == 'SP6') %>% 
  print(n = 86) -> SP6_results

sum(SP6_results$resp_T1 == '0')
sum(SP6_results$resp_T1 == '1')




#Permutation procedure to measure between-speaker consistency and compare w/ within

#used random intercepts model only due to convergence issues. 

#Model for reference. Slope coefficient is what perms are compared to to produce p-value:
int_model <- glmer(resp_T2 ~ resp_T1 + (1|participant) + (1|item), 
                   data = results, family = 'binomial')
summary(int_model)



one_permutation <- function(results){
  
  permuted_data <- results %>% 
    arrange(factor(participant, levels = sample(c(unique(results$participant)))))
  
  results$resp_T1 <- permuted_data$resp_T1
  
  M1 <- glmer(resp_T2 ~ resp_T1 + (1|participant) + (1|item), 
              data = results, family = 'binomial')
  
  return(fixef(M1)['resp_T1'])
}


parallel_permut <- function(results, r, cores){
  
  M2 <- glmer(resp_T2 ~ resp_T1 + (1|participant) + (1|item), 
              data = results, family = 'binomial')
  
  permutation_vec <- pbreplicate(r, one_permutation(results), cl = cores)
  
  perm_result_list <- list(
    'p_value' = mean(as.numeric(permutation_vec >= fixef(M2)['resp_T1'])),
    'coefficients' = permutation_vec)
  
  return(perm_result_list)
}


#perform perm procedure (takes a long time)
#note - technically, results may vary on each run due to random permutations
perm_result_list <- parallel_permut(results, r = 10000, cores = 5)

#get p value
perm_result_list$p_value

#get start of list of coefficients from permutations
perm_result_list$coefficients

#some summary stats from the list of coefficients from permutations
mean(perm_result_list$coefficients)
range(perm_result_list$coefficients)






#--------ANALYSIS 1B--------

#explore if difference in within-speaker consistency for WO vs WC alternations
results %>% 
  mutate(alt_type = if_else(startsWith(item, 'S'), 'WO', 'WC')) -> alt_type_results

#model which INCLUDES 3 most unbalanced variables
alt_type_model <- glmer(resp_T2 ~ resp_T1 * alt_type + 
                          (1 + resp_T1 * alt_type|participant) + 
                          (1 + resp_T1 * alt_type|item),
                        data = alt_type_results, family = 'binomial', control = 
                          glmerControl(optimizer = 'bobyqa'))

#get fixed effects coefficients table
alt_type_fixef_nofilter <- as.data.frame(coef(summary(alt_type_model)))

alt_type_fixef_nofilter %>% 
  mutate(Estimate = round(Estimate, 2),
         `Std. Error` = round(`Std. Error`, 2),
         `z value` = round(`z value`, 2),
         `Pr(>|z|)` = round(`Pr(>|z|)`, 2)) -> tbl_alt_type_fixef_nofilter

tbl_alt_type_fixef_nofilter



#model with 3 unbalanced variables REMOVED
alt_type_results %>% 
  filter(item != "SD2",
         item != "SD5",
         item != "SP6") -> alt_type_results_filtered

alt_type_filtered_model <- glmer(resp_T2 ~ resp_T1 * alt_type + 
                                   (1 + resp_T1 * alt_type|participant) + 
                                   (1 + resp_T1 * alt_type|item),
                                 data = alt_type_results_filtered, family = 'binomial', 
                                 control = glmerControl(optimizer = 'bobyqa'))

#get fixed effects coefficients table
alt_type_filtered_model_fixed_coefs <- as.data.frame(coef(summary(alt_type_filtered_model)))

alt_type_filtered_model_fixed_coefs %>% 
  mutate(Estimate = round(Estimate, 2),
         `Std. Error` = round(`Std. Error`, 2),
         `z value` = round(`z value`, 2),
         `Pr(>|z|)` = round(`Pr(>|z|)`, 2)) -> tbl_alt_type_fixed_effects

tbl_alt_type_fixed_effects




#probability tables (for model with unbalanced variables removed)

alt_type_intercept <- fixef(alt_type_filtered_model)[1]
alt_type_resp_T1_slope <- fixef(alt_type_filtered_model)[2]
alt_type_alt_typeWO_slope <- fixef(alt_type_filtered_model)[3]
alt_type_interaction_slope <- 
  fixef(alt_type_filtered_model)[2] + fixef(alt_type_filtered_model)[4]


#WC probabilities table
WC_topright <- plogis(alt_type_intercept + alt_type_resp_T1_slope * 0)
WC_bottomright <- plogis(alt_type_intercept + alt_type_resp_T1_slope * 1)
WC_topleft <- 1 - plogis(alt_type_intercept + alt_type_resp_T1_slope * 0)
WC_bottomleft <- 1 - plogis(alt_type_intercept + alt_type_resp_T1_slope * 1)

WC_probs <- data.frame(T2_0 = c(WC_topleft, WC_bottomleft),
                       T2_1 = c(WC_topright, WC_bottomright),
                       row.names = c('T1_0', 'T1_1'))

WC_probs %>% 
  mutate(T2_0 = round(T2_0, 2),
         T2_1 = round(T2_1, 2)) -> tbl_WC_probs

tbl_WC_probs


#WO probabilities table
WO_topright <- plogis(alt_type_intercept + alt_type_alt_typeWO_slope * 1)
WO_bottomright <- plogis(alt_type_intercept + alt_type_interaction_slope * 1)
WO_topleft <- 1 - plogis(alt_type_intercept + alt_type_alt_typeWO_slope * 1)
WO_bottomleft <- 1 - plogis(alt_type_intercept + alt_type_interaction_slope * 1)

WO_probs <- data.frame(T2_0 = c(WO_topleft, WO_bottomleft),
                       T2_1 = c(WO_topright, WO_bottomright),
                       row.names = c('T1_0', 'T1_1'))

WO_probs %>% 
  mutate(T2_0 = round(T2_0, 2),
         T2_1 = round(T2_1, 2)) -> tbl_WO_probs

tbl_WO_probs







#--------ANALYSIS 2--------

#Get feedback
data_raw %>% 
  select(participant, feedback_textbox.text) %>% 
  filter(!is.na(feedback_textbox.text)) -> feedback


data_raw %>%
  mutate(expName = str_replace(expName, 'idiolect_t1', 't1'),
         expName = str_replace(expName, 'idiolect_t2', 't2'),
         left_ident_binary = if_else(expName == 't1', 0, 1), 
         right_ident_binary = if_else(expName == 't1', 1, 0),
         response = if_else(response.keys == 'right', 
                            right_ident_binary, left_ident_binary)) %>% 
  select(participant, sent_ident, response.rt, expName, response) %>% 
  filter(!is.na(sent_ident),
         participant %in% valid_parts$x,
         participant != "5fefbee") %>% 
  rename("RT" = "response.rt",
         "Time" = "expName",
         "Sentence" = "sent_ident",
         "Participant" = "participant",
         "Response" = "response") -> RT_data_raw


#Organise data so there is one row per participant per sentence
RT_data_raw %>% 
  pivot_wider(names_from = Time, values_from = Response) %>% 
  rename("resp_T1" = "t1", "resp_T2" = "t2") -> RT_data_na_loss


#Identical RT for participant 60a4b39, sent LT1 (row 4416) which causes loss of row. 
#Manually fix by adding row with correct info and replacing resp_T2 value with NA

RT_data_na_wip <- RT_data_na_loss %>% 
  add_row(Participant = "60a4b39", 
          Sentence = "LT1", 
          RT = 3.9110, 
          resp_T1 = NA,
          resp_T2 = 0,
          .before = 4470)

RT_data_na_wip[4416, 5] = NA


#back on track
RT_data_na_wip %>% 
  mutate(Time = RT_data_raw$Time) %>% 
  pivot_wider(names_from = Time, values_from = RT) %>% 
  rename("RT_T1" = "t1", "RT_T2" = "t2") -> RT_data_na


#remove NAs and fill in resp and RT info for different Times
f <- function(x) {
  x <- na.omit(x)
  if (length(x) > 0) paste(x,collapse='-') else NA
}

RT_data_na %>% 
  group_by(Participant, Sentence) %>% 
  summarise_all(list(f)) -> RT_data_na_removed


#add match info
RT_data_na_removed %>% 
  mutate(Match = if_else(resp_T1 == resp_T2, "same", "different")) -> RT_data_match


#rearrange data so every row is a response, but now tagged with match info
RT_data_match %>% 
  pivot_longer(cols = c("resp_T1", "resp_T2"), 
               names_to = "Time", 
               values_to = "Response") -> RT_long

RT_data_match %>% 
  pivot_longer(cols = c("RT_T1", "RT_T2"),
               names_to = "Time_again", 
               values_to = "RT") -> RT_long_two

RT_long %>% 
  ungroup() %>% 
  mutate(RT = RT_long_two$RT) %>% 
  select(Participant, Sentence, RT, Time, Match) -> RT_results_all

RT_results_all$RT <- as.numeric(RT_results_all$RT)



#log-transformation
RT_results_all$RT <- RT_results_all$RT * 1000

RT_results_all %>% 
  mutate(log_RT = log(RT)) -> RT_results_all_log


#Remove RTs with z-score above 3 or below -3 and corresponding RT from same part, same sent
#and different Time.

scale(RT_results_all_log$log_RT) %>% 
  sort(decreasing = TRUE)

scale(RT_results_all_log$log_RT) %>% 
  sort(decreasing = FALSE)

#There are 71 z-scores above 3 so 71 highest to be removed, along with their 
#corresponding one. There are also 15 z-scores below -3, so 15 lowest to be removed,
#along with their corresponding one

RT_results_all_log <- RT_results_all_log[order(RT_results_all_log$log_RT, decreasing = TRUE),] %>% 
  mutate(part_sent = paste(Participant, Sentence, sep = "_"))

RT_results_all_log %>% 
  slice(1:71) -> high_outliers

RT_results_all_log %>% 
  slice(6866:6880) -> low_outliers

outliers <- rbind(high_outliers, low_outliers)

RT_results_all_log %>% 
  filter(!(part_sent %in% outliers$part_sent)) %>% 
  mutate(Time = str_replace(Time, 'resp_T1', 'T1'),
         Time = str_replace(Time, 'resp_T2', 'T2')) -> RT_results


#Boxplot (RTs of responses which match vs RTs of responses which don't match 
#at T1 and T2)

ggplot(RT_results) +
  aes(x = Time, y = log_RT, fill = Match) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Time", y = "LogRT", fill = "Match")


#RT analysis model
RT_model <- lm(log_RT ~ Time * Match, data = RT_results)

#get coefficients table
tidy(RT_model) %>% 
  mutate(estimate = round(estimate, 2),
         std.error = round(std.error, 2),
         statistic = round(statistic, 2),
         p.value = round(p.value, 2)) -> tbl_RT_model

tbl_RT_model
