library(tidyverse)
library(here)
library(cmdstanr)

# Setup 
## Set Path
path <- here()

## STAN Setup
cmd_path<-paste0("C:\\Coding\\cmdstan-2.32.1")
set_cmdstan_path(path=cmd_path)

# Model COmpilation
m3_CS <- cmdstan_model("Models/M3_ComplexSpan_CS_Choleksy.stan")

# Initialization function for sampling
init_pre <- function()
{
  list(hyper_pars=c(runif(stan.dat$J,10,20)),
       subj_pars=c(runif(stan.dat$N,1,10)))
  
}

# Read in Data ----


# Set Data Set 1 = Visual, 2 = Verbal
dataset <- 1

ifelse(dataset==1, path_data <- paste0(path,"/Data_Visual/New Timing"),
       path_data <- paste0(path,"/Data_Verbal/New Timing"))


# Load Datasets


files_ran <-list.files(path_data,pattern = "m3",full.names = T)
#files_seq <-list.files(path_data,pattern = "Ran",full.names = T)

df_ran <- files_ran %>% map(read_csv) %>% reduce(rbind) 
#df_seq <- files_seq %>% map(read_csv) %>% reduce(rbind)

df_ran <- df_ran %>% mutate(NPLs = str_remove(NPLs,pattern = ",.*"), NPLs=as.numeric(NPLs))# bug for last column is character fixed in next revision
#df_seq <- df_seq %>% mutate(NPLs = str_remove(NPLs,pattern = ",.*"), NPLs=as.numeric(NPLs))

# Sum Up Categorie Choices

choices_ran <- df_ran %>% select(Subject,SumACC_SecTask,IIP:NPLs) %>% filter(!Subject %in% c(21,24)) %>% group_by(Subject) %>% 
  summarise(across(IIP:NPLs, list(sum), .names = "Sum{.col}"), 
            Retrievals = sum(SumIIP,SumIOP,SumDIP,SumDIOP,SumNPLs),
            ACC_Main=SumIIP/Retrievals,
            ACC_Sec = mean(SumACC_SecTask)) %>% mutate(Condition = "Random Recall") %>% 
  relocate(Condition, .after =Subject) %>% #filter(!Subject %in% c(11,19)) %>% 
  mutate(across(SumIIP:SumNPLs, ~.x / Retrievals, .names="p{.col}")) %>% summarise(across(starts_with("p"), mean ))

choices_seq <- df_seq %>% select(Subject,SumACC_SecTask,IIP:NPLs) %>% group_by(Subject) %>% 
  summarise(across(IIP:NPLs, list(sum), .names = "Sum{.col}"), 
            Retrievals = sum(SumIIP,SumIOP,SumDIP,SumDIOP,SumNPLs),
            ACC_Main=SumIIP/Retrievals,
            ACC_Sec = mean(SumACC_SecTask)) %>% mutate(Condition = "Sequential Recall") %>% 
  relocate(Condition, .after =Subject)


df_full <-choices_ran   
df_full %>% group_by(Condition) %>% summarise(meanDIP = mean(SumDIP), meanDIOP = mean(SumDIOP), 
                                              meanIOP=mean(SumIOP),meanNPL = mean(SumNPLs),
                                              meanACC = mean(ACC_Main),
                                              meanACC_secondary=mean(ACC_Sec))


# Plot

df_full %>% group_by(Subject,Condition) %>%
  summarise(meanACC = mean(ACC_Main), meanDIP = mean(SumDIP), meanDIOP = mean(SumDIOP)) %>%
  ggplot(., aes(x=Condition,y=meanACC, fill=Condition)) +
  geom_bar(stat = "identity") + scale_fill_brewer(palette="Dark2")

# M3 - Modeling ----

## Add Weight Matrix for random recall

d_weights <- df_ran %>% select(Subject,n_DIP_Pos1:n_DIP_Pos6, n_DIP_total, DIP, DIOP) %>% 
  mutate(n_DIOP_Pos1 = case_when(
    n_DIP_Pos1 == 0 ~ 0,
    n_DIP_Pos1 == 1 ~ n_DIP_total -1, 
    n_DIP_Pos1 == 2 ~ n_DIP_total -2),
  
  n_DIOP_Pos2 = case_when(   n_DIP_Pos2 == 0 ~ 0,
                             n_DIP_Pos2 == 1 ~ n_DIP_total -1, 
                             n_DIP_Pos2 == 2 ~ n_DIP_total -2),
  n_DIOP_Pos3 = case_when(        n_DIP_Pos3 == 0 ~ 0,
                                  n_DIP_Pos3 == 1 ~ n_DIP_total -1, 
                                  n_DIP_Pos3 == 2 ~ n_DIP_total -2),
  n_DIOP_Pos4 = case_when(n_DIP_Pos4 == 0 ~ 0,
                          n_DIP_Pos4 == 1 ~ n_DIP_total -1, 
                          n_DIP_Pos4 == 2 ~ n_DIP_total -2),
  n_DIOP_Pos5 = case_when(n_DIP_Pos5 == 0 ~ 0,
                          n_DIP_Pos5 == 1 ~ n_DIP_total -1, 
                          n_DIP_Pos5 == 2 ~ n_DIP_total -2),
  n_DIOP_Pos6 = case_when(n_DIP_Pos6 == 0 ~ 0,
                        n_DIP_Pos6 == 1 ~ n_DIP_total -1, 
                        n_DIP_Pos6 == 2 ~ n_DIP_total -2)) %>% group_by(Subject) %>% 
  summarise(DIP_total = sum(n_DIP_total), 
            DIPs = sum(DIP),
            DIOPs = sum(DIOP),
            DIP_weight =  DIP_total/ 114, 
            DIOP_total = sum(across(.cols= n_DIOP_Pos1:n_DIOP_Pos6, sum)), 
            DIOP_weight = DIOP_total /114 ) %>% select(Subject, DIPs,DIP_total, DIP_weight,DIOPs, DIOP_total,DIOP_weight) 


## Weight Matrix DIP and DIOP 
df_full <- cbind(df_full, d_weights$DIP_weight, d_weights$DIOP_weight)
d_weight <- d_weights %>% select(DIP_weight,DIOP_weight) 



## Data creation for STAN ----
N = length(unique(choices_ran$Subject))
respCat = c(1,5,1,5,10)


stan.dat <- list(count = as.matrix(choices_ran[,3:7]), 
                 N = N,
                 K = 5, 
                 J = 3,
                 d_weight = d_weight,
                 retrievals = choices_ran$Retrievals[1],
                 R = respCat,
                 scale_b = 0.1)


## Fit the shit and extract Parameters  ----
fit3_pre <- m3_CS$sample(data = stan.dat,
                         refresh = 100,
                         chains = 4,
                         parallel_chains = 4,
                         iter_warmup = 1500,
                         iter_sampling = 3000,
                         adapt_delta = .99,
                         max_treedepth = 15,
                         init = init_pre,
                         show_messages = F)

## Evaluation of Fit ----

# Extract Parameters for pre cue condition (Multivariate)
M3_hyper <- fit3_pre$summary(c("hyper_pars","mu_f"),mean)
M3_f <- fit3_pre$summary(c("f"), mean,sd)
mean(M3_f$mean)

fit3_pre$summary()

M3_subj <- fit3_pre$summary(c("subj_pars"), mean)
M3_count_rep <- fit3_pre$summary(c("count_rep"),mean)
M3_omega <- fit3_pre$summary("cor_mat_lower_tri",mean)

# Tidy Subject Parameters
subj <- M3_subj  %>% mutate(variable = str_remove_all(variable, "subj_pars")) %>%
  separate(col = variable,into = c("theta","Subject"),sep = ",")  %>%
  mutate(Subject = str_remove(Subject,pattern = "]"), 
         theta = case_when(theta == "[1" ~ "c",
                           theta == "[2" ~ "a",
                           theta == "[3" ~ "logMu_f")) %>%
  pivot_wider(.,names_from = "theta",values_from = c("mean"))%>% mutate(Subject = as.integer(Subject))

f <- M3_f %>% mutate(Subject = seq(1:stan.dat$N), theta = "f") %>%
  relocate(c("Subject","theta"), .before = variable) %>% select(-variable) %>%
  pivot_wider(.,names_from = "theta",values_from = c("mean"))

# One Dataset for all Parameters -- save to File
theta_subject <- left_join(subj,f, by="Subject") %>% 
  relocate(mean_f,.after = mean_a) %>% select(-mean_logMu_f,-rhat_logMu_f,-sd)


## Model Fit - Predictive Modelfit ----
count_rep <- M3_count_rep %>% separate(variable,into = c("Subject","Category"),sep = ",") %>%
  mutate(Subject= str_remove_all(Subject,pattern="count_rep\\["), Category=str_remove_all(Category,"\\]"), mean=round(mean,0)) %>%
  pivot_wider(., names_from = Category, values_from = mean) %>%
  rename("IIP_rep"=`1`,"IOP_rep" = `2`,"DIP_rep" = `3`,"DIOP_rep" = `4`,"NPLs_rep" = `5`) %>% mutate(Subject = as.integer(Subject))


count_fit <- left_join(count_rep,df_clean,by="Subject")  
cor(count_fit$IIP_rep, count_fit$SumIIP)  
cor(count_fit$IOP_rep, count_fit$SumIOP)   
cor(count_fit$DIP_rep, count_fit$SumDIP) 
cor(count_fit$DIOP_rep, count_fit$SumDIOP)
cor(count_fit$NPLs_rep, count_fit$SumNPLs)

rhat()
