# 4b plotting rarified models

temp_full_RDS = "outputs/temp_rarified_full.RDS"

rarified_full <- readRDS(temp_full_RDS)

unique(rarified_full$model)

rarified_full %>%
  select(model, n_presence, pa_AUC,pa_sensitivity,pa_specificity) %>%
  mutate(pa_AUC = as.numeric(pa_AUC)) %>%
  pivot_longer(contains("pa"),
               names_to = "metric") %>%
  #filter(model != "gaussian/none")%>%
  mutate(metric = gsub(pattern = "pa_",replacement = "", x =metric))%>%
  mutate(metric = str_to_title(metric))%>%
  mutate(metric = gsub(pattern = "Auc",replacement = "AUC", x =metric))%>%
  # mutate(model = case_when(model == "gaussian/gaussian" ~ "Gaussian",
  #                           model == "kde/kde" ~ "KDE",
  #                           model == "maxnet" ~ "MaxNet",
  #                           model == "rangebagging/none" ~ "Rangebagging"))%>%
  na.omit() %>%
  rename(Model = model) %>%
  ggplot(mapping = aes(x=n_presence,y=value,color=Model))+
  geom_point(alpha=0.01)+
  # geom_smooth(method = "lm",
  #             se = FALSE)+
  geom_smooth(method = "loess",
              se = FALSE)+
  facet_wrap(~ metric,scales = "free_y")+
  # scale_x_continuous(trans = "log10")+
  xlab("Occurrence Records")+
  ylab(NULL)+
  theme_bw()

# Models indistinguishable from maxnet at large sample sizes

rarified_full %>%
  select(model, n_presence, pa_AUC,pa_sensitivity,pa_specificity) %>%
  mutate(pa_AUC = as.numeric(pa_AUC)) %>%
  pivot_longer(contains("pa"),
               names_to = "metric") %>%
  filter(model %in% c("kde/kde","gaussian/kde","rulsif","maxnet"))%>%
  mutate(metric = gsub(pattern = "pa_",replacement = "", x =metric))%>%
  mutate(metric = str_to_title(metric))%>%
  mutate(metric = gsub(pattern = "Auc",replacement = "AUC", x =metric))%>%
  # mutate(model = case_when(model == "gaussian/gaussian" ~ "Gaussian",
  #                           model == "kde/kde" ~ "KDE",
  #                           model == "maxnet" ~ "MaxNet",
  #                           model == "rangebagging/none" ~ "Rangebagging"))%>%
  na.omit() %>%
  rename(Model = model) %>%
  ggplot(mapping = aes(x=n_presence,y=value,color=Model,lty=Model))+
  geom_point(alpha=0.01)+
  # geom_smooth(method = "lm",
  #             se = FALSE)+
  geom_smooth(method = "loess",
              se = TRUE,
              mapping = aes(fill=Model))+
  facet_wrap(~ metric)+
  # scale_x_continuous(trans = "log10")+
  xlab("Occurrence Records")+
  ylab(NULL)+
  theme_bw()
