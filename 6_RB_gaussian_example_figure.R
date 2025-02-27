# Example rangebagging vs gaussian figure


# Load libraries

library(BIEN)
library(ggplot2)
library(S4DM)
library(sf)
library(terra)
library(tidyterra)

# Make a temporary directory to store climate data

temp <- tempdir()

# Get some occurrence data

#tv <- BIEN_occurrence_species(species = "Trillium vaseyi")
data("sample_points")


# Get environmental data
# To make things a bit faster and easier, we'll limit ourselves to the 2 variables (mean temperature and annual precipitation)


# env <- worldclim_global(var = "bio",
#                          res = 10,
#                          path = temp)
# env <- env[[c(1,12)]]

env <- rast(system.file('ex/sample_env.tif', package="S4DM"))  

env <-env[[1]]

# And we'll rescale the variables as well

env <- scale(env)

# Just to take a look to make sure we didn't mess anything up

plot(env)


# make models

sample_points_bg <- get_env_bg(coords = sample_points[c("longitude","latitude")],
                               env = env,
                               width = 50000) #note that we used a small set of background points to expedite model fitting

sample_points_pres <- get_env_pres(coords = sample_points[c("longitude","latitude")],
                                   env = env)


rbna <- S4DM::fit_plug_and_play(presence = sample_points_pres$env,
                                background = sample_points_bg$env,
                                presence_method = "rangebagging",
                                background_method = "none",
                                d=1)

naga <- S4DM::fit_plug_and_play(presence = sample_points_pres$env,
                                background = sample_points_bg$env,
                                presence_method = "none",
                                background_method = "gaussian")



# make data for plotting


pred_data <- data.frame(wc2.1_10m_bio_1 = seq(from=-3,to=3,by=.1))

rbna_pred <- S4DM::project_plug_and_play(pnp_model = rbna,
                                         data = pred_data)

naga_pred <- S4DM::project_plug_and_play(pnp_model = naga,
                                         data = pred_data)


rbga_pred <- rbna_pred/naga_pred


library(ggplot2)
library(tidyverse)

plot_data <-
  bind_rows(
    data.frame(model = "rangebagging",
               env = pred_data,
               value= rbna_pred),
    
    data.frame(model = "gaussian",
               env = pred_data,
               value = naga_pred),
    
    data.frame(model = "rangebagging/gaussian",
               env = pred_data,
               value = rbga_pred)
    
  )%>%
  rename(env = wc2.1_10m_bio_1)

ggplot(data = plot_data,
       mapping = aes(x=env,
                     y=value,
                     color=model),
       alpha=0.5)+
  geom_line()


##############################

# simpler model to show the overall i`dea

plot_data <-
bind_rows(
  data.frame(env = seq(-3,3,.1),
             model = "Background",
             value= dnorm(mean = 0,sd = 1,
                   x = seq(-3,3,.1))),
  data.frame(env = seq(-3,3,.1),
             model = "Presence",
             value= dunif(min = -1,max = 1,
                   x = seq(-3,3,.1)))
  )%>%
  mutate(value = case_when(model=="Presence" & value > 0 ~1,.default = value))%>%
  pivot_wider(names_from = model,values_from = value) %>%
  mutate("Relative\nOccurrence\nRate" = Presence/Background)%>%
  pivot_longer(cols = 2:4,names_to = "model")

example_plot <-
plot_data %>%
ggplot(mapping = aes(x=env,y=value,color=model))+
  geom_line(linewidth=1,alpha=0.75)+
  theme_bw()+
  scale_color_manual(values = c('pink',"grey","purple"))+
  xlab("environment")+
  ylab("f(environment)")+
  theme(legend.title = element_blank())

ggsave(filename = "figures/example_rb_gaussian.jpeg",
       plot = example_plot,
       width = 5,height = 5,units = "in",dpi = 600)


ggsave(filename = "figures/example_rb_gaussian.png",
       plot = example_plot,
       width = 5,height = 5,units = "in",dpi = 600)
