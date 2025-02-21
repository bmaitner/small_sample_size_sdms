# Figure 2 code

# Load libraries

library(BIEN)
library(geodata)
library(ggplot2)
library(pbsdm)
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
source("C:/Users/Brian Maitner/Desktop/current_projects/pbsdm/R/get_env_pres.R")
source("C:/Users/Brian Maitner/Desktop/current_projects/pbsdm/R/get_env_bg.R")
# Make a temporary directory to store climate data

temp <- tempdir()

# Get some occurrence data

tv <- BIEN_occurrence_species(species = "Trillium vaseyi")
#tv <- BIEN_occurrence_species(species = "Sabal palmetto")

# Get environmental data

env <- worldclim_global(var = "bio",
                        res = 10,
                        path = temp)

# To make things a bit faster and easier, we'll limit ourselves to the 2 variables (mean temperature and annual precipitation)

env <- env[[c(1)]]

# Here, we'll use the same data as before for Trillium vaseyi.

#First, we'll select the background data


library(pbsdm)
tv_bg <- get_env_bg(coords = tv[c("longitude","latitude")],
                    env = env,
                    width = 500000,
                    standardize = FALSE) #note that we used a small set of background points to expedite model fitting

# The returned object 'xs_bg' contains two objects:
# 1) tv_bg$env a matrix of environmental covariates. This is what we need for modeling.
# 2) tv_bg$bg_cells a vector containing the environmental raster cell IDs that are present in tv_bg$env. This is useful for mapping the results.

# Next, we get the presence data:

tv_presence <- get_env_pres(coords = tv[c("longitude","latitude")],
                            env = env)

tv_gaussian <- fit_plug_and_play(presence = tv_presence$env,
                                background = tv_bg$env,
                                method = "gaussian")


tv_rb <- fit_plug_and_play(presence = tv_presence$env[1:6,],
                                 background = tv_bg$env,
                                 presence_method = "rangebagging",
                                 background_method = "none")

pbsdm:::pnp_rangebagging(data = tv_presence$env[1:5,],v = 100,d = 2,p = 0.5,method="fit")



tv_gaussian_predictions <- project_plug_and_play(pnp_model = tv_gaussian,
                                                data = tv_bg$env)


# Figure

  #combine presence, background, and prediction

presence_dist <- pbsdm:::pnp_gaussian(data = tv_bg$env,
                     object = tv_gaussian$f1,
                     method = "predict") %>% exp()

background_dist <- pbsdm:::pnp_gaussian(data = tv_bg$env,
                                        object = tv_gaussian$f0,
                                        method = "predict") %>% exp()

#plot(tv_gaussian_predictions,(presence_dist/background_dist)) #sanity check

  bind_rows(
    data.frame(temp = tv_bg$env,
               dist = "presence",
               values = presence_dist),
    data.frame(temp = tv_bg$env,
               dist = "background",
               values = background_dist),
    data.frame(temp = tv_bg$env,
               dist = "Relative Occurrence Rate",
               values = tv_gaussian_predictions),
    ) %>% `row.names<-`(NULL) %>%
    rename(temp = wc2.1_10m_bio_1) -> plot_data
  
  plot_data$dist <- factor(x = plot_data$dist,levels=c("Relative Occurrence Rate","presence","background"))

    plot_data %>%
      ggplot(mapping = aes(x = temp,
                           y = values,
                           color = dist))+
    geom_line(linewidth=2)+
    xlab(expression('Temperature ('*~degree*C*')'))+
    ylab("Relative Occurrence Rate\n(or Density)")+
    scale_color_manual(values=c("magenta", "#70ffdf", "#9d4dff"),
                         labels=c("Relative\nOccurrence\nRate",
                                  "Presence",
                                  "Background"))+
      geom_vline(xintercept = min(tv_presence$env),lty=2)+
      geom_vline(xintercept = max(tv_presence$env),lty=2)+
      # geom_vline(xintercept = min(tv_bg$env),lty=1)+
      # geom_vline(xintercept = max(tv_bg$env),lty=1)+
    scale_x_continuous(expand=c(0,0))+
    theme_bw()+
    theme(legend.title=element_blank()) -> fig2

    fig2

    ggsave(filename = "figures/Figure2.jpg",
           plot = fig2,
           width = 6,
           height = 4,
           units = "in")    



# version of fig 2 with presence and background points
tv_presence$env
    
pb_points <- data.frame(temp = tv_presence$env$wc2.1_10m_bio_1,
                        Point = "Presence",
                        values = -.2)%>%
  bind_rows(data.frame(temp = tv_bg$env$wc2.1_10m_bio_1,
                       Point = "Background",
                       values = -.5))

fig2_v2 <-
plot_data %>%
  ggplot(mapping = aes(x = temp,
                       y = values,
                       color = dist))+
  
  geom_vline(data = pb_points %>%
               filter(Point == "Background"),
             mapping = aes(xintercept=temp),color="lightgrey")+
  
  geom_vline(data = pb_points %>%
               filter(Point == "Presence"),
             mapping = aes(xintercept=temp),color="black")+

  geom_line(linewidth=2)+
  xlab(expression('Temperature ('*~degree*C*')'))+
  ylab("Relative Occurrence Rate\n(or Density)")+
  scale_color_manual(values=c("magenta", "#70ffdf", "#9d4dff"),
                     labels=c("Relative\nOccurrence\nRate",
                              "Presence",
                              "Background"))+
  # geom_vline(xintercept = min(tv_presence$env),lty=2)+
  # geom_vline(xintercept = max(tv_presence$env),lty=2)+
  # geom_vline(xintercept = min(tv_bg$env),lty=1)+
  # geom_vline(xintercept = max(tv_bg$env),lty=1)+
  scale_x_continuous(expand=c(0,0))+
  theme_bw()+
  theme(legend.title=element_blank())
  #+
  # geom_point(data = pb_points%>%
  #              filter(Point == "Presence"),
  #            mapping = aes(x=temp,y=values),
  #            inherit.aes = FALSE,
  #            alpha=0.1)
  #+


    # geom_jitter(data = pb_points%>%
    #               filter(Point == "Presence"),
    #             mapping = aes(x=temp,y=values),
    #             inherit.aes = FALSE,
    #             height = 0.1)

ggsave(filename = "figures/Figure2v2.jpg",
       plot = fig2_v2,
       width = 6,
       height = 4,
       units = "in")    


######################################


fig2_v3 <-
  plot_data %>%
  ggplot(mapping = aes(x = temp,
                       y = values,
                       color = dist))+
  geom_line(linewidth=2)+
  xlab(expression('Temperature ('*~degree*C*')'))+
  ylab("Relative Occurrence Rate\n(or Density)")+
  scale_color_manual(values=c("magenta", "#70ffdf", "#9d4dff"),
                     labels=c("Relative\nOccurrence\nRate",
                              "Presence",
                              "Background"))+
  # geom_vline(xintercept = min(tv_presence$env),lty=2)+
  # geom_vline(xintercept = max(tv_presence$env),lty=2)+
  # geom_vline(xintercept = min(tv_bg$env),lty=1)+
  # geom_vline(xintercept = max(tv_bg$env),lty=1)+
  scale_x_continuous(expand=c(0,0))+
  theme_bw()+
  theme(legend.title=element_blank())+
  # geom_point(data = pb_points%>%
  #              filter(Point == "Presence"),
  #            mapping = aes(x=temp,y=values),
  #            inherit.aes = FALSE,
  #            alpha=0.1)+
  geom_jitter(data = pb_points%>%
                filter(Point == "Presence"),
              mapping = aes(x=temp,y=values),
              inherit.aes = FALSE,
              height = 0.1)+
  geom_jitter(data = pb_points%>%
                filter(Point == "Background")%>%
                slice_sample(n = 1000),
              mapping = aes(x=temp,y=values),
              inherit.aes = FALSE,
              height = 0.1,
              color="darkgrey")

fig2_v3

ggsave(filename = "figures/Figure2v3.jpg",
       plot = fig2_v3,
       width = 6,
       height = 4,
       units = "in")    
