

library(survival)
library(broom)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(survminer)
library(dplyr)

###################################################
# Extract model estimates
###################################################

forest_dat <-
  broom::tidy(
    cox_source_full,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
  mutate(
    Variable =
      c(
        "Distance to nearest invaded region\n(per 100 km)",
        "Invasion pressure\n(per SD increase)",
        "Neighbourhood transmission intensity\n(per case increase)",
        "Annual mean temperature\n(per 1°C increase)",
        "Relative humidity\n(per 1% increase)",
        "Annual precipitation\n(per SD increase)",
        "Intermediate urbanisation\n(vs urban)",
        "Rural\n(vs urban)"
        
      ),
    
    Domain =
      
      c(
        "Spatial diffusion",
        "Spatial diffusion",
        "Spatial diffusion",
        "Climate",
        "Climate",
        "Climate",
        "Urbanisation",
        "Urbanisation"
      ),
    
    HR_CI =
      
      sprintf(
        "%.2f (%.2f–%.2f)",
        estimate,
        conf.low,
        conf.high
      ),
    
    p_label =
      
      case_when(
        
        p.value < 0.001 ~ "<0.001",
        
        TRUE ~ sprintf("%.3f", p.value)
      )
    
  )

forest_dat$Variable <- factor(
  
  forest_dat$Variable,
  
  levels = rev(
    
    c(
      # Spatial diffusion
      "Distance to nearest invaded region\n(per 100 km)",
      "Invasion pressure\n(per SD increase)",
      "Neighbourhood transmission intensity\n(per case increase)",
      # Climate
      "Annual mean temperature\n(per 1°C increase)",
      "Relative humidity\n(per 1% increase)",
      "Annual precipitation\n(per SD increase)",
      # Land-use
      "Intermediate urbanisation\n(vs urban)",
      "Rural\n(vs urban)"
    )))
forest_dat

library(ggplot2)
library(ggthemes)

Fig_Forest <-
  
  ggplot(
    forest_dat,
    aes(
      y =Variable,
      x = estimate,
      colour =Domain   )
    
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = 2,
    colour = "grey50"
  ) +
  
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high
       ),
    height = 0.15,
    linewidth = 1
  ) +
  geom_point( size = 3 ) +
  geom_text(
    aes(
      x = 2.7,
      label = HR_CI
    ),
   colour = "black",
   hjust = 0,
   size = 3.6
  ) +
  
  geom_text(
    aes(
      x = 5.0,
      label = p_label
    ),
    colour = "black",
    hjust = 0,
    size = 3.6
  #  fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = 2.7,
    y = 8.5,
    label = "HR (95% CI)",
   #fontface = "bold",
   size=4,
    hjust = 0  ) +
  annotate(
    "text",
    x = 5.0,
    y = 8.5,
    label = "P-value",
   # fontface = "bold",
    size=4,
    hjust = 0
  ) +
  
  scale_color_manual(
    values =
      c(
        "Spatial diffusion" = "#CB181D",  # dark red
        "Climate" = "#2171B5",            # blue
        "Urbanisation" = "#238B45"            # green
           )) +
  
  scale_x_log10(
    
    breaks =
      c(
        0.5,
        0.75,
        1,
        1.5,
        2,
        3,
        5
      )
    
  ) +
  
  coord_cartesian(
    
    xlim =
      c(
        0.4,
        6
      )
    
  ) +
  
  labs(
    
    x = "Hazard ratio (HR)",
    
    y = NULL,
    
    colour = NULL
    
  ) +
  
  theme_classic2(
    
    base_size = 13
    
  ) +
  
  theme(
    
    legend.position = "bottom",#c(0.3,0.98),# "top",
    legend.direction = "horizontal",
    axis.text.y =
      element_text(size = 11),
    
    plot.margin =
      margin(
        5,
        90,
        5,
        5
      )
    
  )

Fig_Forest


### Figure 3B Spline and KM analysis:#############



############################################################
# PANEL B FIGURE
############################################################

cox_spline <- coxph(
  
  Surv(
    start,
    stop,
    event
  ) ~
    
    ns(
      distance100,
      df = 4
    ),
  
  data =  analysis_cases
  
)

############################################################
# Prediction grid
############################################################

upper_limit <-
  
  quantile(
    
    first_invaded$nearest_invaded_km,
    
    probs = 0.99
    
  )

newdat <- data.frame(
  
  distance100 =
    
    seq(
      
      min(analysis_cases$distance100),
      
      upper_limit / 100,
      
      length.out = 1000
      
    )
)

############################################################
# Reference distance = 50 km
############################################################

ref_lp <-
  
  predict(
    
    cox_spline,
    
    newdata =
      data.frame(
        distance100 = 1
      ),
    
    type = "lp"
    
  )

############################################################
# Predictions
############################################################

pred <-
  
  predict(
    
    cox_spline,
    
    newdata = newdat,
    
    type = "lp",
    
    se.fit = TRUE
    
  )

newdat <-
  newdat %>%
  mutate(HR =exp(pred$fit - ref_lp ),
    HR_low =
      exp(
        pred$fit -
          ref_lp -
          1.96 *
          pred$se.fit
      ),
    HR_high =
      exp(
        pred$fit -
          ref_lp +
          1.96 *
          pred$se.fit
      ),
    distance_km =
      distance100 * 100
    
  )

FigB <-
  ggplot(
    newdat,
    aes(
      distance_km,
      HR   )
     ) +
  
  geom_ribbon(
    aes(
    ymin = HR_low,
     ymax = HR_high
    ),
    fill = "#FCAE91",
    alpha = 0.35
  ) +
  
  geom_line(
    colour = "#CB181D",
    linewidth = 1.3
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_rug(
    data = first_invaded,
    aes(
      x = nearest_invaded_km
    ),
    
    inherit.aes = FALSE,
    sides = "b",
    
    alpha = 0.4
    
  ) +
  
  labs(
    
    x =
      "Distance to nearest invaded region (km)",
    y =
      "Relative invasion hazard"
    
  ) +
  
  theme_classic2(
    base_size = 13
  )+
  #scale_x_log10()
  coord_cartesian(xlim = c(0,600))+
  annotate(
    
    "text",
    
    x = 300,
    
    y = 1.15,
    
    label = "Reference = 100 km",
    
    size = 3
    
  )



analysis_cases$distance_cat <- cut(
  analysis_cases$nearest_invaded_km,
  breaks = c(0,75,150,300,500,Inf),
  include.lowest = TRUE
  )

fit_km <- survfit(
  
  Surv(
    start,
    stop,
    event
  ) ~ distance_cat,
  data = analysis_cases
  
)


km_inset <- ggsurvplot(
  fit_km,
  conf.int = FALSE,
  risk.table = FALSE,
  censor = FALSE,
  palette = c(
    "#CB181D",
    "#FB6A4A",
    "#FC9272",
    "#FCAE91",
    "#FEE5D9"
  ),
  legend.title = "Distance to \n nearest invaded region",
  legend.labs = c(
    "0–75 km",
    "75–150 km",
    "150–300 km",
    "300–500 km",
    ">500 km"
  ),
  xlab = "Year since first WNV invasion",
  ylab = "Probability a näive region \n survives the WNV insvasion ",
  
  ggtheme =
    theme_classic2(base_size = 7) +
    theme(
      #
      legend.position = c(0.2,0.25),#"right",#
      # legend.direction = "vertical",
      legend.background = element_blank(),
      legend.key = element_blank(),
      #axis.title = element_blank(),
      axis.text = element_text(size = 6)
    )) 




km_inset$plot

km_plot <- km_inset$plot +#scale_x_discrete(limtis=c(0,1.1))
  #scale_x_continuous(limtis=c(0,1.1),
  #  expand = c(0,0)
  #) +
  scale_y_continuous(
    limits = c(0.18,1.02),
    expand = c(0,0)
    
  ) +
  
  theme(
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    plot.margin = margin(     1,1,1,1)
  )
km_df <- surv_summary(km_fit)

last_pts <-
  
  km_df %>%
  #surv_summary(km_fit)%>%
  group_by(strata) %>%
  slice_max(time, n = 1)

km_plot <- km_plot +
  
  geom_text(
    data = last_pts,
    aes(
      x = time - 4,
      y = surv-0.01,
      label = c(
        "0–75 km",
        "75–150 km",
        "150–300 km",
        "300–500 km",
        ">500 km"
      )
      
    ),
    
    hjust = 0,
    size = 2.5
  )

library(cowplot)

Fig2B_withInset <-
  ggdraw() +
  draw_plot(
    
    FigB,
    x = 0,
    y = 0,
    width = 1,
    height = 1
  ) +
  
  draw_plot(
    km_plot,
    x = 0.3,
    y = 0.5,
    width = 0.6,
    height = 0.4
  )



Fig3_Main_AB_Forest_Splines <- ggarrange( Fig_Forest, Fig2B_withInset,
                                          labels = c("A", "B"),ncol=1)


Fig3_Main_AB_Forest_Splines




