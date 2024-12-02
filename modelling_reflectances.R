library(ggplot2)

setwd('D:/Reflectance summary')

raster_files  <-
  list.files(
    path = 'D:/Reflectance summary',
    pattern = 'compiled_spectra',
    recursive = TRUE,
    full.names = TRUE
  ) 
summary<-c('Filename','ControlDev', 'TreatedDev')
for (file in raster_files){
  dataset<-read.csv(file)
  standard_dev<-tapply(dataset$Spectra, dataset$Group, sd)
  row<- list(file, standard_dev[1], standard_dev[2])
  summary<-rbind(summary,row)
}
sum<-summary[1:19,]
setwd('D:/Reflectance summary')
write.csv(sum,'stdev.csv')

reflectance_summary<-read.csv('summaries.txt', sep='\t')

for (i in 1:nrow(reflectance_summary)){
  file_name<-reflectance_summary[i,]$File_name
  games_howell_test<-read.csv(file_name)
  reflectance_summary[i,]$Control.Mean<-games_howell_test$Control
  reflectance_summary[i,]$Treated.Mean<-games_howell_test$Treated
  reflectance_summary[i,]$p.value<-games_howell_test$P_Value
  i=i+1
}
Arsenic<-reflectance_summary[1:9,]
Arsenic$Metal<-'Arsenic'
Selenium<-reflectance_summary[10:18,]
Selenium$Metal<-'Selenium'

labelled_by_metal_grass_treatment<-rbind(Arsenic,Selenium)
write.csv(labelled_by_metal_grass_treatment,'overall_reflectance_summary.txt')

# For plotting overall reflectances
labelled_by_metal_grass_treatment$COMBO<-paste(labelled_by_metal_grass_treatment$Metal, labelled_by_metal_grass_treatment$Grass.Type)

ggplot(labelled_by_metal_grass_treatment, aes(x=Treatment.Concentration..ppm., y=Treated.Mean, color=COMBO)) +
  geom_point(size=3)+
  geom_errorbar(aes(ymin = Treated.Mean - Treated.St.Dev, ymax = Treated.Mean + Treated.St.Dev), 
                                    width = 0.2) +
  labs(x = "Treatment Concentration (ppm)", y = "Reflectance Mean Value") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Optional: remove legend title

# Scatter Plots
# ---Arsenic
ggplot(Arsenic, aes(x=Treatment.Concentration..ppm., y=Treated.Mean, color=Grass.Type)) +
  geom_point(size=3)+
  geom_errorbar(aes(ymin = Treated.Mean - Treated.St.Dev, ymax = Treated.Mean + Treated.St.Dev), 
                  width = 0.2) + 
  labs(x = "Concentration", y = "Mean Value") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Optional: remove legend title

# ---Selenium
ggplot(Selenium, aes(x=Treatment.Concentration..ppm., y=Treated.Mean, color=Grass.Type)) +
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Treated.Mean-Treated.St.Dev, ymax=Treated.Mean+Treated.St.Dev), width=.1) +
  labs(x = "Concentration", y = "Mean Value") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Optional: remove legend title

# Violin Plots (using grey bars only)
arsenic_restructured<-read.csv('SELENIUM.txt', sep='\t')
spectral_long <- arsenic_restructured %>%
  pivot_longer(cols = c(Switch.Grass, Perrenial.Rye.Grass, Bahia),
               names_to = "Grass.Type",
               values_to = "Reflectance")

# ggplot(spectral_long, aes(x = factor(Treatment.Concentration..ppm.), y = Reflectance, fill = Grass.Type)) +
#   geom_violin(trim = FALSE, alpha = 0.4) +  # Violin plot with some transparency
#   geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.6) +  # Boxplot inside violin
#   stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(0.9)) +  # Mean as black dots
#   labs(x = "Treatment Concentration (ppm)", y = "Reflectance") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

ggplot(spectral_long, aes(x = factor(Treatment.Concentration..ppm.), y = Reflectance, fill = Grass.Type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with some transparency
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.8, color = "black") +  # Boxplot inside violin
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(0.9)) +  # Mean as black dots
  scale_fill_grey(start = 0.3, end = 0.9) +  # Shades of gray for fill color
  labs(x = "Treatment Concentration (ppm)", y = "Reflectance") +
  theme_minimal() +
  theme(legend.title = element_blank())


arsenic_restructured<-read.csv('ARSENIC.txt', sep='\t')
spectral_long <- arsenic_restructured %>%
  pivot_longer(cols = c(Switch.Grass, Perrenial.Rye.Grass, Bahia),
               names_to = "Grass.Type",
               values_to = "Reflectance")

# ggplot(spectral_long, aes(x = factor(Treatment.Concentration..ppm.), y = Reflectance, fill = Grass.Type)) +
#   geom_violin(trim = FALSE, alpha = 0.4) +  # Violin plot with some transparency
#   geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.6) +  # Boxplot inside violin
#   stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(0.9)) +  # Mean as black dots
#   labs(x = "Treatment Concentration (ppm)", y = "Reflectance") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

ggplot(spectral_long, aes(x = factor(Treatment.Concentration..ppm.), y = Reflectance, fill = Grass.Type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with some transparency
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.8, color = "black") +  # Boxplot inside violin
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(0.9)) +  # Mean as black dots
  scale_fill_grey(start = 0.3, end = 0.9) +  # Shades of gray for fill color
  labs(x = "Treatment Concentration (ppm)", y = "Reflectance") +
  theme_minimal() +
  theme(legend.title = element_blank())

