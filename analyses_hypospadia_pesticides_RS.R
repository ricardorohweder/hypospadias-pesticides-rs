library(dplyr)
library(sf)
library(spdep)
library(INLA)
library(leaflet)
library(RColorBrewer)

# df with the number of hypospadias cases and live births per municipality
df.hipnv

# calculate the overall ratio of cases per live births
SMR_total <- df.hipnv %>%
  summarize(nascidos = sum(lb), casos = sum(hip), SMR_total = casos/nascidos)

# sum data from all years by municipality
banco2 = df.hipnv %>%
  group_by(mun) %>%
  summarise(across(c(2:3), sum, .names = "soma_{.col}"))

# calculate the expected number of cases per municipality
fator <- as.numeric(SMR_total[3])
banco2 <- banco2 %>%
  mutate(E = soma_lb * fator)

# create a copy for adjusted analysis (metropolitan municipalities with higher than expected observed cases)
banco1 <- banco2

# create a dummy case value for the selected municipalities
banco2_filtrado <- banco2 %>%
  filter(mun %in% c("VIAMAO", "PORTO ALEGRE", "ALVORADA", "CACHOEIRINHA", "GRAVATAI")) %>%
  mutate(soma_hip = round(E, 0))

# update the dataframe banco2 with the filtered and modified rows (adjusted observed numbers)
banco2 <- banco2 %>%
  left_join(banco2_filtrado %>% select(mun, soma_hip), by = "mun", suffix = c("", "_new")) %>%
  mutate(soma_hip = coalesce(soma_hip_new, soma_hip)) %>%
  select(-soma_hip_new)

# create a new dataframe with the codes for each municipality according to the map (main analysis)
d <- as.data.frame(banco2)

# df without adjustment for Viamao, Porto Alegre, Alvorada, Cachoeirinha, gravatai (secondary)
#d <- as.data.frame(banco1)

# calculate SIR (standardized prevalence ratio)
d$SIR = (d$soma_hip/d$E)

# load the map of Rio Grande do Sul (RS)
map = st_read(dsn = file.path("RS_Municipios_2022.shp"), quiet =
                TRUE)

# remove lakes (total municipalities = 497)
map <- map[3:499,]

# standardize municipality names in the column mun
map$mun <- toupper(iconv(map$NM_MUN, to = "ASCII//TRANSLIT"))

# ensure they are in the same order 
d <- d[order(d$mun), ]
map <- map[order(map$mun), ]
all(d$mun == map$mun)

########## BYM2 model

# create a list of neighboring municipalities
nb = poly2nb(map)

# save the list
nb2INLA('maprs', nb)

# load the list of neighboring municipalities
g = inla.read.graph(filename = "maprs")

# indices for the areas
d = d %>% 
  mutate(id=row_number())

# model
f1 = soma_hip ~ f(id, model = "bym2", graph = g)

# fit
m1 = inla(f1, family = "poisson", data = d, E = E,
          control.predictor = list(compute = TRUE))

# model results
summary(m1)

# mapping the relative risk
d$RR <- m1$summary.fitted.values[, "mode"]
d$LL <- m1$summary.fitted.values[, "0.025quant"]
d$UL <- m1$summary.fitted.values[, "0.975quant"]

# mapping the relative risk, only where it is different from 1
if ("0.025quant" %in% colnames(m1$summary.fitted.values) &&
    "0.975quant" %in% colnames(m1$summary.fitted.values)) {
  
  d$RR_sig <- ifelse(m1$summary.fitted.values$`0.025quant` > 1 | m1$summary.fitted.values$`0.975quant` < 1, 
                     m1$summary.fitted.values$mean, 1)
}

# check the relative risk in significant municipalities
d[which(!d$RR_sig==1),c(1,7,8,9)]

# combine relative risk estimates with geographic information to plot the map
map_sf = map %>% 
  left_join(d, by = "mun")

# Transforming the coordinate reference system
map_sf <- st_transform(map_sf, crs = 4326)

# creating a color scale to represent the relative risk
pal <- colorNumeric(palette = "YlOrRd", domain = map_sf$RR)

# creating labels for map
labels <- sprintf("<strong> %s </strong> <br/>
  Observed: %s <br/> Expected: %s <br/>
  SIR: %s <br/> RR: %s (%s, %s)",
                  map_sf$mun, map_sf$soma_hip, round(map_sf$E, 2),
                  round(map_sf$SIR, 2), round(map_sf$RR, 2),
                  round(map_sf$LL, 2), round(map_sf$UL, 2)
) %>% lapply(htmltools::HTML)

#plot the map
lRR <- leaflet(map_sf) %>%
  addTiles() %>%
  addPolygons(
    color = "grey", weight = 1, fillColor = ~ pal(RR),
    fillOpacity = 0.5,
    highlightOptions = highlightOptions(weight = 4),
    label = labels,
    labelOptions = labelOptions(
      style =
        list(
          "font-weight" = "normal",
          padding = "3px 8px"
        ),
      textsize = "15px", direction = "auto"
    )
  )  %>% 
  addScaleBar(position = "bottomleft", options = scaleBarOptions(metric = TRUE, imperial = TRUE)) %>%
  addLegend(
    pal = pal, values = ~RR, opacity = 0.5, title = "Relative Risk",
    position = "bottomleft"
  ) %>%
  # Adicionando o minimapa com camadas de países
  addMiniMap(
    tiles = providers$Esri.WorldStreetMap, 
    position = "bottomright",
    width = 200, height = 250,
    zoomLevelOffset = -5
  ) %>%
  addControl(
    html = '<img src="https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcTHTXjD_bvVt6a1jXaVP0iONWro2rUzHpCxMw&s" style="width:20px;">',
    position = "topright"
  )

# view the map
lRR

#############
# assess the association between exposure (relative agricultural area) and hypospadias

# df with exposure metrics (relative area dedicated to agriculture crops)
crops

# create a vector with exposure variables
agriculture.names <- colnames(crops)[2:48]

# df to store the estimates
df.res.crops <- data.frame(
  agriculture = NA,
  beta = NA,
  betasd = NA,
  beta025 = NA,
  beta975 = NA,
  betamode = NA,
  precis = NA,
  precissd = NA,
  precis025 = NA,
  precis975 = NA,
  precismode = NA,
  phi = NA,
  phisd = NA,
  phi025 = NA,
  phi975 = NA,
  phimode = NA
)

# merge the cases database with exposure database
d_crops <-  merge(d, crops, by = "mun", all.x = TRUE)
all(d_crops$mun == map$mun)

# use a loop to estimate the association for each variable individually
for (var in agriculture.names) {
  # create the formula by substituting the variable
  f1 <- as.formula(paste("soma_hip ~ 1 +", var, "+ f(id, model = 'bym2', graph = g)"))
  
  # fit the model
  m1 <- inla(f1, family = "poisson", data = d_crops, E = E,
             control.predictor = list(compute = TRUE))
  
  # store estimates
  df.res.crops[var,2:6] <- m1$summary.fixed[2,c(1,2,3,5,6)]
  df.res.crops[var,7:11] <- m1$summary.hyperpar[1,c(1,2,3,5,6)]
  df.res.crops[var,12:16] <- m1$summary.hyperpar[2,c(1,2,3,5,6)]
}

# check estimates
df.res.crops

### map the distribution of cultivation accross municipalities
md_crops = map %>% 
  left_join(d_crops, by = "mun")

# Transforming the coordinate reference system
md_crops <- st_transform(md_crops, crs = 4326)

# Calculate the quartiles for the values of tobacco divided by 13years (or other crop)
md_crops$cropyearly <- md_crops$Fumo/13

quartis <- quantile(md_crops$cropyearly[md_crops$cropyearly > 0], probs = c(0.25, 0.5, 0.75, 1))

# Divide the Tobacco variable into 5 ordinal categories
md_crops$cultivation <- cut(
  md_crops$cropyearly,
  breaks = c(-Inf, 0, quartis),
  labels = c("No cultivation", 
             sprintf("%.3f - %.3f", min(md_crops$cropyearly[md_crops$cropyearly > 0]), quartis[1]),
             sprintf("%.3f - %.3f", quartis[1], quartis[2]),
             sprintf("%.3f - %.3f", quartis[2], quartis[3]),
             sprintf("%.3f - %.3f", quartis[3], quartis[4])),
  include.lowest = TRUE
)

# Define the YlOrRd color palette for the 5 categories
pal <- colorFactor(
  palette = brewer.pal(5, "YlOrRd"),
  domain = md_crops$cultivation
)

# Adjust the labels to reflect the categories and other values
labels <- sprintf("<strong> %s </strong> <br/>
  Observed: %s <br/> Expected: %s <br/>
  SIR: %s <br/> Cultivation (relative area): %s <br/> Relative Risk: %s (%s, %s)",
                  md_crops$mun, md_crops$soma_hip, round(md_crops$E, 2),
                  round(md_crops$SIR, 2), round(md_crops$crop, 6),
                  round(md_crops$RR, 2),
                  round(md_crops$LL, 2), round(md_crops$UL, 2)
) %>% lapply(htmltools::HTML)

# Create the map with the five  categories
map_cultivation <- leaflet(md_crops) %>%
  addTiles() %>%
  addPolygons(
    color = "grey", weight = 1, fillColor = ~pal(cultivation),
    fillOpacity = 0.8,
    highlightOptions = highlightOptions(weight = 4),
    label = labels,
    labelOptions = labelOptions(
      style =
        list(
          "font-weight" = "normal",
          padding = "3px 8px"
        ),
      textsize = "15px", direction = "auto"
    )
  )  %>% 
  addScaleBar(position = "bottomleft", options = scaleBarOptions(metric = TRUE, imperial = TRUE)) %>%
  addLegend(
    pal = pal, values = ~cultivation, opacity = 0.7, title = "Crop Cultivation<br/>
    (% area)",
    position = "bottomleft"
  ) %>%
  addControl(
    html = '<img src="https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcTHTXjD_bvVt6a1jXaVP0iONWro2rUzHpCxMw&s" style="width:20px;">',
    position = "topright"
  )

map_cultivation

#############
# assess the association between exposure (pesticides sold) and hypospadias

# df with exposure metrics (amount of pesticides sold per km2)
pest.sold

# create a vector with exposure variables
pest.names <- colnames(pest.sold)[2:353]
pest.names <- colnames(pest.sold)[200]

# df to store the estimates
df.res.pestsold <- data.frame(
  pesticide = NA,
  beta = NA,
  betasd = NA,
  beta025 = NA,
  beta975 = NA,
  betamode = NA,
  precis = NA,
  precissd = NA,
  precis025 = NA,
  precis975 = NA,
  precismode = NA,
  phi = NA,
  phisd = NA,
  phi025 = NA,
  phi975 = NA,
  phimode = NA
)

# merge the cases database with exposure database
d_pestsold <-  merge(d, pest.sold, by = "mun", all.x = TRUE)
all(d_pestsold$mun == map$mun)

# use a loop to estimate the association for each variable individually
for (var in pest.names) {
  # create the formula by substituting the variable
  f1 <- as.formula(paste("soma_hip ~ 1 +", var, "+ f(id, model = 'bym2', graph = g)"))
  
  # fit the model
  m1 <- inla(f1, family = "poisson", data = d_pestsold, E = E,
             control.predictor = list(compute = TRUE))
  
  # store estimates
  df.res.pestsold[var,2:6] <- m1$summary.fixed[2,c(1,2,3,5,6)]
  df.res.pestsold[var,7:11] <- m1$summary.hyperpar[1,c(1,2,3,5,6)]
  df.res.pestsold[var,12:16] <- m1$summary.hyperpar[2,c(1,2,3,5,6)]
}

# check estimates
df.res.pestsold


#############
# assess the association between exposure (probability of detection pesticides) and hypospadias

# df with exposure metrics (amount of pesticides sold per km2)
pest.prob

# create a vector with exposure variables
pestp.names <- colnames(pest.prob)[2:103]

# df to store the estimates
df.res.pestprob <- data.frame(
  pesticide = NA,
  beta = NA,
  betasd = NA,
  beta025 = NA,
  beta975 = NA,
  betamode = NA,
  precis = NA,
  precissd = NA,
  precis025 = NA,
  precis975 = NA,
  precismode = NA,
  phi = NA,
  phisd = NA,
  phi025 = NA,
  phi975 = NA,
  phimode = NA
)

# merge the cases database with exposure database
d_pestprob <-  merge(d, pest.prob, by = "mun", all.x = TRUE)
all(d_pestprob$mun == map$mun)

# use a loop to estimate the association for each variable individually
for (var in pestp.names) {
  # create the formula by substituting the variable
  f1 <- as.formula(paste("soma_hip ~ 1 +", var, "+ f(id, model = 'bym2', graph = g)"))
  
  # fit the model
  m1 <- inla(f1, family = "poisson", data = d_pestprob, E = E,
             control.predictor = list(compute = TRUE))
  
  # store estimates
  df.res.pestprob[var,2:6] <- m1$summary.fixed[2,c(1,2,3,5,6)]
  df.res.pestprob[var,7:11] <- m1$summary.hyperpar[1,c(1,2,3,5,6)]
  df.res.pestprob[var,12:16] <- m1$summary.hyperpar[2,c(1,2,3,5,6)]
}

# check estimates
df.res.pestprob
