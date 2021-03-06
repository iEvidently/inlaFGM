# Geostatistical data. FGM in The Gambia

##1 Load libraries
library(tidyverse)
library(leaflet)
library(raster)
library(INLA)
library(rgdal)


##2 Data preparation
nig <- read_csv("https://raw.githubusercontent.com/iEvidently/inlaFGM/master/nigData.csv") 

d <- nig %>%
  sample_frac(., 0.1)


coo <- cbind(d$lng, d$lat)

nigShp <- getData(name = 'GADM', country = 'NGA', level = 2)

r <- getData(name = 'alt', country = 'NGA', mask = TRUE)
ra <- aggregate(r, fact = 25, fun = mean)
dp <- rasterToPoints(ra)
coop <- dp[, c("x", "y")]




##3 Build mesh

mesh <- inla.mesh.2d(loc = coo, max.edge = 2, cutoff = 0.5)

mesh$n

plot(mesh)
points(coo, col = "red")



##4 Build the SPDE model on the mesh
spde <- inla.spde2.matern(mesh = mesh, alpha = 2)


indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)


A <- inla.spde.make.A(mesh = mesh, loc = coo)

Ap <- inla.spde.make.A(mesh = mesh, loc = coop)


## 5 Stack data for the estimation and prediction
X <- data.frame(b0 = 1,
                age = d[, "age"], 
                muslim = d[, "muslim"],
                yoruba = d[, "yoruba"],
                igbo = d[, "igbo"],
                nc =d[, "nc"], 
                nw =d[, "nw"],
                se =d[, "se"], 
                ss =d[, "ss"],
                sw =d[, "sw"])

stk.e <- inla.stack(tag = "est",
                    data = list(OBS = d$pos),
                    A = list(1, A),
                    effects = list(X, s = indexs))


stk.p <- inla.stack(tag = "pred",
                    data = list(OBS = NA),
                    A = list(Ap),
                    effects = list(s=indexs))



stk.full <- inla.stack(stk.e, stk.p)






##6 Model formula

formula <- OBS ~ 0 + b0 + age + muslim + yoruba + igbo + nc + nw + se + ss + sw +  f(s, model = spde)

 
res = inla(formula,
           data = inla.stack.data(stk.full), 
           family = "binomial", 
           verbose = TRUE,
           control.fixed = list(prec = 0.1, prec.intercept = 0.1), 
           control.predictor = list(
             A = inla.stack.A(stk.full), 
             compute = TRUE, 
             link = 1), 
           control.family = list(link="logit"),
           control.results = list(return.marginals.random = TRUE,
                                  return.marginals.predictor = TRUE),
           control.compute=list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE)) 



##7 Results

summary(res)




##8 Mapping FGM prevalence

index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]


r_prev_mean <- rasterize(x = coop, y = ra, field = prev_mean, fun = mean)

pal <- colorNumeric("inferno", c(0, 1), na.color = "transparent")


leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_mean), title = "Prevalence") %>%
  addScaleBar(position = c("bottomleft"))

