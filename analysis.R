# Geostatistical data. Malaria in The Gambia

# 1 Data
library(tidyverse)
library(leaflet)
library(raster)
library(INLA)
library(rgdal)



nig <- read_csv("nigData.csv") 

d <- nig %>%
  sample_frac(., 0.1)

# d <- nig

coo <- cbind(d$lng, d$lat)

nigShp <- getData(name = 'GADM', country = 'NGA', level = 2)

r <- getData(name = 'alt', country = 'NGA', mask = TRUE)
ra <- aggregate(r, fact = 5, fun = mean)
dp <- rasterToPoints(ra)
coop <- dp[, c("x", "y")]




# 
# # mesh <- inla.mesh.2d(loc = coo, max.edge = c(0.1, 5), cutoff = 0.001)

mesh <- inla.mesh.2d(loc = coo, max.edge = 2, cutoff = 0.5)
# The number of vertices is given by mesh$n and we can plot the mesh with plot(mesh).
# 

mesh$n

plot(mesh)
points(coo, col = "red")


# 3.3 Build the SPDE model on the mesh
# Then, we use the inla.spde2.matern() function to build the SPDE model.


spde <- inla.spde2.matern(mesh = mesh, alpha = 2)


# 3.4 Index set
# Now we generate the index set for the SPDE model. We do this with the function inla.spde.make.index() where we specify the name of the effect (s) and the number of vertices in the SPDE model (spde$n.spde). This creates a list with vector s equal to 1:spde$n.spde, and vectors s.group and s.repl that have all elements equal to 1s and size given by the number of mesh vertices.


indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

# 3.5 Projector matrix
# We need to build a projector matrix A that projects the spatially continuous Gaussian random field at the mesh nodes. The projector matrix A can be built with the inla.spde.make.A() function passing the mesh and the coordinates.


A <- inla.spde.make.A(mesh = mesh, loc = coo)

Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
# 3.6 Prediction data


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


stk.v <- inla.stack(tag = "val",
                    data = list(OBS = NA),
                    A = list(1, A),
                    effects = list(X, s = indexs))

stk.p <- inla.stack(tag = "pred",
                    data = list(OBS = NA),
                    A = list(Ap),
                    effects = list(s=indexs))



stk.full <- inla.stack(stk.e, stk.v, stk.p)








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


# # 3.11 Results
# # We inspect the results by typing summary(res).


summary(res)




# https://rpubs.com/INBOstats/spde


# 
# 4 Mapping malaria prevalence
# Now we map the malaria prevalence predictions in a leaflet map. We can obtain the mean prevalence and lower and upper limits of the 95% credible intervals from res$summary.fitted.values where we need to specify the rows corresponding to the indices of the predictions and the columns "mean", "0.025quant" and "0.975quant". The indices of the stack stk.full that correspond to the predictions are the ones tagged with tag = "pred". We can obtain them by using inla.stack.index() and specifying tag = "pred".


index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]


r_prev_mean <- rasterize(x = coop, y = ra, field = prev_mean, fun = mean)

pal <- colorNumeric("inferno", c(0, 1), na.color = "transparent")

# pal <- colorNumeric(c("green", "blue", "red"), 0:1 , na.color = "transparent")

leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(r_prev_mean), title = "Prevalence") %>%
  addScaleBar(position = c("bottomleft"))


# 
ypred <- exp(res$summary.random$s$mean) / (1 + exp(res$summary.random$s$mean))
index.pred <- inla.stack.index(stack = stk.full, tag = "val")$data
proj.grid <- inla.mesh.projector(mesh, loc = coo, dims = c(300, 300))
projection <- inla.mesh.project(proj.grid, ypred)

ggplot_projection_shapefile(projection, proj.grid, nigShp)

coord.pred <- as_tibble(coo) %>%
  mutate(mean = projection) %>%
  rename(x = V1, y = V2)

coordinates(coord.pred) <- ~ x + y

utmproj <- "+proj=longlat +datum=WGS84"
proj4string(coord.pred) = utmproj
coord.pred2 <- spTransform(coord.pred, crs("+init=epsg:26392"))

rasT <- raster(extent(coord.pred2), ncols = 300, nrows = 300)

rasT <- raster(coord.pred2)

