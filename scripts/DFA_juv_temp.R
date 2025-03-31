## Purpose: To construct a Dynamic Factor Analysis model of juvenile Tanner crab
##          temperature indicators. 
##
## NOTES:
## - unconstrained 1- and 2-trend models are best fit -- go with 1-trend for simplicity?
## - not putting summer surface temperature in there because not bottom, and larval


## Load libraries
library(MARSS)
library(corrplot)
library(oce)


## Read in setup
source("./scripts/setup.R")


## Process data
# Mean summer bottom temperature, summer cold pool extent
env <- read_csv("./outputs/temp_coldpool.csv")

# Summer juvenile Tanner temperature occupancy
occ <- read_csv("./outputs/tanner_temp_occupied.csv")


dat <- full_join(env, occ) %>%
       rename(temp_occ = TEMP_OCC,
              year = YEAR) %>%
       select(-summer_st)


# plot
dat %>%
  pivot_longer(c(2:4), names_to = "name", values_to = "value") %>%
  ggplot(aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab("Value")
ggsave("./figures/juv_temp_indicatorsTS.png", 
       width = 12, height = 5.5)


dfa_dat <- dat %>%
           select(-year) %>%
           t()

colnames(dfa_dat) <- years


# plot correlations
cors <- cor(t(dfa_dat), use = "p")
diag(cors) <- 0

max(cors)
min(cors) 

plot <- as.data.frame(t(dfa_dat))

# plot correlations
corrplot(cors, method = "sq", col.lim = c(-1, 1), 
         col = oceColorsPalette(64), tl.col = "black", 
         cl.cex = 0.7, order = "FPC")


# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# changing convergence criterion to ensure convergence
cntl.list = list(minit = 200, 
                 maxit = 20000, 
                 allow.degen = FALSE, 
                 conv.test.slope.tol = 0.1, 
                 abstol = 0.0001)

# fit models & store results
for(R in levels.R){
  for(m in 1:2){  # considering either 1- or 2-trend model
    
    dfa.model <- list(A = "zero", 
                      R = R, 
                      m = m)
    
    kemz <- MARSS(dfa_dat, 
                  model = dfa.model,
                  form = "dfa", 
                  z.score = TRUE, 
                  control = cntl.list)
    
    model.data <- rbind(model.data,
                        data.frame(R = R,
                                   m = m,
                                   logLik = kemz$logLik,
                                   K = kemz$num.params,
                                   AICc = kemz$AICc,
                                   stringsAsFactors = FALSE))
    
    assign(paste("kemz", m, R, sep = "."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data <- model.data %>%
              mutate(dAICc = AICc - min(AICc)) %>% 
              arrange(dAICc)
model.data # unconstrained is the best model and does converge

# save model selection table
write.csv(model.data, "./outputs/juv_temp_dfa_model_selection.csv",
          row.names = F)


## fit the best model --------------------------------------------------
model.list <- list(A = "zero", 
                   m = 2, 
                   R = "unconstrained") # best model is two-trend

mod <- MARSS(dfa_dat, 
             model = model.list, 
             z.score = TRUE, 
             form = "dfa", 
             control = cntl.list)

# rotate
# get the inverse of the rotation matrix
Z.est <- coef(mod, type = "matrix")$Z

H.inv <- varimax(coef(mod, type = "matrix")$Z)$rotmat

# rotate factor loadings
Z.rot <- Z.est %*% H.inv

# rotate trends
trends.rot <- solve(H.inv) %*% mod$states

# Add CIs to marssMLE object
mod <- MARSSparamCIs(mod)

# Use coef() to get the upper and lower CIs
Z.low <- coef(mod, type = "Z", what = "par.lowCI")
Z.up <- coef(mod, type = "Z", what = "par.upCI")
Z.rot.up <- Z.up %*% H.inv
Z.rot.low <- Z.low %*% H.inv

plot.CI <- data.frame(names = rownames(dfa_dat),
                      mean = as.vector(Z.rot),
                      upCI = as.vector(Z.rot.up),
                      lowCI = as.vector(Z.rot.low))
plot.CI


plot.CI$names <- reorder(plot.CI$names, plot.CI$mean)
plot.CI$trend <- rep(c("T1", "T2"), each = 3)

loadings.plot <- ggplot(plot.CI, aes(x = names, y = mean, fill = trend)) +
                 geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
                 geom_errorbar(aes(ymax = upCI, ymin = lowCI), 
                               position = position_dodge(width = 0.9), width = 0.5) +
                 labs(x = "", y = "Loading") +
                 theme_bw() +
                 theme(legend.title = element_blank(),
                       legend.position = "none") +
                 geom_hline(yintercept = 0)

# plot trend
trend <- data.frame(trend = rep(c("T1", "T2"), each = length(years)),
                    t = years,
                    estimate = as.vector(mod$states),
                    conf.low = as.vector(mod$states) - 1.96*as.vector(mod$states.se),
                    conf.high = as.vector(mod$states) + 1.96*as.vector(mod$states.se))


trend.plot <- ggplot(trend, aes(t, estimate, color = trend, fill = trend)) +
              geom_point() +
              geom_line() +
              geom_ribbon(aes(x = t, ymin = conf.low, ymax = conf.high), 
                          linetype = 0, alpha = 0.1) + 
              geom_hline(yintercept = 0) +
              labs(x = "", y = "Trend") +
              theme_bw()

ggpubr::ggarrange(loadings.plot,
                  trend.plot,
                  ncol = 2,
                  widths = c(0.45, 0.55))
ggsave("./figures/juv_temp_2trend_DFA_loadings_trend.png", 
       width = 12, height = 5.5)


# fit second-best model (1 trend unconstrained)
model.list <- list(A = "zero", 
                   m = 1, 
                   R = "unconstrained") # second-best model - this is the borealization index

mod <- MARSS(dfa_dat, 
             model = model.list, 
             z.score = TRUE, 
             form = "dfa", 
             control = cntl.list)

# # save 
# saveRDS(mod, "./outputs/juv_temp_DFA_model.rds")

# plot fits to data
DFA_pred <- as.data.frame(print(predict(mod))) %>%
            mutate(year = rep(years, 3)) 

DFA_pred <- DFA_pred %>%
            group_by(.rownames) %>%
            # get R^2 for each time series
            summarise(R_sq = cor(y, estimate, use = "pairwise")^2) %>%
            mutate(plot_label = paste(.rownames, " (", round(R_sq, 3), ")", sep = "")) %>%
            right_join(., DFA_pred)

# plot
ggplot(DFA_pred, aes(x = estimate, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~plot_label, ncol = 4, scale = "free") +
  labs(x = "Estimated", y = "Observed") +
  theme_bw()


# process loadings and trend
CI <- MARSSparamCIs(mod)

plot.CI <- data.frame(names = rownames(dfa_dat),
                      mean = CI$par$Z[1:3],
                      upCI = CI$par.upCI$Z[1:3],
                      lowCI = CI$par.lowCI$Z[1:3])
plot.CI$names <- reorder(plot.CI$names, CI$par$Z[1:3])

loadings.plot <- ggplot(plot.CI, aes(x = names, y = mean)) +
                 geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
                 geom_errorbar(aes(ymax = upCI, ymin = lowCI), 
                               position = position_dodge(width = 0.9), width = 0.5) +
                 labs(x = "", y = "Loading") +
                 theme_bw() +
                 theme(legend.title = element_blank(),
                       legend.position = "none") +
                 geom_hline(yintercept = 0)



# plot trend
trend <- data.frame(t = years,
                    estimate = as.vector(mod$states),
                    conf.low = as.vector(mod$states)-1.96*as.vector(mod$states.se),
                    conf.high = as.vector(mod$states)+1.96*as.vector(mod$states.se))

trend.plot <- ggplot(trend, aes(t, estimate)) +
              geom_point() +
              geom_line() +
              geom_ribbon(aes(x = t, ymin = conf.low, ymax = conf.high), 
                          linetype = 0, alpha = 0.1) + 
              geom_hline(yintercept = 0) +
              labs(x = "", y = "Trend") +
              theme_bw()



ggpubr::ggarrange(loadings.plot,
                  trend.plot,
                  ncol = 2,
                  widths = c(0.45, 0.55))
ggsave("./figures/juv_temp_1trend_DFA_loadings_trend.png", 
       width = 12, height = 5.5)


# and save loadings and trend
write.csv(plot.CI, "./outputs/juv_temp_dfa_loadings.csv", row.names = F)

trend %>%
  rename(year = t,
         dfa_temp = estimate) %>%
  write.csv("./outputs/juv_temp_dfa_trend.csv", row.names = F)




