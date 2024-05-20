#Visualizations/summaries of model results

#after running PrepData.R and Model.R
load("boston-res.RData")
load("boston-mcmc.RData")


# [ Paper Figure 6]

#finding tau

#count
rank <- c(rep(K, ncol(preds) - 2), 1, 10)
jpos <- 0
data <- data.frame(intercept = 1:200)
for (j in seq_along(rank)) {
  data[,j+1] <- cumsum(exp(res_em$logptau[jpos + 1:rank[j]]))
  jpos <- jpos + rank[j]
}

#indicator
jpos <- 0
data2 <- data.frame(intercept = 1:200)
for (j in seq_along(rank)) {
  data2[,j+1] <- cumsum(exp(res_em$logptau_hz[jpos + 1:rank[j]]))
  jpos <- jpos + rank[j]
}

data[11:200, 6] <- 1 #college binary indcator
data2[11:200, 6] <- 1

df <- data %>%
  select(intercept, V2, V3,V4, V6) %>%
  gather(key = "covariates", value = "value", -intercept)

df2 <- data2 %>%
  select(intercept, V2, V3,V4, V6) %>%
  gather(key = "covariates", value = "value", -intercept)

p1 <- ggplot(df, aes(x = intercept, y = value)) + 
  geom_line(aes(color = covariates), show.legend = FALSE) +
  xlab("Basis Rank Expansion") + 
  ylab("Cumulative Probability") +
  labs(title = "Count Level") +
  geom_hline(yintercept=0.5)

p2 <- ggplot(df2, aes(x = intercept, y = value)) + 
  geom_line(aes(color = covariates)) +
  xlab("Basis Rank Expansion") + 
  ylab("") +
  scale_color_hue(labels = c("intercept", "tax", "distance to \n police", "college \n housing")) +
  labs(title = "Bernoulli Level")+
  geom_hline(yintercept = 0.5)

grid.arrange(p1, p2, ncol = 2)

# [ Paper Figure 7]

# binned deviance residual plot
p_bg <- plogis(X_hz %*% fit_em$coef_hz)
mu <- drop(p_bg * exp(X_bg %*% fit_em$coef_bg) + (1 - p_bg) * exp(X %*% fit_em$coef))

# cluster by predictors
nc <- 500
ii <- factor(cutree(hclust(dist(X[, -1]), method = "ward.D"), nc))

wf <- tapply(y, ii, length)
yf <- tapply(y, ii, sum) / wf # average response

p_bg <- drop((1 / as.vector(wf)) * rowsum(fit_em$p_bg, ii)) # average pi

xf <- (1 / as.vector(wf)) * rowsum(X_bg, ii) # average preds (bg)
mf0 <- drop(exp(xf %*% fit_em$coef_bg))

xf <- (1 / as.vector(wf)) * rowsum(X, ii) # average preds (main)
mf1 <- drop(exp(xf %*% fit_em$coef))

mf <- p_bg * mf0 + (1 - p_bg) * mf1 # E[y]
sf <- sqrt(mf + p_bg * (1 - p_bg) * (mf0 - mf1) ^ 2) # SD[y]
rf <- sqrt(wf) * (yf - mf) / sf # Pearson residuals

rfn <- data_frame(mf, rf)
op <- par(mfrow = c(1, 2))
pbg <- data_frame(y, fit_em$p_bg)

g0 <- ggplot(data = pbg, aes(x = 1-fit_em$p_bg, y = y)) +
  geom_jitter(width = 0) +
  labs(x = "Indicator Probability", y = "Observed Counts, slightly jittered")+
  theme(axis.title.x = element_text(size = 14),  # Adjust the size as needed
        axis.title.y = element_text(size = 14))

g <- ggplot(data = rfn, aes(x = mf, y = rf))+
  geom_point() +
  labs(x = "Pooled Fitted Values", y = "Pooled Pearson Residuals")+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_smooth(se = FALSE) +
  theme(axis.title.x = element_text(size = 14),  # Adjust the size as needed
        axis.title.y = element_text(size = 14))


g1 <- ggplot(data = rfn, aes(sample=rf)) +
  geom_qq() +
  geom_qq_line() + 
  labs(x = "Theoretical Quanitles", y = "Sample Quantiles") + 
  theme(axis.title.x = element_text(size = 14),  # Adjust the size as needed
        axis.title.y = element_text(size = 14))

g2 <- grid.arrange(g0, g, g1, ncol = 3)

g2 <- arrangeGrob(g0, g, g1, ncol=3) #generates g

# [ Paper Figure 8: Need shape files from PrepData ]

# [ Which streets should we marginalize out of the network? ] 
keep <- street_buffer %>%
  filter(lengths(st_intersects(., res_parcels)) > 0)
rem <- dplyr::symdiff(V(gl)$id, keep$OBJECTID) # These are the ids, not the location
rem2 <- which(V(gl)$id %in% rem == TRUE)

shape2 <- shape %>% filter(OBJECTID %in% V(gl)$id)
shape2 <- shape2 %>% filter(!OBJECTID %in% rem)

ii <- c(1:tau[1]) # intercept 
shape2$intercept <- drop(e.vectors[, ii] %*% fit_em$coef[ii])
shape2$logtax_effect <- drop(e.vectors[, 1:tau[2], drop = FALSE] %*%
                               fit_em$coef[tau[1] + 1:tau[2]])
shape2$pol_effect <- drop(e.vectors[, 1:tau[3], drop = FALSE] %*%
                            fit_em$coef[1+tau[1]+tau[2] + 1:tau[3]])
shape2$hz_prob <- 1-drop(plogis(-X_hz %*% fit_em$coef_hz))
shape2$indicator_prob <- 1 - res_em$p_bg

breaks_hz <- c(0,0.25, .5, .75, 1)
breaks_po <- c(-0.038, -0.02, -.01, 0, .04)
breaks_hz <- c(0,0.05, .09, .3,.8, 1)


my_pal1 <- c('gray80',"#C7E9C0","#74C476" ,"#238B45","#00441B")
my_pal2 <- c('darkgray',"lightgray",'yellow','darkorange','red')

t <- tmap_arrange(
  tm_shape(shape2) + 
    tm_lines(col = "intercept", palette = "-RdYlBu", size = 0.05, style = "quantile", midpoint = NA, lwd = 1.5),
  tm_shape(shape2) + 
    tm_lines(col = "logtax_effect", palette = my_pal1, size = 0.05, style = "quantile", lwd = 1.5),
  tm_shape(shape2) +
    tm_lines(col = "pol_effect", palette = "-Blues", size = 0.05, style = "quantile", lwd = 1.5) +
    tm_shape(police) + tm_symbols(col = "yellow", size = .25),
  tm_shape(shape2) + 
    tm_lines(col = "indicator_prob", palette = my_pal2, size = 0.05, breaks = breaks_hz, lwd = 1.5),
  nrow = 2
)

shape3 <- shape2 %>%
  filter(indicator_prob > .5)

tmap_mode("plot")
  
tm_shape(shape3) + 
    tm_lines(col = "intercept", palette = "-RdYlBu", size = 0.05, style = "quantile", midpoint = NA, lwd = 5) +
    tm_layout(legend.show = FALSE, frame.lwd = 0)
  
tm_shape(shape3) + 
    tm_lines(col = "logtax_effect", palette = my_pal1, size = 0.05, style = "quantile", lwd = 5) +
    tm_layout(legend.show = FALSE, frame.lwd = 0)
  
tm_shape(shape3) +
    tm_lines(col = "pol_effect", palette = "-Blues", size = 0.05, style = "quantile", lwd = 5) +
    tm_layout(legend.show = FALSE, frame.lwd = 0) +
    tm_shape(police) + tm_symbols(col = "yellow", size = .25)




# [ Paper Figure 9 ]

iii <- c(1:tau_hz[1]) # intercept 
shape2$intercept_hz <- drop(e.vectors[, iii] %*% fit_em$coef_hz[iii])
shape2$logtax_effect_hz <- drop(e.vectors[, 1:tau_hz[2], drop = FALSE] %*%
                                  fit_em$coef_hz[tau_hz[1] + 1:tau_hz[2]])
shape2$pol_effect_hz <- drop(e.vectors[, 1:tau_hz[3], drop = FALSE] %*%
                               fit_em$coef_hz[1+tau_hz[1]+tau_hz[2] + 1:tau_hz[3]])

tmap_arrange(
  tm_shape(shape2) + 
    tm_lines(col = "intercept_hz", palette = "RdYlBu", size = 0.05, style = "quantile", midpoint = NA, lwd = 1.5),
  tm_shape(shape2) + 
    tm_lines(col = "logtax_effect_hz", palette = rev(my_pal1), size = 0.05, style = "quantile", lwd = 1.5),
  tm_shape(shape2) + 
    tm_lines(col = "pol_effect_hz", palette = rev("-Blues"), size = 0.05, style = "quantile", lwd = 1.5) +
    tm_shape(police) + tm_symbols(col = "yellow", size = .25),
  nrow = 2
)


# [ MCMC exploration ]

burn_in <- .2; thin <- 10
n <- nrow(X) 
nt <- ncol(X); nz <- ncol(X_bg); nw <- ncol(X_hz)
valid <- (control$nsamples * burn_in):control$nsamples # discard burn-in
nsamp <- 12500 
valid <- (nsamp*burn_in):nsamp
valid <- valid[valid %% thin == 0] # thin


for (chain in 1:control$nchain) {
  theta_s[[chain]] <- tcrossprod(X, sims[valid, chain, 1:nt])
  theta_f <- X %*% fit_em$coef[1:nt]
  out[[chain]] <- ggplot(data_frame(as.vector(theta_s[[chain]])), aes(as.vector(theta_s[[chain]])))+
    geom_histogram(bins = 25, fill = "gray", color = "black")+
    geom_vline(xintercept = mean(theta_f), lwd = 1.5)+
    labs(title = paste("chain",chain), y ="", x = "theta")
  
  zeta_s[[chain]] <- tcrossprod(X_bg, sims[valid, chain, nt + 1:nz])
  zeta_f <- X_bg %*% fit_em$coef_bg
  out1[[chain]] <- ggplot(data_frame(as.vector(zeta_s[[chain]])), aes(as.vector(zeta_s[[chain]])))+
    geom_histogram(bins = 25, fill = "gray", color = "black")+
    geom_vline(xintercept = mean(zeta_f), lwd = 1.5)+
    labs(title = paste("chain",chain), y ="", x = "zeta")
  
  omega_s[[chain]] <- tcrossprod(X_hz, sims[valid, chain, nt + nz + 1:nw])
  omega_f <- X_hz %*% fit_em$coef_hz
  out2[[chain]] <- ggplot(data_frame(as.vector(omega_s[[chain]])), aes(as.vector(omega_s[[chain]])))+
    geom_histogram(bins = 25, fill = "gray", color = "black")+
    geom_vline(xintercept = mean(omega_f), lwd = 1.5)+
    labs(title = paste("chain",chain), y ="", x = "omega")
}

#four chains
grid.arrange(out[[1]], out1[[1]],out2[[1]],
             out[[2]], out1[[2]], out2[[2]],
             out[[3]],out1[[3]], out2[[3]],
             out[[4]],out1[[4]],out2[[4]],ncol=3)



