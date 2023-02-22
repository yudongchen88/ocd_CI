# This R file provides the codes to generate Tables 1, 2, S1 and Figure 2 of the accompanied paper
source('ocd_CI.R')


# Generate Table 1 (main simulation table)
niter <- 2000
df1 <- expand.grid(p=c(100), s=c(2,10,100), vartheta=c(2,1), ratio=c(2,1,0.5))
df2 <- expand.grid(p=c(500), s=c(2,22,500), vartheta=c(2,1), ratio=c(2,1,0.5))
df <- rbind(df1, df2)
df['beta'] <- df['vartheta'] * df['ratio']
df <- df[, -4]
df <- df[with(df, order(p, s, -vartheta, -beta)), ]
df['kaul_CI_len'] <- df['kaul_CI_cov'] <- df['ocd_CI_len'] <- df['ocd_CI_cov'] <- df['ocd_delay'] <- rep(0, nrow(df))
se_df <- data.frame(df)

for (i in nrow(df)){
  ocd_delay <- ocd_CI_cov<- ocd_CI_len <- kaul_CI_cov <- kaul_CI_len <- rep(0, niter)
  for (iter in 1:niter){
    ret <- with(df, ocd_CI(p=p[i], s=s[i], z=1000, vartheta=vartheta[i], signal_shape='random', spatial='identity', beta=beta[i],
                    ocd_thresholds='MC', l_choice='zero', d1_choice='default', alpha=0.05, kaul_comp = TRUE, print_type ='muted'))
    
    ocd_CI_cov[iter] <- ifelse(ret$CI[1] <= 1000 & ret$CI[2] >= 1000, 100, 0)
    ocd_delay[iter] <- max(ret$CI[2] - 1000, 0)
    ocd_CI_len[iter] <- ret$CI[2] - ret$CI[1]
    
    kaul_CI_cov[iter] <- ifelse(ret$kaul_CI[1] <= 1000 & ret$kaul_CI[2] >= 1000, 100, 0)
    kaul_CI_len[iter] <- ret$kaul_CI[2] - ret$kaul_CI[1]
  }
  df$ocd_delay[i] <- mean(ocd_delay)
  df$ocd_CI_cov[i] <- mean(ocd_CI_cov)
  df$ocd_CI_len[i] <- mean(ocd_CI_len)
  df$kaul_CI_cov[i] <- mean(kaul_CI_cov)
  df$kaul_CI_len[i] <- mean(kaul_CI_len)
  
  se_df$ocd_delay[i] <- sd(ocd_delay)/sqrt(niter)
  se_df$ocd_CI_cov[i] <- sd(ocd_CI_cov)/sqrt(niter)
  se_df$ocd_CI_len[i] <- sd(ocd_CI_len)/sqrt(niter)
  se_df$kaul_CI_cov[i] <- sd(kaul_CI_cov)/sqrt(niter)
  se_df$kaul_CI_len[i] <- sd(kaul_CI_len)/sqrt(niter)
}
# Table 1 (main)
df
# Table 1 (standard errors)
se_df
# End of generating Table 1




# Generate Table 2 (support recovery)
niter <- 2000
df <- expand.grid(s=c(5, 50), vartheta=c(2,1), signal_shape=c('uniform', 'inv_sqrt', 'harmonic'))
df <- df[with(df, order(signal_shape, s, -vartheta)), ]
df$signal_shape <- as.character(df$signal_shape)
df['supp_aug_includes'] <- df['supp_included'] <- rep(0, nrow(df))
se_df <- data.frame(df)

for (i in nrow(df)){
  supp_included <- supp_aug_includes <- rep(0, niter)
  for (iter in 1:niter){
    ret <- with(df, ocd_CI(p=100, s=s[i], z=1000, vartheta=vartheta[i], signal_shape=signal_shape[i], spatial='identity', beta=vartheta[i],
                  ocd_thresholds='MC', l_choice='given', d1_choice='variant1', alpha=0.05, print_type='muted'))
    supp <- ret$support
    supp_aug <- ret$support_aug
    
    S <- support(p=100, s=df$s[i], vartheta=df$vartheta[i], beta=df$vartheta[i], signal_shape=df$signal_shape[i])
    supp_included[iter] <- ifelse((length(supp) == 0) | (max(supp) <= S[1]) , 100, 0)
    supp_aug_includes[iter] <- ifelse((length(supp_aug) >= S[2]) & (supp_aug[S[2]]==S[2]) , 100, 0)
  }
  df$supp_included[i] <- mean(supp_included)
  df$supp_aug_includes[i] <- mean(supp_aug_includes)
  se_df$supp_included[i] <- sd(supp_included)/sqrt(niter)
  se_df$supp_aug_includes[i] <- sd(supp_aug_includes)/sqrt(niter)
}
# Table 2 (main)
df
# Table 2 (standard errors)
se_df
# End of generating Table 2




# Generate Table S1 (spatial dependence)
niter <- 2000
df <- expand.grid(rho=c(0.5, 0.75), s=c(2,10,100), vartheta=c(2,1,0.5))
df <- df[with(df, order(rho, s, -vartheta)), ]
df['ocd_CI_len'] <- df['ocd_CI_cov'] <- df['ocd_delay'] <- rep(0, nrow(df))
se_df <- data.frame(df)

for (i in nrow(df)){
  ocd_delay <- ocd_CI_cov<- ocd_CI_len <- rep(0, niter)
  for (iter in 1:niter){
    ret <- with(df, ocd_CI(p=100, s=s[i], z=1000, vartheta=vartheta[i], signal_shape='random', spatial='toeplitz', rho=rho[i], beta=vartheta[i],
                           ocd_thresholds='MC', l_choice='zero', d1_choice='default', alpha=0.05, print_type ='muted'))
    
    ocd_CI_cov[iter] <- ifelse(ret$CI[1] <= 1000 & ret$CI[2] >= 1000, 100, 0)
    ocd_delay[iter] <- max(ret$CI[2] - 1000, 0)
    ocd_CI_len[iter] <- ret$CI[2] - ret$CI[1]
  }
  df$ocd_delay[i] <- mean(ocd_delay)
  df$ocd_CI_cov[i] <- mean(ocd_CI_cov)
  df$ocd_CI_len[i] <- mean(ocd_CI_len)
  
  se_df$ocd_delay[i] <- sd(ocd_delay)/sqrt(niter)
  se_df$ocd_CI_cov[i] <- sd(ocd_CI_cov)/sqrt(niter)
  se_df$ocd_CI_len[i] <- sd(ocd_CI_len)/sqrt(niter)
}
# Table S1 (main)
df
# Table S1 (standard errors)
se_df
# End of generating Table S1




# Generate Figure 1 left panel
niter <- 500
df <- data.frame(p=c(100,100,100,500,500,500), s=c(2,10,100,2,22,500))
d1_coef <- seq(0, 0.795, by=0.005)
cov <- matrix(0, 160, 6)
for (j in 1:6){
  d1_all <- d1_coef * sqrt(log(df$p[j]/0.05))
  for (i in 1:160){
    temp_cov <- rep(0, niter)
    for (iter in 1:niter){
      ret <- ocd_CI(p=df$p[j], s=df$s[j], z=1000, vartheta=1, signal_shape='random', spatial='identity', beta=1,
                    ocd_thresholds='MC', l_choice='zero', d1_choice='other', d1=d1_all[i], alpha=0.05, print_type='muted')
      temp_cov[iter] <- ifelse(ret$CI[1] <= 1000 & ret$CI[2] >= 1000, 1, 0)
    }
    cov[i, j] <- mean(temp_cov)
  }
}
par(mfrow=c(1, 1), mar=c(2.5, 2.5, 0.7, 0.7), mgp=c(1.5,0.5,0), tcl=-0.2)
palet <- matplotlib_palette(10)
plot(d1_coef, cov[, 1], type='l', col=palet[1], xlab='tuning parameter c', ylab='Coverage probability', lwd=2)
lines(d1_coef, cov[, 2], lty=2, col=palet[1], lwd=2)
lines(d1_coef, cov[, 3], lty=3, col=palet[1], lwd=2)
lines(d1_coef, cov[, 4], lty=1, col=palet[2], lwd=2)
lines(d1_coef, cov[, 5], lty=2, col=palet[2], lwd=2)
lines(d1_coef, cov[, 6], lty=3, col=palet[2], lwd=2)
abline(h=0.95, lty=3, lwd=1)
legend.col = c(palet[1],  palet[2], 'black', 'black', 'black')
legend.lty = c(1,1,1,2,3)
legend.txt = c('p=100', 'p=500', 's=2', 's=sqrt(p)', 's=p') 
legend('bottomright', legend=legend.txt, lty=legend.lty, col=legend.col, lwd=2, cex=0.7, pt.cex=1.1)
# dev.off()


# Generate Figure 1 right panel
niter <- 500
df <- data.frame(vartheta=rep(c(2,1,0.5), c(3,3,3)), beta=c(8,2,0.5,4,1,0.25,2,0.5,0.125))
d1_coef <- seq(0, 0.795, by=0.005)
d1_all <- d1_coef * sqrt(log(100/0.05))
cov <- matrix(0, 160, 9)
for (j in 1:9){
  for (i in 1:160){
    temp_cov <- rep(0, niter)
    for (iter in 1:niter){
      ret <- ocd_CI(p=100, s=10, z=1000, vartheta=df$vartheta[j], signal_shape='random', spatial='identity', beta=df$beta[j],
                    ocd_thresholds='MC', l_choice='zero', d1_choice='other', d1=d1_all[i], alpha=0.05, print_type='muted')
      temp_cov[iter] <- ifelse(ret$CI[1] <= 1000 & ret$CI[2] >= 1000, 1, 0)
    }
    cov[i, j] <- mean(temp_cov)
  }
}
par(mfrow=c(1, 1), mar=c(2.5, 2.5, 0.7, 0.7), mgp=c(1.5,0.5,0), tcl=-0.2)
palet <- matplotlib_palette(10)
plot(d1_coef, cov[, 1], type='l', col=palet[1], xlab='tuning parameter c', ylab='Coverage probability', lwd=2)
lines(d1_coef, cov[, 2], lty=2, col=palet[1], lwd=2)
lines(d1_coef, cov[, 3], lty=3, col=palet[1], lwd=2)
lines(d1_coef, cov[, 4], lty=1, col=palet[2], lwd=2)
lines(d1_coef, cov[, 5], lty=2, col=palet[2], lwd=2)
lines(d1_coef, cov[, 6], lty=3, col=palet[2], lwd=2)
lines(d1_coef, cov[, 7], lty=1, col=palet[3], lwd=2)
lines(d1_coef, cov[, 8], lty=2, col=palet[3], lwd=2)
lines(d1_coef, cov[, 9], lty=3, col=palet[3], lwd=2)
abline(h=0.95, lty=3, lwd=1)
legend.col = c(palet[1],  palet[2], palet[3], 'black', 'black', 'black')
legend.lty = c(1,1,1,1,2,3)
legend.txt = c('vartheta=2', 'vartheta=1', 'vartheta=0.5', 'beta=vartheta', 'beta=4*vartheta', 'beta=vartheta/4') 
legend('bottomright', legend=legend.txt, lty=legend.lty, col=legend.col, lwd=2, cex=0.7, pt.cex=1.1)
# dev.off()
# End of generating Figure 1



# Generate Figure 2 left panel
niter <- 500
df <- data.frame(s=rep(c(5,50),3), signal_shape=rep(c('uniform', 'inv_sqrt', 'harmonic'), c(2,2,2)))
df$signal_shape <- as.character(df$signal_shape)
d1_all <- 0.46 + 0.04*(1:125)
supp_included <- supp_aug_includes <- matrix(0, 125, 6)
for (i in 1:125){
  for (j in 1:6){
    temp_supp_included <-temp_supp_aug_includes <- rep(0, niter)
    for (iter in 1:niter){
      ret <- ocd_CI(p=100, s=df$s[j], z=1000, vartheta=2, signal_shape=df$signal_shape[j], spatial='identity', beta=2,
                             ocd_thresholds='MC', l_choice='zero', d1_choice='other', d1=d1_all[i], alpha=0.05, print_type='muted')
      supp <- ret$support
      supp_aug <- ret$support_aug
      
      S <- support(p=100, s=df$s[j], vartheta=2, beta=2, signal_shape=df$signal_shape[j])
      temp_supp_included[iter] <- ifelse((length(supp) == 0) | (max(supp) <= S[1]) , 1, 0)
      temp_supp_aug_includes[iter] <- ifelse((length(supp_aug) >= S[2]) & (supp_aug[S[2]]==S[2]) , 1, 0)
    }
    supp_included[i, j] <- mean(temp_supp_included)
    supp_aug_includes[i, j] <- mean(temp_supp_aug_includes)
  }
}

par(mfrow=c(1, 1), mar=c(2.5, 2.5, 0.7, 0.7), mgp=c(1.5,0.5,0), tcl=-0.2)
palet <- matplotlib_palette(10)
plot(1-supp_included[, 1], supp_aug_includes[, 1], type ='l', col = palet[1], 
     xlab = '1-belong', ylab = 'contain', lwd = 2)
lines(1-supp_included[, 2], supp_aug_includes[, 2], lty = 2, col = palet[1], lwd = 2)
lines(1-supp_included[, 3], supp_aug_includes[, 3], col = palet[2], lwd = 2)
lines(1-supp_included[, 4], supp_aug_includes[, 4], lty = 2, col = palet[2], lwd = 2)
lines(1-supp_included[, 5], supp_aug_includes[, 5], col = palet[3], lwd = 2)
lines(1-supp_included[, 6], supp_aug_includes[, 6], lty = 2, col = palet[3], lwd = 2)
abline(v = 0.05, lty = 3, lwd = 2)

# sqrt(2*log(p)) = 3.03, closet to d1_all[64]; sqrt(2*log(p/alpha)) = 3.89, closest to d1_all[86]
points(1-supp_included[64, 1], supp_aug_includes[64, 1], pch = 19, col = palet[1], cex = 1.1)
points(1-supp_included[86, 1], supp_aug_includes[86, 1], pch = 17, col = palet[1], cex = 1.1)
points(1-supp_included[64, 2], supp_aug_includes[64, 2], pch = 19, col = palet[1], cex = 1.1)
points(1-supp_included[86, 2], supp_aug_includes[86, 2], pch = 17, col = palet[1], cex = 1.1)
points(1-supp_included[64, 3], supp_aug_includes[64, 3], pch = 19, col = palet[2], cex = 1.1)
points(1-supp_included[86, 3], supp_aug_includes[86, 3], pch = 17, col = palet[2], cex = 1.1)
points(1-supp_included[64, 4], supp_aug_includes[64, 4], pch = 19, col = palet[2], cex = 1.1)
points(1-supp_included[86, 4], supp_aug_includes[86, 4], pch = 17, col = palet[2], cex = 1.1)
points(1-supp_included[64, 5], supp_aug_includes[64, 5], pch = 19, col = palet[3], cex = 1.1)
points(1-supp_included[86, 5], supp_aug_includes[86, 5], pch = 17, col = palet[3], cex = 1.1)
points(1-supp_included[64, 6], supp_aug_includes[64, 6], pch = 19, col = palet[3], cex = 1.1)
points(1-supp_included[86, 6], supp_aug_includes[86, 6], pch = 17, col = palet[3], cex = 1.1)
legend.col = c(palet[1],  palet[2], palet[3], 'black', 'black')
legend.lty = c(1,1,1,1,2)
legend.txt = c('uniform', 'inverse sqrt', 'harmonic', 's=5', 's=50') 
legend('right', legend=legend.txt, lty=legend.lty, col=legend.col, lwd=2, cex=0.7, pt.cex=1.1)
dev.off()
# End of generating Figure 2 left panel


# Generate Figure 2 right panel
niter <- 500
count <- rep(0, 100)
for (iter in 1:niter){
  ret <- ocd_CI(p=100, s=20, z=1000, vartheta=2, signal_shape='inv_sqrt', spatial='identity', beta=2,
                ocd_thresholds='MC', l_choice='zero', d1_choice='variant2', alpha=0.05, print_type='muted')
  supp_aug <- ret$support_aug
  count[supp_aug] = count[supp_aug] + 1
}
coor_freq <- count/niter

par(mfrow=c(1, 1), mar=c(2.5, 3, 2.5, 0.2), oma=c(1, 1, 1, 1), mgp=c(1.2,0.2,0), tcl=-0.2)
my_bar <- barplot(coor_freq, col=c(rep('red', 20), rep('black', 80)), xaxt='n',
                  xlab = 'coordinate', ylab = 'frequency', ylim = c(0,1))
axis(side=1,at=my_bar[c(1,seq(10,100,10))],labels=c(1,seq(10,100,10)))

# dev.off()
# End of generating Figure 2 right panel

