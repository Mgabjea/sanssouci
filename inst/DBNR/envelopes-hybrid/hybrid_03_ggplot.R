resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- file.path(resPath,Sys.Date()) #ligne ajoutée pour trouvé le bon fichié 
resPath <- R.utils::Arguments$getReadablePath(resPath)


#figPath <- file.path(resPath, "fig/DBNR")  # lignes rajouté, utilité non comprise encore
#figPath <- R.utils::Arguments$getReadablePath(figPath)

figPath <- "fig/DBNR"  #ligne ajoutée
figPath <- R.utils::Arguments$getWritablePath(figPath) #ligne ajoutée

figPath <- R.utils::Arguments$getReadablePath("fig/DBNR") #ligne originale


configs <- expand.grid(
    grouped = groupeds,
    setting = settings,
    stringsAsFactors = FALSE)

# grouped <- TRUE
# setting <- "rgauss"

for (grouped in unique(configs$grouped)) {
    for (setting in unique (configs$setting)) {
  
      
        #filename <- sprintf("conf-env_m=12800_s=100_K1=8_d=0.75_barmu=2_grouped=TRUE_setting=%s.rds", setting) #ligne test ajoutée
        filename <- sprintf("all-conf-env-%s.rds", setting)
        # filename <- sprintf("conf-env-alpha_grouped=%s_setting=%s.rds", grouped, setting) #ligne originale
        pathname <- file.path(resPath, filename)
        dat <- readRDS(pathname)
        rm(pathname)
        dim(dat)
        
        names(dat)
        order_ <- "p.value"
        sdat <- subset(dat, s == 100 & K1 == 8 & d > 0.5 & order == order_)
        
        dim(sdat)
        
        library("magrittr")
        library("dplyr")
        #gdat <- sdat %>% group_by(idxs, method, d, barmu, alpha) #line removed due to barmu column removed
        gdat <- sdat %>% group_by(idxs, method, d, barmu, alpha)
        pdat <- gdat %>% summarise(value = mean(value))
        
        library("ggplot2")
        m1max <-  floor(unique(sdat$s)*max(sdat$d))*unique(sdat$K1)
        xymax <- 4/3*m1max
        
        ff <- gsub("\\.rds$", ".pdf", filename)
        pathname <- file.path(figPath, ff)

        pdat <- subset(pdat, alpha %in% alphas)
        pdat$alpha <- as.factor(pdat$alpha)
        ggplot(pdat, aes(idxs, value, colour = method, linetype = alpha)) + 
            geom_line() + facet_grid(d ~ barmu, labeller = label_both) +
            xlim(1, xymax) + ylim(0, xymax) +
            ylab("Upper bound on the number of false positives") +
            xlab(sprintf("Hypotheses sorted by %s", order_))
        ggsave(filename = pathname)
    }
}
#egze
