library("future")
source("inst/DBNR/envelopes-hybrid/hybrid_simu.hulk.R")
source('~/Documents/stage_Stats/code/sanssouci/inst/DBNR/envelopes-hybrid/hybrid_00_setup.R')

# plan(multisession, workers = 100) # multiprocess est déprécié

resPath <- "resData/DBNR/confidenceEnvelopes"
resPath <- file.path(resPath, Sys.Date())
#resPath <- file.path(resPath, "2018-05-31")
resPath <- R.utils::Arguments$getWritablePath(resPath)

nb <- nrow(configs)
#nb <- 40
#cc <- 220 ## (a nice one)

#alphas <- c(0.001, 0.01, 0.05, 0.2, 0.5)
#alphas <- c(0.001, 0.05, 0.05-0.001)
alphas <- c(0.0001, 0.001, 0.05, 0.05-0.0001, 0.05-0.001)[3] #alphas == 0.05

for (cc in 1:nb) {
    print(cc)
    conf <- configs[cc, ]
    # > conf
    # m   s K1   d barmu grouped setting
    # 220 10000 100  8 0.9     4    TRUE   const
    stag <- paste("m=", conf[["m"]], "_",
                  "s=", conf[["s"]], "_",
                  "K1=", conf[["K1"]], "_",
                  "d=", conf[["d"]], "_",
                  "barmu=", conf[["barmu"]], "_",
                  "grouped=", conf[["grouped"]], "_",
                  "setting=", conf[["setting"]], sep = "")
    filename <- sprintf("conf-env_%s.rds", stag)
    print(filename)
    
        sdatList <- list()
        for (rr in 1:repl) {
            res <- simu.hulk(m = conf[["m"]], 
                         s = conf[["s"]], 
                         K1 = conf[["K1"]], 
                         d = conf[["d"]], 
                         barmu = conf[["barmu"]],
                         grouped = conf[["grouped"]], 
                         setting = conf[["setting"]],
                         # methods_ = c("tree", "part", "Simes", "Oracle", "hyb1", "hyb2", "hyb3"),
                         # methods_ = c("hyb1", "part", "Simes", "tree"),
                         methods_ = c("tree", "Oracle", "hyb2", "hyb3"),
                         alpha = alphas, verbose = FALSE)
            sdat <- Reduce(rbind, res)
            sdat$replication <- rr
            sdatList[[rr]] <- sdat
        
        dat <- Reduce(rbind, sdatList)
        rownames(conf) <- NULL
        dat <- cbind(dat, conf)
        pathname <- file.path(resPath, filename)
        saveRDS(dat, pathname)
        dat <- readRDS(pathname)
    }
}

