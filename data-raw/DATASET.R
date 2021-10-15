## code to prepare `DATASET` dataset goes here

#usethis::use_data("DATASET")

library(mice)
df_obs <- create_multilevel()
patt0 <- mice::ampute(df_obs, prop = 0.3, mech = "MCAR")$patterns
patt1 <- patt0[-1,]

data_example <- mice::ampute(df_obs, patterns = patt1, prop = 0.3, mech = "MCAR")$amp

usethis::use_data(data_example, overwrite = TRUE)
