## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

popularity <- read.csv("/Users/maireadshaw/My Cloud/r2MLM Wrapper/hox2000.csv") %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(mean_extraversion = mean(extrav)) %>%
  dplyr::mutate(mean_sex = mean(sex)) %>%
  dplyr::mutate(extravCWC = extrav - mean_extraversion) %>%
  dplyr::mutate(sexCWC = sex - mean_sex) %>%
  select(-c(popteach, Zextrav, Zsex, Ztexp, Zpopular, Zpopteach, Cextrav, Ctexp, Csex))
