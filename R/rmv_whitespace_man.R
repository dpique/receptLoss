#remove leading whitespace from files in "man"

fls <- list.files(here::here("man"))

for (f in fls){
    fn <- here::here("man", f) #here::here("man/nhrs.Rd")
    test <- readLines(fn)
    test2 <- trimws(x = test)
    fc<-file(fn)
    writeLines(test2, fc)
    close(fc)
}
