BiocVersion <- "3.14"
fs <- list.files(pattern="*metadata*.R")
for (f in fs) {
  source(f)
}