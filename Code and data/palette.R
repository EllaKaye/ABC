palette <- function(start, end, length) {
  r_1 <- col2rgb(start)[1]/255
  r_2 <- col2rgb(end)[1]/255
  g_1 <- col2rgb(start)[2]/255
  g_2 <- col2rgb(end)[2]/255
  b_1 <- col2rgb(start)[3]/255
  b_2 <- col2rgb(end)[3]/255
  out <- cbind(seq(r_1, r_2, length.out = length),seq(g_1, g_2, length.out = length),seq(b_1, b_2, length.out =  length))
  rgb(out)
}