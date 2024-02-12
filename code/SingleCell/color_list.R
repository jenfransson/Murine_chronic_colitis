
# get_color_list = function(n){
#   indeces = c(1,11,5,7,14,13,15,12,6,9,16,2,10,8,3,4)
#   
#   if(n>16) stop("n must be 16 or smaller")
#   c("#9E1531","#C66A27","#F3953A","#C9971D","#F9C027","#DDE48A",
#     "#97C13A","#64A638","#174E3F","#43ADA8","#3C87EA","#7863A9",
#     "#1C124A","#6D1C8C","#CF2561","#EC8483")[sort(indeces[1:n])]
# }


get_color_list = function(n){
  indeces = c(1,11,5,7,15,14,16,13,6,9,17,2,10,8,3,4,12,18)
  
  if(n>18) stop("n must be 18 or smaller")
  c("#9E1531","#C66A27","#F3953A","#C9971D","#F9C027","#DDE48A",
    "#97C13A","#64A638","#174E3F","#43ADA8","#3C87EA","#7863A9",
    "#C9CDFF","#1C124A","#6D1C8C","#CF2561","#EC8483","#ECAFC6")[sort(indeces[1:n])]
}
