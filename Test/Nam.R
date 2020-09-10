library(reticulate)
library(ggplot2)
use_python(python = "C:\\Users\\cagri\\AppData\\Local\\Programs\\Python\\Python38\\python.exe", required = TRUE)
py_config()
source_python("./Test/nam.py")

file <- "Alihoca"
params <- c(6.96780205e+00, 4.86098809e+02, 6.66247792e-01, 5.42601108e+02
            , 2.43815545e+01, 8.21285865e-01, 1.00000000e-02, 1.00000000e-02
            , 7.37979357e+02, 9.64180895e-01, 2.06295770e+00)

df <- run(97,params,T,file)

df <- data.frame(df[[1]])
df$date <- as.Date(row.names(df), "%Y-%m-%d") 

ggplot(df, aes(x=date)) + 
  geom_line(aes(y = Q), color = "darkred") + 
  geom_line(aes(y = Qsim), color="steelblue", linetype="twodash") 
