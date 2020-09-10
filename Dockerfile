# get shiny serves plus tidyverse packages image
FROM rocker/shiny-verse:latest

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libgdal-dev \
    libproj-dev \
    libudunits2-0 \
    libudunits2-dev

  

# install R packages required 
# (change it dependeing on the packages you need)
RUN R -e 'install.packages(c("shiny","RSQLite","stringr","SPEI","ggplot2","zoo","sf","rgdal","tidyverse","raster","leaflet.extras","shinydashboard","plotly"), repos="http://cran.rstudio.com/")'

# copy the app to the image


COPY Dashboard.Rproj /srv/shiny-server/covid/
COPY app.R /srv/shiny-server/covid/
COPY run.R /srv/shiny-server/covid/

# select port
EXPOSE 3838

# allow permission
RUN sudo chown -R shiny:shiny /srv/shiny-server
# COPY shiny-server.sh /usr/bin/shiny-server.sh

# run app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/covid', host = '0.0.0.0', port = 3838)"]

#CMD ["/usr/bin/shiny-server.sh"]
