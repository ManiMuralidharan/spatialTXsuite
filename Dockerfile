FROM rocker/r-ver:4.4.2
RUN apt-get update && apt-get install -y --no-install-recommends \
 libcurl4-openssl-dev libssl-dev libxml2-dev libhdf5-dev libgeos-dev \
 libgdal-dev libproj-dev libudunits2-dev libglpk-dev \
 libfontconfig1-dev libharfbuzz-dev libfribidi-dev && rm -rf /var/lib/apt/lists/*
RUN R -e 'install.packages(c("BiocManager","shiny","ggplot2","Matrix","igraph","harmony"), repos="https://cloud.r-project.org")'
RUN R -e 'BiocManager::install(c("SpatialExperiment","scater","scran","scuttle","edgeR","spdep","BiocNeighbors","DropletUtils","SpatialExperimentIO"),ask=FALSE,update=FALSE)'
WORKDIR /app
COPY app.R /app/app.R
COPY R/ /app/R/
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
