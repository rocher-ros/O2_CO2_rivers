options(timeout = 500)


url <- "https://zenodo.org/records/13990674/files/prepared%20data.zip?download=1"


download.file(url, destfile = "prepared data.zip")

unzip("prepared data.zip")
