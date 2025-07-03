#Simple script to retrieve the datasets from the repository in Zenodo.

options(timeout = 500)

#url of the download link in Zenodo
url <- "https://zenodo.org/records/15798750/files/prepared_data.zip?download=1"

#download the files
download.file(url, destfile = "prepared data.zip")

#unzip the data
unzip("prepared data.zip")


#remove the zip file
file.remove("prepared data.zip")
