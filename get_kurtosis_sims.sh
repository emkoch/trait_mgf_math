#!/bin/bash

# based off of script to download large files from google drive found on stack exchange:
# https://stackoverflow.com/questions/48133080/how-to-download-a-google-drive-url-via-curl-or-wget/48133859

fileid_norm="1rAhe3xilUlqUJZ0n4VNqCHaXyKE14v8e"
fileid_lap="1cbzjAbEtpuuwS8Y3-jxW-m9KqrEanOrz"
fileid_norm_growth="1cZdNAoUqiB4bcUs2lXgH55Qrn70Ges6Q"
fileid_lap_growth="1gtP5J8fb5XUvLPG-oSES76FHqKxQfz3n"

fname_norm="sim_phens_norm.pyc"
fname_lap="sim_phens_lap.pyc"
fname_norm_growth="sim_phens_norm_growth.pyc"
fname_lap_growth="sim_phens_lap_growth.pyc"

curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid_norm}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid_norm}" -o ${fname_norm}

curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid_lap}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid_lap}" -o ${fname_lap}

curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid_norm_growth}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid_norm_growth}" -o ${fname_norm_growth}

curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid_lap_growth}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid_lap_growth}" -o ${fname_lap_growth}
