# hash:sha256:14445418e992381bf8f4103fb9c2ebb57ac253afc256d5a324b0e9ac3077e896
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip3 install -U --no-cache-dir \
    pandas==2.0.2 \
    pillow==9.5.0 \
    plotly==5.18.0 \
    scipy==1.11.2 \
    statsmodels==0.14.0 \
    streamlit==1.27.0
