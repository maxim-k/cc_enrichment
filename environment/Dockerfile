# hash:sha256:9ce9b4493db94a9eb39ec5171a8f60067a337330a6097a7541345b283089d6f2
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip3 install -U --no-cache-dir \
    pandas==2.0.3 \
    pillow==9.5.0 \
    plotly==5.18.0 \
    scipy==1.11.2 \
    statsmodels==0.14.0 \
    streamlit==1.27.0
