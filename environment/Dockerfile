# hash:sha256:49116b54c919f18a9e8f2015edfa85c69d2a083201879b72fc1342c30953310e
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip3 install -U --no-cache-dir \
    pandas==2.0.3 \
    pillow==9.5.0 \
    plotly \
    scipy==1.11.2 \
    statsmodels==0.14.0 \
    streamlit==1.27.0
