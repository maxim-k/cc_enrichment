# hash:sha256:32b031f1b326c7248aba58757ed181e4cecbe8c6b0bd2c50c5591444ca6ad4f9
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip3 install -U --no-cache-dir \
    pandas==2.0.3 \
    plotly==5.16.1 \
    python-dotenv==1.0.0 \
    scipy==1.11.2 \
    statsmodels==0.14.0 \
    streamlit==1.27.0
