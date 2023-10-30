# hash:sha256:c5cc8db686bc3896ef19a194241a2c2b540184f02f5581c57d49b72b8b34d0a7
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip3 install -U --no-cache-dir \
    pandas==2.0.3 \
    pillow==9.5.0 \
    plotly==5.16.1 \
    python-dotenv==1.0.0 \
    scipy==1.11.2 \
    statsmodels==0.14.0 \
    streamlit==1.27.0
