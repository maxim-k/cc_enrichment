# hash:sha256:54cee994e097331101c3d063bfff8a454273f426c9aa8598ade866b433aa7c01
FROM registry.apps.codeocean.com/codeocean/mambaforge3:22.11.1-4-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip3 install -U --no-cache-dir \
    pandas==2.0.3 \
    pillow==10.0.1 \
    plotly==5.16.1 \
    python-dotenv==1.0.0 \
    scipy==1.11.2 \
    statsmodels==0.14.0 \
    streamlit==1.27.0
