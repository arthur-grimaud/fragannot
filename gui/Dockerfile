# Dockerfile for Fragannot GUI
# author: Micha Birklbauer
# version: 1.0.0

FROM ubuntu:22.04

LABEL maintainer="micha.birklbauer@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    git \
    python3-distutils \
    python3-lxml \
    python3-pip

RUN git clone https://github.com/arthur-grimaud/fragannot.git
WORKDIR fragannot
RUN git fetch -all
RUN git checkout standalone-gui 

RUN pip3 install -r requirements.txt
RUN pip3 install -r developer_requirements.txt

WORKDIR publish
RUN python3 setup.py install

WORKDIR ../gui
RUN pip3 install -r requirements.txt

CMD  ["streamlit", "run", "streamlit_app.py"]
