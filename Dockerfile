FROM nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04

RUN apt-get update \
    && apt-get install -y software-properties-common \
    && add-apt-repository ppa:jonathonf/python-3.6 \
    && apt-get update \
    && apt-get install -y \
        build-essential cmake \
        python3.6-dev python3-pip python3-tk \
        libpoppler-cpp-dev libmagic-dev libxrender-dev \
        libsm6 libxext6 libglib2.0-0 \
        libreoffice poppler-utils \
    && ln -s /usr/bin/python3.6 /usr/local/bin/python \
    && python -m pip install -U pip==18.1

RUN pip install torch==0.4.1

ADD CRF++-0.58.tar.gz /opt/

RUN cd /opt/CRF++-0.58 \
    && ./configure && make && make install && cd python \
    && cp /opt/CRF++-0.58/crfpp.h . \
    && python setup.py build && ldconfig \
    && python setup.py install

COPY requirements.txt /

RUN pip install -r requirements.txt \
    && python -c "import nltk; nltk.download('punkt'); nltk.download('averaged_perceptron_tagger')" \
    && rm -f requirements.txt

WORKDIR /app

CMD ["sleep", "infinity"]
