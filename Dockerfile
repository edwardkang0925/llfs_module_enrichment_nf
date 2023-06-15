FROM python:latest

WORKDIR /app

RUN pip3 install pandas numpy

COPY preprocess.py .

ENV PATH /app:$PATH