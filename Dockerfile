FROM python:latest

WORKDIR /app

RUN pip3 install pandas numpy statsmodels

ENV PATH /app:$PATH