FROM ubuntu:18.04

RUN apt-get update && apt-get -y upgrade && \
    apt-get install -y sudo && \
    apt-get install -y build-essential && \
    apt-get install -y emacs && \
    apt-get install -y wget && \
    apt-get install -y zlib1g && \
    apt-get install -y unzip && \
    apt-get install -y git && \
    apt-get install -y cmake && \
    apt-get install -y libssl-dev

# Install compilers.
RUN apt-get install -y gcc && \
    apt-get install -y g++

# SET path to compilers.
# https://stackoverflow.com/questions/17275348/how-to-specify-new-gcc-path-for-cmake
ENV CC=/usr/bin/gcc \
    CXX=/usr/bin/g++

# OpenBlas, Lapack, eigen
RUN apt-get install -y libeigen3-dev

WORKDIR $HOME/usr/src

CMD bash
