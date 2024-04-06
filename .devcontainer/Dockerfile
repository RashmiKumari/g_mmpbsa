FROM ubuntu:20.04

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Europe/London apt-get -y install tzdata
RUN apt-get install -y gcc g++ curl libssl-dev build-essential pkg-config git python3 python3-dev python3-pip python3-pkgconfig

RUN mkdir /build
RUN mkdir /build/external
WORKDIR /build/external

RUN curl -O https://cmake.org/files/v3.27/cmake-3.27.7.tar.gz && tar -zxvf cmake-3.27.7.tar.gz
RUN curl -L -O https://ftp.gromacs.org/gromacs/gromacs-2023.2.tar.gz && tar -zxvf gromacs-2023.2.tar.gz

WORKDIR /build/external/cmake-3.27.7
RUN ./bootstrap && make && make install

WORKDIR /build/external/gromacs-2023.2
RUN mkdir build installed
WORKDIR /build/external/gromacs-2023.2/build
RUN /usr/local/bin/cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=on -DGMX_FFT_LIBRARY=fftpack \
    -DCMAKE_INSTALL_PREFIX=/build/external/gromacs-2023.2/installed && \
    make && \
    make install

ENV GMX_PATH=/build/external/gromacs-2023.2/installed
ENV GMX_SRC=/build/external/gromacs-2023.2

CMD ["echo", "testing"]
