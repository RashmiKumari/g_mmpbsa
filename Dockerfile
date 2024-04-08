FROM ubuntu:22.04

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Europe/London apt-get -y install tzdata
RUN apt-get install -y cmake gcc g++ curl libssl-dev build-essential pkg-config git python3 python3-dev python3-pip python3-pkgconfig

RUN mkdir /build
RUN mkdir /build/external
WORKDIR /build/external

RUN curl -L -O https://ftp.gromacs.org/gromacs/gromacs-2023.2.tar.gz && tar -zxvf gromacs-2023.2.tar.gz

WORKDIR /build/external/gromacs-2023.2
RUN mkdir build installed
WORKDIR /build/external/gromacs-2023.2/build
RUN cmake .. -DGMX_SIMD=SSE2 -DGMX_GPU=off -DGMXAPI=OFF -DGMX_INSTALL_LEGACY_API=ON -DGMX_FFT_LIBRARY=fftpack \
            -DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=ON \
            -DGMX_SIMD=SSE2 \
            -DCMAKE_INSTALL_PREFIX=/build/external/gromacs-2023.2/installed && \
    make && \
    make install

ENV GMX_SRC=/build/external/gromacs-2023.2
ENV GMX_PATH=/build/external/gromacs-2023.2/installed

WORKDIR /workspace
RUN mkdir build
WORKDIR /workspace/build
RUN cmake .. -DGMX_SRC=$GMX_SRC -DGMX_PATH=$GMX_PATH
RUN make
    
CMD ["/workspace/build/src/g_mmpbsa", "-h"]