FROM quay.io/pypa/manylinux_2_28_x86_64


COPY external/gromacs /app-src/external/gromacs
COPY external/apbs /app-src/external/apbs
COPY src /app-src/src

RUN yum install -y epel-release
RUN dnf -y install https://repo.almalinux.org/almalinux/8/devel/x86_64/os/Packages/suitesparse-static-4.4.6-11.el8.x86_64.rpm

RUN dnf -y install \
unzip \
wget \
arpack-devel \
arpack-static \
f2c \
eigen3-devel \
boost-devel \
openblas-serial64 \
openblas-static \
openblas-devel \
lapack-devel \
lapack-static \
suitesparse-devel 


RUN bash /io/scripts/static_dep_build_docker.sh
    
CMD ["/workspace/build/src/g_mmpbsa", "-h"]