FROM quay.io/pypa/manylinux2010_x86_64

LABEL maintainer="Martino Pilia <martino.pilia@gmail.com>"

# Install CMake
RUN curl -L "https://github.com/Kitware/CMake/releases/download/v3.15.4/cmake-3.15.4.tar.gz" --output cmake-3.15.4.tar.gz \
&&  tar xzfv cmake-3.15.4.tar.gz \
&&  cd cmake-3.15.4 \
&&  ./bootstrap -- -DCMAKE_INSTALL_PREFIX:PATH="/usr" -DCMAKE_BUILD_TYPE:STRING="Release" \
&&  make \
&&  make install \
&&  cd .. \
&&  rm -rf "cmake-3.15.4" "cmake-3.15.4.tar.gz"

# Install CUDA
RUN yum -y update \
&&  yum -y install yum-utils \
&&  yum-config-manager --add-repo http://developer.download.nvidia.com/compute/cuda/repos/rhel6/x86_64/cuda-rhel6.repo \
&&  yum -y update \
&&  yum -y install \
        cuda-compiler-10-1 \
        cuda-cudart-dev-10-1 \
&&  yum -y remove yum-utils \
&&  yum clean all

ENV PATH=/usr/local/cuda-10.1/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/local/cuda-10.1/lib64:$LD_LIBRARY_PATH
