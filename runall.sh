#!/usr/bin/env sh
#
# Copyright (c) 2019-2023, ARM Limited. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
#

# Build if necessary.
if [ ! -d build ]; then
    git submodule update --init --recursive
    cmake \
        -Bbuild \
        -DCMAKE_BUILD_TYPE=Release \
        -DPBRT_FLOAT_AS_DOUBLE=1 \
        -DPBRT_SAMPLED_SPECTRUM=1 \
        -DCMAKE_INSTALL_PREFIX=${local_prefix_dir} &&
        cmake --build build --config Release --parallel 4 &&
        cmake --install build --config Release
fi

# Run all scenes.
./build/pbrt --quick ./scenes/killeroo-simple.pbrt --outfile killeroo-simple.exr
./build/pbrt --quick ./scenes/eg2020/cornell/cornell-irrad.pbrt --outfile cornell-irrad.exr
./build/pbrt --quick ./scenes/eg2020/cornell-ri/cornell-ri.pbrt --outfile cornell-ri.exr
./build/pbrt --quick ./scenes/eg2020/cornell-ri/cornell-ri-conventional.pbrt --outfile cornell-ri-conventional.exr
#./build/pbrt --quick ./scenes/eg2020/cornell-adv/cornell-adv.pbrt --outfile cornell-adv.exr
