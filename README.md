# Astrophysical Hashed Oct Tree-based Smoothed Particle Hydrodynamics code submodule: 'orig-tree16'

## Description

This is the original SPH + tree code that I received many moons ago in grad school. It is based on the old version of the hashed oct tree code. It will eventually be replaced by what's in the '2hot-apps' repo (once that's working).

## Building

This example assumes you want to build the code on a Haswell node on Darwin (this has been tested and works):

Clone the repo and cd into it:
```
git clone git@gitlab.lanl.gov:sph/orig-tree16.git
cd orig-tree16
```

Get interactive allocation (not necessary for machines where front and back ends have the same architecture and building on the front end is allowed).
```
salloc -t 2:00:00 -N 1 -A all_users -p general -C "cpu_family:haswell"
```

Load appropriate GCC and MPI compilers:
```
module load gcc/12.2.0 openmpi/4.1.5-gcc_12.2.0
```

build with:
```
make ARCH=darwin-hsw PAROS=mpi
```

n.b.: On Rocinante, use `ARCH=x86_64`
Check that you have an executable:
```
ls sph+nln/sph.darwin-hsw-mpi
```

To remove the build, use both (notice the 'C' vs 'c' in 'clean', they do different things)
```
make ARCH=<my-arch> PAROS=mpi clean
make ARCH=<my-arch> PAROS=mpi Clean
```


## Contributing
tbd

## License and copyright


BSD 3-Clause License

Copyright (c) 2026, Los Alamos National Laboratory. O5196.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
