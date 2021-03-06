## Synopsis

FastPaSS is a spectral library search alogrithm for peptide-spectrum matching based on SpectraST 4.0. The SpectraST algorithm uses inner (dot) products to compute similarity between query and library spectra. In FastPaSS, these are replaced by matrix multiplications on the GPU for increased performance. The options now support more features from SpectraST 5.0.

## Installation

Requirements:
- CUDA toolkit
- CMAKE build automation
- Pugixml (included)
- Stringencoders (included)

Use the available macros in the Makefile or use the cmake itself:

$ mkdir build; cd build; cmake -DCMAKE_BUILD_TYPE=Release ..; make install

## API Reference

The FastPaSS executable is located in the 'bin' folder. Use FastPaSS -h to view the available options and description.

## Third party licenses

### Pugixml v1.5

Repository: http://pugixml.org/

Copyright (C) 2006-2014, by Arseny Kapoulkine (arseny.kapoulkine@gmail.com)
Report bugs and download new versions at http://pugixml.org/

This library is distributed under the MIT License:

Copyright (c) 2006-2014 Arseny Kapoulkine

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

### Stringencoders v3.10.3

Repository: https://code.google.com/p/stringencoders/

  Copyright 2005, 2006, 2007
  Nick Galbreath -- nickg [at] modp [dot] com
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:
 
    Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 
    Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 
    Neither the name of the modp.com nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.
 
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
  This is the standard "new" BSD license:
  http://www.opensource.org/licenses/bsd-license.php
 
