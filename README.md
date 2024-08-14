# iowahills_dsp
Repo for Code Kit from http://www.iowahills.com/

Site is reportedly down since September, 2021. Wayback machine has a backup: https://web.archive.org/web/20171110201917/http://iowahills.com/

A platform-independent C/C++ library with many **DSP** (digital signal processing) functions, amongst also *FIR* and *IIR* filter design - but also *FFT*, *DFT*, *Goertzel* and *Windowing* functions. Find a detailed description at http://www.iowahills.com/A7ExampleCodePage.html


## Build and Install

This library (and the test programs) are built with CMake, a cross-platform and open-source build system.
Typically the program is not built directly inside the source directory.
My preference is to build the program in a directory named ``build``, located in the root of the source directory or besides the source directory.

```
git clone https://github.com/hayguen/iowahills_dsp.git
cmake -S iowahills_dsp -B build_iowahills_dsp -DCMAKE_BUILD_TYPE=Release
cmake --build build_iowahills_dsp
sudo cmake --build build_iowahills_dsp --target install
```


## License
Iowa Hills Software, LLC, has put several sources online on their site http://www.iowahills.com/.
There is also a *Code Kit Download(zip)* provided at http://www.iowahills.com/A7ExampleCodePage.html

Theres is (or was) no license information on the website or inside the provided files. After requesting clarification and permission to publish on github, Daniel Klostermann (Iowa Hills Software, LLC) clarified:

> No license is required. Do whatever you want with it

He also invited me to publish on github.

Despite his very permissive words, i interpret as *Public Domain*, i put this repository under **MIT License** for having legal protection - still allowing everyone free use.


## Development

The sources contained global variables; thus, the library wasn't safe for multithreaded use. Hope, these are completely eliminated now.

I've slightly modified the sources, to get them compile and link - mostly without warnings:
- removed some unused variables
- added cases and return (code) for missing/uncatched enum values in a switch/case
- removed/renamed duplicate Sinc() function
- restructured include and src files: modified include directives
..

Despite above changes, i've also added minimal cmake support with support for the install and uninstall targets.
Thus, the example(s) directory has an own CMakeLists.txt entry point, which requires the library, to be installed.

## Related

* BSL, C++: https://github.com/hayguen/spuce
* MIT, C++: https://github.com/electro-smith/DaisySP
* MIT, C/C++: https://github.com/jgaeddert/liquid-dsp
* MIT, C++: https://github.com/signalsmith-audio/dsp
* BSL, C++: https://github.com/hbe72/cdsp
* LGPL, C++: https://github.com/AlexandreRouma/dsp
* ???, C++: https://github.com/AlexandreRouma/cdsp
* LGPL: https://github.com/lsp-plugins/lsp-dsp-lib
* LGPL: https://github.com/dac1976/dsp
* GPL: https://github.com/kfrlib/kfr
* GPL: https://github.com/mohabouje/eDSP
* Matlab: https://github.com/dario-pilori/dsp-library
