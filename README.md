# iowahills_dsp
Repo for Code Kit from http://www.iowahills.com/

A platform-independent C/C++ library with many **DSP** (digital signal processing) functions, amongst also *FIR* and *IIR* filter design - but also *FFT*, *DFT*, *Goertzel* and *Windowing* functions. Find a detailed description at http://www.iowahills.com/A7ExampleCodePage.html


## License
Iowa Hills Software, LLC, has put several sources online on their site http://www.iowahills.com/.
There is also a *Code Kit Download(zip)* provided at http://www.iowahills.com/A7ExampleCodePage.html

Theres is (or was) no license information on the website or inside the provided files. After requesting clarification and permission to publish on github, Daniel Klostermann (Iowa Hills Software, LLC) clarified:

> No license is required. Do whatever you want with it

He also invited me to publish on github.

Despite his very permissive words, i interpret as *Public Domain*, i put this repository under **MIT License** for having legal protection - still allowing everyone free use.


## Development

Note: The source contains global variables; thus, the library isn't safe for multithreaded use!

I've slightly modified the sources, to get them compile and link - mostly without warnings:
- removed some unused variables
- added cases and return (code) for missing/uncatched enum values in a switch/case
- removed/renamed duplicate Sinc() function
- restructured include and src files: modified include directives
..

Despite above changes, i've also added minimal cmake support with support for the install and uninstall targets.
Thus, the example(s) directory has an own CMakeLists.txt entry point, which requires the library, to be installed.
