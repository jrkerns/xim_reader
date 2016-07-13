Readme
======

XimReader is an open source tool for TrueBeam DeveloperMode provided by Varian Medical Systems, Palo Alto. 
HND compression algorithm is used to compress xim files. For a brief description of  HND compression algorithm 
please refer to the xim_readme.txt file.

XimReader is licensed under the Veritas Open Source (VOS) 1.0 License.
You may obtain a copy of the License at:

    website: http://radiotherapyresearchtools.com/license/

For questions, please send us an email at: TrueBeamDeveloper@varian.com                   

Features:
=========

* XimReader reads binary xim image file.
* HND decompression algorithm is used to decompress XIM image.

Getting Started:
================
usage: ximReader.py [-h] -f FNAME [-s SHOWIMAGE] [-v]

TrueBeam(TM) xim image reader.

optional arguments:
  -h, --help            show this help message and exit
  -f FNAME, --filename FNAME
                     Please enter name of the binary xim file
  -s SHOWIMAGE, --showImage SHOWIMAGE
                     Show xim image (Optional, 0 or 1)
  -v, --version         show program's version number and exit


**Developer Mode is intended for non-clinical use only and is NOT cleared for use on humans**


----------------------------------


Matlab Version
==============
There is also a separate xim image reader for Matlab in the Downloads section. Because it is based on the .NET framework, it can only be used on Windows machines. See the included readme.txt for more details.