# Copyright (c) 2014, Varian Medical Systems, Inc. (VMS)
# All rights reserved.
#
# ximReader is an open source tool for reading .xim file (both compressed and uncompressed)
# HND compression algorithm is used to compress xim files. For a brief description of  HND
# compression algorithm please refer to the xim_readme.txt file.
#
# ximReader is licensed under the VarianVeritas License.
# You may obtain a copy of the License at:
#
#       website: http://radiotherapyresearchtools.com/license/
#
# For questions, please send us an email at: TrueBeamDeveloper@varian.com
#
# Developer Mode is intended for non-clinical use only and is NOT cleared for use on humans.
#
# Created on: 12:04:06 PM, Sept. 26, 2014
# Authors: Pankaj Mishra and Thanos Etmektzoglou
#
# Modified on : 10:58:045 AM, Jul. 24, 2015
# Modified by Nilesh Gorle

import textwrap
import os
import struct, sys, numpy as np
from matplotlib import pyplot as plt
from argparse import ArgumentParser

LINE_SPACE = 150
XIMREADER_FILENAME = "XimReaderData.txt"
XIMREADER_IMG_NAME = "XimReaderImage.png"

class XimFileInfo(object):
    '''
    XimFileInfo is the main class to store header data, histogram data
    and property data in text format and saving the plot image.  
    '''

    def __init__(self, **kwargs):
        '''
        XimFileInfo Constructor
        '''
        self.headerDataDict = kwargs.get('headerDataDict')
        self.histogramDataDict = kwargs.get('histogramDataDict')
        self.propertyDataList = kwargs.get('propertyDataList')

        self.outputFolder = os.path.join(os.path.dirname(__file__), 'XimData')
        if not os.path.exists(self.outputFolder):
            os.mkdir(self.outputFolder)

        outputPath = os.path.join(self.outputFolder, XIMREADER_FILENAME)
        print "%s is created on %s" % (XIMREADER_FILENAME, self.outputFolder)

        # Open file to write data.
        self.ximFile = open(outputPath, "w")

        # Flushing all contents from existing text file and appending new contents
        self.ximFile.seek(0)
        self.ximFile.truncate()

    def saveHeaderInfo(self):
        '''
        saving header information
        '''
        title = "Header Data".center(LINE_SPACE)

        self.ximFile.writelines("="*LINE_SPACE + "\n" + title + "\n" + "="*LINE_SPACE + "\n")

        thead = "FormatIdentifier".center(20), "|", "FormatVersion".center(20), \
                "|", "Width".center(10), "|", "Height".center(10), "|", \
                "BitsPerPixel".center(20), "|", "BytesPerPixel".center(20), "|", \
                "CompressionIndicator".center(20)
        thead = "".join(thead)
        self.ximFile.writelines(thead + "\n" + "="*LINE_SPACE)

        FormatIdentifier = (self.headerDataDict["FormatIdentifier"]).replace("\x00", "")
        tbody = FormatIdentifier.center(20), "|", \
                str(self.headerDataDict.get("FormatVersion")).center(20), "|", \
                str(self.headerDataDict.get("Width")).center(10), "|", \
                str(self.headerDataDict.get("Height")).center(10), "|", \
                str(self.headerDataDict.get("BitsPerPixel")).center(20), "|", \
                str(self.headerDataDict.get("BytesPerPixel")).center(20), "|", \
                str(self.headerDataDict.get("CompressionIndicator")).center(20)
        tbody = "".join(tbody)
        self.ximFile.writelines("\n" + tbody + "\n" + "="*LINE_SPACE)
        print "Header data stored into file successfully."

    def saveHistogramInfo(self):
        '''
        saving histogram information
        '''
        title = "Histogram Data".center(LINE_SPACE)
        self.ximFile.writelines("\n"*4 + "="*LINE_SPACE + "\n" + title + "\n" + "="*LINE_SPACE + "\n")
        thead = "NumberOfBins".center(20), "|", \
                "Value".center(20)

        thead = "".join(thead)
        self.ximFile.writelines(thead + "\n" + "="*LINE_SPACE + "\n")

        valList = textwrap.wrap(str(self.histogramDataDict.get("Value")), width=130)


        tbody = str(self.histogramDataDict.get("NumberOfBins")).center(20) + "|\n"

        for i in xrange(len(valList)):
            tbody += (" "*20 + "|  " + (valList[i]).center(20) + "\n")

        tbody = "".join(tbody)
        self.ximFile.writelines(tbody + "\n" + "="*LINE_SPACE)
        print "Histogram data stored into file successfully."

    def savePropertyInfo(self):
        '''
        saving property information
        '''
        DATATYPE_DICT = {0: "Integer", 1: "Double", 2: "String",
                         4: "Double Array", 5: "Integer Array"}

        title = "Property Data".center(LINE_SPACE)
        self.ximFile.writelines("\n"*4 + "="*LINE_SPACE + "\n" + title + "\n" + "="*LINE_SPACE + "\n")

        thead = "Length".center(6), "|", "Name".center(50), "|", \
                "Type".center(15), "|", "Value".center(25)
        thead = "".join(thead)
        self.ximFile.writelines(thead + "\n" + "="*LINE_SPACE + "\n")

        tbody = ""
        for length, name, ptype, value in self.propertyDataList:
            if isinstance(value, str):
                value = value.replace('\n', ' ').replace('\r', ',')
            value = str(value)

            if len(value) > 54:
                if value.startswith('<') :
                    import xml.dom.minidom
                    xml_string = xml.dom.minidom.parseString(value)
                    pretty_xml_as_string = xml_string.toprettyxml()
                    valList = pretty_xml_as_string.split("\n")

                    tbody = (str(length).center(6) + "|" + str(name).center(50) + "|" + \
                      (DATATYPE_DICT[ptype]).center(15) + "|\n")

                    for i in xrange(len(valList)):
                        tbody += (" "*73 + "|" + (valList[i]).center(25) + "\n")
                else:
                    valList = textwrap.wrap(str(value), width=80)

                    tbody = (str(length).center(6) + "|" + str(name).center(50) + "|" + \
                             (DATATYPE_DICT[ptype]).center(15) + "|\n")

                    for i in xrange(len(valList)):
                            tbody += (" "*73 + "|  " + (valList[i]).center(25) + "\n")

            else:
                tbody = (str(length).center(6) + "|" + str(name).center(50) + "|" + \
                        (DATATYPE_DICT[ptype]).center(15) + "|" + str(value).center(25))

            self.ximFile.writelines(tbody + "\n" + "-"*LINE_SPACE + "\n")

        self.ximFile.writelines("="*LINE_SPACE)
        print "Property data stored into file successfully."

    def closeFile(self):
        '''
        Closing ximData.txt file
        '''
        self.ximFile.close()

    def saveXimInfo(self):
        '''
        Saving headerData, histogramData and propertyData in text file
        '''
        # storing headerData
        self.saveHeaderInfo()

        # storing histogramData
        self.saveHistogramInfo()

        # storing propertyData
        self.savePropertyInfo()

        # Close file
        self.closeFile()

class XimReader():
    '''
    XimReader is the main class for converting an xim file to a two 
    dimensional image. This class reads header, pixel data, histogram 
    and properties of a given xim file. If the xim image is compressed 
    then HND decompression algorithm is used to for decompression. 
    Note: HND is a lossless compression algorithm
    '''

    def __init__(self, filename=None):
        '''
        Open the given file 
        :param filename:
        '''
    self.filename = filename
        self.openFile()

    def openFile(self):
        '''
        Check for the existence of the xim file
        and open a file handler
        '''
        try:
            # Open the binary xim file for reading
            self.f = open(self.filename, 'rb')
        except IOError:
            # No xim file by the given name exists
            print "xim file doesn't exist"

    def headerData(self):
        '''
        Header has a fixed length of 32 bytes. 
        Integers and floats are stored in little-endian format        
        '''
        self.ximHeader = dict()  # Dictionary of header values
        self.ximHeader['FormatIdentifier'] = self.f.read(8)
        self.ximHeader['FormatVersion'] = struct.unpack('<i', self.f.read(4))[0]
        self.ximHeader['Width'] = struct.unpack('<i', self.f.read(4))[0]
        self.ximHeader['Height'] = struct.unpack('<i', self.f.read(4))[0]
        self.ximHeader['BitsPerPixel'] = struct.unpack('<i', self.f.read(4))[0]
        self.ximHeader['BytesPerPixel'] = struct.unpack('<i', self.f.read(4))[0]
        self.ximHeader['CompressionIndicator'] = struct.unpack('<i', self.f.read(4))[0]


    def pixelData(self):
        '''
        Pixel values in an HND image is stored in pixelData field. Pixel data are either 
        compressed or uncompressed which can be determined by the "Compression indicator" 
        field in the header data.
        '''
        w = self.ximHeader['Width']
        h = self.ximHeader['Height']
        bpp = self.ximHeader['BytesPerPixel']

        # Image pixels are stored uncompressed in the xim image file.
        if not self.ximHeader['CompressionIndicator']:
            # Read in int4 (32 bit) image pixe values
            uncompressedPixelBufferSize = struct.unpack('<%i', self.f.read(4))[0]
            # Read in pixel values in 1D array
            uncompressedPixelBuffer = np.asarray(struct.unpack('<%ii' % (uncompressedPixelBufferSize / 4), \
                                                                    self.f.read(uncompressedPixelBufferSize)))

        # Decompress the pixelData using HND decompression algorithm.
        else:
            self.LUTSize = struct.unpack('<i', self.f.read(4))[0]  # Lookup table size
            LUT = np.asarray(struct.unpack('<%iB' % self.LUTSize, self.f.read(self.LUTSize)))  # Lookup table
            compressedBufferSize = struct.unpack('<i', self.f.read(4))[0]  # Compressed pixel buffer size
            uncompressedPixelBuffer = self.uncompressHnd(w, h, bpp, LUT)  # Uncompress the pixel data
            uncompressedBufferSize = struct.unpack('<i', self.f.read(4))[0]  # Uncompressed pixel image size

        # Reshape uncompressed image into 2D array
        self.uncompressedImage = np.reshape(uncompressedPixelBuffer, (h, w))


    def histogramData(self):

        self.histogram = dict()
        self.histogram['NumberOfBins'] = struct.unpack('<i', self.f.read(4))[0]
        self.histogram['Value'] = struct.unpack('<%ii' % self.histogram['NumberOfBins'], \
                                           self.f.read(4 * self.histogram['NumberOfBins']))


    def propertiesData(self):
        """
        Get property data for images    
        """
        self.propertyDataList = []

        value = None
        propertyValList = []

        propertyCount = struct.unpack('<i', self.f.read(4))[0]

        PROPERTY_TYPE_DICT = {0 : ('<i', 4),
                              1 : ('<d', 8),
                              2 : ('<i', 4),
                              }

        def get_value(fmt, fmt_length):
            """
            Getting integer or double value and appending it into property value list
            """
            try:
                value = struct.unpack(fmt, self.f.read(fmt_length))[0]
                propertyValList.append(value)
                return value
            except:
                return None

        if propertyCount:
            for i in xrange(propertyCount):
                if not value:
                    propertyNameLength = struct.unpack('<i', self.f.read(4))[0]
                else:
                    propertyNameLength = value
                    value = None
                    propertyValList = []

                propertyName = struct.unpack('<%is' % propertyNameLength , self.f.read(propertyNameLength))[0]
                propertyType = struct.unpack('<i', self.f.read(4))[0]


                if propertyType in PROPERTY_TYPE_DICT.keys():
                    propertyValue = struct.unpack(PROPERTY_TYPE_DICT[propertyType][0],
                                                   self.f.read(PROPERTY_TYPE_DICT[propertyType][1]))[0]

                    if propertyType == 2:
                        propertyValue = struct.unpack('<%is' % propertyValue, self.f.read(propertyValue))[0]

                    rstTpl = propertyNameLength, propertyName, propertyType, propertyValue
                    self.propertyDataList.append(rstTpl)

                elif propertyType in [4, 5]:
                        if propertyType == 4:
                            fmt, fmt_length = '<d', 8
                        else:
                            fmt, fmt_length = '<i', 4

                        value = get_value(fmt, fmt_length)
                        while value != None:
                            if len(str(value)) == 2:
                                break

                            value = get_value(fmt, fmt_length)

                        rstTpl = propertyNameLength, propertyName, propertyType, propertyValList

                        self.propertyDataList.append(rstTpl)

                        if not value:
                            break
                else:
                    print "Format Type not valid"

        else:
            print "Property not exist"

    def saveInfo(self):
        """
        Storing headerData, histogramData and propertyData into txt file and saving plot image.
        """
        kwargs = {"headerDataDict" : self.ximHeader,
                  "histogramDataDict" : self.histogram,
                  "propertyDataList" : self.propertyDataList
                  }

        self.ximFileInfoObj = XimFileInfo(**kwargs)
        self.ximFileInfoObj.saveXimInfo()


    def uncompressHnd(self, w, h, bpp, lut):
        '''
        Uncompress the xim file based on HND algorithm. The first row and the 
        first pixel of the second row are stored uncompressed. The remainders 
        of the pixels are compressed by storing only the difference between 
        neighboring pixels.
        
        E.g. consider the following hypothetical 12 pixel image:
                R11    R12    R13    R14
                R21    R22    R23    R24
                R31    R32    R33    R34
        Pixels R11 through R14 and R21 are stored uncompressed, while pixels 
        R22 through R34 are compressed by storing only the difference: 
        
        diff = R11 + R22 - R21 - R12
        
        Exploiting the fact that most images exhibit similarity in neighboring 
        pixel values, the above difference can be stored using fewer bytes, 
        e.g. 1, 2 or 4 bytes.
         
        For decompression, the algorithm needs to know the byte size of each 
        stored difference. To accomplish this, a lookup table is placed at the 
        beginning of the image. The lookup table contains a 2-bit flag for each 
        pixel which defines the byte size for each compressed pixel difference. 
        So a flag value of 0 means the difference fits into one byte while 
        1 and 2 mean a two and four byte difference respectively.
          
        :param w: Uncompressed image width
        :param h: Uncompressed image height
        :param bpp: byte per pixel
        :param lut: look up table
        '''

        # Initialize uncompressed image variable
        imagePix = np.zeros((h * w), dtype='int32')

        # Read in the first row
        ind = 0  # Index variable
        for i in xrange(w):
            imagePix[ind] = struct.unpack('<i', self.f.read(4))[0]
            ind += 1

        # ... and the first pixel of the second row
        imagePix[ind] = struct.unpack('<i', self.f.read(4))[0]
        ind += 1

        # lookup table 'bit' flag to byte conversion
        byteConversion = {'00':1, '01':2, '02':4}

        # Determine the number of unused 2-bit flag fields
        # in the last byte of the look up table
        completeBytes, partialByte = divmod((w * (h - 1) - 1), bpp)

        # Calculate current pixel value based  on "diff"
        # and adjacent pixel values as following:
        # R22 (current pixel) = diff + R21 + R12 - R11
        for i in xrange(completeBytes):

            # Convert the lookup table byte to bits
            bitFlags = '{0:08b}'.format(lut[i])
            for j in reversed(xrange(4)):

                # Determine the byte size of "diff" pixel
                byteSize = byteConversion[bitFlags[2 * j: 2 * (j + 1)]]
                # Calculate "diff" value based on the byte size
                diff = self.char2Int(byteSize)
                # R22 (current pixel) = diff + R21 + R12 - R11
                imagePix[ind] = diff + imagePix[ind - 1] + imagePix[ind - w] - imagePix[ind - w - 1]
                ind += 1

        # Calculate pixel corresponding to the last byte (partilaByte)
        # of the look up table
        bitFlags = '{0:08b}'.format(lut[i + 1])

        for j in xrange(3, 3 - partialByte, -1):

            # Determine the byte size of "diff" pixel
            byteSize = byteConversion[bitFlags[2 * j: 2 * (j + 1)]]
            # Calculate "diff" value based on the byte size
            diff = self.char2Int(byteSize)
            # R22 (current pixel) = diff + R21 + R12 - R11
            imagePix[ind] = diff + imagePix[ind - 1] + imagePix[ind - w] - imagePix[ind - w - 1]
            ind += 1

        return imagePix

    def char2Int(self, sz):
        '''
        Convert little-endian chars to a 32 bit integer
        Character size can be 1 byte: signed char 
                              2 bytes : short
                              4 bytes : int4 
        :param sz:
        '''
        if sz == 1:
            value = struct.unpack('<b', self.f.read(1))[0]  # b: signed char
        elif sz == 2:
            value = struct.unpack('<h', self.f.read(2))[0]  # h: short
        elif sz == 4:
            value = struct.unpack('<i', self.f.read(4))[0]  # i: int4

        return value

def process_arguments(args):

    # Construct the parser
    parser = ArgumentParser(description='TrueBeam(TM) xim image reader.')

    # Add expected arguments
    # Name of the xim image file
    parser.add_argument('-f', '--filename', dest='fname', type=str, required=True, \
                        help="Please enter name of the binary xim file")

    # Add image display option (optional)
    parser.add_argument('-s', '--showImage', dest='showImage', type=int, default=1,
                        help="Show xim image (Optional, 0 or 1)")
    # Version number
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    # Apply the parser to the argument list
    options = parser.parse_args(args)

    return vars(options)

def main():

    # Read the command line argument
    options = process_arguments(sys.argv[1:])
    # Create a file object
    fp = XimReader(options['fname'])
    # Read header data
    fp.headerData()
    # Read xim image, decompress if needed
    fp.pixelData()
    # Histogram data
    fp.histogramData()

    # Properties data
    fp.propertiesData()

    # Saving headerData, histogramData and propertiesData.
    fp.saveInfo()

    # Now show image
    if (options['showImage']):
        m = np.mean(fp.uncompressedImage.flatten())
        s = np.mean(fp.uncompressedImage.flatten())
        plt.imshow(fp.uncompressedImage, vmin=0, vmax=m + 0.1 * s, cmap=plt.gray())
        plt.savefig(os.path.join(fp.ximFileInfoObj.outputFolder, XIMREADER_IMG_NAME))
        print "%s is stored on %s" % (XIMREADER_IMG_NAME, fp.ximFileInfoObj.outputFolder)
        plt.show()

if __name__ == "__main__":
    # Let's get started
    main()
