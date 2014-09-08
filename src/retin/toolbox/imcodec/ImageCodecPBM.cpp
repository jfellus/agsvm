/*
Copyright © CNRS 2012. 
Authors: Philippe-Henri Gosselin, David Picard, Romain Négrel
Contact: philippe-henri.gosselin@ensea.fr

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/
/**
 * \file ImageCodecPBM.cpp
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#include "ImageCodecPBM.h"

#ifdef RETIN_ENABLE_PBM

#include <string.h>

namespace retin {

void    ImageDecoderPBM::skipComment(char* c)
{
    while(*c == '#')
    {
        while(fread(c,1,1,file) == 1 && *c != '\n');
        if (fread(c,1,1,file) != 1)
        throwException ("Read error");
    }
}

void    ImageDecoderPBM::readLine(char* buf)
{
    int i = 1;
    while( (fread(buf+i,1,1,file) == 1 ) && (buf[i] != '\n') && (buf[i] != '\r') && (i<79) ) i++;
    buf[i] = 0;
}

void	ImageDecoderPBM::read_open (const std::string& filename)
{
    strFilename = filename;
    openFile ("rb");
}

void	ImageDecoderPBM::read_header ()
{
    char header[4];
    header[3] = '\0';
    if (fread(header,1,3,file) != 3)
        throwException ("Read error: header");
    if (strcmp(header,"P5\n") == 0) {
        mode = 5;
        nChannels = 1;
    }
    else {
        throwException ("Unsupported format");
    }
    bHasPalette = false;
    
    char buf[0x100];
    if (fread(buf,1,1,file) != 1)
        throwException ("Read error: size");
    skipComment(buf);
    readLine(buf);
    int x,y;
    sscanf(buf, "%d %d", &x,&y);
    nWidth = x;
    nHeight = y;
    if (fread(buf,1,1,file) != 1)
        throwException ("Read error: max");
    skipComment(buf);
    readLine(buf);
}

void	ImageDecoderPBM::read_data (unsigned char* pixels,unsigned char* palette)
{
    if (mode == 5) {
        if (fread(pixels,1,nWidth*nHeight,file) != nWidth*nHeight)
            throwException ("Read error");
    }
    else {
        throwException ("Unsupported mode");
    }
}

void	ImageDecoderPBM::read_close ()
{
    closeFile ();
}


}

#endif

