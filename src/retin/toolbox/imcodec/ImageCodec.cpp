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
 * \file ImageCodec.cpp
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#include "ImageCodec.h"

#include <iostream>
#include <cctype>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>

#ifdef RETIN_ENABLE_PNG
#include "ImageCodecPNG.h"
#endif

#ifdef RETIN_ENABLE_JPEG
#include "ImageCodecJPEG.h"
#endif

#ifdef RETIN_ENABLE_PBM
#include "ImageCodecPBM.h"
#endif

namespace retin {
    
struct tolower
{ 
    char operator()(char c) const 
    {
        return std::tolower(static_cast<unsigned char>(c));
    } 
};

std::string extension (const std::string& str) {
    size_t idx = str.rfind('.');
    std::string s;
    if (idx == std::string::npos)
        s = std::string();
    s = str.substr(idx);
    std::transform(s.begin(), s.end(), s.begin(), tolower());
    return s;
}

void	ImageCodec::throwException (const char* msg)
{
        std::string fmt;
        fmt  =  "ImageCodec(";
        fmt += strFilename;
        fmt += ") : ";
        fmt += msg;
	throw std::runtime_error (fmt);
}

void	ImageCodec::openFile (const char* pcAccess)
{
	closeFile ();
	file = fopen (strFilename.c_str(),pcAccess);
	if (!file)
		throwException ("Erreur lors de l'ouverture du fichier");
}

void	ImageCodec::closeFile ()
{
	if (file)
	{
		fclose(file);
		file = NULL;
	}
}

ImageEncoder*	getImageEncoder (const std::string& strFile)
{
	std::string strExt = extension(strFile);
#ifdef RETIN_ENABLE_PNG
	if (strExt == ".png")
		return new ImageEncoderPNG();
#endif
	return NULL;
}


ImageDecoder*	getImageDecoder (const std::string& strFile)
{
	std::string strExt = extension(strFile);
#ifdef RETIN_ENABLE_PNG
	if (strExt == ".png")
		return new ImageDecoderPNG();
#endif
#ifdef RETIN_ENABLE_JPEG
	if (strExt == ".jpeg" || strExt == ".jpg")
		return new ImageDecoderJPEG();
#endif
#ifdef RETIN_ENABLE_PBM
        if (strExt == ".pgm")
                return new ImageDecoderPBM();
#endif
	return NULL;

}

void		loadImage (unsigned char*& pixels,unsigned char*& palette,size_t& width,size_t& height,size_t& channels,const std::string& strFile)
{
    std::auto_ptr<ImageDecoder> decoder ( getImageDecoder (strFile) );
    if (!decoder.get())
            throw std::runtime_error("Unsupported file format "+strFile);
    decoder->read_open(strFile);
    decoder->read_header();
    width = decoder->width();
    height = decoder->height();
    channels = decoder->channels();
    if (decoder->palette())
            palette = new unsigned char[256*3];
    pixels = new unsigned char[width*height*channels];
    decoder->read_data (pixels,palette);
    decoder->read_close ();
}

void		loadImageGrayscale(unsigned char*& pixels, size_t& width, size_t& height, const std::string& strFile) {
	unsigned char *_pixels, *_palette, index, r, g, b;
	size_t _channels;
	
	_palette = NULL;
	loadImage(_pixels, _palette, width, height, _channels, strFile);

	// indexed
	if(_palette) {
		pixels = new unsigned char[width*height];

		for(size_t i = 0 ; i < width*height ; i++) {
			index = _pixels[i];
			r = _palette[index*3+0];
			g = _palette[index*3+1];
			b = _palette[index*3+2];

			pixels[i] = round((r+g+b)/3.0);
		}
	}
	// color
	else if(_channels == 3) {
		pixels = new unsigned char[width*height];
		rgb2int(pixels, _pixels, width*height);
	}
	// already grayscale
	else if(_channels == 1) {
		pixels = _pixels;
	}
	// whatever else
	else {
		pixels = NULL;
		width = 0;
		height = 0;
	}
}

void		loadImageRGB(unsigned char*& pixels, size_t& width, size_t& height, const std::string& strFile) {
	unsigned char *_pixels, *_palette, index, r, g, b;
	size_t _channels;
	
	_palette = NULL;
	loadImage(_pixels, _palette, width, height, _channels, strFile);

	// indexed
	if(_palette) {
		pixels = new unsigned char[3*width*height];

		for(size_t i = 0 ; i < width*height ; i++) {
			index = _pixels[i];
			r = _palette[index*3+0];
			g = _palette[index*3+1];
			b = _palette[index*3+2];

			pixels[3*i+0] = r;
			pixels[3*i+1] = g;
			pixels[3*i+2] = b;			
		}
	}
	// already color
	else if(_channels == 3) {
		pixels = _pixels;
	}
	// grayscale
	else if(_channels == 1) {
		pixels = new unsigned char[3*width*height];
		for(size_t i = 0 ; i < width*height ; i++) {
			pixels[3*i+0] = _pixels[i];
			pixels[3*i+1] = _pixels[i];
			pixels[3*i+2] = _pixels[i];
		}
	}
	// whatever else
	else {
		pixels = NULL;
		width = 0;
		height = 0;
	}
}
void		saveImage (const unsigned char* pixels,const unsigned char* palette,size_t width,size_t height,size_t channels,const std::string& strFile)
{
    std::auto_ptr<ImageEncoder> encoder ( getImageEncoder (strFile) );
    if (!encoder.get())
            throw std::runtime_error("Unsupported file format "+strFile);

    encoder->setWidth (width);
    encoder->setHeight (height);
    encoder->setChannels (channels);
    encoder->setPalette (palette);

    encoder->write_open(strFile);
    encoder->write_header();
    encoder->write_data (pixels,palette);
    encoder->write_close ();
}
 
void		saveImage (const float* v,size_t width,size_t height,const std::string& strFile)
{
    size_t i,n=width*height;
    float vMin = v[0];
    float vMax = v[0];
    for (i=1;i<n;i++)
    {
        if (vMin > v[i]) vMin = v[i];
        if (vMax < v[i]) vMax = v[i];
    }
    if (vMax == vMin) vMax = vMin+1;
    std::vector<unsigned char> image(width*height);
    for (i=0;i<n;i++) {
        float x = 255*(v[i] - vMin)/(vMax - vMin);
        if (x < 0) x = 0;
        else if (x > 255) x = 255;
        image[i] = (unsigned char)x;
    }
    saveImage(&image[0],NULL,width,height,1,strFile);
}

void            checkImageCodec (const std::string& inputFile,const std::string& outputFile)
{
    unsigned char *pixels = NULL, *palette = NULL;
    size_t  width,height,channels;
    loadImage(pixels,palette,width,height,channels,inputFile);
    saveImage(pixels,palette,width,height,channels,outputFile);
    if (pixels) delete[] pixels;
    if (palette) delete[] palette;
}


}
