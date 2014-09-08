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
 * \file ImageCodecPNG.cpp
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#include "ImageCodecPNG.h"

#ifdef RETIN_ENABLE_PNG

namespace retin {

void	ImageDecoderPNG::read_open (const std::string& filename)
{
	png_ptr = NULL;
	info_ptr = NULL;

	strFilename = filename;
	openFile ("rb");

}

void	ImageDecoderPNG::read_header ()
{
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
	if (png_ptr == NULL)
		throwException ("Erreur d'allocation");

	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL)
		throwException ("Erreur d'allocation");

	if (setjmp(png_jmpbuf(png_ptr)))
		throwException ("Erreur d'allocation");

	png_init_io(png_ptr, file);

	png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING, png_voidp_NULL);

	if (info_ptr->bit_depth != 8)
		throwException ("Seules les images en précision 8-bit sont supportées");

	nWidth = info_ptr->width;
	nHeight = info_ptr->height;

	if (png_ptr->color_type == PNG_COLOR_TYPE_RGB)
	{
		nChannels = 3;
		bHasPalette = false;
	}
	else if (png_ptr->color_type == PNG_COLOR_TYPE_RGB_ALPHA)
	{
		nChannels = 4;
		bHasPalette = false;
	}
	else if (png_ptr->color_type == PNG_COLOR_TYPE_GRAY)
	{
		nChannels = 1;
		bHasPalette = false;
	}
	else if (png_ptr->color_type == PNG_COLOR_TYPE_PALETTE)
	{
		nChannels = 1;
		bHasPalette = true;
	}
	else
		throwException ("Format des pixels non supporté");
}

void	ImageDecoderPNG::read_data (unsigned char* pixels,unsigned char* palette)
{
	png_bytep* row_pointers = png_get_rows(png_ptr, info_ptr);

	for (size_t j=0;j<nHeight;j++)
	{
		memcpy (pixels+j*nChannels*nWidth,row_pointers[j],nChannels*nWidth);
	}
	if (bHasPalette && palette)
	{
		int n;
		png_colorp png_matricePalette;
		png_get_PLTE (png_ptr,info_ptr,&png_matricePalette,&n);
		for (int j=0;j<n;j++)
		{
			palette[j*3+0] = png_matricePalette[j].red;
			palette[j*3+1] = png_matricePalette[j].green;
			palette[j*3+2] = png_matricePalette[j].blue;
		}
	}
}

void	ImageDecoderPNG::read_close ()
{
	if (png_ptr)
	{
		if (info_ptr)
			png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
		else
			png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
		png_ptr = NULL;
		info_ptr = NULL;
	}
	closeFile ();
}

// ***********************************************************************************

void	ImageEncoderPNG::write_open (const std::string& filename)
{
	png_ptr = NULL;
	info_ptr = NULL;

	strFilename = filename;
	openFile ("wb");
}


void	ImageEncoderPNG::write_header ()
{
	png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, 
		NULL,NULL,NULL);
    	if (!png_ptr)
		throwException ("Erreur d'allocation");

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
		throwException ("Erreur d'allocation");

	if (setjmp(png_jmpbuf(png_ptr)))
		throwException ("Erreur");
	
	info_ptr->width = nWidth;
	info_ptr->height = nHeight;
	info_ptr->bit_depth = 8;
	info_ptr->interlace_type = PNG_INTERLACE_NONE;
	info_ptr->compression_type = PNG_COMPRESSION_TYPE_DEFAULT;

	if (bHasPalette)
		info_ptr->color_type = PNG_COLOR_TYPE_PALETTE;
	else if (nChannels == 1)
		info_ptr->color_type = PNG_COLOR_TYPE_GRAY;
	else if (nChannels == 3)
		info_ptr->color_type = PNG_COLOR_TYPE_RGB;
	else if (nChannels == 4)
		info_ptr->color_type = PNG_COLOR_TYPE_RGB_ALPHA;
	else
		throwException ("Format non supporté");
}

void	ImageEncoderPNG::write_data (const unsigned char* pixels,const unsigned char* palette)
{
	png_bytep* row_pointers = new png_bytep[nHeight];
	for (size_t j=0;j<nHeight;j++)
	{
		row_pointers[j] = new png_byte[nWidth*nChannels];
		memcpy (row_pointers[j],pixels+j*nChannels*nWidth,nChannels*nWidth);
	}
	png_set_rows(png_ptr, info_ptr, row_pointers);

	if (bHasPalette && palette)
	{
		int n = 256;
		png_color* png_matricePalette = new png_color[n];
		for (int j=0;j<n;j++)
		{
			png_matricePalette[j].red = palette[j*3+0];
			png_matricePalette[j].green = palette[j*3+1];
			png_matricePalette[j].blue = palette[j*3+2];
		}
		png_set_PLTE (png_ptr,info_ptr,png_matricePalette,256);
                delete[] png_matricePalette;
	}

	png_init_io(png_ptr, file);
	png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_INVERT_ALPHA, NULL);

	for (size_t j=0;j<nHeight;j++)
		delete [] row_pointers[j];
        delete[] row_pointers;

}

void	ImageEncoderPNG::write_close ()
{
	if (png_ptr)
	{
		if (info_ptr)
			png_destroy_write_struct(&png_ptr, &info_ptr);
		else
			png_destroy_write_struct(&png_ptr, png_infopp_NULL);
		png_ptr = NULL;
		info_ptr = NULL;
	}
	closeFile ();
}

}

#endif
