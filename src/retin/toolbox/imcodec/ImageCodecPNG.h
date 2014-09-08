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
 * \file ImageCodecPNG.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __ImageCodecPNG_h__
#define __ImageCodecPNG_h__

#include "ImageCodec.h"

#ifdef RETIN_ENABLE_PNG

extern "C" {
#include <png.h>
}

namespace retin {

//! Encodeur PNG.
/*! \ingroup ImageCodec
*/
class	ImageEncoderPNG : public ImageEncoder
{
	png_structp png_ptr;
	png_infop info_ptr;
public:
	~ImageEncoderPNG () { write_close(); }

virtual	void		write_open (const std::string& filename);
virtual	void		write_header ();
virtual	void		write_data (const unsigned char* pixels,const unsigned char* palette);
virtual	void		write_close ();

};

//! Decodeur PNG.
/*! \ingroup ImageCodec
*/
class	ImageDecoderPNG : public ImageDecoder
{
	png_structp png_ptr;
	png_infop info_ptr;
public:
	~ImageDecoderPNG () { read_close(); }

virtual	void		read_open (const std::string& filename);
virtual	void		read_header ();
virtual	void		read_data (unsigned char* pixels,unsigned char* palette);
virtual	void		read_close ();

};

}

#endif

#endif
