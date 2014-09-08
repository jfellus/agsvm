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
 * \file ImageCodecJPEG.cpp
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#include "ImageCodecJPEG.h"

#ifdef RETIN_ENABLE_JPEG

#include <setjmp.h>
#include <string.h>

#ifdef RETIN_ENABLE_BOOST
#include <boost/thread.hpp>
#endif

namespace retin {

#ifdef RETIN_ENABLE_BOOST
boost::mutex	jpegMutex;
#endif
    
struct my_error_mgr {
        struct jpeg_error_mgr pub;
        jmp_buf setjmp_buffer;
};

void my_error_exit (j_common_ptr cinfo)
{
        my_error_mgr* myerr = (my_error_mgr*) cinfo->err;
        (*cinfo->err->output_message) (cinfo);
        longjmp(myerr->setjmp_buffer, 1);
}
    
void	ImageDecoderJPEG::read_open (const std::string& filename)
{
	strFilename = filename;
	openFile ("rb");
}

void	ImageDecoderJPEG::read_header ()
{
        #ifdef RETIN_ENABLE_BOOST
        boost::mutex::scoped_lock lock(jpegMutex);
        #endif
    
        struct jpeg_decompress_struct cinfo;
        struct my_error_mgr jerr;
        
 	cinfo.err = jpeg_std_error(&jerr.pub);
        jpeg_create_decompress(&cinfo);

        cinfo.err = jpeg_std_error(&jerr.pub);
        jerr.pub.error_exit = my_error_exit;
        if (setjmp(jerr.setjmp_buffer)) {
            jpeg_destroy_decompress(&cinfo);
            throwException ("Error");
        }
        
        jpeg_stdio_src(&cinfo, file);        
        jpeg_read_header(&cinfo, TRUE);        
        jpeg_start_decompress(&cinfo);
        
        bHasPalette = false;
        nWidth = cinfo.output_width;
        nHeight = cinfo.output_height;
        nChannels = cinfo.output_components;
        
        if (image) delete[] image;
        image = new unsigned char[nWidth*nHeight*nChannels];
        if (!image) {
            jpeg_destroy_decompress(&cinfo);
            throwException ("Out of memory");
        }
        int row_stride = cinfo.output_width * cinfo.output_components;
        JSAMPARRAY buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
        unsigned char* ptr = image;
        while (cinfo.output_scanline < cinfo.output_height) {
            jpeg_read_scanlines(&cinfo, buffer, 1);
            memcpy ( ptr, buffer[0], row_stride );
            ptr += row_stride;
        }
        
        jpeg_finish_decompress(&cinfo);
        jpeg_destroy_decompress(&cinfo);
}

void	ImageDecoderJPEG::read_data (unsigned char* pixels,unsigned char* palette)
{
    if (!image)
        throwException ("no image");
    memcpy(pixels,image,nWidth*nHeight*nChannels);
}

void	ImageDecoderJPEG::read_close ()
{
    if (image) {
        delete[] image;
        image = NULL;
    }
    closeFile ();
}


}

#endif
