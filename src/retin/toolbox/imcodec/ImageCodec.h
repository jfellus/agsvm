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
 * \file ImageCodec.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __ImageCodec_h__
#define __ImageCodec_h__

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <stdexcept>

#include "retin/feature/color/ImageProcessing.h"

namespace retin {

//! Classe mère pour les Encodeur/Decodeur d'images.
/*! \ingroup ImageCodec
*/
class	ImageCodec
{
protected:
	std::string	strFilename;
	FILE*		file;
	bool		bHasPalette;
	size_t		nWidth,nHeight,nChannels;

	void		throwException (const char* msg);
	void		openFile (const char* pcAccess);
	void		closeFile ();

public:
	ImageCodec () : file(NULL),bHasPalette(false),nWidth(0),nHeight(0),nChannels(0) { }
virtual	~ImageCodec () { if (file) fclose(file); }

	bool		palette () const { return bHasPalette; }
	size_t		width () const { return nWidth; }
	size_t		height () const { return nHeight; }
	size_t		channels () const { return nChannels; }

	void		setPalette (bool b) { bHasPalette = b; }
	void		setWidth (size_t x) { nWidth = x; }
	void		setHeight (size_t x) { nHeight = x; }
	void		setChannels (size_t x) { nChannels = x; }

};

class	ImageEncoder : public ImageCodec
{
public:

	//! Ouvre le fichier en lecture et prépare l'encodage.
virtual	void		write_open (const std::string& filename) = 0;
	//! Ecrit l'entête du fichier (type, taille, etc.).
virtual	void		write_header () = 0;
	//! Ecrit les données.
virtual	void		write_data (const unsigned char* pixels,const unsigned char* palette=NULL) = 0;
	//! Ferme le fichier et détruit les données temporaires.
virtual	void		write_close () = 0;

};

class	ImageDecoder : public ImageCodec
{
public:

	//! Ouvre le fichier en lecture et prépare le décodage.
virtual	void		read_open (const std::string& filename) = 0;
	//! Lit l'entête du fichier (type, taille, etc.).
virtual	void		read_header () = 0;
	//! Lit les données.
virtual	void		read_data (unsigned char* pixels,unsigned char* palette=NULL) = 0;
	//! Ferme le fichier et détruit les données temporaires.
virtual	void		read_close () = 0;

};

//! Renvoie le codec correspondant à une extention de fichier.
ImageEncoder*	getImageEncoder (const std::string& strFile);
//! Renvoie le codec correspondant à une extention de fichier.
ImageDecoder*	getImageDecoder (const std::string& strFile);

//! Charge une image.
void		loadImage (unsigned char*& pixels,unsigned char*& palette,size_t& width,size_t& height,size_t& channels,const std::string& strFile);
void		loadImageGrayscale(unsigned char*& pixels, size_t& width, size_t& height, const std::string& strFile);
void		loadImageRGB(unsigned char*& pixels, size_t& width, size_t& height, const std::string& strFile);
//! Sauve une image.
void		saveImage (const unsigned char* pixels,const unsigned char* palette,size_t width,size_t height,size_t channels,const std::string& strFile);
//! Sauve une image en noir et blanc, avec un recalage automatique des valeurs entre 0 et 255.
void		saveImage (const float* pixels,size_t width,size_t height,const std::string& strFile);


//! Vérifie que le codec fonctionne bien.
void            checkImageCodec (const std::string& inputFile,const std::string& outputFile);

}

#endif
