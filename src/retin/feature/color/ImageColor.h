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
#ifndef __ImageColor_h__
#define __ImageColor_h__

#include "ImageProcessing.h"

#include "retin/toolbox/document/Parameters.h"
#include "retin/toolbox/document/BytesMatrix.h"
#include "retin/toolbox/document/FeatureMatrix.h"

namespace retin {

template<class Matrix>
void		convertColors (boost::shared_ptr<Matrix>& output,boost::shared_ptr<Matrix> input,ParametersRef params)
{
	std::string inputSpace = params->get("inputSpace");
	if (inputSpace == "rgb") {
	    if (input->dim() == 1)
		inputSpace = "int";
	    else if (input->dim() != 3)
		retinThrowException1("Invalid input dimension (%d)",input->dim());
	}
	else if (inputSpace == "lab") {
	    if (input->dim() != 3)
		retinThrowException1("Invalid input dimension (%d)",input->dim());
	}
	else if (inputSpace == "int") {
	    if (input->dim() != 1)
		retinThrowException1("Invalid input dimension (%d)",input->dim());
	}
	else
	    retinThrowException1("Invalid input space %s",inputSpace.c_str());

	std::string outputSpace = params->get("outputSpace");

	int outputDim;
	if (outputSpace == "lab" || outputSpace == "rgb")
		outputDim = 3;
	else if (outputSpace == "ab")
		outputDim = 2;
	else if (outputSpace == "lum" || outputSpace == "int")
		outputDim = 1;
	else
		retinThrowException1("Invalid output space %s",outputSpace.c_str());

 	output = boost::make_shared<Matrix>(input->width(),input->height(),outputDim);

	if (inputSpace == outputSpace)
	{
		output = input;
	}
	else if (inputSpace == "rgb" && outputSpace == "lab")
	{
		rgb2lab (output->data(),input->data(),input->width()*input->height());
	}
	else if (inputSpace == "rgb" && outputSpace == "lum")
	{
		rgb2lum (output->data(),input->data(),input->width()*input->height());
	}
	else if (inputSpace == "rgb" && outputSpace == "ab")
	{
		rgb2ab (output->data(),input->data(),input->width()*input->height());
	}
	else if (inputSpace == "rgb" && outputSpace == "int")
	{
		rgb2int (output->data(),input->data(),input->width()*input->height());
	}
	else if (inputSpace == "int" && outputSpace == "lab")
	{
		int2lab (output->data(),input->data(),input->width()*input->height());
	}
	else if (inputSpace == "int" && outputSpace == "lum")
	{
		int2lum (output->data(),input->data(),input->width()*input->height());
	}
	else if (inputSpace == "lab" && outputSpace == "rgb")
	{
		lab2rgb (output->data(),input->data(),input->width()*input->height());
	}
	else
		retinThrowException2("Conversion %s => %s unsupported",inputSpace.c_str(),outputSpace.c_str());

	if (params->has("equalize"))
		equalize(output->data(),output->width()*output->height(),output->dim());
}

}

#include <jni.h>

extern "C" {

JNIEXPORT jobject JNICALL Java_retin_feature_color_ImageColor_runBytesMatrix
  (JNIEnv *, jobject, jobject);

JNIEXPORT jobject JNICALL Java_retin_feature_color_ImageColor_runFeatureMatrix
  (JNIEnv *, jobject, jobject);

}

#endif
