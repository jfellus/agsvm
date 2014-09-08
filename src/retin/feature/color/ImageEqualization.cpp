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
#include "ImageEqualization.h"
#include "ImageProcessing.h"
#include "retin/toolbox/document/serialization/Java.h"
#include "retin/toolbox/document/Parameters.h"
#include "retin/toolbox/document/BytesMatrix.h"
#include "retin/toolbox/document/FeatureMatrix.h"
#include <math.h>

using namespace std;
using namespace boost;
using namespace retin;

template<class Matrix>
void		equalizeColors (boost::shared_ptr<Matrix>& output,boost::shared_ptr<Matrix> input,ParametersRef params)
{
        
	if (input->dim() != 3 && input->dim() != 1)
		retinThrowException1("Invalid input dimension (%d)",input->dim());

	string inputSpace = params->get("inputSpace");
	string workSpace = params->get("workSpace");

	shared_ptr<FeatureMatrix> temp;
	if (inputSpace == "rgb" && workSpace == "lab") {
	    temp = make_shared<FeatureMatrix>(input->width(),input->height(),3);
            if(input->dim() == 3)
                rgb2lab (temp->data(), input->data(), input->width()*input->height());
            else
                int2lab (temp->data(), input->data(), input->width()*input->height());
	    equalize1n (temp);
	    output = make_shared<Matrix>(input->width(),input->height(),3);
            lab2rgb (output->data(),temp->data(),temp->width()*temp->height());	    
	}
	else
	    retinThrowException2("Invalid combination or input spase %s and work space %s",inputSpace.c_str(),workSpace.c_str());
}



JNIEXPORT jobject JNICALL Java_retin_feature_color_ImageEqualization_runFeatureMatrix
  (JNIEnv * env, jobject thisObject, jobject inputObj)
{
	if (!inputObj)
		return NULL;

	RETIN_JAVA_METHOD_BEGIN

	JavaEnv je(env);
	ParametersRef params = je.getDocumentField<Parameters> (thisObject,"params");
	FeatureMatrixRef input = je.createDocument<FeatureMatrix>(inputObj);
 	FeatureMatrixRef output;
	equalizeColors<FeatureMatrix> (output,input,params);
	return je.createObject (output);

	RETIN_JAVA_METHOD_END
	return NULL;
}

JNIEXPORT jobject JNICALL Java_retin_feature_color_ImageEqualization_runBytesMatrix
  (JNIEnv * env, jobject thisObject, jobject inputObj)
{
	if (!inputObj)
		return NULL;

	RETIN_JAVA_METHOD_BEGIN

	JavaEnv je(env);
	ParametersRef params = je.getDocumentField<Parameters> (thisObject,"params");
	BytesMatrixRef input = je.createDocument<BytesMatrix>(inputObj);
 	BytesMatrixRef output;
	equalizeColors<BytesMatrix> (output,input,params);
	return je.createObject (output);

	RETIN_JAVA_METHOD_END
	return NULL;
}
