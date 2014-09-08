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

#include "ImageProcessing.h"
#include "retin/toolbox/document/Document.h"
#include <math.h>

using namespace std;
using namespace retin;

template<class Matrix>
void		equalize1n (boost::shared_ptr<Matrix> input)
{
	typedef typename Matrix::value_type Value;

	double min = 1E9;
	double max = 0;
	double mean = 0;
	double sqr = 0;
	Value* data = input->data();
	uint n = input->width()*input->height();
	uint p = input->dim();
	for (uint i=0;i<n;i++) {
		double x = read<Value>(data[i*p]);
		mean += x;
		sqr += x*x;
		if (x > max) max = x;
		if (x < min) min = x;
	}

	mean /= n;
	sqr = sqrt(sqr/n - mean*mean);
	/*double min = mean - 2*sqr;
	double max = mean + sqr;*/
	//printf ("%f < %f +/- %f < %f\n",min,mean,sqr,max);

	for (uint i=0;i<n;i++) {
	    double x = read<Value>(data[i*p]);
	    data[i*p] = convert<Value>( (x-min)/(max-min) );
	}

}


#include <jni.h>

extern "C" {

JNIEXPORT jobject JNICALL Java_retin_feature_color_ImageEqualization_runBytesMatrix
  (JNIEnv *, jobject, jobject);

JNIEXPORT jobject JNICALL Java_retin_feature_color_ImageEqualization_runFeatureMatrix
  (JNIEnv *, jobject, jobject);

}
