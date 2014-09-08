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
#ifndef __retin_ImageResample_h__
#define __retin_ImageResample_h__

#include <stddef.h>
#include <math.h>

namespace retin {
template<class Value>
static inline Value Lanczos(Value x,float radius=3.0)
{
  if (x == 0.0) return 1.0;
  float pi_x = x * M_PI;
  return radius * sin(pi_x) * sin(pi_x / radius) / (pi_x * pi_x);
}

template<class Value>
void	resampleRowsTrans (Value* dest,uint outputWidth,const Value* source,uint inputWidth,uint inputHeight,uint channels,const uint radius=8)
{
	float blur = 1.0;
	float factor = float(outputWidth) / float(inputWidth);
 	float scale   = std::min(factor, 1.0f) / blur;
	float support = radius / scale;
	if (support <= 0.5) { support = 0.5 + 1E-12; scale = 1.0; }

	uint filterSize = std::min(inputWidth, 5+uint(2*support));
	std::vector<Value> filter (filterSize);
	std::vector<Value> temp(channels);

	for (uint x=0;x<outputWidth;x++)
	{
		float center = (x+0.5) / factor;
		uint start = (uint)std::max(center-support+0.5f, (float)0);
		uint stop = (uint)std::min(center+support+0.5f, (float)inputWidth);
		uint steps = stop - start;
		float s = start - center + 0.5;
		float density = 0.0;
		for (uint n=0;n<steps;n++,s++)
		{
			float w = Lanczos (s*scale,radius);
			filter[n] = w;
			density += w;
		}
		if (fabs(density) > 1E-7)
			density = 1.0 / density;
		else
			density = 1.0;
		for (uint n=0;n<steps;n++)
			filter[n] *= density;

		for (uint y=0;y<inputHeight;y++)
		{
			for (uint k=0;k<channels;k++)
				temp[k] = 0;
			for (uint n=0;n<steps;n++)
			{
				float w = filter[n];
				for (uint k=0;k<channels;k++)
				  temp[k] += w*source[ k+(start+n+y*inputWidth)*channels ];
			}
			for (uint k=0;k<channels;k++)
				dest[ k+(y+x*inputHeight)*channels ] = temp[k];
		}
	}
}


template<class Value>
void resample_image (Value*& out, Value* in, size_t& h, size_t& w, int channels, int maxSize)
{
	int nh, nw;

	if (w <= maxSize && h <= maxSize)
	{
                out = new Value [w*h*channels];
                memcpy(out, in, sizeof(Value) * w * h * channels);
                return;
	}

	if (w >= h) {
		nh = maxSize*h/w;
		nw = maxSize;
	}
	else {
		nw = maxSize*w/h;
		nh = maxSize;
	}

	Value* temp = new Value [nw*h*channels];
	resampleRowsTrans (temp,nw,in,w,h,channels);
	out = new Value [nw*nh*channels];
	resampleRowsTrans (out,nh,temp,h,nw,channels);
	h = nh;
	w = nw;
	delete [] temp;
}

}

#endif
