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
#ifndef __retin_ImageProcessing_h__
#define __retin_ImageProcessing_h__

#include <stddef.h>
#include <math.h>

namespace retin {

inline double flab (double t)
{
	return (t > 0.008856)?pow(t,1.0/3.0):(7.787*t+16.0/116.0);
}

template<class Value>
inline double	read(Value x) { return x; }
template<>
inline double	read(unsigned char x) { return x / 255.0; }

template<class Value>
inline double	read2(Value x) { return x; }
template<>
inline double	read2(unsigned char x) { return x / 255.0 - 0.5; }

template<class Value>
inline Value	convert(double x) { return x; }
template<>
inline unsigned char	convert(double x) {
	x *= 255.0;
	if (x < 0.0) x = 0.0;
	if (x > 255.0) x = 255.0;
	return round(x);
}

template<class Value>
inline Value	convert2(double x) { return x; }
template<>
inline unsigned char	convert2(double x) {
	x = (x+0.5) * 255.0;
	if (x < 0.0) x = 0.0;
	if (x > 255.0) x = 255.0;
	return round(x);
}



template<class Value>
void		equalize (Value* data,size_t n,size_t dim)
{
        double mean = 0;
	double sqr = 0;
	for (size_t i=0;i<n*dim;i++) {
		double x = read<Value>(data[i]);
		mean += x;
		sqr += x*x;
	}

	mean /= dim*n;
	sqr = 2*sqrt(sqr/(dim*n) - mean*mean);
	double min = mean - sqr;
	double max = mean + sqr;
//	printf ("%f < %f +/- %f < %f\n",min,mean,sqr,max);

	for (size_t i=0;i<n*dim;i++) {
		double x = read<Value>(data[i]);
		data[i] = convert<Value>( (x-min)/(max-min) );
	}

}


template<class ValueOut,class ValueIn>
void	rgb2lab (ValueOut* out,const ValueIn* in,int n)
{
	for (int i=0;i<n;i++)
	{
		double r,g,b,X,Y,Z,L,as,bs;
		r = read<ValueIn>(in[3*i+0]);
		g = read<ValueIn>(in[3*i+1]);
		b = read<ValueIn>(in[3*i+2]);
                
		X = 0.412453 * r + 0.357580 * g + 0.189423 * b;
		Y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
		Z = 0.019334 * r + 0.119193 * g + 0.950227 * b;

		L = (Y>0.008856)?1.16*pow(Y,1.0/3.0)-0.16:9.033*Y;
		as = 5*(flab(X)-flab(Y));
		bs = 2*(flab(Y)-flab(Z));

		out[3*i+0] = convert<ValueOut>(L);
		out[3*i+1] = convert2<ValueOut>(as);
		out[3*i+2] = convert2<ValueOut>(bs);
	}
}

template<class ValueOut,class ValueIn>
void	int2lab (ValueOut* out,const ValueIn* in,int n)
{
	for (int i=0;i<n;i++)
	{
		double g,X,Y,Z,L,as,bs;
		g = read<ValueIn>(in[i]);

		X = 0.412453 * g + 0.357580 * g + 0.189423 * g;
		Y = 0.212671 * g + 0.715160 * g + 0.072169 * g;
		Z = 0.019334 * g + 0.119193 * g + 0.950227 * g;

		L = (Y>0.008856)?1.16*pow(Y,1.0/3.0)-0.16:9.033*Y;
		as = 5*(flab(X)-flab(Y));
		bs = 2*(flab(Y)-flab(Z));

		out[3*i+0] = convert<ValueOut>(L);
		out[3*i+1] = convert2<ValueOut>(as);
		out[3*i+2] = convert2<ValueOut>(bs);
	}
}

template<class ValueOut,class ValueIn>
void	int2lum (ValueOut* out,const ValueIn* in,int n)
{
	for (int i=0;i<n;i++)
	{
		double g,Y,L;
		g = read<ValueIn>(in[i]);

		Y = 0.212671 * g + 0.715160 * g + 0.072169 * g;

		L = (Y>0.008856)?1.16*pow(Y,1.0/3.0)-0.16:9.033*Y;

		out[i] = convert<ValueOut>(L);
	}
}

template<class ValueOut,class ValueIn>
void	lab2rgb (ValueOut* out,const ValueIn* in,int n)
{
	for (int i=0;i<n;i++)
	{
		double L,as,bs,X,Y,Z,P,r,g,b;
		L = read<ValueIn>(in[3*i+0]) * 100.0;
		as = read2<ValueIn>(in[3*i+1]) * 100.0;
		bs = read2<ValueIn>(in[3*i+2]) * 100.0;

		P = (L+16)/116;
		X = pow(P+as/500,3);
		Y = pow(P,3);
		Z = pow(P-bs/200,3);

		r =  3.2388593*X - 1.5312038*Y - 0.5293567*Z;
		g = -0.9687698*X + 1.874211*Y  + 0.0507745*Z;
		b =  0.0556188*X - 0.2039392*Y + 1.0567818*Z;

		out[3*i+0] = convert<ValueOut>(r);
		out[3*i+1] = convert<ValueOut>(g);
		out[3*i+2] = convert<ValueOut>(b);
	}
}

template<class ValueOut,class ValueIn>
void	rgb2lum (ValueOut* out,const ValueIn* in,int n)
{
	for (int i=0;i<n;i++)
	{
		double r,g,b,Y,L;
		r = read<ValueIn>(in[3*i+0]);
		g = read<ValueIn>(in[3*i+1]);
		b = read<ValueIn>(in[3*i+2]);

		Y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
		L = (Y>0.008856)?1.16*pow(Y,1.0/3.0)-0.16:9.033*Y;

		out[i] = convert<ValueOut>(L);
	}
}

template<class Value>
void	rgb2ab (Value* out,const Value* in,int n)
{
	for (int i=0;i<n;i++)
	{
		double r,g,b,X,Y,Z,a;
		r = read<Value>(in[3*i+0]);
		g = read<Value>(in[3*i+1]);
		b = read<Value>(in[3*i+2]);

		X = 0.412453 * r + 0.357580 * g + 0.189423 * b;
		Y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
		Z = 0.019334 * r + 0.119193 * g + 0.950227 * b;

		a = 5*(flab(X)-flab(Y));
		b = 2*(flab(Y)-flab(Z));

		out[2*i+0] = convert2<Value>(a);
		out[2*i+1] = convert2<Value>(b);
	}
}

template<class Value>
void	rgb2int (Value* out,const Value* in,int n)
{
	for (int i=0;i<n;i++)
	{
		double r,g,b;
		r = read<Value>(in[3*i+0]);
		g = read<Value>(in[3*i+1]);
		b = read<Value>(in[3*i+2]);
		out[i] = convert<Value>( (r+g+b) / 3.0 );
	}
}

}

#endif
