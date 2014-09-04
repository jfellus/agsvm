/*
 * utils.cpp
 *
 *  Created on: 22 nov. 2013
 *      Author: jfellus
 */



#include "utils.h"
#include <string.h>

string stringprintf(const char* fmt, ...) {
	char* s = new char[strlen(fmt)+512];
	va_list vl;
	va_start(vl, fmt);
	vsprintf(s, fmt, vl);
	va_end(vl);
	string ss(s);
	delete s;
	return ss;
}

void fappend(const string& filename, const string& line) {
	FILE* f = fopen(filename.c_str(), "a");
		fputs(line.c_str(), f);
	fclose(f);
}
