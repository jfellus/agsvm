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

void foverwrite(const string& filename, const string& line) {
	FILE* f = fopen(filename.c_str(), "w");
		fputs(line.c_str(), f);
	fclose(f);
}

string get_config_str(const char* what, const char* default_val) {
	FILE* f = fopen("config.properties", "r");
	char line[512];
	char* c = 0;
	string v = default_val;
	while ( fgets (line , 512 , f) != NULL ) {
	      if(!strncmp(what, line, strlen(what))) {
	    	  c = line + strlen(what);
	    	  if(*c!=' ' && *c!='=') continue;
	    	  while(*c==' ' || *c=='=') c++;
	    	  v = c;
	    	  break;
	      }
	}
	fclose(f);
	v.erase(v.find_last_not_of(" \n\r\t")+1);
	return v;
}


void DBGVECTOR(float* x, int D) {
	printf("[");
	for(int i=0; i<D; i++) printf(" %f", x[i]);
	printf(" ]\n");
}

inline void shell(const string& s) {if(system(s.c_str())) throw "error";}

void setenv(const char* name, int val) {
	shell(fmt("export %s=%u", name, val));
}
void setenv(const char* name, double val) {
	shell(fmt("export %s=%f", name, val));
}
