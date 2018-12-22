#ifndef VERSION_H
#define VERSION_H

	//Date Version Types
	static const char DATE[] = "02";
	static const char MONTH[] = "10";
	static const char YEAR[] = "2018";
	static const char UBUNTU_VERSION_STYLE[] =  "18.10";
	
	//Software Status
	static const char STATUS[] =  "Alpha";
	static const char STATUS_SHORT[] =  "a";
	
	//Standard Version Type
	static const long MAJOR  = 1;
	static const long MINOR  = 0;
	static const long BUILD  = 4;
	static const long REVISION  = 16;
	
	//Miscellaneous Version Types
	static const long BUILDS_COUNT  = 5;
	#define RC_FILEVERSION 1,0,4,16
	#define RC_FILEVERSION_STRING "1, 0, 4, 16\0"
	static const char FULLVERSION_STRING [] = "1.0.4.16";
	
	//These values are to keep track of your versioning state, don't modify them.
	static const long BUILD_HISTORY  = 4;
	

#endif //VERSION_H
