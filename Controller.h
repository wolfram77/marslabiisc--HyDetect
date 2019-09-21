#ifndef CONTROLLER_H_
#define CONTROLLER_H_


typedef enum{
	MST = 0,
	COLORING,
	COMM_DETECT,
	TC
}application;

typedef enum{
	HOST = 0,
	DEVICE,
	HYBRID
}arch;

//to store the input attributes
typedef struct option{
	application apps;
	char *fileName;
	arch mode;
}options;


//to read and parse input parameters	
options parse_cmdline(int, char*);

#endif
