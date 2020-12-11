#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

void Quit_With_Error_Msg(char szMsg[]);

CMol Mol;

FILE *fFile_Run_Log;

int main(int argc, char *argv[])
{
	fFile_Run_Log = fopen("log-pdb-to-crd.txt", "w");

	Mol.ReadPSF("mol.xpsf", 0);	//the first parameter is the xpsf file name; the second parameter is 0 (only one molecule in xpsf) or 1 (two molecules in xpsf). 
	Mol.ReadPDB("mol-opt.pdb");
	Mol.WriteCRDFile("mol.crd");
	Mol.WriteInpFile("mol.inp");

	fclose(fFile_Run_Log);

	return 0;
}


void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("../error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in update-torsion-para.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}
