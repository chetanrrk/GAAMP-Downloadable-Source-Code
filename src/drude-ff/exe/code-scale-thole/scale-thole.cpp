#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_N	(8192)
#define MAX_LEN	(256)

#ifndef min(a,b)
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

char szRec[MAX_N][MAX_LEN];
double thole[MAX_N];

int nAtom=0, nLine=0, Line_List[MAX_N];

void Read_xpsf(char szName[]);
void Scale_thole(double scale);
void Update_thole(int Idx, double thole);
void Save_xpsf(char szName[]);

int main(int argc, char *argv[])
{
	if(argc != 3)	{
		printf("Usage: scale-thole drude-mol.xpsf scale-factor\n");
		exit(1);
	}

	Read_xpsf(argv[1]);

	Scale_thole(atof(argv[2]));

	Save_xpsf(argv[1]);

	return 0;
}

void Save_xpsf(char szName[])
{
	FILE *fOut;
	int i;

	fOut = fopen(szName, "w");

	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szRec[i]);
	}
	
	fclose(fOut);
}

void Update_thole(int Idx, double thole)
{
	char szBuff[256];
	
	sprintf(szBuff, "%10.5lf", thole);
	memcpy(szRec[Idx]+105, szBuff, 10);
}

void Scale_thole(double scale)
{
	int i, Idx;

	for(i=0; i<nAtom; i++)	{
		if(thole[i] > 0.001)	{	// a heavy atom
			Idx = Line_List[i];	// the line of the heavy atom
			Update_thole(Idx, scale*thole[i]);
		}
	}
}

void Read_xpsf(char szName[])
{
	int ReadItem, iTmp, nDrude=0;
	FILE *fIn;
	char szLine[256], *ReadLine, szResName[256], szAtom[256], szChem[256], szTmp[256];
	double chg=0.0, polarizability, mass, thole_Local;

	nAtom = 0;
	nLine = 0;

	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit.\n", szName);
		exit(1);
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		ReadItem = sscanf(szLine, "%d%s%d%s%s%s%lf%lf%d%lf%lf", 
			&iTmp, &szTmp, &iTmp, szResName, szAtom, szChem, &chg, &mass, &iTmp, &polarizability, &thole_Local);

		if(ReadItem == 11)	{
			thole[nAtom] = thole_Local;
			Line_List[nAtom] = nLine;
			
			nAtom++;
		}

		strcpy(szRec[nLine], szLine);
		nLine++;
	}

	fclose(fIn);

}
