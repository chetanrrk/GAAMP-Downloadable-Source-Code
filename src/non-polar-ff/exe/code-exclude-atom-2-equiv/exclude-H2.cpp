#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//start	data for one line
int nItem_Line;
char szItemLine[64][128];
//end	data for one line

int Split_Str_One_Line(char szBuff[]);
int Is_Atom_2_In_List(int nItem);

int Idx_H_excluded;

int main(int argc, char *argv[])
{
	FILE *fIn, *fOut;
	char *ReadLine, szLine[1024], szTmp[256];
	int nItem, Idx, i;

	if(argc != 2)	{
		printf("Usage: exlude_H index\n");
		exit(1);
	}
	Idx_H_excluded = atoi(argv[1]);


	fIn = fopen("equiv-org.txt", "r");
	fOut = fopen("equiv-mod.txt", "w");

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 1024, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		nItem = Split_Str_One_Line(szLine);
		Idx = Is_Atom_2_In_List(nItem);

		if(Idx < 0)	{
			fprintf(fOut, "%s", szLine);
		}
		else	{	// to modify and output
			if(nItem > 3)	{	// more than two equivalent atoms
				sprintf(szLine, "%s ", szItemLine[0]);
				for(i=1; i<nItem; i++)	{
					if( i!=Idx )	{
						sprintf(szTmp, "%4s", szItemLine[i]);
						strcat(szLine, szTmp);
					}
				}
				fprintf(fOut, "%s\n", szLine);
			}
		}

		
	}

	fclose(fIn);
	fclose(fOut);

	return 0;
}


int Split_Str_One_Line(char szBuff[])
{
	int iLen, iPos, iPos_End, CountChange=0, iLen_Atom_Name;
	char c, szBak[256];

	nItem_Line = 0;

	iLen = strlen(szBuff);
	strcpy(szBak, szBuff);

	iPos = 0;
	while(1)	{
		while(1)	{	// to find the beginning of a string
			c = szBuff[iPos];
			if( (c != 0xA) && (c != 0xD) && (c != 0x20)  && (c != '\t') && (c != 0) )	{
				break;
			}
			iPos++;
			if(iPos >= iLen)	{
				break;
			}
		}
		if(iPos >= iLen)	{
			break;
		}

		iPos_End = iPos;
		while(1)	{	// to find the end of a string
			c = szBuff[iPos_End];
			if( (c == 0xA) || (c == 0xD) || (c == 0x20)  || (c == '\t') || (c == 0) )	{
				break;
			}
			iPos_End++;
			if(iPos >= iLen)	{
				break;
			}
		}

		iLen_Atom_Name = iPos_End-iPos;
		memcpy(szItemLine[nItem_Line], szBuff+iPos, iLen_Atom_Name);
		szItemLine[nItem_Line][iLen_Atom_Name] = 0;	// end
		nItem_Line++;


		iPos = iPos_End;
	}

	return nItem_Line;
}

int Is_Atom_2_In_List(int nItem)
{
	int i;
	char szIdx[16];

	sprintf(szIdx, "%d", Idx_H_excluded);

	for(i=1; i<nItem; i++)	{
		if(strcmp(szItemLine[i],szIdx)==0)	{
			return i;
		}
	}
	return (-1);
}
