#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_LINE_PRM	(2048)

char szPrmTxt[MAX_LINE_PRM][256];

void To_Upper_Case(char szBuff[]);
void Remove_LJ14_Para(char szName[]);
void Quit_With_Error_Msg(char szMsg[]);
void Remove_Redundant_Entries(char szName[]);

int main(int argc, char *argv[])
{
	char ErrorMsg[256];

	if(argc != 2)	{
		sprintf(ErrorMsg, "Usage: cut-lj14 mol.prm\n");
		Quit_With_Error_Msg(ErrorMsg);
	}

//	Remove_LJ14_Para(argv[1]);
	Remove_Redundant_Entries(argv[1]);


	return 0;
}

void Remove_LJ14_Para(char szName[])
{
	FILE *fIn, *fOut;
	int nLine=0, i, ReadItem;
	char *ReadLine, Line_Nonbonded, ChemName[64];
	double fTmp=0.0, Emin, Rmin;

	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open mol.prm for read in Rewrite_Prm_File().\nQuit\n");
	}
	nLine=0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szPrmTxt[nLine], 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			nLine++;
			if(nLine >= MAX_LINE_PRM)	{
				fclose(fIn);
				Quit_With_Error_Msg("nLine >= MAX_LINE_PRM in Remove_LJ14_Para().\nQuit\n");
			}
		}
	}
	fclose(fIn);

	Line_Nonbonded=-1;
	for(i=0; i<nLine; i++)	{
		if(strncmp(szPrmTxt[i], "NONBONDED", 9)==0)	{
			Line_Nonbonded = i;
			break;
		}
	}
	if(Line_Nonbonded < 0)	{
		Quit_With_Error_Msg("Fail to find the entry for NONBONDED in mol.prm\n");
	}

	fOut = fopen(szName, "w");
	for(i=0; i<nLine; i++)	{
		if(i > Line_Nonbonded)	{
			ReadItem = sscanf(szPrmTxt[i], "%s%lf%lf%lf%lf%lf%lf", ChemName, &fTmp, &Emin, &Rmin, &fTmp, &fTmp, &fTmp);
			if(ReadItem == 7)	{	// containing LJ 14 parameters
				fprintf(fOut, "%-5s   0.00%10.4lf%10.4lf\n", ChemName, Emin, Rmin);
			}
			else	{
				fprintf(fOut, "%s", szPrmTxt[i]);
			}
		}
		else	{
			fprintf(fOut, "%s", szPrmTxt[i]);
		}
	}
	fclose(fOut);

}

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("../error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in cut-LJ14.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}

#define MAX_LINE	(2048)
#define MAX_LEN		(256)

char szTxtPrm[MAX_LINE][MAX_LEN];

void Remove_Redundant_Entries(char szName[])
{
	int n_Line, Line_Bond_Begin=-1, Line_Bond_End=-1, Line_Angle_Begin=-1, Line_Angle_End=-1, Line_Dihedral_Begin=-1, Line_Dihedral_End=-1, Line_Improper_Begin=-1, Line_Improper_End=-1, Line_LJ_Begin=-1, Line_LJ_End=-1;
	int Found_Improper=0, ToOutput[MAX_LINE], i, j, ReadItem;
	FILE *fOut;
	char ErrorMsg[256], *ReadLine;
	char szAtom_I_1[256], szAtom_I_2[256], szAtom_I_3[256], szAtom_I_4[256], szAtom_J_1[256], szAtom_J_2[256], szAtom_J_3[256], szAtom_J_4[256];
	double Para_1, Para_2, Para_3;
	double Para_1_j, Para_2_j, Para_3_j;

	fOut = fopen(szName, "r");
	if(fOut == NULL)	{
		sprintf(ErrorMsg, "Error in open file: %s\nQuit\n", szName);
		Quit_With_Error_Msg(ErrorMsg);
	}
	else	{
		n_Line = 0;

		while(1)	{
			if(feof(fOut))	{
				break;
			}
			ReadLine = fgets(szTxtPrm[n_Line], MAX_LEN, fOut);
			To_Upper_Case(szTxtPrm[n_Line]);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				if(strncmp(szTxtPrm[n_Line], "BOND", 4)==0)	{
//					sprintf(szTxtPrm[n_Line], "BONDS\n");			// rename "BOND" to "BONDS"
					Line_Bond_Begin = n_Line+1;
				}
				else if(strncmp(szTxtPrm[n_Line], "ANGLE", 5)==0)	{
//					sprintf(szTxtPrm[n_Line], "ANGLES\n");			// rename "ANGLE" to "ANGLES"
					Line_Bond_End = n_Line-1;
					Line_Angle_Begin = n_Line+1;
				}
				else if(strncmp(szTxtPrm[n_Line], "DIHEDRAL", 8)==0)	{
//					sprintf(szTxtPrm[n_Line], "DIHEDRALS\n");			// rename "DIHEDRAL" to "DIHEDRALS"
					Line_Angle_End = n_Line-1;
					Line_Dihedral_Begin = n_Line+1;
				}
				else if(strncmp(szTxtPrm[n_Line], "IMPHI", 5)==0)	{
					Found_Improper = 1;
					Line_Dihedral_End = n_Line-1;
					Line_Improper_Begin = n_Line+1;
				}
				else if(strncmp(szTxtPrm[n_Line], "NONBONDED", 9)==0)	{
					Line_Improper_End = n_Line-1;
					Line_LJ_Begin = n_Line+1;

					if(Found_Improper == 0)	{	// not exist
						Line_Dihedral_End = n_Line-1;
						Line_Improper_Begin = n_Line;
					}
				}


				ToOutput[n_Line] = 1;
				n_Line++;
				if(n_Line >= MAX_LINE)	{
					sprintf(ErrorMsg, "n_Line >= MAX_LINE\nQuit\n");
					fclose(fOut);
					Quit_With_Error_Msg(ErrorMsg);
				}
			}
		}
		fclose(fOut);
		Line_LJ_End = n_Line-1;
	}


	//start	to check bond parameters
	for(i=Line_Bond_Begin; i<=Line_Bond_End; i++)	{
		ReadItem = sscanf(szTxtPrm[i], "%s%s%lf%lf", szAtom_I_1, szAtom_I_2, &Para_1, &Para_2);
		if(ReadItem == 4)	{
			for(j=i+1; j<=Line_Bond_End; j++)	{
				ReadItem = sscanf(szTxtPrm[j], "%s%s%lf%lf", szAtom_J_1, szAtom_J_2, &Para_1_j, &Para_2_j);
				if(ReadItem == 4)	{
					if( ((strcmp(szAtom_I_1, szAtom_J_1)==0) && (strcmp(szAtom_I_2, szAtom_J_2)==0)) || ((strcmp(szAtom_I_1, szAtom_J_2)==0) && (strcmp(szAtom_I_2, szAtom_J_1)==0)) )	{	// same entry
						ToOutput[j] = 0;
						printf("To delete: %s", szTxtPrm[j]);
					}
				}
			}
		}
	}
	//end	to check bond parameters

	//start	to check angle parameters
	for(i=Line_Angle_Begin; i<=Line_Angle_End; i++)	{
		ReadItem = sscanf(szTxtPrm[i], "%s%s%s%lf%lf", szAtom_I_1, szAtom_I_2, szAtom_I_3, &Para_1, &Para_2);
		if(ReadItem == 5)	{
			for(j=i+1; j<=Line_Angle_End; j++)	{
				ReadItem = sscanf(szTxtPrm[j], "%s%s%s%lf%lf", szAtom_J_1, szAtom_J_2, szAtom_J_3, &Para_1_j, &Para_2_j);
				if(ReadItem == 5)	{
					if( ((strcmp(szAtom_I_1, szAtom_J_1)==0) && (strcmp(szAtom_I_2, szAtom_J_2)==0) && (strcmp(szAtom_I_3, szAtom_J_3)==0)) 
						|| ((strcmp(szAtom_I_1, szAtom_J_3)==0) && (strcmp(szAtom_I_2, szAtom_J_2)==0) && (strcmp(szAtom_I_3, szAtom_J_1)==0)) )	{	// same entry
						ToOutput[j] = 0;
						printf("To delete: %s", szTxtPrm[j]);
					}
				}
			}
		}
	}
	//end	to check angle parameters

	//start	to check dihedral parameters
	for(i=Line_Dihedral_Begin; i<=Line_Dihedral_End; i++)	{
		ReadItem = sscanf(szTxtPrm[i], "%s%s%s%s%lf%lf%lf", szAtom_I_1, szAtom_I_2, szAtom_I_3, szAtom_I_4, &Para_1, &Para_2, &Para_3);
		if(ReadItem == 7)	{
			for(j=i+1; j<=Line_Dihedral_End; j++)	{
				ReadItem = sscanf(szTxtPrm[j], "%s%s%s%s%lf%lf%lf", szAtom_J_1, szAtom_J_2, szAtom_J_3, szAtom_J_4, &Para_1_j, &Para_2_j, &Para_3_j);
				if(ReadItem == 7)	{
					if( ( ((strcmp(szAtom_I_1, szAtom_J_1)==0) && (strcmp(szAtom_I_2, szAtom_J_2)==0) && (strcmp(szAtom_I_3, szAtom_J_3)==0) && (strcmp(szAtom_I_4, szAtom_J_4)==0)) 
						|| ((strcmp(szAtom_I_1, szAtom_J_4)==0) && (strcmp(szAtom_I_2, szAtom_J_3)==0) && (strcmp(szAtom_I_3, szAtom_J_2)==0) && (strcmp(szAtom_I_4, szAtom_J_1)==0)) )
						&& (fabs(Para_2_j-Para_2) < 1.0E-5) )	{	// same entry
						ToOutput[j] = 0;
						printf("To delete: %s", szTxtPrm[j]);
					}
				}
			}
		}
	}
	//end	to check dihedral parameters

	//start	to check improper parameters
	for(i=Line_Improper_Begin; i<=Line_Improper_End; i++)	{
		ReadItem = sscanf(szTxtPrm[i], "%s%s%s%s%lf%lf%lf", szAtom_I_1, szAtom_I_2, szAtom_I_3, szAtom_I_4, &Para_1, &Para_2, &Para_3);
		if(ReadItem == 7)	{
			for(j=i+1; j<=Line_Improper_End; j++)	{
				ReadItem = sscanf(szTxtPrm[j], "%s%s%s%s%lf%lf%lf", szAtom_J_1, szAtom_J_2, szAtom_J_3, szAtom_J_4, &Para_1_j, &Para_2_j, &Para_3_j);
				if(ReadItem == 7)	{
					if( ( ((strcmp(szAtom_I_1, szAtom_J_1)==0) && (strcmp(szAtom_I_2, szAtom_J_2)==0) && (strcmp(szAtom_I_3, szAtom_J_3)==0) && (strcmp(szAtom_I_4, szAtom_J_4)==0)) 
						|| ((strcmp(szAtom_I_1, szAtom_J_4)==0) && (strcmp(szAtom_I_2, szAtom_J_3)==0) && (strcmp(szAtom_I_3, szAtom_J_2)==0) && (strcmp(szAtom_I_4, szAtom_J_1)==0)) )
						&& (fabs(Para_2_j-Para_2) < 1.0E-5) )	{	// same entry
						ToOutput[j] = 0;
						printf("To delete: %s", szTxtPrm[j]);
					}
				}
			}
		}
	}
	//end	to check improper parameters

	//start	to check LJ parameters
	for(i=Line_LJ_Begin; i<=Line_LJ_End; i++)	{
		ReadItem = sscanf(szTxtPrm[i], "%s%lf%lf%lf", szAtom_I_1, &Para_1, &Para_2, &Para_3);
		if(ReadItem == 4)	{
			for(j=i+1; j<=Line_LJ_End; j++)	{
				ReadItem = sscanf(szTxtPrm[j], "%s%lf%lf%lf", szAtom_J_1, &Para_1, &Para_2, &Para_3);
				if(ReadItem == 4)	{
					if( strcmp(szAtom_I_1, szAtom_J_1)==0 )	{	// same entry
						ToOutput[j] = 0;
						printf("To delete: %s", szTxtPrm[j]);
					}
				}
			}
		}
	}
	//end	to check LJ parameters

	fOut = fopen(szName, "w");
	for(i=0; i<n_Line; i++)	{
		if(ToOutput[i])	{
			fprintf(fOut, "%s", szTxtPrm[i]);
		}
	}
	fclose(fOut);
}


void To_Upper_Case(char szBuff[])
{
	char gap;
	int i=0;

	gap = 'A'-'a';

	while(1)	{
		if(szBuff[i]==0)	{
			break;
		}
		if( (szBuff[i]>='a') && (szBuff[i]<='z') )	{
			szBuff[i] += gap;
		}
		i++;
	}
}
