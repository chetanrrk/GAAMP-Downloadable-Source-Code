#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

#define MAX_LEN	(256)
#define MAX_ATOM	(512)
#define MAX_LINE_RTF	(1024)
#define MAX_BOND	(2048)
#define MAX_LINE_INPUT	(2048)

#define MAX_ATOM_TYPE	(2048)
#define LEN_TXT	(131072)

int nAtom, nLineRec, nBond, nImproper, nDihedral, nAtomType, UsedAtomType[MAX_ATOM_TYPE];
char szAtomName[MAX_ATOM][16];
char szChemName[MAX_ATOM][16];
char szRtf[MAX_LINE_RTF][MAX_LEN];
char szAtom_Type_Txt[MAX_ATOM_TYPE][MAX_LEN], szAtom_Type_List[MAX_ATOM_TYPE][16];
char DistMatrix[MAX_ATOM][MAX_ATOM];
char szMassTxt[LEN_TXT];
char szAddBond[LEN_TXT], szAddAngle[LEN_TXT], szAddDihedral[LEN_TXT], szAddImprop[LEN_TXT];

int BondList[MAX_BOND][2];
int ImproperList[MAX_BOND][4];
int DihedralList[MAX_BOND][4];
int AtomBond[MAX_ATOM], Bond_Array[MAX_ATOM][6];

//start	data for one line
int nItem_Line;
char szItemLine[64][128];
//end	data for one line


int FindTags_As_First_String(char szBuff[], char szTag[]);
int Extract_Atom_Info();
int Extract_Bond_Info();
int Extract_ImproperDihedral_Info();
int Query_Atom_Name(char szAtom[]);
int Split_Str_One_Line(char szBuff[]);
void Setup_Dist_matrix(void);
void Quit_With_Error_Msg(char szMsg[]);
void Setup_FF_Parameters(CForceField *ff);
void Save_Mini_FF(void);
void Read_Input(char szName_Input[]);
void Extract_Atom_Type_Info(char szName_Full_Rtf[]);
int Is_A_Used_Chem_Name(char szChemName_Check[]);
void Output_Rtf(void);
void Find_String_in_File(FILE *fRead, char szTag[]);

CForceField ForceField;
CForceField ForceFieldCheck;
FILE *fFile_Run_Log;	// will be shared by other source code
int net_charge;

int main(int argc, char *argv[])
{
	char ErrorMsg[256];

	if(argc != 5)	{
		printf("Usage: mini_ff input your.rtf your.prm net_charge\nQuit\n");
		exit(1);
	}

	fFile_Run_Log = fopen("log-extract-ff.txt", "w");
	if(fFile_Run_Log==NULL)	{
		sprintf(ErrorMsg, "Fail to create the log file.\n\n");
		Quit_With_Error_Msg(ErrorMsg);
	}

	net_charge = atoi(argv[4]);

	Read_Input(argv[1]);
	Extract_Atom_Info();
	Extract_Bond_Info();
	Extract_ImproperDihedral_Info();

	Extract_Atom_Type_Info(argv[2]);

	Output_Rtf();

	ForceField.ReadForceField(argv[3]);


	Setup_Dist_matrix();

	Setup_FF_Parameters(&ForceField);

	Save_Mini_FF();

	ForceFieldCheck.ReadForceField("mol.prm");
	ForceFieldCheck.Quit_on_Error = 1;
	Setup_FF_Parameters(&ForceFieldCheck);	// to check whether prm file is complete or not


	fclose(fFile_Run_Log);

	return 0;
}

void Output_Rtf(void)
{
	FILE *fOut;
	int i;

	fOut = fopen("mol.rtf", "w");
	fprintf(fOut, "%s\nAUTO ANGLES DIHE\n\n", szMassTxt);

	for(i=0; i<nLineRec; i++)	{
		fprintf(fOut, "%s", szRtf[i]);
	}

	fclose(fOut);
}

void Extract_Atom_Type_Info(char szName_Full_Rtf[])
{
	FILE *fIn;
	char szError[256], szLine[256], *ReadLine, szBuff[256], szTmp[256];
	int ReadItem, ChemID;
	double mass=0.0;

	nAtomType = 0;

	fIn = fopen(szName_Full_Rtf, "r");

	if(fIn == NULL)	{
		sprintf(szError, "Fail to open file %s\nQuit\n", szName_Full_Rtf);
		Quit_With_Error_Msg(szError);
	}

	szMassTxt[0] = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine != NULL)	{
			if(strncmp(szLine, "MASS", 4)==0)	{
				ReadItem = sscanf(szLine, "%s%d%s%lf", szTmp, &ChemID, szBuff, &mass);
				if(ReadItem == 4)	{
					strcpy(szAtom_Type_Txt[nAtomType], szLine);
					strcpy(szAtom_Type_List[nAtomType], szBuff);
					nAtomType++;

					if(Is_A_Used_Chem_Name(szBuff))	{
						strcat(szMassTxt, szLine);
					}
				}
			}
		}
	}

	fclose(fIn);
}

int Is_A_Used_Chem_Name(char szChemName_Check[])
{
	int i;

	for(i=0; i<nAtom; i++)	{
		if(strcmp(szChemName_Check, szChemName[i])==0)	{	// used
			return 1;
		}
	}

	return 0;
}

void Find_String_in_File(FILE *fRead, char szTag[])
{
	char szLine[256], *ReadLine, szError[256];
	int nLen;

	nLen = strlen(szTag);
	while(1)	{
		if(feof(fRead))	{
			sprintf(szError, "Fail to find tag: %s\nQuit\n", szTag);
			fclose(fRead);
			Quit_With_Error_Msg(szError);
		}
		ReadLine = fgets(szLine, 256, fRead);
		if(ReadLine == NULL)	{
			sprintf(szError, "Fail to find tag: %s\nQuit\n", szTag);
			fclose(fRead);
			Quit_With_Error_Msg(szError);
		}
		else	{
			if(strncmp(szLine, szTag, nLen)==0)	{
				return;
			}
		}
	}
}

void Read_Input(char szName_Input[])
{
	FILE *fIn;
	char szError[256], szLine[256], *ReadLine, szResName[24]="RESI MOL              ", szTmp[256];
	int ToSave, ReadItem;
	double net_chg_Local=0.0;

	fIn = fopen(szName_Input, "r");
	if(fIn == NULL)	{
		sprintf(szError, "Fail to open file %s\nQuit\n", szName_Input);
		Quit_With_Error_Msg(szError);
	}

	// to find "RESI "
	while(1)	{
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			fclose(fIn);
			Quit_With_Error_Msg("Fail to find the residue information.\nQuit\n");
		}
		else	{
			if(FindTags_As_First_String(szLine, "RESI")==1)	{
				break;
			}
		}
	}

	nLineRec = 1;
	memcpy(szLine, szResName, 16);
	strcpy(szRtf[0], szLine);

	//start	to validate the net charge
	ReadItem = sscanf(szLine, "%s%s%lf", szTmp, szTmp, &net_chg_Local);
	if(ReadItem != 3)	{
		fclose(fIn);
		Quit_With_Error_Msg("Fail to extract the net charge info.\nQuit\n");
	}
	if( fabs(net_chg_Local-1.0*net_charge) > 1.0E-4 )	{
		fclose(fIn);
		Quit_With_Error_Msg("The net charge from CGenFF is NOT consistent with the charge provided by user.\nQuit\n");
	}
	//end	to validate the net charge

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		ToSave = 0;

		if(FindTags_As_First_String(szLine, "ATOM")==1)	{	// to 
			ToSave = 1;
		}
		else if(FindTags_As_First_String(szLine, "BOND")==1)	{
			ToSave = 1;
		}
		else if(FindTags_As_First_String(szLine, "DOUBLE")==1)	{
			ToSave = 1;
		}
		else if(FindTags_As_First_String(szLine, "IMPR")==1)	{
			ToSave = 1;
		}

		// auto angle and dihedral

		if(ToSave)	{
			strcpy(szRtf[nLineRec], szLine);
			nLineRec++;
		}

		if(FindTags_As_First_String(szLine, "END")==1)	{
			strcpy(szRtf[nLineRec], "\nEND\n");
			nLineRec++;
			break;
		}

	}

	Find_String_in_File(fIn, "BONDS");

	szAddBond[0] = szAddAngle[0] = szAddDihedral[0] = szAddImprop[0] = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		if(strncmp(szLine, "ANGLE", 5)==0)	{
			break;
		}
		else	{
			strcat(szAddBond, szLine);
		}
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		if(strncmp(szLine, "DIHEDRALS", 9)==0)	{
			break;
		}
		else	{
			strcat(szAddAngle, szLine);
		}
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		if(strncmp(szLine, "IMPROPER", 8)==0)	{
			break;
		}
		else	{
			strcat(szAddDihedral, szLine);
		}
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		if(strncmp(szLine, "END", 3)==0)	{
			break;
		}
		else	{
			strcat(szAddImprop, szLine);
		}
	}


	fclose(fIn);
}

void Save_Mini_FF(void)
{
	FILE *fOut;
	int i;
	double fTmp=0.0;

	fOut = fopen("mol.prm", "w");

	fprintf(fOut, "* Mini Force Field Parameter File.\n*\n\nBONDS\n");

	for(i=0; i<ForceField.n_Rec_Bond; i++)	{
		if(ForceField.Active_Bond[i])	{
			fprintf(fOut, "%-8s%-8s %8.3lf   %8.4lf\n", 
				ForceField.Bond_Rec[i].Chem[0], ForceField.Bond_Rec[i].Chem[1], ForceField.Bond_Rec[i].para[0], ForceField.Bond_Rec[i].para[1]);
		}
	}
	fprintf(fOut, "\n%s\n", szAddBond);

	fprintf(fOut, "\nANGLES\n");

	for(i=0; i<ForceField.n_Rec_Angle; i++)	{
		if(ForceField.Active_Angle[i])	{
			fprintf(fOut, "%-8s%-8s%-8s %8.3lf   %8.4lf", 
				ForceField.Angle_Rec[i].Chem[0], ForceField.Angle_Rec[i].Chem[1], ForceField.Angle_Rec[i].Chem[2], ForceField.Angle_Rec[i].para[0], ForceField.Angle_Rec[i].para[1]);
			if( (ForceField.Angle_Rec[i].para[2] > 1.0E-3) && (ForceField.Angle_Rec[i].para[3] > 1.0E-3) )	{	// UREY-BRADLEY
				fprintf(fOut, "  %8.3lf  %9.5lf", ForceField.Angle_Rec[i].para[2], ForceField.Angle_Rec[i].para[3]);
			}
			fprintf(fOut, "\n");
		}
	}
	fprintf(fOut, "\n%s\n", szAddAngle);

	fprintf(fOut, "\nDIHEDRALS\n");

	for(i=0; i<ForceField.n_Rec_Dihedral; i++)	{
		if(ForceField.Active_Dihedral[i])	{
			fprintf(fOut, "%-8s%-8s%-8s%-8s %8.4lf  %8.0lf   %8.1lf\n", 
				ForceField.Dihedral_Rec[i].Chem[0], ForceField.Dihedral_Rec[i].Chem[1], ForceField.Dihedral_Rec[i].Chem[2], ForceField.Dihedral_Rec[i].Chem[3], 
				ForceField.Dihedral_Rec[i].para[0], ForceField.Dihedral_Rec[i].para[1], ForceField.Dihedral_Rec[i].para[2]);
		}
	}
	fprintf(fOut, "\n%s\n", szAddDihedral);


//	fprintf(fOut, "\nIMPHI\n");
	fprintf(fOut, "\nIMPROPERS\n");

	for(i=0; i<ForceField.n_Rec_ImproDihedral; i++)	{
		if(ForceField.Active_ImproDihedral[i])	{
			fprintf(fOut, "%-8s%-8s%-8s%-8s %8.4lf  %8.0lf   %8.1lf\n", 
				ForceField.ImproDihedral_Rec[i].Chem[0], ForceField.ImproDihedral_Rec[i].Chem[1], ForceField.ImproDihedral_Rec[i].Chem[2], ForceField.ImproDihedral_Rec[i].Chem[3], 
				ForceField.ImproDihedral_Rec[i].para[0], ForceField.ImproDihedral_Rec[i].para[1], ForceField.ImproDihedral_Rec[i].para[2]);
		}
	}
	fprintf(fOut, "\n%s\n", szAddImprop);

	fprintf(fOut, "\nNONBONDED\n");

	for(i=0; i<ForceField.n_Rec_LJ; i++)	{
		if(ForceField.Active_LJ[i])	{
			if( (fabs(ForceField.LJ_Rec[i].para[4]-ForceField.LJ_Rec[i].para[1]) > 1.0E-6) || (fabs(ForceField.LJ_Rec[i].para[5]-ForceField.LJ_Rec[i].para[2]) > 1.0E-6) )	{	// LJ 1-4
				fprintf(fOut, "%-8s 0.00   %10.6lf  %10.6lf   0.00   %10.6lf  %10.6lf\n", 
					ForceField.LJ_Rec[i].Chem, ForceField.LJ_Rec[i].para[1], ForceField.LJ_Rec[i].para[2], ForceField.LJ_Rec[i].para[4], ForceField.LJ_Rec[i].para[5]);
			}
			else	{
				fprintf(fOut, "%-8s 0.00   %10.6lf  %10.6lf\n", 
					ForceField.LJ_Rec[i].Chem, ForceField.LJ_Rec[i].para[1], ForceField.LJ_Rec[i].para[2]);
			}
		}
	}

	fprintf(fOut, "\n\nEND\n");

	fclose(fOut);

}

#define MAX_DIH_ITEM	(7)	// 1, ..., 6
void Setup_FF_Parameters(CForceField *ff)
{
	int i, j, k, l;
	char szChemName_Local[8][N_LEN_CHEM_NAME];
	double Para_List[MAX_DIH_ITEM*3];

	memset(ff->Active_Bond, 0, sizeof(int)*ff->n_Rec_Bond);
	memset(ff->Active_Angle, 0, sizeof(int)*ff->n_Rec_Angle);
	memset(ff->Active_Dihedral, 0, sizeof(int)*ff->n_Rec_Dihedral);
	memset(ff->Active_ImproDihedral, 0, sizeof(int)*ff->n_Rec_ImproDihedral);
	memset(ff->Active_LJ, 0, sizeof(int)*ff->n_Rec_LJ);
	memset(ff->Active_NBFix, 0, sizeof(int)*ff->n_Rec_NBFix);

	for(i=0; i<nAtom; i++)	{	// bond
		for(j=i+1; j<nAtom; j++)	{
			if(DistMatrix[i][j] == 1)	{	// bond
				strcpy(szChemName_Local[0], szChemName[i]);
				strcpy(szChemName_Local[1], szChemName[j]);
				ff->GetPara_Bond(szChemName_Local, Para_List);
			}
		}
	}

	for(i=0; i<nAtom; i++)	{	// angle
		for(j=0; j<nAtom; j++)	{
			if(DistMatrix[i][j] == 1)	{	// a bond
				for(k=0; k<nAtom; k++)	{
					if( (k!=i) && (DistMatrix[j][k] == 1) )	{	// the third atom, now an angle is formed
						strcpy(szChemName_Local[0], szChemName[i]);
						strcpy(szChemName_Local[1], szChemName[j]);
						strcpy(szChemName_Local[2], szChemName[k]);
						ff->GetPara_Angle(szChemName_Local, Para_List);

						for(l=0; l<nAtom; l++)	{
							if( (l!=j) && (l!=i) && (DistMatrix[k][l] == 1) )	{	// the fourth atom, form a dihedral
								strcpy(szChemName_Local[3], szChemName[l]);
								ff->GetPara_Dihedral(szChemName_Local, Para_List);
							}
						}
					}
				}
			}
		}
	}

	for(i=0; i<nAtom; i++)	{	// LJ parameter
		strcpy(szChemName_Local[0], szChemName[i]);
		ff->GetPara_LJ(szChemName_Local, Para_List);
	}

	if(ff->n_Rec_NBFix > 0)	{
		for(i=0; i<nAtom; i++)	{
			for(j=i+1; j<nAtom; i++)	{
				if(DistMatrix[i][j] >= 3)	{	// nonbonded pair
					ff->QueryNBFix(szChemName[i], szChemName[j]);
				}
			}
		}
	}

	for(i=0; i<nImproper; i++)	{
		strcpy(szChemName_Local[0], szChemName[ImproperList[i][0]]);
		strcpy(szChemName_Local[1], szChemName[ImproperList[i][1]]);
		strcpy(szChemName_Local[2], szChemName[ImproperList[i][2]]);
		strcpy(szChemName_Local[3], szChemName[ImproperList[i][3]]);

		ff->GetPara_ImproDIhedral(szChemName_Local, Para_List);
	}
}

int Query_Atom_Name(char szAtom[])
{
	int i;

	for(i=0; i<nAtom; i++)	{
		if(strcmp(szAtom, szAtomName[i])==0)	{
			return i;
		}
	}

	printf("Fail to find atom: %s\nQuit\n", szAtom);
	exit(1);

	return -1;
}

int Extract_ImproperDihedral_Info()
{
	int i, j, nItem, iAtom_1, iAtom_2, iAtom_3, iAtom_4;

	nImproper = 0;

	memset(ImproperList, 0, sizeof(int)*MAX_BOND*4);

	for(i=0; i<nLineRec; i++)	{
		if( FindTags_As_First_String(szRtf[i], "IMPR")==1 )	{
			nItem = Split_Str_One_Line(szRtf[i]);
			
			for(j=1; j<nItem; j+=4)	{
				iAtom_1 = Query_Atom_Name(szItemLine[j]);
				iAtom_2 = Query_Atom_Name(szItemLine[j+1]);
				iAtom_3 = Query_Atom_Name(szItemLine[j+2]);
				iAtom_4 = Query_Atom_Name(szItemLine[j+3]);

				ImproperList[nImproper][0] = iAtom_1;
				ImproperList[nImproper][1] = iAtom_2;
				ImproperList[nImproper][2] = iAtom_3;
				ImproperList[nImproper][3] = iAtom_4;

				nImproper++;
			}
		}
	}

	return nImproper;
}



int Extract_Bond_Info()
{
	int i, j, nItem, iAtom_1, iAtom_2;

	nBond = 0;

	memset(AtomBond, 0, sizeof(int)*nAtom);

	for(i=0; i<nLineRec; i++)	{
		if( (FindTags_As_First_String(szRtf[i], "BOND")==1) || (FindTags_As_First_String(szRtf[i], "DOUBLE")==1) )	{
			nItem = Split_Str_One_Line(szRtf[i]);
			
			for(j=1; j<nItem; j+=2)	{
				if(strcmp(szItemLine[j], "!")==0)	{
					continue;
				}

				iAtom_1 = Query_Atom_Name(szItemLine[j]);
				iAtom_2 = Query_Atom_Name(szItemLine[j+1]);

				BondList[nBond][0] = iAtom_1;
				BondList[nBond][1] = iAtom_2;

				Bond_Array[iAtom_1][AtomBond[iAtom_1]] = iAtom_2;
				AtomBond[iAtom_1] ++;

				Bond_Array[iAtom_2][AtomBond[iAtom_2]] = iAtom_1;
				AtomBond[iAtom_2] ++;

				nBond++;
			}
		}
	}

	return nBond;
}

int Extract_Atom_Info()
{
	int i, ReadItem;
	double chg=0.0;
	char szTmp[256];

	nAtom = 0;
	for(i=0; i<nLineRec; i++)	{
		if(FindTags_As_First_String(szRtf[i], "ATOM")==1)	{
			ReadItem = sscanf(szRtf[i], "%s%s%s%lf", szTmp, szAtomName[nAtom], szChemName[nAtom], &chg);
			if(ReadItem == 4)	{
				nAtom++;
			}
			else	{
				printf("Error in line: %s\nQuit\n", szRtf[i]);
				exit(1);
			}
		}
	}

	return nAtom;
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

int FindTags_As_First_String(char szBuff[], char szTag[])
{
	int ReadItem;
	char szFirstStr[256];

	ReadItem = sscanf(szBuff, "%s", szFirstStr);;
	if(ReadItem == 1)	{
		if(strcmp(szFirstStr, szTag) == 0)	{
			return 1; 
		}
	}
	return 0;
}

void Setup_Dist_matrix(void)
{
	int i, j, k, l, Atom_i, Atom_j, Atom_k, Atom_l;

	nDihedral = 0;
	
	//start	to constrcut distance matrix
	memset(DistMatrix, 99, sizeof(char)*MAX_ATOM*MAX_ATOM);	//99, an arbitrary large numer
	for(i=0; i<nAtom; i++)	{	// i-i, itself. |d| = zero
		DistMatrix[i][i] = 0;
	}

	for(i=0; i<nAtom; i++)	{	//start	enumeration
		Atom_i = i;

		for(j=0; j<AtomBond[Atom_i]; j++)	{	//i->j
			Atom_j = Bond_Array[Atom_i][j];

			if(DistMatrix[Atom_i][Atom_j] > 1)	{
				DistMatrix[Atom_i][Atom_j] = 1;
				DistMatrix[Atom_j][Atom_i] = 1;
			}

			for(k=0; k<AtomBond[Atom_j]; k++)	{	//i->j->k
				Atom_k = Bond_Array[Atom_j][k];

				if(DistMatrix[Atom_i][Atom_k] > 2)	{
					DistMatrix[Atom_i][Atom_k] = 2;
					DistMatrix[Atom_k][Atom_i] = 2;
				}
				
				for(l=0; l<AtomBond[Atom_k]; l++)	{	//i->j->k->l
					Atom_l = Bond_Array[Atom_k][l];
					
					if(DistMatrix[Atom_i][Atom_l] > 3)	{	// all dihedral pair
						DistMatrix[Atom_i][Atom_l] = 3;
						DistMatrix[Atom_l][Atom_i] = 3;

						DihedralList[nDihedral][0] = Atom_i;
						DihedralList[nDihedral][1] = Atom_j;
						DihedralList[nDihedral][2] = Atom_k;
						DihedralList[nDihedral][3] = Atom_l;
						nDihedral ++;
					}
				}
			}
		}
	}
	//end	to constrcut distance matrix

}


void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("./error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in extract_mini_ff.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}
