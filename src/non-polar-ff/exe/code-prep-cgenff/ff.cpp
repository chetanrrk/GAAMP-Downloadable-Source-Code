#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ff.h"

#define MAX_LEN_LINE	(196)
#define ITEM_IN_LINE_8	(8)
#define ITEM_IN_LINE_9	(9)

extern FILE *fFile_Run_Log;	// will be shared by other source code
extern void Quit_With_Error_Msg(char szMsg[]);


CForceField::CForceField()
{
	n_Rec_Bond = n_Rec_Angle = n_Rec_Dihedral = n_Rec_ImproDihedral = n_Rec_LJ = n_Rec_NBFix = 0;
	Quit_on_Error = 0;
}

CForceField::~CForceField()
{
}

void CForceField::ReadForceField(char szNamePara[])
{
	FILE *fIn;
	char szLine[256], *ReadLine;
//	int n_Len, ReadItem, i, j, iTmp, iCount;
	int n_Len, ReadItem;

	fIn = fopen(szNamePara, "r");
	if(fIn == NULL)	{
//		printf("CForceField::ReadForceField> Fail to open file %s\nQuit\n", szNamePara);
		fprintf(fFile_Run_Log, "CForceField::ReadForceField> Fail to open file %s\nQuit\n", szNamePara);
		fflush(fFile_Run_Log);
		exit(1);
	}

	while(1)	{	//bond parameters
		if(feof(fIn))	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for BOND is found.\n");
			fflush(fFile_Run_Log);
			break;
		}
		ReadLine = fgets(szLine, MAX_LEN_LINE, fIn);
		if(ReadLine == NULL)	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for BOND is found.\n");
			fflush(fFile_Run_Log);
			break;
		}

		if( strncmp(szLine, "BONDS", 5) == 0 )	{
			break;
		}
	}

	while(1)	{
		if(feof(fIn))	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for ANGLE is found.\n");
			fflush(fFile_Run_Log);
			break;
		}
		ReadLine = fgets(szLine, MAX_LEN_LINE, fIn);
		if(ReadLine == NULL)	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for ANGLE is found.\n");
			fflush(fFile_Run_Log);
			break;
		}

		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		if(strncmp(szLine, "ANGLES", 6) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %s %lf %lf", 
			Bond_Rec[n_Rec_Bond].Chem[0], Bond_Rec[n_Rec_Bond].Chem[1], 
			 &(Bond_Rec[n_Rec_Bond].para[0]), &(Bond_Rec[n_Rec_Bond].para[1]));
		if( ReadItem == 4 )	{
			n_Rec_Bond++;
		}
	}
	if(n_Rec_Bond > N_REC_MAX)	{
//		printf("n_Rec_Bond > N_REC_MAX. n_Rec_Bond = %d\nQuit\n", n_Rec_Bond);
		fprintf(fFile_Run_Log, "n_Rec_Bond > N_REC_MAX. n_Rec_Bond = %d\nQuit\n", n_Rec_Bond);
		fflush(fFile_Run_Log);
		exit(1);
	}

	while(1)	{
		if(feof(fIn))	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for DIHEDRALS is found.\n");
			fflush(fFile_Run_Log);
			break;
		}
		ReadLine = fgets(szLine, MAX_LEN_LINE, fIn);
		if(ReadLine == NULL)	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for DIHEDRALS is found.\n");
			fflush(fFile_Run_Log);
			break;
		}

		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		if(strncmp(szLine, "DIHEDRALS", 9) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %s %s %lf %lf %lf %lf", 
			Angle_Rec[n_Rec_Angle].Chem[0], Angle_Rec[n_Rec_Angle].Chem[1], Angle_Rec[n_Rec_Angle].Chem[2], 
			 &(Angle_Rec[n_Rec_Angle].para[0]), &(Angle_Rec[n_Rec_Angle].para[1]), 
			 &(Angle_Rec[n_Rec_Angle].para[2]), &(Angle_Rec[n_Rec_Angle].para[3]));
		if( ReadItem == 5 )	{
			Angle_Rec[n_Rec_Angle].para[2] = 0.0;
			Angle_Rec[n_Rec_Angle].para[3] = 0.0;

			n_Rec_Angle++;
		}
		else if( ReadItem == 7 )	{
			n_Rec_Angle++;
		}
	}
	if(n_Rec_Angle > N_REC_MAX)	{
//		printf("n_Rec_Angle > N_REC_MAX. n_Rec_Angle = %d\nQuit\n", n_Rec_Angle);
		fprintf(fFile_Run_Log, "n_Rec_Angle > N_REC_MAX. n_Rec_Angle = %d\nQuit\n", n_Rec_Angle);
		fflush(fFile_Run_Log);
		exit(1);
	}


	while(1)	{
		if(feof(fIn))	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for IMPROPER is found.\n");
			fflush(fFile_Run_Log);
			break;
		}
		ReadLine = fgets(szLine, MAX_LEN_LINE, fIn);
		if(ReadLine == NULL)	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for IMPROPER is found.\n");
			fflush(fFile_Run_Log);
			break;
		}
		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		if(strncmp(szLine, "IMPROPER", 8) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %s %s %s %lf %lf %lf", 
			Dihedral_Rec[n_Rec_Dihedral].Chem[0], Dihedral_Rec[n_Rec_Dihedral].Chem[1], Dihedral_Rec[n_Rec_Dihedral].Chem[2], Dihedral_Rec[n_Rec_Dihedral].Chem[3], 
			 &(Dihedral_Rec[n_Rec_Dihedral].para[0]), &(Dihedral_Rec[n_Rec_Dihedral].para[1]), &(Dihedral_Rec[n_Rec_Dihedral].para[2]));
		if(ReadItem != 7)	{
			continue;
		}
		n_Rec_Dihedral++;
	}
	if(n_Rec_Dihedral > N_REC_MAX)	{
//		printf("n_Rec_Dihedral > N_REC_MAX. n_Rec_Dihedral = %d\nQuit\n", n_Rec_Dihedral);
		fprintf(fFile_Run_Log, "n_Rec_Dihedral > N_REC_MAX. n_Rec_Dihedral = %d\nQuit\n", n_Rec_Dihedral);
		fflush(fFile_Run_Log);
		exit(1);
	}


	while(1)	{
		if(feof(fIn))	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for NONBONDED is found.\n");
			fflush(fFile_Run_Log);
			break;
		}
		ReadLine = fgets(szLine, MAX_LEN_LINE, fIn);
		if(ReadLine == NULL)	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for NONBONDED is found.\n");
			fflush(fFile_Run_Log);
			break;
		}


		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		if(strncmp(szLine, "NONBONDED", 9) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %s %s %s %lf %lf %lf", 
			ImproDihedral_Rec[n_Rec_ImproDihedral].Chem[0], ImproDihedral_Rec[n_Rec_ImproDihedral].Chem[1], ImproDihedral_Rec[n_Rec_ImproDihedral].Chem[2], ImproDihedral_Rec[n_Rec_ImproDihedral].Chem[3], 
			 &(ImproDihedral_Rec[n_Rec_ImproDihedral].para[0]), &(ImproDihedral_Rec[n_Rec_ImproDihedral].para[1]), &(ImproDihedral_Rec[n_Rec_ImproDihedral].para[2]));
		if(ReadItem == 7)	{
			n_Rec_ImproDihedral++;
		}
	}
	if(n_Rec_ImproDihedral > N_REC_MAX)	{
		fprintf(fFile_Run_Log, "n_Rec_ImproDihedral > N_REC_MAX. n_Rec_ImproDihedral = %d\nQuit\n", n_Rec_ImproDihedral);
		fflush(fFile_Run_Log);
		exit(1);
	}

/*
	// only read the first entry for CMAP
	//start	to read CMap parameters 24X24 from topar file
	double CMap_data[5];	//five items per line
	for(i=0; i<CMAP_DIM; i++)	{
		while(1)	{
			fgets(szLine, MAX_LEN_LINE, fIn);
			ReadItem = sscanf(szLine+1, "%d", &iTmp);
			if( (ReadItem == 1) && (szLine[0] == '!') )	{	//the entry of CMap data
				break;
			}
		}

		iCount = 0;
		while(1)	{
			fgets(szLine, MAX_LEN_LINE, fIn);
			ReadItem = sscanf(szLine, "%lf %lf %lf %lf %lf", 
				&(CMap_data[0]), &(CMap_data[1]), &(CMap_data[2]), &(CMap_data[3]), &(CMap_data[4]));
			if( ReadItem == 5 )	{
				for(j=0; j<5; j++)	{
					mctp[0][i][iCount+j] = CMap_data[j];
				}
				iCount += 5;
			}
			else if( (iCount == 20) && (ReadItem == 4) )	{
				for(j=0; j<5; j++)	{
					mctp[0][i][iCount+j] = CMap_data[j];
				}
				iCount += 4;
				break;
			}
			else	{
				fprintf(fFile_Run_Log, "Please check your parameter file around the entry of CMap.\n");
				fflush(fFile_Run_Log);
				exit(1);
			}
		}
	}
//	setcmap(24, 12, mctp, 360.0/24.0);
	setcmap(24, 12, mctp, 15.0);
	//end	to read CMap parameters 24X24 from topar file
*/

/*	
	while(1)	{
		if(feof(fIn))	{
			fprintf(fFile_Run_Log, "Fail to find the entry for NONBONDED parameters.\nQuit\n");
			fflush(fFile_Run_Log);
			exit(1);
		}
		fgets(szLine, MAX_LEN_LINE, fIn);
		
		if(strncmp(szLine, "NONBONDED", 9) == 0 )	{
			break;
		}
	}
*/
	
	while(1)	{
		if(feof(fIn))	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for NBFIX is found.\n");
			fflush(fFile_Run_Log);
			break;
		}
		ReadLine = fgets(szLine, MAX_LEN_LINE, fIn);
		if(ReadLine == NULL)	{
			fprintf(fFile_Run_Log, "FitCharge> Warning: No entry for NBFIX is found.\n");
			fflush(fFile_Run_Log);
			break;
		}

		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}
		if(strncmp(szLine, "NBFIX", 5) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %lf %lf %lf %lf %lf %lf", 
			LJ_Rec[n_Rec_LJ].Chem, 
			&(LJ_Rec[n_Rec_LJ].para[0]), &(LJ_Rec[n_Rec_LJ].para[1]), &(LJ_Rec[n_Rec_LJ].para[2]), 	//LJ parameters
			&(LJ_Rec[n_Rec_LJ].para[3]), &(LJ_Rec[n_Rec_LJ].para[4]), &(LJ_Rec[n_Rec_LJ].para[5]));	//LJ 1-4 parameters
		if(ReadItem == 4)	{
			LJ_Rec[n_Rec_LJ].para[3] = LJ_Rec[n_Rec_LJ].para[0];	//assign 1-4 parameters same as non-1-4 parameters defaultly
			LJ_Rec[n_Rec_LJ].para[4] = LJ_Rec[n_Rec_LJ].para[1];
			LJ_Rec[n_Rec_LJ].para[5] = LJ_Rec[n_Rec_LJ].para[2];

			n_Rec_LJ++;
		}
		else if(ReadItem == 7)	{	//with full parameters
			n_Rec_LJ++;
		}
		if(n_Rec_LJ > N_REC_MAX)	{
			fprintf(fFile_Run_Log, "n_Rec_LJ > N_REC_MAX. n_Rec_LJ = %d\nQuit\n", n_Rec_LJ);
			fflush(fFile_Run_Log);
			exit(1);
		}
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szLine, MAX_LEN_LINE, fIn);
		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		ReadItem = sscanf(szLine, "%s %s %lf %lf", 
			NBFix_Rec[n_Rec_NBFix].Chem_1, NBFix_Rec[n_Rec_NBFix].Chem_2, 
			&(NBFix_Rec[n_Rec_NBFix].para[0]), &(NBFix_Rec[n_Rec_NBFix].para[1]));
		if(ReadItem == 4)	{
			n_Rec_NBFix++;
			if(n_Rec_NBFix > N_REC_MAX)	{
				fprintf(fFile_Run_Log, "n_Rec_NBFix > N_REC_MAX. n_Rec_NBFix = %d\nQuit\n", n_Rec_NBFix);
				fflush(fFile_Run_Log);
				exit(1);
			}
		}
	}


	fclose(fIn);
}

int Compare_Chem_Type(char Type_1[], char Type_2[])
{
	if( (strcmp(Type_1, "X")==0) || (strcmp(Type_2, "X")==0) )	{	// same type
		return 0;
	}
	return (strcmp(Type_1, Type_2));
}


void CForceField::GetPara_Bond(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i;

	for(i=0; i<n_Rec_Bond; i++)	{
		if( (strcmp(szChemName[0], Bond_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], Bond_Rec[i].Chem[1])==0) )	{
			Para[0] = Bond_Rec[i].para[0];
			Para[1] = Bond_Rec[i].para[1];

			Active_Bond[i] = 1;
			return;
		}
	}

	for(i=0; i<n_Rec_Bond; i++)	{
		if( (strcmp(szChemName[1], Bond_Rec[i].Chem[0])==0) && (strcmp(szChemName[0], Bond_Rec[i].Chem[1])==0) )	{
			Para[0] = Bond_Rec[i].para[0];
			Para[1] = Bond_Rec[i].para[1];

			Active_Bond[i] = 1;
			return;
		}
	}

//	printf("\n\nFatal error!\nCan't find the parameters for bond %s - %s\nQuit.\n", szChemName[0], szChemName[1]);
	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for bond %s - %s\nQuit.\n", szChemName[0], szChemName[1]);
	fflush(fFile_Run_Log);

	if(Quit_on_Error)	exit(1);


	return;
}

void CForceField::GetPara_Angle(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i;

	for(i=0; i<n_Rec_Angle; i++)	{
		if( (strcmp(szChemName[0], Angle_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], Angle_Rec[i].Chem[1])==0) && (strcmp(szChemName[2], Angle_Rec[i].Chem[2])==0) )	{
			Para[0] = Angle_Rec[i].para[0];
			Para[1] = Angle_Rec[i].para[1];
			Para[2] = Angle_Rec[i].para[2];
			Para[3] = Angle_Rec[i].para[3];

			Active_Angle[i] = 1;
			return;
		}
	}

	for(i=0; i<n_Rec_Angle; i++)	{
		if( (strcmp(szChemName[2], Angle_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], Angle_Rec[i].Chem[1])==0) && (strcmp(szChemName[0], Angle_Rec[i].Chem[2])==0) )	{
			Para[0] = Angle_Rec[i].para[0];
			Para[1] = Angle_Rec[i].para[1];
			Para[2] = Angle_Rec[i].para[2];
			Para[3] = Angle_Rec[i].para[3];

			Active_Angle[i] = 1;
			return;
		}
	}

	for(i=0; i<n_Rec_Angle; i++)	{
		if( (Compare_Chem_Type(szChemName[0], Angle_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[1], Angle_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[2], Angle_Rec[i].Chem[2])==0) )	{
			Para[0] = Angle_Rec[i].para[0];
			Para[1] = Angle_Rec[i].para[1];
			Para[2] = Angle_Rec[i].para[2];
			Para[3] = Angle_Rec[i].para[3];

			Active_Angle[i] = 1;
			return;
		}
	}

	for(i=0; i<n_Rec_Angle; i++)	{
		if( (Compare_Chem_Type(szChemName[2], Angle_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[1], Angle_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[0], Angle_Rec[i].Chem[2])==0) )	{
			Para[0] = Angle_Rec[i].para[0];
			Para[1] = Angle_Rec[i].para[1];
			Para[2] = Angle_Rec[i].para[2];
			Para[3] = Angle_Rec[i].para[3];

			Active_Angle[i] = 1;
			return;
		}
	}

//	printf("\n\nFatal error!\nCan't find the parameters for angle %s, %s, %s\nQuit.\n", 
	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for angle %s, %s, %s\nQuit.\n", 
		szChemName[0], szChemName[1], szChemName[2]);
	fflush(fFile_Run_Log);
	if(Quit_on_Error)	exit(1);

	return;
}

void CForceField::GetPara_Dihedral(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i, Idx, Count, iPos;

	Count = 0;
	iPos = 0;

	for(i=0; i<n_Rec_Dihedral; i++)	{
		if( (strcmp(szChemName[0], Dihedral_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], Dihedral_Rec[i].Chem[1])==0) && (strcmp(szChemName[2], Dihedral_Rec[i].Chem[2])==0) && (strcmp(szChemName[3], Dihedral_Rec[i].Chem[3])==0))	{
			Idx = (int)(Dihedral_Rec[i].para[1] + 1.0E-6);
			iPos = Idx*3;
			Para[iPos  ] = Dihedral_Rec[i].para[0];
			Para[iPos+1] = Dihedral_Rec[i].para[1];
			Para[iPos+2] = Dihedral_Rec[i].para[2];

			Active_Dihedral[i] = 1;
			Count++;
		}
	}
	if(Count > 0)	{
		return;
	}
	for(i=0; i<n_Rec_Dihedral; i++)	{
		if( (strcmp(szChemName[3], Dihedral_Rec[i].Chem[0])==0) && (strcmp(szChemName[2], Dihedral_Rec[i].Chem[1])==0) && (strcmp(szChemName[1], Dihedral_Rec[i].Chem[2])==0) && (strcmp(szChemName[0], Dihedral_Rec[i].Chem[3])==0))	{
			Idx = (int)(Dihedral_Rec[i].para[1] + 1.0E-6);
			iPos = Idx*3;
			Para[iPos  ] = Dihedral_Rec[i].para[0];
			Para[iPos+1] = Dihedral_Rec[i].para[1];
			Para[iPos+2] = Dihedral_Rec[i].para[2];

			Active_Dihedral[i] = 1;
			Count++;
		}
	}
	if(Count > 0)	{
		return;
	}


	for(i=0; i<n_Rec_Dihedral; i++)	{
		if( (Compare_Chem_Type(szChemName[0], Dihedral_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[1], Dihedral_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[2], Dihedral_Rec[i].Chem[2])==0) && (Compare_Chem_Type(szChemName[3], Dihedral_Rec[i].Chem[3])==0))	{
			Idx = (int)(Dihedral_Rec[i].para[1] + 1.0E-6);
			iPos = Idx*3;
			Para[iPos  ] = Dihedral_Rec[i].para[0];
			Para[iPos+1] = Dihedral_Rec[i].para[1];
			Para[iPos+2] = Dihedral_Rec[i].para[2];

			Active_Dihedral[i] = 1;
			Count++;
		}
	}
	if(Count > 0)	{
		return;
	}
	for(i=0; i<n_Rec_Dihedral; i++)	{
		if( (Compare_Chem_Type(szChemName[3], Dihedral_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[2], Dihedral_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[1], Dihedral_Rec[i].Chem[2])==0) && (Compare_Chem_Type(szChemName[0], Dihedral_Rec[i].Chem[3])==0))	{
			Idx = (int)(Dihedral_Rec[i].para[1] + 1.0E-6);
			iPos = Idx*3;
			Para[iPos  ] = Dihedral_Rec[i].para[0];
			Para[iPos+1] = Dihedral_Rec[i].para[1];
			Para[iPos+2] = Dihedral_Rec[i].para[2];

			Active_Dihedral[i] = 1;
			Count++;
		}
	}
	if(Count > 0)	{
		return;
	}


	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for dihedral %s, %s, %s, %s\nQuit.\n", 
		szChemName[0], szChemName[1], szChemName[2], szChemName[3]);
	fflush(fFile_Run_Log);
	if(Quit_on_Error)	exit(1);

	return;
}

void CForceField::GetPara_ImproDIhedral(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i;

	for(i=0; i<n_Rec_ImproDihedral; i++)	{
		if( (strcmp(szChemName[0], ImproDihedral_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], ImproDihedral_Rec[i].Chem[1])==0) && (strcmp(szChemName[2], ImproDihedral_Rec[i].Chem[2])==0) && (strcmp(szChemName[3], ImproDihedral_Rec[i].Chem[3])==0))	{
			Para[0] = ImproDihedral_Rec[i].para[0];
			Para[1] = ImproDihedral_Rec[i].para[2];
			Para[2] = ImproDihedral_Rec[i].para[1];	// potential type

			Active_ImproDihedral[i] = 1;
			return;
		}
	}

	for(i=0; i<n_Rec_ImproDihedral; i++)	{
		if( (strcmp(szChemName[3], ImproDihedral_Rec[i].Chem[0])==0) && (strcmp(szChemName[2], ImproDihedral_Rec[i].Chem[1])==0) && (strcmp(szChemName[1], ImproDihedral_Rec[i].Chem[2])==0) && (strcmp(szChemName[0], ImproDihedral_Rec[i].Chem[3])==0))	{
			Para[0] = ImproDihedral_Rec[i].para[0];
			Para[1] = ImproDihedral_Rec[i].para[2];
			Para[2] = ImproDihedral_Rec[i].para[1];	// potential type

			Active_ImproDihedral[i] = 1;
			return;
		}
	}

	for(i=0; i<n_Rec_ImproDihedral; i++)	{
		if( (Compare_Chem_Type(szChemName[0], ImproDihedral_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[1], ImproDihedral_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[2], ImproDihedral_Rec[i].Chem[2])==0) && (Compare_Chem_Type(szChemName[3], ImproDihedral_Rec[i].Chem[3])==0))	{
			Para[0] = ImproDihedral_Rec[i].para[0];
			Para[1] = ImproDihedral_Rec[i].para[2];
			Para[2] = ImproDihedral_Rec[i].para[1];	// potential type

			Active_ImproDihedral[i] = 1;
			return;
		}
	}

	for(i=0; i<n_Rec_ImproDihedral; i++)	{
		if( (Compare_Chem_Type(szChemName[3], ImproDihedral_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[2], ImproDihedral_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[1], ImproDihedral_Rec[i].Chem[2])==0) && (Compare_Chem_Type(szChemName[0], ImproDihedral_Rec[i].Chem[3])==0))	{
			Para[0] = ImproDihedral_Rec[i].para[0];
			Para[1] = ImproDihedral_Rec[i].para[2];
			Para[2] = ImproDihedral_Rec[i].para[1];	// potential type

			Active_ImproDihedral[i] = 1;
			return;
		}
	}

	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for improper dihedral %s, %s, %s, %s\nQuit.\n", 
		szChemName[0], szChemName[1], szChemName[2], szChemName[3]);
	fflush(fFile_Run_Log);
	if(Quit_on_Error)	exit(1);

	return;
}

int CForceField::GetPara_LJ(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i;

	if(strcmp(szChemName[0], "-----") == 0)	{	//drude particle
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return (-1);
	}

	if(strcmp(szChemName[0], "DOH2") == 0)	{	//drude particle
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return (-1);
	}
	if(strcmp(szChemName[0], "DRUD") == 0)	{	//drude particle
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return (-1);
	}

	if(strcmp(szChemName[0], "LP") == 0)	{	//drude particle
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return (-1);
	}

	for(i=0; i<n_Rec_LJ; i++)	{
		if( strcmp(szChemName[0], LJ_Rec[i].Chem)==0 )	{
			Para[0] = LJ_Rec[i].para[0];
			Para[1] = LJ_Rec[i].para[1];
			Para[2] = LJ_Rec[i].para[2];
			Para[3] = LJ_Rec[i].para[3];
			Para[4] = LJ_Rec[i].para[4];
			Para[5] = LJ_Rec[i].para[5];

			Active_LJ[i] = 1;
			return i;
		}
	}



//	printf("\n\nFatal error!\nCan't find the parameters for LJ, %s\nQuit.\n", 
	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for LJ, %s\nQuit.\n", 
		szChemName[0]);
	fflush(fFile_Run_Log);
	exit(1);

	return (-1);
}

int  CForceField::QueryNBFix(char szChem_1[], char szChem_2[])
{
	int i;

	for(i=0; i<n_Rec_NBFix; i++)	{
		if( ((strcmp(szChem_1, NBFix_Rec[i].Chem_1)==0) && (strcmp(szChem_2, NBFix_Rec[i].Chem_2)==0)) || ((strcmp(szChem_2, NBFix_Rec[i].Chem_1)==0) && (strcmp(szChem_1, NBFix_Rec[i].Chem_2)==0)) )	{

			Active_NBFix[i] = 1;
			return i;
		}
	}
	return -1;
}
