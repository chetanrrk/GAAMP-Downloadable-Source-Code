#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#else
#include <unistd.h>	//for linux
#endif

//#include <mpi.h>
#include "ff.h"

char szForceFiled[]="mol.prm";
char szXpsfFile[]="mol.xpsf";
char szCrdFile[]="mol-opt.crd";
char szPhiToScan[]="qm-1d-states.dat";

char szExe_G09[256];
char szQMParaFile[]="../QM-para.txt";
char szMyPath[]="../mypath.txt";
char szFileElem[]="elem-list.txt";

char szQM_Level_Rotamer_Cmd[512]="%mem=1600MB\n%nproc=4\n# HF/6-31G* nosymm opt\n\nMol rotamer\n\n%d 1\n";	// kbt
const char szGAUSS_SCRDIR[]="GAUSS_SCRDIR";
char szGAUSS_SCRDIR_Base[256];


#define N_MAX_DIH	(16)
#define MAX_N_STATE		(6)
#define MAX_ATOM	(256)
#define MAX_BIN_PHI	(128)
#define BIN_SIZE	(5)
#define BIN_SIZE_QM	(10)
#define N_MAX_QM_PART	(6)
#define MAX_CONF	(65536)
#define N_TOR_PARA	(11)

#define QM_LEVEL_HF		(1)
#define QM_LEVEL_MP2	(2)


CMol Mol;
CForceField ForceField;
int nDihedral, QM_Level;

FILE *fFile_Run_Log;
void Quit_With_Error_Msg(char szMsg[]);

int n_Phi=0, Counter=0;
int DihList[N_MAX_DIH][4], IdxDihSelect[N_MAX_DIH], n_State_List[N_MAX_DIH], n_State_List_Save[N_MAX_DIH];
double Phi_Set[N_MAX_DIH][MAX_N_STATE], Phi_Set_Save[N_MAX_DIH][MAX_N_STATE], Phi_To_Set[N_MAX_DIH], Phi_To_Set_Scan[N_MAX_DIH];
double x_Save[MAX_ATOM], y_Save[MAX_ATOM], z_Save[MAX_ATOM];
double E_Phi[MAX_BIN_PHI], E_Scan[MAX_BIN_PHI];

double Para_Tor_List[N_MAX_DIH][N_TOR_PARA];

int n_QM_Part=0, Phi_QM_n_Conf[N_MAX_QM_PART];
char szQM_Scan[N_MAX_QM_PART][256];
double Phi_QM_Start[N_MAX_QM_PART];
double E_Conf[MAX_CONF], Phi_Conf[MAX_CONF][N_MAX_DIH];
char szElemName[MAX_ATOM][8];

int	ProgID=0, nProc=1;

extern int FindString(char szBuff[], char szTag[]);
void Replace_D_E(char szStr[]);
int To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[], int Quit);
int Skip_N_Line(FILE *fIn, char szFileName[], int n);
double Get_Energy(char szLine[]);
void ReadElementList(void);

int Read_Soft_DihedralList(void);
int Get_Dihedral_Index(int ia, int ib, int ic, int id);
void Enumerate_Dihedrals(int IdxDih, int iState);
void BackupCoordinates(void);
void RestoreCoordinates(void);
void SaveOptPdb(char szName[]);
void Output_Gaussian_File(void);
void Get_Netcharge_From_Xpsf(void);
double Extract_Rotamer_Coord_E(char szName[]);

void Get_EXE_Path(char szExeFile[], char szExePath[]);
void Setup_QM_Level(void);
int Split_Into_Two_String(char szBuff[], char str1[], char str2[]);
void Replace_NewLine(char szBuff[]);

void Read_Tor_Para_1D_Fitting(void);
void Assign_Torsion_Parameters(void);
void Determine_QM_Level(char szQMLevel[]);

int netcharge_mol;
int n_Conf=0, n_Rotamer=0;
double E_Min=1.0E100;




//start	data and functions related with job scheduler
#define MAX_JOB	(1024)

enum { IDLE, BUSY };
enum { Job_ToRun, Job_Running, Job_Done };

class CWORKER
{
public:
	int JobID, Status, nCore;
	char szInput[256], szOutput[256];

	void Start_Job(int CoreNum, int ID, char szGjf[], char szOut[]);	// start running a job
	int Get_Job_Status(void);
	void Update_Core_Number_Input_File(void);
};

CWORKER *pWorker;
int nJob=0, nJobDone=0, nWorker=0, nCore_Per_Node=1, n_QM_Part_List[N_MAX_DIH];
int JobStatus[MAX_JOB];
char szDirList[MAX_JOB][256], szInputList[MAX_JOB][256], szOutputList[MAX_JOB][256];

int Count_Finished_Job();
void RunAllJobs(void);
int Get_First_Available_Worker(void);
void Exatract_All_Info_QM_1D_Scan(void);
//end	data and functions related with job scheduler


#define E_RANGE	(5.0)
void Enumerate_Dihedrals(int IdxDih)	// depth
{
	int iState, i;
//	char szName[256];

	if(IdxDih < n_Phi)	{
		for(iState=0; iState<n_State_List[IdxDih]; iState++)	{
			Phi_To_Set[IdxDih] = Phi_Set[IdxDih][iState];
			Enumerate_Dihedrals(IdxDih+1);
		}
	}
	else	{
		RestoreCoordinates();
		for(i=0; i<n_Phi; i++)	{
			Mol.QueryDihedral(IdxDihSelect[i]);
			Mol.Edit_Dihedral(IdxDihSelect[i], Phi_To_Set[i]);
		}
		E_Conf[n_Conf] = Mol.Cal_E(0);
		memcpy(Phi_Conf[n_Conf], Phi_To_Set, sizeof(double)*n_Phi);
		if(E_Conf[n_Conf] < E_Min)	{
			E_Min = E_Conf[n_Conf];
		}
		n_Conf++;

//		sprintf(szName, "ini-%d.pdb", n_Conf);
//		Mol.SavePdb(szName);

		return;
	}
}

char szCurDir[256], szNewDir[256], szGjfFile[]="mol-qm-rotamer.gjf";
int main(int argc, char **argv)
{
	int i, j, nAtom;
	double E_Rotamer, E_Rotamer_Min=1.0E100;
	char szCmd[256], szOutput[]="output.out", *szEnv;
	FILE *fOut, *fIn;

	fFile_Run_Log = fopen("qm-rotamer-scan.log", "w");

	szEnv = getenv(szGAUSS_SCRDIR);
	if (szEnv == NULL) {
		Quit_With_Error_Msg("Environment variable $GAUSS_SCRDIR is NOT set.\nQuit\n");
	}
	strcpy(szGAUSS_SCRDIR_Base, szEnv);

	Get_EXE_Path("G09_EXE_PATH", szExe_G09);

	ForceField.ReadForceField(szForceFiled);
	Mol.ReadPSF(szXpsfFile, 0);
	Get_Netcharge_From_Xpsf();

	Setup_QM_Level();
	ReadElementList();

	Mol.AssignForceFieldParameters(&ForceField);
	Mol.ReadCRD(szCrdFile);
	BackupCoordinates();
	nAtom = Mol.nAtom;

	Read_Soft_DihedralList();

	Enumerate_Dihedrals(0);

	getcwd(szCurDir, 256);

	n_Rotamer = 0;
	nJob = 0;
	for(i=0; i<n_Conf; i++)	{
		if( n_Conf > 200 )	{
			if(E_Conf[i] > (E_Min + 100.0) )	{
				continue;
			}
		}
		else	{
			if(E_Conf[i] > (E_Min + 600.0) )	{
				continue;
			}
		}
		chdir(szCurDir);
		RestoreCoordinates();
		for(j=0; j<n_Phi; j++)	{
			Mol.QueryDihedral(IdxDihSelect[j]);
			Mol.Edit_Dihedral(IdxDihSelect[j], Phi_Conf[i][j]);
		}
		
		sprintf(szNewDir, "%s/rotamer-%d", szCurDir, n_Rotamer+1);
		sprintf(szCmd, "mkdir rotamer-%d", n_Rotamer+1);
		system(szCmd);
		chdir(szNewDir);

		Output_Gaussian_File();

		strcpy(szDirList[nJob], szNewDir);
		strcpy(szInputList[nJob], szGjfFile);
		sprintf(szOutputList[nJob], "%s/%s", szNewDir, szOutput);
		
		n_Rotamer++;
		nJob++;
	}

	RunAllJobs();
	chdir(szCurDir);
		
	fOut = fopen("all-rotamer.dat", "w");

	int Finished;
	char szLine[256];
	for(i=0; i<nJob; i++)	{
		fIn = fopen(szOutputList[i], "r");
		Finished = To_Find_Tag(fIn, szOutput, " Optimization completed", szLine, 0);
		fclose(fIn);
		if(Finished < 0)	{	// un-successful optimization
			continue;
		}

		E_Rotamer = Extract_Rotamer_Coord_E(szOutputList[i]);
		if(E_Rotamer != 0.0)	{	// a valid energy
			fprintf(fOut, "E_Rotamer %.13E Conf %4d\nCoordinate\n", E_Rotamer, i+1);
			for(j=0; j<nAtom; j++)	{
				fprintf(fOut, "%14.6lf%14.6lf%14.6lf\n", Mol.x[j], Mol.y[j], Mol.z[j]);
			}
		}

		if(E_Rotamer < E_Rotamer_Min)	{
			Mol.SavePdb("QM-min.pdb");
			E_Rotamer_Min = E_Rotamer;
		}
	}

	fclose(fOut);


	fclose(fFile_Run_Log);

	return 0;
}



double Extract_Rotamer_Coord_E(char szName[])
{
	FILE *fIn;
	int Finished, j, iTmp, nAtom, ReadItem;
	char szLine[256], ErrorMsg[256], szStr_Energy[256], *ReadLine;
	double E_Rotamer;

	nAtom = Mol.nAtom;
	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Error to open file: %s\nQuit\n", szName);
		Quit_With_Error_Msg(ErrorMsg);
	}
	
	Finished = To_Find_Tag(fIn, szName, " Optimization completed", szLine, 0);
	if(Finished < 0)	{	// skip unfinished jobs
		printf(ErrorMsg, "Warning> Can't find the tag of Optimization completed in file: %s\n", szName);
		fclose(fIn);
		return 0.0;	// an invalid energy
	}

	
	To_Find_Tag(fIn, szName, " orientation:", szLine, 1);
	Skip_N_Line(fIn, szName, 4);	// including the coordinates of X1 and X2.
	
	for(j=0; j<nAtom; j++)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szLine, 256, fIn);
		ReadItem = sscanf(szLine, "%d %d %d %lf %lf %lf", &iTmp, &iTmp, &iTmp, &(Mol.x[j]), &(Mol.y[j]), &(Mol.z[j]));
		if(ReadItem!=6)	{
			break;
		}
	}
	
	if(j==nAtom)	{	// correct number of atoms as expected
		fseek(fIn, 0, SEEK_SET);
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				if( (FindString(szLine, "SCF Done:  E(RHF)") >= 0) && (QM_Level == QM_LEVEL_HF) )	{
					sscanf(szLine+20, "%s", szStr_Energy);
					Replace_D_E(szStr_Energy);
					ReadItem = sscanf(szStr_Energy, "%lf", &E_Rotamer);
					if(ReadItem != 1)	{
						sprintf(ErrorMsg, "Fail to extract energy from : %s\nQuit\n", szStr_Energy);
						Quit_With_Error_Msg(ErrorMsg);
					}
				}
				else if( (FindString(szLine, "EUMP2 =") >= 0) && (QM_Level == QM_LEVEL_MP2) )	{
					sscanf(szLine+34, "%s", szStr_Energy);
					Replace_D_E(szStr_Energy);
					ReadItem = sscanf(szStr_Energy, "%lf", &E_Rotamer);
					if(ReadItem != 1)	{
						sprintf(ErrorMsg, "Fail to extract energy from : %s\nQuit\n", szStr_Energy);
						Quit_With_Error_Msg(ErrorMsg);
					}
				}
			}
			
		}
	}
	
	fclose(fIn);

	return E_Rotamer;
}

void Output_Gaussian_File(void)
{
	FILE *fOut;
	int i, nAtom;

	nAtom = Mol.nAtom;
	
	fOut = fopen(szGjfFile, "w");
	
	fprintf(fOut, "%s", szQM_Level_Rotamer_Cmd);
	for(i=0; i<nAtom; i++)	{
//		fprintf(fOut, "%c  %9.5lf %9.5lf %9.5lf\n", Mol.AtomName[i][0], Mol.x[i], Mol.y[i], Mol.z[i]);
		fprintf(fOut, "%s  %9.5lf %9.5lf %9.5lf\n", szElemName[i], Mol.x[i], Mol.y[i], Mol.z[i]);
	}
	fprintf(fOut, "\n");
	
	fclose(fOut);

//	char szName[256];
//	sprintf(szName, "ini-%d.pdb", n_Rotamer);
//	Mol.SavePdb(szName);

}

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("../error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in QM_rotamer_scan.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}

int Read_Soft_DihedralList(void)
{
	FILE *fIn;
	int ReadItem;
	char szLine[256], *ReadLine;

	n_Phi = 0;

	fIn = fopen(szPhiToScan, "r");

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 128, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		ReadItem = sscanf(szLine, "%d %d %d %d %lf %lf %lf %lf %lf %lf", 
			&(DihList[n_Phi][0]), &(DihList[n_Phi][1]), &(DihList[n_Phi][2]), &(DihList[n_Phi][3]), 
			&(Phi_Set[n_Phi][0]), &(Phi_Set[n_Phi][1]), &(Phi_Set[n_Phi][2]), &(Phi_Set[n_Phi][3]), &(Phi_Set[n_Phi][4]), &(Phi_Set[n_Phi][5]));

		if(ReadItem > 4)	{
			DihList[n_Phi][0]--;
			DihList[n_Phi][1]--;
			DihList[n_Phi][2]--;
			DihList[n_Phi][3]--;
			IdxDihSelect[n_Phi] = Mol.Query_Dihedral_Index(DihList[n_Phi][0], DihList[n_Phi][1], DihList[n_Phi][2], DihList[n_Phi][3]);

			if(IdxDihSelect[n_Phi] < 0)	{
				Quit_With_Error_Msg("Fail to identify the index of one soft dihedral.\n");
			}
			Mol.BuildSegmentList_Dihedrals(IdxDihSelect[n_Phi]);

			n_State_List[n_Phi] = ReadItem - 4;
			n_Phi++;
		}
		else	{
			break;
		}
	}

	fclose(fIn);

	memcpy(Phi_Set_Save, Phi_Set, sizeof(double)*N_MAX_DIH*MAX_N_STATE);
	memcpy(n_State_List_Save, n_State_List, sizeof(int)*N_MAX_DIH);

	Read_Tor_Para_1D_Fitting();	// to assign the fitted torsion parameters

	return n_Phi;
}

void BackupCoordinates(void)
{
	int nAtom;
	nAtom = Mol.nAtom;

	memcpy(x_Save, Mol.x, sizeof(double)*nAtom);
	memcpy(y_Save, Mol.y, sizeof(double)*nAtom);
	memcpy(z_Save, Mol.z, sizeof(double)*nAtom);
}

void RestoreCoordinates(void)
{
	int nAtom;
	nAtom = Mol.nAtom;

	memcpy(Mol.x, x_Save, sizeof(double)*nAtom);
	memcpy(Mol.y, y_Save, sizeof(double)*nAtom);
	memcpy(Mol.z, z_Save, sizeof(double)*nAtom);
}


void SaveOptPdb(char szName[])
{
	int i;
	for(i=0; i<n_Phi; i++)	{
		Mol.QueryDihedral(IdxDihSelect[i]);
		Mol.Edit_Dihedral(IdxDihSelect[i], Phi_To_Set_Scan[i]);
	}
	Mol.SavePdb(szName);
}


void Get_Netcharge_From_Xpsf(void)
{
	int i, nAtom;
	double cg_sum=0.0;
	char szCharge[16];

	nAtom = Mol.nAtom;
	for(i=0; i<nAtom; i++)	{
		cg_sum += Mol.CG[i];
	}
	sprintf(szCharge, "%.0lf", cg_sum);

	netcharge_mol = atoi(szCharge);	// assuming the netcharge is an integer !!!
}


void Replace_D_E(char szStr[])
{
	int nLen, i;

	nLen = strlen(szStr);
	for(i=0; i<nLen; i++)	{
		if(szStr[i] == 'D')	{
			szStr[i] = 'E';
			return;
		}
	}
	return;
}

double Get_Energy(char szLine[])
{
	char szEnergy[256], ErrorMsg[256];
	int i=0, Count=0;

	while(1)	{
		if(szLine[i]=='=')	{
			i++;
			break;
		}
		else if(szLine[i]==0x0)	{
			sprintf(ErrorMsg, "Fail to extract the energy from, \n%s\nQuit\n", szLine);
			Quit_With_Error_Msg(ErrorMsg);
		}
		else	{
			i++;
		}
	}

	while(1)	{
		if(szLine[i] != ' ')	{
			break;
		}
		else if(szLine[i]==0x0)	{
			sprintf(ErrorMsg, "Fail to extract the energy from, \n%s\nQuit\n", szLine);
			Quit_With_Error_Msg(ErrorMsg);
		}
		else	{
			i++;
		}
	}

	while(1)	{
		if(szLine[i] == ' ')	{
			break;
		}
		else if(szLine[i]==0x0)	{
			if(Count == 0)	{
				sprintf(ErrorMsg, "Error in extracting the energy from, \n%s\nQuit\n", szLine);
				Quit_With_Error_Msg(ErrorMsg);
			}
			else	{
				break;
			}
		}
		else	{
			szEnergy[Count] = szLine[i];
			Count++;
			i++;
		}
	}
	szEnergy[Count] = 0;

	Replace_D_E(szEnergy);

	return atof(szEnergy);
}

int To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[], int Quit)
{
	char *ReadLine, ErrorMsg[256];

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
			if(Quit)	{
				fclose(fIn);
				Quit_With_Error_Msg(ErrorMsg);
			}
			else	{
				break;
			}
		}

		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
			if(Quit)	{
				fclose(fIn);
				Quit_With_Error_Msg(ErrorMsg);
			}
			else	{
				break;
			}
		}
		else	{
			if(FindString(szLine, szTag) >= 0)	{
				return 1;
			}
		}
	}

	return -1;
}

int Skip_N_Line(FILE *fIn, char szFileName[], int n)
{
	int i;
	char szLine[256], *ReadLine, ErrorMsg[256];

	for(i=0; i<n; i++)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail in Skip_N_Line(%s, %d).\nQuit\n", szFileName, n);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Fail in Skip_N_Line(%s, %d).\nQuit\n", szFileName, n);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}
	}
	return 1;
}


void Get_EXE_Path(char szExeFile[], char szExePath[])
{
	FILE *fIn;
	char szLine[256], szBuff[256], *ReadLine, ErrorMsg[256];
	int nLen;

	fIn = fopen(szMyPath, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szMyPath);
		Quit_With_Error_Msg(ErrorMsg);
	}

	nLen = strlen(szExeFile);
	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to find the path of %s in %s\nQuit\n", szExeFile, szMyPath);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		ReadLine = fgets(szLine, 256, fIn);

		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Fail to find the path of %s in %s\nQuit\n", szExeFile, szMyPath);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}
		else	{
			if(strncmp(szLine, szExeFile, nLen)==0)	{
				sscanf(szLine, "%s %s", szBuff, szExePath);
				break;
			}
		}
	}
	fclose(fIn);

	fIn = fopen(szExePath, "rb");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "The file %s for %s doesn't exist.\nQuit\n", szExePath, szExeFile);
		Quit_With_Error_Msg(ErrorMsg);
	}
	else	{
		fclose(fIn);
	}
}

void Setup_QM_Level(void)
{
	FILE *fIn;
	char *ReadLine, szLine[256], szSubStr_1[256], szSubStr_2[256], ErrorMsg[256];
	char szKey_QM_MEM[256], szKey_QM_NPROC[256], szQM_Level_Opt[256], szQM_Level_Dimer_Opt[256], szQM_Level_ESP[256], szQM_Level_E_Monomer[256], szQM_Level_1D_Scan[256], szQM_Level_Rotamer[256];
	int nStr, KeyCount;

	fIn = fopen(szQMParaFile, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file %s\nQuit\n", szQMParaFile);
		Quit_With_Error_Msg(ErrorMsg);
	}

	KeyCount = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		Replace_NewLine(szLine);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			nStr = Split_Into_Two_String(szLine, szSubStr_1, szSubStr_2);
			if(nStr == 2)	{
				if(strcmp(szSubStr_1,"QM_MEM")==0)	{
					strcpy(szKey_QM_MEM, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_NPROC")==0)	{
					strcpy(szKey_QM_NPROC, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_OPT")==0)	{
					strcpy(szQM_Level_Opt, szSubStr_2);
					KeyCount++;
					nCore_Per_Node = atoi(szKey_QM_NPROC);
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_DIMER_OPT")==0)	{
					strcpy(szQM_Level_Dimer_Opt, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_ESP")==0)	{
					strcpy(szQM_Level_ESP, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_E_MONOMER")==0)	{
					strcpy(szQM_Level_E_Monomer, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_1D_SCAN")==0)	{
					strcpy(szQM_Level_1D_Scan, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_ROTAMER")==0)	{
					strcpy(szQM_Level_Rotamer, szSubStr_2);
					KeyCount++;
				}
			}
		}
	}

	fclose(fIn);

	Determine_QM_Level(szQM_Level_1D_Scan);

	if(KeyCount != 8)	{	// !!!
		sprintf(ErrorMsg, "Setup_QM_Level> Error: incomplete entries in %s\nQuit\n", szQMParaFile);
		Quit_With_Error_Msg(ErrorMsg);
	}

	sprintf(szQM_Level_Rotamer_Cmd, "%%mem=%s%%nproc=%s%sMol rotamer scan\n\n%d,1\n",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_Rotamer, netcharge_mol);
//	printf("%s", szQM_Level_Rotamer_Cmd);
}


void Determine_QM_Level(char szQMLevel[])
{
	char ErrorMsg[256];

	if(FindString(szQMLevel, "HF") >= 0)	{
		QM_Level = QM_LEVEL_HF;
	}
	else if(FindString(szQMLevel, "MP2") >= 0)	{
		QM_Level = QM_LEVEL_MP2;
	}
	else	{
		sprintf(ErrorMsg, "Fail to determine the QM level.\n%s\nQuit\n", szQMLevel);
		Quit_With_Error_Msg(ErrorMsg);
	}
	return;
}

int Split_Into_Two_String(char szBuff[], char str1[], char str2[])
{
	int nLen, i, iBegin_First, Count, iBegin_Second, nLen_1, nLen_2;

	nLen = strlen(szBuff);
	str1[0] = 0;
	str2[0] = 0;

	for(i=0; i<nLen; i++)	{
		if( (szBuff[i] != ' ') && (szBuff[i] != '\t') ) 	{	// To find the first character of the first string
			break;
		}
	}

	iBegin_First = i;

	Count = 0;
	for(i=iBegin_First; i<nLen; i++)	{
		if( (szBuff[i] == ' ') || (szBuff[i] == '\t') ) 	{	// To find the last character of the first string
			break;
		}
		else	{
			str1[Count] = szBuff[i];
			Count++;
		}
	}
	str1[Count] = 0;
	nLen_1 = Count;

	for(; i<nLen; i++)	{
		if( (szBuff[i] != ' ') && (szBuff[i] != '\t')  && (szBuff[i] != 0x22) ) 	{	// To find the first character of the second string
			break;
		}
	}
	iBegin_Second = i;

	Count = 0;
	for(i=iBegin_Second; i<nLen; i++)	{
		if( (szBuff[i] == 0x0) || (szBuff[i] == 0x22) ) 	{	// To find the last character of the second string
//		if( (szBuff[i] == 0x0) || (szBuff[i] == 0x0D) || (szBuff[i] == 0x0A) || (szBuff[i] == 0x22) ) 	{	// To find the last character of the second string
			break;
		}
		else	{
			str2[Count] = szBuff[i];
			Count++;
		}
	}
	str2[Count] = 0;
	nLen_2 = Count;

	if(nLen_2 > 0)	{
		return 2;
	}
	else if(nLen_1 > 0)	{
		return 1;
	}
	else	{
		return 0;
	}
}

void Replace_NewLine(char szBuff[])
{
	int nLen, i;

	nLen =strlen(szBuff);

	for(i=1; i<nLen; i++)	{
		if( (szBuff[i-1]==0x5C) && (szBuff[i]==0x6E) )	{
			szBuff[i-1] = 0x20;
			szBuff[i] = 0x0A;
		}
	}
}


void Read_Tor_Para_1D_Fitting(void)
{
	FILE *fIn;
	char szName[256], ErrorMsg[256];
	int j, Idx_Phi, ReadItem;

	for(Idx_Phi=0; Idx_Phi<n_Phi; Idx_Phi++)	{
		sprintf(szName, "torsion-para-%d.dat", Idx_Phi+1);
		fIn = fopen(szName, "r");
		if(fIn == NULL)	{
			sprintf(ErrorMsg, "Fail to open file: %s.\nQuit\n", szName);
			Quit_With_Error_Msg(ErrorMsg);
		}

		for(j=0; j<(N_TOR_PARA-1); j++)	{
			ReadItem = fscanf(fIn, "%lf", &(Para_Tor_List[Idx_Phi][j]));
			if(ReadItem != 1)	{
				break;
			}
		}
		if(j!=(N_TOR_PARA-1))	{	// not a valid record
			fclose(fIn);
			sprintf(ErrorMsg, "Error in reading file: %s.\nQuit\n", szName);
			Quit_With_Error_Msg(ErrorMsg);
		}

		fclose(fIn);
	}

	Assign_Torsion_Parameters();
}

void Assign_Torsion_Parameters(void)
{
	int i_Phi, Idx;
	double *pPara_k_Dih, *pPara_phi;

	for(i_Phi=0; i_Phi<n_Phi; i_Phi++)	{
		Idx = IdxDihSelect[i_Phi];
		pPara_k_Dih = Mol.Para_k_Dih[Idx];
		pPara_phi = Mol.Para_phi[Idx];

		pPara_k_Dih[1] = Para_Tor_List[i_Phi][0];	// k1
		pPara_phi[1] = Para_Tor_List[i_Phi][1];		// sigma_1
		pPara_k_Dih[2] = Para_Tor_List[i_Phi][2];	// k2
		pPara_phi[2] = Para_Tor_List[i_Phi][3];		// sigma_2
		pPara_k_Dih[3] = Para_Tor_List[i_Phi][4];	// k3
		pPara_phi[3] = Para_Tor_List[i_Phi][5];		// sigma_3
		pPara_k_Dih[4] = Para_Tor_List[i_Phi][6];	// k4
		pPara_phi[4] = Para_Tor_List[i_Phi][7];		// sigma_4
		pPara_k_Dih[5] = 0.0;						// k5
		pPara_phi[5] = 0.0;							// sigma_5

		pPara_k_Dih[6] = Para_Tor_List[i_Phi][8];	// k6
		pPara_phi[6] = Para_Tor_List[i_Phi][9];		// sigma_6
	}
}


void ReadElementList(void)
{
	int n_Atom_In_Mol, i;
	FILE *fIn;
	char ErrorMsg[256];

	n_Atom_In_Mol = Mol.nAtom;

	fIn = fopen(szFileElem, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file %s\nQuit\n", szFileElem);
		Quit_With_Error_Msg(ErrorMsg);
	}

	for(i=0; i<n_Atom_In_Mol; i++)	{
		fscanf(fIn, "%s", szElemName[i]);
	}
	fclose(fIn);

}



void CWORKER::Start_Job(int CoreNum, int ID, char szGjf[], char szOut[])
{
	char szCmd[256], szEnv[256], szErrorMsg[256];
	FILE *fIn;

	nCore = CoreNum;
	Status = BUSY;
	
	JobID = ID;
	strcpy(szInput, szGjf);
	strcpy(szOutput, szOut);

	sprintf(szEnv, "%s/%d", szGAUSS_SCRDIR_Base, ID);
#ifdef WIN32
	if( SetEnvironmentVariable(szGAUSS_SCRDIR, szEnv) == 0) {
#else
	if (setenv(szGAUSS_SCRDIR, szEnv, 1) < 0) {
#endif
		sprintf(szErrorMsg, "Error setting env for %s with ID = %d.\n", szGAUSS_SCRDIR, ID);
		Quit_With_Error_Msg(szErrorMsg);
	}

	sprintf(szCmd, "mkdir %s", szEnv);
	system(szCmd);

	chdir(szDirList[ID]);

	fIn = fopen(szOutput, "r");
	if(fIn == NULL)	{
		sprintf(szCmd, "%s < %s > %s &", szExe_G09, szInput, szOutput);
		Update_Core_Number_Input_File();
		system(szCmd);
	}
	else	{
		fclose(fIn);
		Get_Job_Status();
	}
}

#define SIZE_CHECK	(512)
int CWORKER::Get_Job_Status(void)
{
	FILE *fIn;
	char szErrorMsg[256], szBuff[1024];
	int i, nLen;

	if(Status == IDLE)	{
		return Status;
	}

	fIn = fopen(szOutput, "rb");

	if(fIn == NULL)	{
		sprintf(szErrorMsg, "Error in open file %s for reading.\nQuit\n", szOutput);
		Quit_With_Error_Msg(szErrorMsg);
	}

	fseek(fIn, -SIZE_CHECK, SEEK_END);

	nLen = fread(szBuff, 1, SIZE_CHECK, fIn);

	for(i=0; i<nLen; i++)	{
		if(strncmp(szBuff+i, "Job cpu time:", 13) == 0)	{
			fclose(fIn);

			JobStatus[JobID] = Job_Done;
			Status = IDLE;
			nJobDone++;
			return Status;
		}
	}

	fclose(fIn);

	Status = BUSY;
	return Status;
}
#undef SIZE_CHECK

void CWORKER::Update_Core_Number_Input_File(void)
{
	FILE *fIn, *fOut;
	char szErrorMsg[256], szNewName[]="new.gjf", szLine[256], *ReadLine, szCmd[256];;

	fIn = fopen(szInput, "r");

	if(fIn == NULL)	{
		sprintf(szErrorMsg, "Fail to open input file %s\nQuit\n", szInput);
		Quit_With_Error_Msg(szErrorMsg);
	}

	fOut = fopen(szNewName, "w");

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			if(strncmp(szLine, "%nproc", 6)==0)	{
				sprintf(szLine, "%s%d\n", "%nproc=", nCore);
			}

			fprintf(fOut, "%s", szLine);
		}
	}


	fclose(fIn);
	fclose(fOut);

	sprintf(szCmd, "mv %s %s", szNewName, szInput);
	system(szCmd);
}

int Count_Finished_Job()
{
	int i, Count=0;

	for(i=0; i<nJob; i++)	{
		if(JobStatus[i] == Job_Done)	{
			Count++;
		}
	}

	return Count;
}

int Get_First_Available_Worker(void)	// to check the status for all workers
{
	int iFirst=-1;

	for(int i=0; i<nWorker; i++)	{
		if(pWorker[i].Get_Job_Status() == IDLE)	{
			if(iFirst < 0)	{
				iFirst = i;
			}
		}
	}

	return iFirst;
}

void RunAllJobs(void)
{
	int i, Idx_Worker, CoreNum, nJob_Submitted=0;
	int CoreNum_List[]={12, 6, 4, 3, 4, 2, 3, 3, 4, 4, 4, 2, 4, 4, 4, 3};	// for 12 cores
//	int CoreNum_List[]={24, 16, 10, 8, 6, 5, 4, 4, 3, 3, 5, 5, 4, 4, 4, 4};	// for 32 cores

	if(nJob > 16)	{
		CoreNum = 4;
	}
	else	{
		CoreNum = CoreNum_List[nJob-1];
	}
	nWorker = nCore_Per_Node / CoreNum;
	if(nWorker == 0)	{
		nWorker = 1;
	}
	pWorker = new CWORKER[nWorker];

	for(i=0; i<nJob; i++)	{
		JobStatus[i] = Job_ToRun;
	}

	for(i=0; i<nWorker; i++)	{
		pWorker[i].Status = IDLE;
	}

	nJob_Submitted = 0;
	while( nJobDone < nJob)	{
		if(Count_Finished_Job() == nJob)	{
			break;
		}
		Idx_Worker = Get_First_Available_Worker();
		if( (Idx_Worker >= 0) && (nJob_Submitted < nJob) )	{
			pWorker[Idx_Worker].Start_Job(CoreNum, nJob_Submitted, szInputList[nJob_Submitted], szOutputList[nJob_Submitted]);
			nJob_Submitted++;
			system("sleep 3");
		}
		else	{
			system("sleep 5");
		}
	}

	delete []pWorker;
}


