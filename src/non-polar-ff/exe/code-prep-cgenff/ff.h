#define N_LEN_CHEM_NAME	(16)
#define N_PARA_BOND		(2)
#define N_PARA_ANGLE	(4)	//k_theta, theta_0, k_UB, S_0
#define N_PARA_DIHEDRAL	(3)
#define N_PARA_IMPRDIHEDRAL	(3)
#define N_PARA_LJ	(6)	// non-1-4 and 1-4 parameters
#define N_PARA_NBFIX	(2)	// E_Min and r_Min




typedef struct	{
	char Chem[2][N_LEN_CHEM_NAME];
	double para[N_PARA_BOND];
}BOND_REC;

typedef struct	{
	char Chem[3][N_LEN_CHEM_NAME];
	double para[N_PARA_ANGLE];
}ANGLE_REC;

typedef struct	{
	char Chem[4][N_LEN_CHEM_NAME];
	double para[N_PARA_DIHEDRAL];
}DIHEDRAL_REC;

typedef struct	{
	char Chem[4][N_LEN_CHEM_NAME];
	double para[N_PARA_IMPRDIHEDRAL];
}IMPRODIHEDRAL_REC;

typedef struct	{
	char Chem[N_LEN_CHEM_NAME];
	double para[N_PARA_LJ];
}LJ_REC;

typedef struct	{
	char Chem_1[N_LEN_CHEM_NAME];
	char Chem_2[N_LEN_CHEM_NAME];
	double para[N_PARA_NBFIX];
}NBFix_REC;

#define N_REC_MAX	(8192)

class CForceField
{
public:
	int n_Rec_NBFix, n_Rec_LJ;
	NBFix_REC NBFix_Rec[N_REC_MAX];
	LJ_REC LJ_Rec[N_REC_MAX];

	int Quit_on_Error;
	int n_Rec_Bond, n_Rec_Angle, n_Rec_Dihedral, n_Rec_ImproDihedral;
	int Active_Bond[N_REC_MAX], Active_Angle[N_REC_MAX], Active_Dihedral[N_REC_MAX], Active_ImproDihedral[N_REC_MAX], Active_LJ[N_REC_MAX], Active_NBFix[N_REC_MAX];

	BOND_REC Bond_Rec[N_REC_MAX];
	ANGLE_REC Angle_Rec[N_REC_MAX];
	DIHEDRAL_REC Dihedral_Rec[N_REC_MAX];
	IMPRODIHEDRAL_REC ImproDihedral_Rec[N_REC_MAX];


public:

	CForceField();
	~CForceField();
	void ReadForceField(char szNamePara[]);
	void ReadUpdatedCMap(void);
	void GetPara_Bond(char szChemName[][N_LEN_CHEM_NAME], double Para[]);
	void GetPara_Angle(char szChemName[][N_LEN_CHEM_NAME], double Para[]);
	void GetPara_Dihedral(char szChemName[][N_LEN_CHEM_NAME], double Para[]);
	void GetPara_ImproDIhedral(char szChemName[][N_LEN_CHEM_NAME], double Para[]);
	int  GetPara_LJ(char szChemName[][N_LEN_CHEM_NAME], double Para[]);	// return the index in LJ record list
	int  QueryNBFix(char szChem_1[], char szChem_2[]);
};
