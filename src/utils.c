#include "utils.h"

// Useful stuff :

string latestSavedImageName = NULL;

string getLatestSavedImageName(void) {
	return latestSavedImageName;
}

int strEndsWith(const char *str, const char *suffix) {
	if (!str || !suffix)
		return 0;
	int lenstr = strlen(str);
	int lensuffix = strlen(suffix);
	if (lensuffix > lenstr)
		return 0;
	return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

int array_min_idx(int* arr, int size) {
	int idx = 0, i;
	for (i = 0; i < size; i++)
	if (arr[i] < arr[idx]) idx = i;
	return idx;
}

int array_max_idx(int* arr, int size) {
	int idx = 0, i;
	for (i = 0; i < size; i++)
	if (arr[i] > arr[idx]) idx = i;
	return idx;
}

u16** new_u16_mat(uint colonnes, uint lignes)
{
	uint i;
	u16** mat = (u16**)malloc(lignes * sizeof(u16*));
	if (!mat) return NULL;

	for (i = 0; i < lignes; i++) {
		mat[i] = (u16*)malloc(colonnes * sizeof(u16));
		if (!mat[i]) return NULL;
	}

	return mat;
}

void free_u16_mat(u16** mat, uint lignes)
{
	uint i;
	for (i = 0; i < lignes; i++)
		if (mat) secure_free(mat[i]);
	secure_free(mat);
}

DonneesImageRGB* new_ImageRGB(uint colonnes, uint lignes)
{
	DonneesImageRGB* imageRGB = (DonneesImageRGB*)malloc(sizeof(DonneesImageRGB));
	if (!imageRGB) return NULL;
	imageRGB->donneesRGB = (unsigned char*)malloc(colonnes * lignes * 3 * sizeof(unsigned char));
	if (!imageRGB->donneesRGB) return NULL;
	imageRGB->hauteurImage = lignes;
	imageRGB->largeurImage = colonnes;

	return imageRGB;
}

DonneesImageRGB* clone_imageRGB(DonneesImageRGB * image_orig)
{
	DonneesImageRGB* img = NULL;
	img = (DonneesImageRGB*)malloc(sizeof(DonneesImageRGB));
	if (!img) return NULL;
	img->donneesRGB = (uchar*)calloc(image_orig->hauteurImage * image_orig->largeurImage * 3, sizeof(uchar));
	if (!img->donneesRGB) return NULL;
	memcpy(img->donneesRGB, image_orig->donneesRGB, image_orig->hauteurImage * image_orig->largeurImage * 3);
	img->hauteurImage = image_orig->hauteurImage;
	img->largeurImage = image_orig->largeurImage;

	return img;
}

void do_saveBMPwithName(DonneesImageRGB * image, const string name, const string suffix)
{
	if (latestSavedImageName) secure_free(latestSavedImageName);

	uint len = strlen(name) + strlen(suffix) + 2;
	latestSavedImageName = (string)calloc(len, sizeof(char));
	if (!latestSavedImageName) exit(-1);
	strcpy(latestSavedImageName, name);
	strcat(latestSavedImageName, "_");
	strcat(latestSavedImageName, suffix);

	ecrisBMPRGB_Dans(image, latestSavedImageName);
}

void saveBMPwithCurrentName(DonneesImageRGB * image, const string name)
{
	do_saveBMPwithName(image, nomFichier, name);
}

void createDirectory(const string name)
{
#ifdef __linux__		
	mkdir(name, 0700);
#else
	(void)_mkdir(name);
#endif
}

char* getCurrentDirectory(void)
{
	const int size = 1024;
	string cwd = (string)calloc(size, sizeof(char));
#ifdef __linux__
	if (getcwd(cwd, size))
#else
	if (_getcwd(cwd, size))
#endif
		return cwd;
	else
		perror("getcwd() error");
	return "(error)";
}

void changeDirectory(const string name)
{
#ifdef __linux__		
	chdir(name);
#else
	(void)_chdir(name);
#endif
	char* cwd = getCurrentDirectory();
	//debugPrint("changed directory : %s\n", cwd);
	secure_free(cwd);
}

uint* readHistoFromFile(const string filename)
{
	FILE* fp;

	uint* histo = NULL;

	if (!(fp = fopen(filename, "rb"))) {
		error("Cannot open file for reading...\n");
		return NULL;
	} else {
		histo = (uint*)calloc(GRAYLEVELS, sizeof(uint));
		if (!histo) { fclose(fp);  return NULL; }

		fread(histo, sizeof(uint), GRAYLEVELS, fp);

		fclose(fp);

		return histo;
	}
}

int writeHistoToFile(const string filename, uint* histo)
{
	FILE* fp;
	int retVal = 0;

	if (!(fp = fopen(filename, "wb"))) {
		error("Cannot open file for writing...\n");
		retVal = -1;
	} else {
		fwrite(histo, sizeof(uint), GRAYLEVELS, fp);

		fclose(fp);
	}
	
	return retVal;
}
