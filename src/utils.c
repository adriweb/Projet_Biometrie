#include "utils.h"

// Useful stuff :

char* latestSavedImageName = NULL;

char* getLatestSavedImageName() {
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

void do_saveBMPwithName(DonneesImageRGB * image, const char* name, const char* suffix)
{
	secure_free(latestSavedImageName);

	uint len = strlen(name) + strlen(suffix) + 2;
	latestSavedImageName = (char*)calloc(len, sizeof(char));
	if (!latestSavedImageName) exit(-1);
	strcpy(latestSavedImageName, name);
	strcat(latestSavedImageName, "_");
	strcat(latestSavedImageName, suffix);

	ecrisBMPRGB_Dans(image, latestSavedImageName);
}

void saveBMPwithCurrentName(DonneesImageRGB * image, const char* name)
{
	do_saveBMPwithName(image, nomFichier, name);
}

void createDirectory(const char* name)
{
#ifdef __linux__		
		mkdir(name, 0700);
#else
		_mkdir(name);
#endif
}

void changeDirectory(const char* name)
{
#ifdef __linux__		
		chdir(name);
#else
		_chdir(name);
#endif
}
