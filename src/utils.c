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

void saveBMPwithCurrentName(DonneesImageRGB * image, const char* name)
{
	secure_free(latestSavedImageName);

	uint len = strlen(nomFichier) + strlen(name) + 2;
	latestSavedImageName = (char*)calloc(len, sizeof(char));
	if (!latestSavedImageName) exit(-1);
	strcpy(latestSavedImageName, nomFichier);
	strcat(latestSavedImageName, "_");
	strcat(latestSavedImageName, name);

	ecrisBMPRGB_Dans(image, latestSavedImageName);
}
