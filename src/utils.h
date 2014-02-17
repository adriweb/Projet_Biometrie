#ifndef __UTILS_H__
#define __UTILS_H__

#include "common.h"

/**
* \brief	Getter pour le wrapper C++
*/
char* getLatestSavedImageName();

int strEndsWith(const char *str, const char *suffix);

int array_min_idx(int* arr, int size);
int array_max_idx(int* arr, int size);

DonneesImageRGB* new_ImageRGB(uint colonnes, uint lignes);

void saveBMPwithCurrentName(DonneesImageRGB * image, const char* name);


#endif
