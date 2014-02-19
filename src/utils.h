#ifndef __UTILS_H__
#define __UTILS_H__

#include "common.h"
#include <sys/types.h>

#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

/**
* \brief	Getter pour le wrapper C++
*/
char* getLatestSavedImageName();

int strEndsWith(const char *str, const char *suffix);

int array_min_idx(int* arr, int size);
int array_max_idx(int* arr, int size);

DonneesImageRGB* new_ImageRGB(uint colonnes, uint lignes);

void do_saveBMPwithName(DonneesImageRGB * image, const char* name, const char* suffix);
void saveBMPwithCurrentName(DonneesImageRGB * image, const char* name);

void createDirectory(const char* name);
void changeDirectory(const char* name);

int readHistoFromFile(uint* histo, const char* filename);
int writeHistoToFile(const char* filename, uint* histo);

#endif
