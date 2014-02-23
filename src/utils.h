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
string getLatestSavedImageName();

int strEndsWith(const char *str, const char *suffix);

int array_min_idx(int* arr, int size);
int array_max_idx(int* arr, int size);

DonneesImageRGB* new_ImageRGB(uint colonnes, uint lignes);

void do_saveBMPwithName(DonneesImageRGB * image, const string name, const string suffix);
void saveBMPwithCurrentName(DonneesImageRGB * image, const string name);

void createDirectory(const string name);
void changeDirectory(const string name);

uint* readHistoFromFile(const string filename);
int writeHistoToFile(const string filename, uint* histo);

#endif
