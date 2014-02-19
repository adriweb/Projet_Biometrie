// Adrien Bertrand
// Biométrie - LBP
// v1.15 - 19/02/2014

#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <locale.h>
#include <time.h>
#include <errno.h>

#include "BmpLib.h"


#define CONSOLE
#ifdef CONSOLE
#define MAIN_NAME main
#else
#define MAIN_NAME main_console
#endif

typedef unsigned char		u8;
typedef unsigned short int	u16;
typedef unsigned int		uint;

typedef u16**				img_gris;

// just in case (Mac OS X doesn't have those by default, on QT Creator at least)
#ifndef INT_MIN
#define INT_MIN     (-2147483647 - 1)   /* minimum (signed) int value */
#define INT_MAX       2147483647        /* maximum (signed) int value */
#endif

/**
* \brief	no-VS compatibility stuff ("secure" functions)
*/
#ifndef _MSC_VER
#define scanf_s			scanf
#define gets_s(a,b)		gets(a)
#define strcpy_s(a,b,c)	strncpy(a,c,b)
#define MAIN_NAME		main
#endif

#define NUMARGS(...)  (sizeof((int[]){__VA_ARGS__})/sizeof(int))

#define GRAYLEVELS	256

/**
* \brief	free(NULL) shouldn't be an issue on decent compilers, but on others... Also, sets to NULL the freed pointer.
*/
#define secure_free(x)	do { if (x) { free(x); x = NULL; } } while(0)

/**
* \brief	fprintf vers le flux d'erreur + flush
*/
#define error(...) do { fprintf(stderr, __VA_ARGS__); fflush(stderr); } while(0)

/**
* \brief	Sachant que R=G=B pour les niveaux de gris, on utilise la fonction générale avec les 3 mêmes paramètres.
*/
#define sauveImageNG(img, img_ng)	creeImage(img, img_ng, img_ng, img_ng)


extern int img_w, img_h;
extern char* nomFichier;
extern char* latestSavedImageName;

#endif
