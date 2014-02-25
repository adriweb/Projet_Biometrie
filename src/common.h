// Adrien Bertrand
// Biométrie - LBP
// v1.20 - 25/02/2014

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

#define IS_DEBUG	1

#define CONSOLE
#ifdef CONSOLE
#define MAIN_NAME	main
#else
#define MAIN_NAME	main_console
#endif

typedef unsigned char		u8;
typedef unsigned short int	u16;
typedef unsigned int		uint;

typedef u8		uchar;
typedef u16**	img_gris;

typedef char*	string; // for the lulz

// just in case (Mac OS X doesn't have those by default, on QT Creator at least)
#ifndef INT_MIN
#define INT_MIN     (-2147483647 - 1)   /* minimum (signed) int value */
#define INT_MAX       2147483647        /* maximum (signed) int value */
#endif

// float.h
#define FLT_EPSILON     1.192092896e-07F

/**
* \brief	no-VS compatibility stuff ("secure" functions)
*/
#ifndef _MSC_VER
#define scanf_s				scanf
#define gets_s(a,b)			gets((a))
#define strcpy_s(a,b,c)		strncpy((a),(c),(b))
#define sprintf_s(a,b,c)	snprintf((a),(b),(c))
#define MAIN_NAME			main
#endif

#define NUMARGS(...)  (sizeof((int[]){__VA_ARGS__})/sizeof(int))

// from chromium
#define array_count(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))

#define GRAYLEVELS	256

/**
* \brief	fprintf vers le flux d'erreur + flush (normalement inutile car non bufferisé, mais bon...)
*/
#define error(...) do { fprintf(stderr, __VA_ARGS__); fflush(stderr); } while(0)

#ifdef IS_DEBUG
#define debugPrint(...) do { fprintf(stdout, __VA_ARGS__); fflush(stdout); } while(0)
#else
#define debugPrint(...)
#endif

/**
* \brief	free(NULL) shouldn't be an issue on decent compilers, but on others... Also, sets to NULL the freed pointer.
*/
#define secure_free(x)	do { if ((x)) { free((x)); (x) = NULL; } else { error("Trying to free NULL (%s) at line %d (%s)\n", #x, __LINE__, __FUNCTION__); } } while(0)

/**
* \brief	Sachant que R=G=B pour les niveaux de gris, on utilise la fonction générale avec les 3 mêmes paramètres.
*/
#define sauveImageNG(img, img_ng)	creeImage(img, img_ng, img_ng, img_ng)


extern int img_w, img_h;
extern string nomFichier;
extern string latestSavedImageName;

#endif
