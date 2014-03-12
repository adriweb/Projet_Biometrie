// Adrien Bertrand
// Biométrie - LBP
// v1.5 - 12/03/2014

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
#include <limits.h>
#include "BmpLib.h"

//#define IS_DEBUG

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


/**
* \brief	no-VS compatibility stuff ("secure" functions)
*/
#ifndef _MSC_VER
#define scanf_s				scanf
#define gets_s(a,b)			gets((a))
#define strcpy_s(a,b,c)		strncpy((a),(c),(b))
#define sprintf_s(...)		snprintf(__VA_ARGS__)
#define MAIN_NAME			main
#endif

#define NUMARGS(...)  (sizeof((int[]){__VA_ARGS__})/sizeof(int))

#define array_count(x) (sizeof(x)/sizeof(0[x]))

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

#define ABS(x)	(((x)<0) ? (-x) : (x))

/**
* \brief	free(NULL) shouldn't be an issue on decent compilers, but on others... Also, sets to NULL the freed pointer.
*/
#define secure_free(x)	do { if ((x)) { free((x)); (x) = NULL; } else { debugPrint("Warning ! Tried to free NULL (%s) at line %d (%s)\n", #x, __LINE__, __FUNCTION__); } } while(0)

/**
* \brief	Sachant que R=G=B pour les niveaux de gris, on utilise la fonction générale avec les 3 mêmes paramètres.
*/
#define sauveImageNG(img, img_ng)	creeImage(img, img_ng, img_ng, img_ng)




/**
* \brief	Simple rectangle défini par son sommet en haut à gauche, et ses dimensions
*/
typedef struct _rect_t {
	uint x, y, w, h;
} rect_t;
typedef rect_t face_rect_t;

typedef struct _point_t {
	uint x, y;
} point_t;




extern int img_w, img_h;
extern string nomFichier;
extern string latestSavedImageName;

#endif
