// Adrien Bertrand
// Biométrie - LBP
// v1.4 - 07/03/2014

#ifndef __FILTERING_H__
#define __FILTERING_H__

#include "common.h"


#define NBR_FILTRES	2

// typedef void *(*void_func_ptr)(void);
// void_func_ptr filter_methods[NBR_FILTRES];

/**
* \brief	Liste des filtres possibles
*/
typedef enum _filters_types {
	flt_LBP,
	flt_Median
} filters_types;

/**
* \brief	Méthodes disponibles pour l'implémentation des filtres
*/
typedef enum _filter_method_t {
	flt_m_LBP,
	flt_m_Median
} filter_method_t;

/**
* \brief	Propriétés d'un filtre.
*/
typedef struct _filter_t {
	filter_method_t method;
	int** mask;
	uint size;
	uint div;
	bool needNormalization;
} filter_t;

extern filter_t** filters;

void freeFilter(filter_t* flt);

/**
* \brief	Applique un filtre sur une image en niveau de gris (matricielle)
* \param    src			l'image (matricielle) en niveau de gris
* \param    filter		2ème couleur
* \return   Pointeur vers la nouvelle image (matricielle) en niveau de gris, traitée (avec le filtre)
*/
u16** apply_filter(u16** src, filter_t* filter);
u16** do_apply_filter(u16** src, filter_t* filter, int imgw, int imgh);

void mask_copy_3(int** dest, int src[3][3]);
void mask_copy_5(int** dest, int src[5][5]);
void matrix_copy(int** dest, int** src, uint cols, uint rows);

void do_reNormalize(int** imageNG, int imgw, int imgh);

/**
* \brief	Crée le tableau des filtres de base.
* \details	Crée le tableau de pointeur de filtres définis dans @ref _filters_types
* \return   Pointeur vers le tableau
*/
filter_t** createFilters(void);
#endif
