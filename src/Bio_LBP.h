// Adrien Bertrand
// Biométrie - LBP
// v1.0 - 12/02/2014

#ifndef __BIO_LBP_H__
#define __BIO_LBP_H__

#include "common.h"


/*** Globals ***/

extern DonneesImageRGB * image_orig;
extern DonneesImageRGB * image;
extern DonneesImageRGB * histo_img;
extern DonneesImageRGB * tmp_img;
extern u16 ** tmp_ng1;
extern u16 ** tmp_ng2;
extern u16 ** tmp_ng3;
extern u16 ** matrice_bleue;
extern u16 ** matrice_rouge;
extern u16 ** matrice_verte;
extern u16 ** image_ng;
extern uint * histo;
extern uint * histoCumule;

extern int man_seuil;
extern int taille_gaussien;
extern int nbrCouleursReduc;

/**
* \brief	Simple rectangle défini par son sommet en haut à gauche, et ses dimensions
*/
typedef struct _rect_t {
	uint x, y, w, h;
} rect_t;

/*** Prototypes ***/

/**
* \brief	Sépare les composantes RGB d'une image en 3 matrices.
* \details	Remplit les 3 matrices mat_bleue, mat_rouge, et mat_verte avec la valeur des composantes RGB de l'image source
* \param	mat_bleue	Pointeur vers la matrice 2D de la composante bleue (destination)
* \param	mat_rouge	Pointeur vers la matrice 2D de la composante rouge (destination)
* \param	mat_verte	Pointeur vers la matrice 2D de la composante verte (destination)
* \param	image		L'image source
*/
void	cree3matrices(u16** mat_bleue, u16** mat_rouge, u16** mat_verte, DonneesImageRGB* image);

/**
* \brief	Crée dans 'image', l'image formée par les 3 composantes couleurs RGB
* \param	image		Pointeur vers l'image de destination
* \param	mat_bleue	Pointeur vers la matrice de la composante bleue de l'image (source)
* \param	mat_rouge	Pointeur vers la matrice de la composante rouge de l'image (source)
* \param	mat_verte	Pointeur vers la matrice de la composante verte de l'image (source)
*/
void	creeImage(DonneesImageRGB* image, u16** mat_bleue, u16** mat_rouge, u16** mat_verte);

/**
* \brief	Inverse in-place les couleurs à partir des matrices des composantes RGB.
* \details	Remplace les valeurs de chaque composantes couleurs en paramètre par leur inverse (255-x)
* \param	mat_bleue	Pointeur vers la matrice de la composante bleue de l'image (source)
* \param	mat_rouge	Pointeur vers la matrice de la composante rouge de l'image (source)
* \param	mat_verte	Pointeur vers la matrice de la composante verte de l'image (source)
*/
void	negatifImage(u16** mat_bleue, u16** mat_rouge, u16** mat_verte);

/**
* \brief	Crée une image en niveau-de-gris, à partir des matrices de composantes RGB de l'image couleur source
* \param	mat_bleue	Pointeur vers la matrice de la composante bleue de l'image (source)
* \param	mat_rouge	Pointeur vers la matrice de la composante rouge de l'image (source)
* \param	mat_verte	Pointeur vers la matrice de la composante verte de l'image (source)
* \param	perceptive	Indique si la fonction doit produire un niveau de gris être perceptif ou mathématique.
* \return	Matrice 2D représentant l'image en niveaux de gris
*/
u16**	couleur2NG(u16** mat_bleue, u16** mat_rouge, u16** mat_verte, bool perceptive);

/**
* \brief	Réalise le seuillage d'une image en niveau de gris (matricielle)
* \param	imageNG		Pointeur vers l'image en niveau de gris (matricielle)
* \param	seuil		Le niveau (seuil), nombre compris entre 0 et 255
*/
void	seuilleImageNG(u16** imageNG, uint seuil);

/**
* \brief	Construit l'histogramme (sous forme de tableau d'entier) à partir de l'image en niveau de gris (matricielle)
* \param	imageNG		Pointeur vers l'image en niveau de gris (matricielle)
* \return	Histogramme sous forme de tableau d'entiers
*/
uint*	histogramme(u16** imageNG);

/**
* \brief	Renvoie l'indice de l'élément du tableau d'entiers dont la valeur est minimale
* \param    arr			le tableau d'entier
* \param    size		la taille du tableau
* \return   L'indice (entier positif) de l'élément min. du tableau
*/
int		array_min_idx(int* arr, int size);

/**
* \brief	Renvoie l'indice de l'élément du tableau d'entiers dont la valeur est maximale
* \param    arr			le tableau d'entier
* \param    size		la taille du tableau
* \return   L'indice (entier positif) de l'élément max. du tableau
*/
int		array_max_idx(int* arr, int size);

/**
* \brief	Crée une représentation graphique (image) de l'histogramme
* \param    histo		Histogramme sous forme de tableau d'entiers
* \return   Pointeur vers l'image représentant l'histogramme.
*/
DonneesImageRGB* imageHistogramme(uint* histo);

/**
* \brief	Teste si les couleurs RGB en paramètres sont égales.
* \details	Détermination automatique du seuil par utilisation de la méthode d'Otsu ("maximization of inter-class variance"). Algorithme basé sur http://www.labbookpages.co.uk/software/imgProc/otsuThreshold.html
* \param    src			L'image
* \return   La valeur (entier [0-255]) de seuil optimal.
*/
int		get_seuil_otsu(u16** src);

void do_Seuil(int seuil);

void do_PaletteReduction(int level);

u16** get_subimage(u16** src, uint src_w, uint src_h, uint x, uint y, uint w, uint h);

uint extract_subimages_and_save(u16** image_ng, uint img_w, uint img_h);

/**
* \brief	Initialise les données de base
*/
void initData(int argc, char *argv[]);

/**
* \brief	Relatif au menu
*/
void choixAction(int choix);

/**
* \brief	Libère la mémoire utilisée par les données allouées
*/
void freeStuff();

/**
* \brief	main au nom variable selon la compilation en console ou en GUI
*/
int MAIN_NAME(int argc, char *argv[]);


#endif
