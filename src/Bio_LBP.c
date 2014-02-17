// Adrien Bertrand
// Biométrie - LBP
// v1.0 - 12/02/2014

#include "Bio_LBP.h"
#include "filtering.h"
#include "utils.h"


// Petites variables globales et pointeurs globaux

int img_w, img_h;
char* nomFichier = NULL;

int man_seuil = 100;
int taille_gaussien = 0;
int nbrCouleursReduc = 0;

DonneesImageRGB * image_orig;
DonneesImageRGB * image;
DonneesImageRGB * histo_img;
DonneesImageRGB * tmp_img;
u16 ** tmp_ng1;
u16 ** tmp_ng2;
u16 ** tmp_ng3;
u16 ** matrice_bleue;
u16 ** matrice_rouge;
u16 ** matrice_verte;
u16 ** image_ng;
uint * histo;
uint * histoCumule;



void cree3matrices(u16** mat_bleue, u16** mat_rouge, u16** mat_verte, DonneesImageRGB* image) {
	int i, j, k = 0;
	for (j = 0; j < image->hauteurImage; j++) {
		for (i = 0; i < image->largeurImage; i++) {
			mat_bleue[j][i] = (u16)image->donneesRGB[k++];
			mat_verte[j][i] = (u16)image->donneesRGB[k++];
			mat_rouge[j][i] = (u16)image->donneesRGB[k++];
		}
	}
}

void creeImage(DonneesImageRGB* image, u16** mat_bleue, u16** mat_rouge, u16** mat_verte) {
	int i, j, k = 0;
	for (j = 0; j < image->hauteurImage; j++) {
		for (i = 0; i < image->largeurImage; i++) {
			image->donneesRGB[k++] = (u8)mat_bleue[j][i];
			image->donneesRGB[k++] = (u8)mat_verte[j][i];
			image->donneesRGB[k++] = (u8)mat_rouge[j][i];
		}
	}
}

void negatifImage(u16** mat_bleue, u16** mat_rouge, u16** mat_verte) {
	int i, j;
#pragma omp parallel for
	for (j = 0; j < img_h; j++) {
		for (i = 0; i < img_w; i++) {
			mat_rouge[j][i] = 255 - mat_rouge[j][i];
			mat_verte[j][i] = 255 - mat_verte[j][i];
			mat_bleue[j][i] = 255 - mat_bleue[j][i];
		}
	}
}

u16** couleur2NG(u16** mat_bleue, u16** mat_rouge, u16** mat_verte, bool perceptive) {
	int i, j;
	double coeff_r = perceptive ? 0.2125 : 1,
		coeff_g = perceptive ? 0.7154 : 1,
		coeff_b = perceptive ? 0.0721 : 1;
	uint divi = perceptive ? 1 : 3;
	u16** imageNG = (u16**)malloc(img_h * sizeof(u16*));
	if (!imageNG) return NULL;
	for (j = 0; j < img_h; j++) {
		imageNG[j] = (u16*)calloc(img_w, sizeof(u16));
		if (!imageNG[j]) return NULL;
		for (i = 0; i < img_w; i++)
			imageNG[j][i] = (u16)((mat_rouge[j][i] * coeff_r + mat_verte[j][i] * coeff_g + mat_bleue[j][i] * coeff_b) / divi);
	}
	return imageNG;
}

void seuilleImageNG(u16** imageNG, uint seuil) {
	int i, j = 0;
#pragma omp parallel for
	for (j = 0; j < img_h; j++)
	for (i = 0; i < img_w; i++)
		imageNG[j][i] = (imageNG[j][i] > seuil) ? 255 : 0;
}

uint* histogramme(u16** imageNG) {
	histo = (uint*)calloc(GRAYLEVELS, sizeof(uint));
	if (!histo) return NULL;

	int i, j;
	for (j = 0; j < img_h; j++)
	for (i = 0; i < img_w; i++)
		histo[imageNG[j][i]]++;

	return histo;
}

uint* histogrammeCumule(u16** imageNG) {
	histoCumule = histogramme(imageNG);

	int i;
	for (i = 1; i < GRAYLEVELS; i++)
		histoCumule[i] += histoCumule[i - 1];

	return histoCumule;
}

void histo_egalisation(u16** imageNG) {
	histoCumule = histogrammeCumule(imageNG);

	int i, j = 0;
	float ratio = 255.0f / (float)(img_h*img_w);
#pragma omp parallel for
	for (j = 0; j < img_h; j++)
	for (i = 0; i < img_w; i++)
		imageNG[j][i] = (u16)((float)histoCumule[imageNG[j][i]] * ratio);

	secure_free(histoCumule);
}

DonneesImageRGB* imageHistogramme(uint* histo)
{
	int i, j;
	uint colonnes = GRAYLEVELS, lignes = 180;

	DonneesImageRGB* imageHisto = new_ImageRGB(colonnes, lignes);
	memset(imageHisto->donneesRGB, 255, colonnes * lignes * 3 * sizeof(unsigned char)); // fond blanc

	float ratio = histo[array_max_idx((int*)histo, GRAYLEVELS)] / (float)lignes;

	uint* histoNorm = (uint*)malloc(GRAYLEVELS*sizeof(uint));
	if (!histoNorm) return NULL;

#pragma omp parallel for
	for (i = 0; i < GRAYLEVELS; i++)
		histoNorm[i] = (uint)(histo[i] / ratio); // normalisation

#pragma omp parallel for
	for (i = 0; i < GRAYLEVELS; i++)
	for (j = 0; j < (int)(histoNorm[i]); j++)
		memset(&(imageHisto->donneesRGB[(3 * (i + j * colonnes))]), (j >= (int)(histoNorm[i-1]) || j >= (int)(histoNorm[i+1])) ? 0 : j % 235, 3);

	secure_free(histoNorm);

	return imageHisto;
}

int get_seuil_otsu(u16** src)
{
	int* histo = (int*)histogramme(src);

	int i, threshold = 0, taille = img_w * img_h;
	int wB = 0, wF = 0;
	float mB, mF, varBetween;
	float sum = 0, sumB = 0, varMax = 0;

	for (i = 0; i < GRAYLEVELS; i++)
		sum += (float)i * histo[i];

	for (i = 0; i < GRAYLEVELS; i++) {
		wB += histo[i];					// Weight Background
		if (wB == 0) continue;

		wF = taille - wB;				// Weight Foreground
		if (wF == 0) break;

		sumB += (float)(i * histo[i]);

		mB = sumB / wB;					// Mean Background
		mF = (sum - sumB) / wF;			// Mean Foreground

		// Calculate Between Class Variance
		varBetween = (float)wB * (float)wF * (mB - mF) * (mB - mF);

		// Check if new maximum found
		if (varBetween > varMax) {
			varMax = varBetween;
			threshold = i;
		}
	}

	secure_free(histo);

	return threshold;
}

void do_Seuil(int seuil)
{
#ifdef CONSOLE
	seuil = -1;
	while (seuil < 0 || seuil > 255) {
		printf("Seuil ? (0-255)\n");
		scanf_s("%d", &seuil);
		printf("\n");
	}
#endif
	if (seuil > 0 && seuil < 256) {
		seuilleImageNG(image_ng, seuil);
		sauveImageNG(image, image_ng);
		saveBMPwithCurrentName(image, "seuil.bmp");
	}
}

u16 ** paletteReduction(u16 ** src, int levelsAmount)
{
	levelsAmount--;
	int i, j, k;
	float interval = 256 / (float)levelsAmount;
	float half_interval = interval / 2;

	// palette init
	u16* levels = (u16*)malloc((levelsAmount + 1) * sizeof(u16));
	if (!levels) return NULL;
#pragma omp parallel for
	for (i = 0; i < levelsAmount + 1; i++)
		levels[i] = (u16)(i * interval);
	levels[levelsAmount] = 255;

	//	for (i = 0; i < levelsAmount + 1; i++)
	//		printf("levels[%d] = %d\n", i, levels[i]);

	// reduced img init + processing
	u16** reduced = (u16**)malloc(img_h * sizeof(u16*));
	if (!reduced) return NULL;
	for (j = 0; j < img_h; j++) {
		reduced[j] = (u16*)calloc(img_w, sizeof(u16)); // will make everything mm_default (0)
		if (!reduced[j]) return NULL;
#pragma omp parallel for
		for (i = 0; i < img_w; i++)
		for (k = 0; k < levelsAmount + 1; k++)
		if ((src[j][i] >= levels[k] - (k>0 ? half_interval + 0.5 : 0)) && (src[j][i] < levels[k] + half_interval + 0.5))
			reduced[j][i] = levels[k];
	}

	secure_free(levels);

	return reduced;
}

void do_PaletteReduction(int level)
{
#ifdef CONSOLE
	level = -1;
	while (level < 2 || level > 254) {
		printf("Combien de niveaux ? (2-254)\n");
		scanf_s("%d", &level);
		printf("\n");
	}
#endif
	if (level > 1 && level < 255) {
		u16** tmp = paletteReduction(image_ng, level);
		sauveImageNG(image, tmp);
		saveBMPwithCurrentName(image, "reduction_palette.bmp");
		int i;
		for (i = 0; i < img_h; i++) secure_free(tmp[i]);
		secure_free(tmp);
	}
}


u16** get_subimage(u16** src, uint src_w, uint src_h, uint x, uint y, uint w, uint h)
{
	uint i, j;
	u16** sub;

	sub = (u16**)malloc(h * sizeof(u16*));
	if (!sub) return NULL;

	for (j = y; j < y+h; j++) {
		sub[j-y] = (u16*)calloc(w, sizeof(u16));
		if (!sub[j-y]) return NULL;
		for (i = x; i < x + w; i++)
			sub[j-y][i-x] = (j < src_h && i < src_w) ? src[j][i] : 0;
	}
	return sub;
}

// retourne le nombre de sous-images extraites
uint extract_subimages_and_save(u16** image_ng, uint img_w, uint img_h)
{
	uint i, j, tmp, counter = 0;
	uint nbr_sub_x = 5;
	uint sub_size = (uint)roundf(((float)img_w / (float)nbr_sub_x)); // == width == height (square)
	uint loop_step = sub_size >> 1; // génération d'une sous-image à chaque taille/2.
	uint nbr_sub_y = (uint)roundf(((float)img_h / (float)sub_size));
/*
	img_gris** subimages_mat; // actual type : u16**** (!)

	subimages_mat = (img_gris**)malloc((nbr_sub_y-1) * sizeof(img_gris*)); // -1 because last is incomplete (thus useless)
	if (!subimages_mat) return 0;
	for (j = 0; j < nbr_sub_y; j++) {
		subimages_mat[j] = (img_gris*)calloc((nbr_sub_x - 1), sizeof(img_gris));
		if (!subimages_mat[j]) return 0;
	}
*/
	u16** sub = NULL;
	char* sub_filename = NULL;
	sub_filename = calloc(strlen(nomFichier) + 25, sizeof(char)); // +15 => "_sub_x_y_w.bmp" + extra safety
	DonneesImageRGB* dest_imgRGB = new_ImageRGB(sub_size, sub_size);

	for (j = 0; j < img_h - loop_step; j += loop_step) {
		for (i = 0; i < img_w - loop_step; i += loop_step) {
			error("processing sub : %d ; %d\n", j, i);
			sub = get_subimage(image_ng, img_w, img_h, i, j, sub_size, sub_size);
			sprintf(sub_filename, "%s_sub_%d_%d_%d.bmp", nomFichier, j, i, sub_size);
			sauveImageNG(dest_imgRGB, sub);
			ecrisBMPRGB_Dans(dest_imgRGB, sub_filename);

			for (tmp = 0; tmp < sub_size; tmp++)
				secure_free(sub[tmp]);
			secure_free(sub);
			
			counter++;
		}
	}

	secure_free(sub_filename);

	return counter;
}

void choixAction(int choix)
{
	int i;

	static bool isDoingAll = false;

	bool end = (choix == 0);

	while (!end) {

		cree3matrices(matrice_bleue, matrice_rouge, matrice_verte, image_orig);
		image_ng = couleur2NG(matrice_bleue, matrice_rouge, matrice_verte, true);

#ifdef CONSOLE

		if (!isDoingAll) {

			printf("******************************\n");
			printf("**** Bertrand - Debournoux ***\n");
			printf("******* Biometrie - LBP ******\n");
			printf("*****  v1.0 - 12/02/2014 *****\n");
			printf("******************************\n\n");
			printf("Image en cours : %s\n\n", nomFichier);
			printf("* 1) LBP avec mediane \n");
			printf("* 10) LBP sans mediane \n");
			printf("* 2) Détection de visage(s) \n");
			printf("* 3) Mediane \n");
			printf("* 4) Histogramme \n");
			printf("* 5) Extractions sous-images \n");
			printf("* 0) Quitter \n\n");
			printf("******************************\n");
			printf("Choix ? \n");
			scanf_s("%d", &choix);
			printf("\n");

		}
#endif

		switch (choix) {
		
		case 1:
			tmp_ng1 = apply_filter(image_ng, filters[flt_Median]);
			tmp_ng2 = apply_filter(tmp_ng1, filters[flt_LBP]);

			sauveImageNG(image, tmp_ng2);
			saveBMPwithCurrentName(image, "lbp-with-median.bmp");
			break;
		case 10:
			tmp_ng1 = apply_filter(image_ng, filters[flt_LBP]);

			sauveImageNG(image, tmp_ng1);
			saveBMPwithCurrentName(image, "lbp.bmp");
			break;
		case 2:
			
			error("Unimplemented !\n");
			break;
		case 3: // median test
			tmp_ng1 = apply_filter(image_ng, filters[flt_Median]);

			sauveImageNG(image, tmp_ng1);
			saveBMPwithCurrentName(image, "median.bmp");
			break;
		case 4:
			cree3matrices(matrice_bleue, matrice_rouge, matrice_verte, image_orig);
			image_ng = couleur2NG(matrice_bleue, matrice_rouge, matrice_verte, false);
			histo = histogramme(image_ng);
			histo_img = imageHistogramme(histo);
			saveBMPwithCurrentName(histo_img, "histogramme.bmp");
			secure_free(histo);
			break;
		case 5:
			cree3matrices(matrice_bleue, matrice_rouge, matrice_verte, image_orig);
			image_ng = couleur2NG(matrice_bleue, matrice_rouge, matrice_verte, false);
			int counter = extract_subimages_and_save(image_ng, img_w, img_h);
			error("counter : %d\n", counter);
			break;
		case 0:
			end = true;
			break;
		default:
			printf("Mauvais choix !\n\n");
			break;
		}

#pragma omp parallel for
		for (i = 0; i < img_h; i++) {
			if (image_ng) secure_free(image_ng[i]);
			if (tmp_ng1) secure_free(tmp_ng1[i]);
			if (tmp_ng2) secure_free(tmp_ng2[i]);
			if (tmp_ng3) secure_free(tmp_ng3[i]);
		}
		secure_free(image_ng);
		secure_free(tmp_ng1);
		secure_free(tmp_ng2);
		secure_free(tmp_ng3);
		libereDonneesImageRGB(&histo_img);
		libereDonneesImageRGB(&tmp_img);


#ifdef CONSOLE
		if (!isDoingAll) {
			if (!end) system("pause");
			system("cls");
		}
#else
		end = true;
#endif

		if (isDoingAll) break;
	}
}


void initData(int argc, char *argv[])
{
	nomFichier = (char*)calloc(350, sizeof(char));
	if (!nomFichier) return;

#ifndef CONSOLE
	error("loaded : %s\n", argv[1]);
#endif
	strcpy_s(nomFichier, (argc > 1) ? strlen(argv[1]) + 1 : 14, (argc > 1) ? argv[1] : "D:\\image.bmp");

	if (!strEndsWith(nomFichier, ".bmp")) {
		printf("Image name doesn't contain the extension, adding it...\n");
		strcat(nomFichier, ".bmp");
	}

	if (!(image = lisBMPRGB(nomFichier))) {
		error("Erreur de lecture, fermeture... \n");
		secure_free(nomFichier);
		exit(-1);
	}
	if (!(image_orig = lisBMPRGB(nomFichier))) {
		error("Erreur de mémoire \n");
		secure_free(nomFichier);
		exit(-1);
	}

	filters = createFilters();

	int i;

	img_w = image->largeurImage;
	img_h = image->hauteurImage;

	matrice_rouge = (u16**)malloc(img_h * sizeof(u16*));
	matrice_verte = (u16**)malloc(img_h * sizeof(u16*));
	matrice_bleue = (u16**)malloc(img_h * sizeof(u16*));
	if (!(matrice_rouge && matrice_verte && matrice_bleue)) exit(-1);

	for (i = 0; i < img_h; i++) {
		matrice_rouge[i] = (u16*)malloc(img_w * sizeof(u16));
		matrice_verte[i] = (u16*)malloc(img_w * sizeof(u16));
		matrice_bleue[i] = (u16*)malloc(img_w * sizeof(u16));
	}
	if (!(matrice_rouge[0] && matrice_verte[0] && matrice_bleue[0])) exit(-1);
}

void freeStuff()
{
	int i;
	for (i = 0; i < NBR_FILTRES; i++)
		freeFilter(filters[i]);
	secure_free(filters);

	secure_free(nomFichier);

#pragma omp parallel for
	for (i = 0; i < img_h; i++) {
		if (matrice_bleue) secure_free(matrice_bleue[i]);
		if (matrice_rouge) secure_free(matrice_rouge[i]);
		if (matrice_verte) secure_free(matrice_verte[i]);
	}
	secure_free(matrice_bleue);
	secure_free(matrice_rouge);
	secure_free(matrice_verte);
	libereDonneesImageRGB(&image);
	libereDonneesImageRGB(&image_orig);
}

int MAIN_NAME(int argc, char *argv[])
{
#ifdef CONSOLE
	setlocale(LC_ALL, ""); // support unicode

	initData(argc, argv);

	choixAction(-1);

	freeStuff();
#else
	error("Why calling the console launcher from the GUI ... ?\n");
#endif

	return 0;
}
