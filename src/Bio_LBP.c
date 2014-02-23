// Adrien Bertrand
// Biométrie - LBP
// v1.17 - 22/02/2014

#include "Bio_LBP.h"
#include "filtering.h"
#include "utils.h"


/*  --------- Petites variables globales et pointeurs globaux ---------  */

int img_w, img_h;
string nomFichier = NULL;

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

histo_model_t* histo_db;
uint histo_db_size = 0;


/*  --------- Fonctions ---------  */

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

uint* do_histogramme(u16** imageNG, uint w, uint h)
{
	uint i, j;

	histo = (uint*)calloc(GRAYLEVELS, sizeof(uint));
	if (!histo) return NULL;

	for (j = 0; j < h; j++)
	for (i = 0; i < w; i++)
		histo[imageNG[j][i]]++;

	return histo;
}

uint* histogramme(u16** imageNG) {
	return do_histogramme(imageNG, img_w, img_h);
}

uint* histogrammeCumule(u16** imageNG) {
	histoCumule = histogramme(imageNG);
	if (!histoCumule) return NULL;

	int i;
	for (i = 1; i < GRAYLEVELS; i++)
		histoCumule[i] += histoCumule[i - 1];

	return histoCumule;
}

void histo_egalisation(u16** imageNG) {
	histoCumule = histogrammeCumule(imageNG);
	if (!histoCumule) return;

	int i, j = 0;
	float ratio = 255.0f / (float)(img_h*img_w);

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

// returns the number of elements in DB
uint make_histo_db()
{
	uint i, sample, idx = 0;
	
	const unsigned int nbr_features = 4;
	string feature_dirs[] = { "bouche", "nez", "oeild", "oeilg" }; // size : nbr_features
	// same order as the _feature_type_t enum

	const unsigned int nbr_samples = 1; // for each directory (feature).
	// todo : find out this number by reading the directories.

	// array of _histo_model_t
	histo_db = (histo_model_t*)malloc(nbr_samples * nbr_features * sizeof(histo_model_t));
	if (!histo_db) { error("error allocating histo_db\n"); return 0; };

	string feature_dirName = calloc(20, sizeof(char));
	string feature_fileName = calloc(30, sizeof(char));
	
	changeDirectory("modeles");
	for (i = 0; i < nbr_features; i++) {
		strcpy_s(feature_dirName, strlen(feature_dirs[i])+1, feature_dirs[i]);
		changeDirectory(feature_dirName);
		for (sample = 1; sample <= nbr_samples; sample++) {
			sprintf_s(feature_fileName, strlen(feature_dirName)+10, "%s_%02d.bmp", feature_dirName, sample);
			histo_db[idx].histo = readHistoFromFile(feature_fileName);
			histo_db[idx].type = (feature_type_t)i;
			debugPrint("read %s and saved it into : histo_db[%u] (%p)\n", feature_fileName, idx, &(histo_db[i]));
			if (!(histo_db[idx].histo)) {
				error("*** histo_db[%u].histo is NULL ! Going to the next element overwriting histo_db[%d] ... ***\n", idx, idx);
			} else {
				idx++;
			}
		}
		changeDirectory("..");
	}
	changeDirectory("..");

	secure_free(feature_fileName);
	secure_free(feature_dirName);

	return idx; // count
}

uint compare_two_histograms(uint* histo1, uint* histo2)
{
	uint val; // ressemblance value;

	//debug
	val = 0;

	return val;
}

// todo : change type
void compare_histo_with_models(uint* histo)
{
	/* En résumé :
	 * new histo_modele;
	 * int tmp, bestFit
	 * 
	 * bestFit = 0
	 * bestType = ""
	 * 
	 *     tmp = compare_two_histograms(histo_modele, histo)
	 *     if (tmp>bestFit) { 
	 *         bestFit = tmp;
	 *         bestType = type;
	 *     }
	 * }
	 * 
	 * printf stuff.
	 *
	 **/
}

u16** get_subimage(u16** src, int src_w, int src_h, int x, int y, int w, int h)
{
	int i, j;
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
uint extract_subimages_and_save(u16** image_ng, int width, int height)
{
	int i, j;
	uint tmp, counter = 0;
	uint nbr_sub_x = 5;
	uint sub_size = (uint)roundf(((float)width / (float)nbr_sub_x)); // == width == height (square)
	int loop_step = sub_size >> 1; // génération d'une sous-image à chaque taille/2.
	uint nbr_sub_y = (uint)roundf(((float)height / (float)sub_size));

	u16** sub = NULL;
	uint* sub_histo = NULL;
	string sub_filename = NULL;
	sub_filename = calloc(strlen(nomFichier) + 35, sizeof(char)); // +15 => "_sub_x_y_w.bmp" + "histo" + extra safety

	sprintf(sub_filename, "%s_sub_images", nomFichier);
	createDirectory(sub_filename);
	changeDirectory(sub_filename);

	DonneesImageRGB* dest_imgRGB = new_ImageRGB(sub_size, sub_size);
	DonneesImageRGB* sub_histo_img = NULL;
	
	for (j = 0; j < height - 2*loop_step; j += loop_step) {
		for (i = 0; i < width - 2*loop_step; i += loop_step) {

			sub = get_subimage(image_ng, width, height, i, j, sub_size, sub_size);
			if (!sub) return counter;
			sprintf(sub_filename, "%s_sub_%d_%d_%u.bmp", nomFichier, j, i, sub_size);
			sauveImageNG(dest_imgRGB, sub);
			ecrisBMPRGB_Dans(dest_imgRGB, sub_filename);

			sub_histo = do_histogramme(sub, sub_size, sub_size);

			// Todo : compare stuff.

			//sprintf(sub_filename, "%s_sub_%d_%d_%u_histo.bin", nomFichier, j, i, sub_size);
			//writeHistoToFile(sub_filename, sub_histo);
			
			//sub_histo_img = imageHistogramme(sub_histo);
			//do_saveBMPwithName(sub_histo_img, sub_filename, "histogramme.bmp");

			secure_free(sub_histo);
			libereDonneesImageRGB(&sub_histo_img);

			for (tmp = 0; tmp < sub_size; tmp++)
				secure_free(sub[tmp]);
			secure_free(sub);
			
			counter++;
		}
	}

	changeDirectory("..");

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
		image_ng = couleur2NG(matrice_bleue, matrice_rouge, matrice_verte, false);

		tmp_ng1 = apply_filter(image_ng, filters[flt_Median]);
		for (i = 0; i < img_h; i++) {
			if (image_ng) secure_free(image_ng[i]);
		}
		secure_free(image_ng);
		image_ng = tmp_ng1;
		tmp_ng1 = NULL;


#ifdef CONSOLE

		if (!isDoingAll) {

			printf("******************************\n");
			printf("**** Bertrand - Debournoux ***\n");
			printf("******* Biometrie - LBP ******\n");
			printf("***** v1.17 - 22/02/2014 *****\n");
			printf("******************************\n\n");
			printf("Image en cours : %s\n\n", nomFichier);
			printf("******************************\n\n");
			printf("(Un filtre median est fait avant toute chose)\n");
			printf("* 1) LBP \n");
			printf("* 2) Filtre median seul \n");
			printf("* 3) Histogramme (image + bin)\n");
			printf("* 4) Détection de visage(s) \n");
			printf("* 5) LBP + Extractions sous-images + Histo \n");
			printf("* 10) Init Models DB \n");
			printf("* 0) Quitter \n\n");
			printf("******************************\n");
			printf("Choix ? \n");
			scanf_s("%d", &choix);
			printf("\n");

		}
#endif

		switch (choix) {
		
		case 1:
			tmp_ng1 = apply_filter(image_ng, filters[flt_LBP]);
			sauveImageNG(image, tmp_ng1);
			saveBMPwithCurrentName(image, "lbp-with-median.bmp");
			break;
		case 2: // median test
			sauveImageNG(image, image_ng);
			saveBMPwithCurrentName(image, "median.bmp");
			break;
		case 3:
			histo = histogramme(image_ng);
			histo_img = imageHistogramme(histo);
			
			saveBMPwithCurrentName(histo_img, "histogramme.bmp");

			string binhisto_name = (string)calloc(strlen(nomFichier) + 15, sizeof(char));
			if (!binhisto_name) { secure_free(histo); break; }
			sprintf(binhisto_name, "%s_histo.bin", nomFichier);
			writeHistoToFile(binhisto_name, histo);

			secure_free(histo);
			secure_free(binhisto_name);
			break;
		case 4:
			debugPrint("Detection de visage :  niveau-de-gris > mediane > lbp > sous-images > histogrammes > comparaison avec modeles > deductions \n");
			break;
		case 5:
			tmp_ng1 = apply_filter(image_ng, filters[flt_LBP]);
			sauveImageNG(image, tmp_ng1);
			saveBMPwithCurrentName(image, "lbp.bmp");
			int counter = extract_subimages_and_save(tmp_ng1, img_w, img_h);
			debugPrint("counter : %d\n", counter);
			break;
		case 10:
			printf("Making Models DB ...\n");
			histo_db_size = make_histo_db();
			printf("DB OK. Number of elements : %d\n", histo_db_size);
			break;
		case 0:
			end = true;
			break;
		default:
			printf("Mauvais choix !\n\n");
			break;
		}

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
	nomFichier = (string)calloc(350, sizeof(char));
	if (!nomFichier) return;

#ifndef CONSOLE
	debugPrint("loaded : %s\n", argv[1]);
#endif
	strcpy_s(nomFichier, (argc > 1) ? strlen(argv[1]) + 1 : 14, (argc > 1) ? argv[1] : "image.bmp");

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
	if (!(matrice_rouge[0] && matrice_verte[0] && matrice_bleue[0]))
		exit(-1);
}

void freeStuff()
{
	int i;
	for (i = 0; i < NBR_FILTRES; i++)
		freeFilter(filters[i]);
	secure_free(filters);

	secure_free(nomFichier);

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
