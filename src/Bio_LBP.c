// Adrien Bertrand
// Biométrie - LBP
// v1.3 - 28/02/2014

#include "Bio_LBP.h"
#include "filtering.h"
#include "utils.h"


// TODO : check distance between face features to make sure it's a face.


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

histo_model_t* histo_models_db;
uint histo_db_size = 0;

histo_ownImage_t* histo_ownImage_db;

/*  --------- Fonctions ---------  */

unsigned long int getIfromXYinImage(DonneesImageRGB* img, uint x, uint y) {
	return (x < (uint)image->largeurImage && y < (uint)image->hauteurImage) ? (3 * (x + y * image->largeurImage)) : (unsigned long int)(-1);
}

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
//#pragma omp parallel for
	for (j = 0; j < img_h; j++) {
		for (i = 0; i < img_w; i++) {
			mat_rouge[j][i] = 255 - mat_rouge[j][i];
			mat_verte[j][i] = 255 - mat_verte[j][i];
			mat_bleue[j][i] = 255 - mat_bleue[j][i];
		}
	}
}

// real function with widht and height parameters
u16** do_couleur2NG(u16** mat_bleue, u16** mat_rouge, u16** mat_verte, bool perceptive, uint w, uint h)
{
	uint i, j;
	double	coeff_r = perceptive ? 0.2125 : 1,
			coeff_g = perceptive ? 0.7154 : 1,
			coeff_b = perceptive ? 0.0721 : 1;
	uint divi = perceptive ? 1 : 3;

	u16** imageNG = new_u16_mat(w, h);
	if (!imageNG) return NULL;

	for (j = 0; j < h; j++)
		for (i = 0; i < w; i++)
			imageNG[j][i] = (u16)((mat_rouge[j][i] * coeff_r + mat_verte[j][i] * coeff_g + mat_bleue[j][i] * coeff_b) / divi);

	return imageNG;
}

// wrapper for global image with known size.
u16** couleur2NG(u16** mat_bleue, u16** mat_rouge, u16** mat_verte, bool perceptive) {
	return do_couleur2NG(mat_bleue, mat_rouge, mat_verte, perceptive, img_w, img_h);
}

void seuilleImageNG(u16** imageNG, uint seuil)
{
	int i, j = 0;
//#pragma omp parallel for
	for (j = 0; j < img_h; j++)
	for (i = 0; i < img_w; i++)
		imageNG[j][i] = (imageNG[j][i] > seuil) ? 255 : 0;
}

uint* new_histo(void)
{
	return (uint*)calloc(GRAYLEVELS, sizeof(uint));
}

uint* do_histogramme(u16** imageNG, uint w, uint h)
{
	uint i, j;

	histo = new_histo();
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

	uint* histoNorm = new_histo();
	if (!histoNorm) return NULL;

//#pragma omp parallel for
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
//#pragma omp parallel for
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
//#pragma omp parallel for
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
uint make_histo_db(void)
{
	uint i, sample, idx = 0;

	DonneesImageRGB* tmp_img_rgb = NULL;
	uint modele_w, modele_h;

	u16** tmp_img_ng = NULL;
	u16** tmp_img_ng_lbp = NULL;
	u16** tmp_mat_rouge = NULL;
	u16** tmp_mat_vert = NULL;
	u16** tmp_mat_bleu = NULL;

	string feature_dirs[] = { "bouche", "nez", "oeild", "oeilg" };
	// same order as the _feature_type_t enum
	const uint nbr_features = array_count(feature_dirs);

	const uint nbr_samples = 8; // max models for each feature directory
	// todo : find out this number by reading the directories.

	// array of _histo_model_t
	histo_models_db = (histo_model_t*)malloc(nbr_samples * nbr_features * sizeof(histo_model_t));
	if (!histo_models_db) { error("error allocating histo_db\n"); return 0; };

	string feature_dirName = calloc(20, sizeof(char));
	string feature_fileName = calloc(30, sizeof(char));
	if (!feature_fileName || !feature_dirName) return 0;
	
	changeDirectory("modeles");

	for (i = 0; i < nbr_features; i++) {
		strcpy_s(feature_dirName, strlen(feature_dirs[i])+1, feature_dirs[i]);
		changeDirectory(feature_dirName);
		for (sample = 1; sample <= nbr_samples; sample++) {
			sprintf_s(feature_fileName, strlen(feature_dirName)+10, "%s_%02d.bmp", feature_dirName, sample);

			if (tmp_img_rgb) libereDonneesImageRGB(&tmp_img_rgb);
			tmp_img_rgb = lisBMPRGB(feature_fileName);
			if (!tmp_img_rgb) continue;

			modele_w = tmp_img_rgb->largeurImage;
			modele_h = tmp_img_rgb->hauteurImage;

			tmp_mat_bleu = new_u16_mat(modele_w, modele_h);
			tmp_mat_rouge = new_u16_mat(modele_w, modele_h);
			tmp_mat_vert = new_u16_mat(modele_w, modele_h);
			if (!(tmp_mat_vert[modele_h - 1])) return 0;
			
			cree3matrices(tmp_mat_bleu, tmp_mat_rouge, tmp_mat_vert, tmp_img_rgb);
			tmp_img_ng = do_couleur2NG(tmp_mat_bleu, tmp_mat_rouge, tmp_mat_vert, false, modele_w, modele_h);

			free_u16_mat(tmp_mat_bleu, modele_h);
			free_u16_mat(tmp_mat_rouge, modele_h);
			free_u16_mat(tmp_mat_vert, modele_h);

			histo_models_db[idx].histo = do_histogramme(tmp_img_ng, modele_w, modele_h);
			histo_models_db[idx].type = (feature_type_t)i;

			if (tmp_img_ng) secure_free(tmp_img_ng);

			debugPrint("- %s\t histo saved it into : histo_db[%u] (%p)\n", feature_fileName, idx, &(histo_models_db[i]));
			if (!(histo_models_db[idx].histo)) {
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
	
	libereDonneesImageRGB(&tmp_img_rgb);

	return idx; // count
}

uint compare_two_histograms(uint* h1, uint* h2)
{
	int i;
	double tmp, result = 0; // distance value;

	uint max_h1 = h1[array_max_idx((int*)h1, GRAYLEVELS)];
	uint max_h2 = h2[array_max_idx((int*)h2, GRAYLEVELS)];
	double ratio1, ratio2;
	ratio1 = 1.0; ratio2 = 1.0;
	if (max_h1 > max_h2) {
		ratio2 = (double)max_h1 / (double)max_h2;
	} else {
		ratio1 = (double)max_h2 / (double)max_h1;
	}

	// Chi-Squared comparison method.
	for (i = 0; i < GRAYLEVELS; i++) {
		tmp = (double)(h1[i] * ratio1) - (double)(h2[i] * ratio2);
		if ((h1[i] * ratio1) > 0)
			result += ((double)(tmp*tmp) / (double)(h1[i] * ratio1));
	}

	return (uint)(result);
}


histo_ownImage_t* compare_histo_with_models(uint* histo)
{
	uint i, tmp;
	
	// todo : set minimum similarity value ?
	
	histo_ownImage_t* match = (histo_ownImage_t*)malloc(sizeof(histo_ownImage_t));
	if (!match) return NULL;

	match->histo = histo;
	match->distance = UINT_MAX;
	match->type = feat_VOID;
	
	for (i = 0; i < histo_db_size; i++) {
		tmp = compare_two_histograms(histo_models_db[i].histo, histo);
		if (tmp < match->distance) {
			match->distance = tmp;
			match->type = histo_models_db[i].type;
		}
	}

	return match;
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

face_feat_t* extract_subimages_and_compare(u16** image_ng, int width, int height)
{
	int i, j;
	uint k;

	face_feat_t* face_features = NULL;

	histo_ownImage_t* tmp = NULL;
	uint nbr_sub = 5;
	uint sub_size = (uint)roundf(((float)(width<height ? width : height) / (float)nbr_sub));
	int loop_step = sub_size >> 2; // génération d'une sous-image à chaque taille/2.

	u16** sub = NULL;
	uint* sub_histo = NULL;
	string sub_filename = NULL;
	sub_filename = calloc(strlen(nomFichier) + 35, sizeof(char)); // +15 => "_sub_x_y_w.bmp" + "histo" + extra safety

	DonneesImageRGB* dest_imgRGB = new_ImageRGB(sub_size, sub_size);
	DonneesImageRGB* sub_histo_img = NULL;

	// calloc so that every field is set to 0
	histo_ownImage_db = (histo_ownImage_t*)calloc(histo_ownImage_db_size, sizeof(histo_ownImage_t));
	if (!histo_ownImage_db) return 0;
	for (k = 0; k < histo_ownImage_db_size; k++) {
		if (!(histo_ownImage_db[k].histo = new_histo()))
			return 0;
		histo_ownImage_db[k].distance = UINT_MAX;
	}

	for (j = 0; j < height - 2*loop_step; j += loop_step) {
		for (i = 0; i < width - 2*loop_step; i += loop_step) {

			sub = get_subimage(image_ng, width, height, i, j, sub_size, sub_size);
			if (!sub) return NULL;

			//sprintf(sub_filename, "%s_sub_%d_%d_%u.bmp", nomFichier, j, i, sub_size);
			//sauveImageNG(dest_imgRGB, sub);
			//ecrisBMPRGB_Dans(dest_imgRGB, sub_filename);

			sub_histo = do_histogramme(sub, sub_size, sub_size);
			
			tmp = compare_histo_with_models(sub_histo);

			if (tmp->type != feat_VOID) {
				if (tmp->distance < histo_ownImage_db[(int)tmp->type].distance) {
					//debugPrint("found better match of type %d (score : %u)!\n", tmp->type, tmp->reliability);
					memcpy(histo_ownImage_db[(int)tmp->type].histo, tmp->histo, GRAYLEVELS*sizeof(uint));
					histo_ownImage_db[(int)tmp->type].distance = tmp->distance;
					histo_ownImage_db[(int)tmp->type].type = tmp->type;
					histo_ownImage_db[(int)tmp->type].x = i;
					histo_ownImage_db[(int)tmp->type].y = j;
				}
			}

			secure_free(sub_histo);
			libereDonneesImageRGB(&sub_histo_img);

			free_u16_mat(sub, sub_size);
		}
	}

	string features_tmp[] = { "bouche", "nez", "oeild", "oeilg" };

	printf("Best matches saved.\n");
	for (k = 0; k < histo_ownImage_db_size; k++) {
		i = histo_ownImage_db[k].x;
		j = histo_ownImage_db[k].y;
		debugPrint("At (%u,%u) : type %u (score : %u)\n", i, j, histo_ownImage_db[k].type, histo_ownImage_db[k].distance);
		sub = get_subimage(image_ng, width, height, i, j, sub_size, sub_size);
		sprintf(sub_filename, "%s_%s.bmp", nomFichier, features_tmp[histo_ownImage_db[k].type]);
		sauveImageNG(dest_imgRGB, sub);
		ecrisBMPRGB_Dans(dest_imgRGB, sub_filename);
	}

	face_features = (face_feat_t*)calloc(histo_ownImage_db_size, sizeof(face_feat_t));
	for (k = 0; k < histo_ownImage_db_size; k++) {
		face_features[k].size = sub_size;
		face_features[k].type = histo_ownImage_db[k].type;
		face_features[k].x = histo_ownImage_db[k].x;
		face_features[k].y = histo_ownImage_db[k].y;
	}

	for (k = 0; k < histo_ownImage_db_size; k++)
		secure_free(histo_ownImage_db[k].histo);
	secure_free(histo_ownImage_db);

	secure_free(sub_filename);

	libereDonneesImageRGB(&dest_imgRGB);

	return face_features;
}

// will handle colors later
void drawPixelOnImage(uint x, uint y, DonneesImageRGB* img)
{
	unsigned long int i = getIfromXYinImage(img, x, y);
	if (i > (unsigned long int)(3 * img->largeurImage*img->hauteurImage)) return;
	img->donneesRGB[i] = 0;		// b
	img->donneesRGB[i+1] = 0;	// g
	img->donneesRGB[i+2] = 255; // r
}

// Bressenham line drawing algorithm
// adapted from http://cboard.cprogramming.com/game-programming/67832-line-drawing-algorithm.html#post485086
void drawLineOnImage(uint x1, uint y1, uint x2, uint y2, DonneesImageRGB* img)
{
	int dx, dy, inx, iny, e;

	dx = x2 - x1;
	dy = y2 - y1;
	inx = dx > 0 ? 1 : -1;
	iny = dy > 0 ? 1 : -1;

	dx = ABS(dx);
	dy = ABS(dy);

	if (dx >= dy) {
		dy <<= 1;
		e = dy - dx;
		dx <<= 1;
		while (x1 != x2) {
			drawPixelOnImage(x1, y1, img);
			if (e >= 0) {
				y1 += iny;
				e -= dx;
			}
			e += dy; x1 += inx;
		}
	} else {
		dx <<= 1;
		e = dx - dy;
		dy <<= 1;
		while (y1 != y2) {
			drawPixelOnImage(x1, y1, img);
			if (e >= 0) {
				x1 += inx;
				e -= dy;
			}
			e += dx; y1 += iny;
		}
	}
	drawPixelOnImage(x1, y1, img);
}

void drawRectangleOnImage(const rect_t* rect, DonneesImageRGB* img)
{
	drawLineOnImage(rect->x, rect->y, rect->x + rect->w, rect->y, img); // top
	drawLineOnImage(rect->x, rect->y, rect->x, rect->y + rect->h, img); // left
	drawLineOnImage(rect->x, rect->y + rect->h, rect->x + rect->w, rect->y + rect->h, img); // bottom
	drawLineOnImage(rect->x + rect->w, rect->y, rect->x + rect->w, rect->y + rect->h, img); // right
}

void mark_face_features(DonneesImageRGB* image_orig, face_feat_t* face_features)
{
	uint k;
	uint xmin, ymin, xmax, ymax;

	rect_t tmp_rect, face_rect;
	string fileName = NULL;
	DonneesImageRGB* img = NULL;
	
	img = clone_imageRGB(image_orig);
	if (!img) return;

	xmin = ymin = UINT_MAX;
	xmax = ymax = 0;

	for (k = 0; k < histo_ownImage_db_size; k++) {
		if (face_features[k].x == 0 && face_features[k].y == 0) {
			debugPrint("ignoring a probably badly-detected feature (%d)\n", k);
			continue;
		}
		tmp_rect = (rect_t){ face_features[k].x, face_features[k].y, face_features[k].size, face_features[k].size };
		drawRectangleOnImage(&tmp_rect, img);
		if (face_features[k].x < xmin) xmin = face_features[k].x;
		if (face_features[k].y < ymin) ymin = face_features[k].y;
		if (face_features[k].x > xmax) xmax = face_features[k].x;
		if (face_features[k].y > ymax) ymax = face_features[k].y;
	}

	int size = face_features[0].size;
	face_rect.x = xmin - (size >> 1);
	face_rect.y = ymin - (uint)((double)size/1.5);
	face_rect.w = (xmax - xmin) + size * 2;
	face_rect.h = (ymax - ymin) + (uint)(size * 2.5);

	drawRectangleOnImage(&face_rect, img);

	fileName = (string)calloc(strlen(nomFichier) + 10, sizeof(char));
	if (!fileName) return;

	sprintf_s(fileName, strlen(nomFichier)+10, "%s_face.bmp", nomFichier);
	ecrisBMPRGB_Dans(img, fileName);

	secure_free(fileName);
	libereDonneesImageRGB(&img);
}

void choixAction(int choix)
{
	static bool isDoingAll = false;
	bool end = (choix == 0);

	face_feat_t* face_features = NULL;

	while (!end) {

		cree3matrices(matrice_bleue, matrice_rouge, matrice_verte, image_orig);
		image_ng = couleur2NG(matrice_bleue, matrice_rouge, matrice_verte, false);

		tmp_ng1 = apply_filter(image_ng, filters[flt_Median]);
		
		free_u16_mat(image_ng, img_h);

		image_ng = tmp_ng1;
		tmp_ng1 = NULL;


#ifdef CONSOLE

		if (!isDoingAll) {

			printf("******************************\n");
			printf("**** Bertrand - Debournoux ***\n");
			printf("******* Biometrie - LBP ******\n");
			printf("***** v1.3 - 28/02/2014 *****\n");
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
			if (histo_db_size == 0) histo_db_size = make_histo_db();
			face_features = extract_subimages_and_compare(tmp_ng1, img_w, img_h);
			debugPrint("Marking detected face features on image.\n");
			mark_face_features(image_orig, face_features);
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

		free_u16_mat(image_ng, img_h);
		// to avoid useless calls
		if (tmp_ng1) free_u16_mat(tmp_ng1, img_h);
		if (tmp_ng2) free_u16_mat(tmp_ng2, img_h);
		if (tmp_ng3) free_u16_mat(tmp_ng3, img_h);
		if (face_features) secure_free(face_features);

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

	img_w = image->largeurImage;
	img_h = image->hauteurImage;

	matrice_bleue = new_u16_mat(img_w, img_h);
	matrice_rouge = new_u16_mat(img_w, img_h);
	matrice_verte = new_u16_mat(img_w, img_h);
	if (!(matrice_verte[img_h - 1])) exit(-1); // check last alloc
	
}

void freeStuff(void)
{

	// todo : free histo_models_db et histo_ownImage_db

	int i;
	for (i = 0; i < NBR_FILTRES; i++)
		freeFilter(filters[i]);
	secure_free(filters);

	secure_free(nomFichier);

	free_u16_mat(matrice_bleue, img_h);
	free_u16_mat(matrice_rouge, img_h);
	free_u16_mat(matrice_verte, img_h);

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
