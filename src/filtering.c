// Adrien Bertrand
// Biométrie - LBP
// v1.0 - 12/02/2014

#include "filtering.h"

filter_t** filters;

void freeFilter(filter_t* flt)
{
	if (flt && flt->mask) {
		uint i;
		for (i = 0; i < flt->size; i++)
			secure_free(flt->mask[i]);
		secure_free(flt->mask);
	}
	secure_free(flt);
}

#define INT_SWAP(a,b) { register int t=(a);(a)=(b);(b)=t; }

int getMedian(int* arr, uint n)
{
	int low, high;
	int median;
	int middle, ll, hh;

	low = 0; high = n - 1; median = (low + high) >> 1;
	for (;;) {

		if (high <= low) // One element only
			return arr[median];

		if (high == low + 1) {  // Two elements only
			if (arr[low] > arr[high])
				INT_SWAP(arr[low], arr[high]);
			return arr[median];
		}

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) >> 1;
		if (arr[middle] > arr[high])    INT_SWAP(arr[middle], arr[high]);
		if (arr[low] > arr[high])       INT_SWAP(arr[low], arr[high]);
		if (arr[middle] > arr[low])     INT_SWAP(arr[middle], arr[low]);

		/* Swap low item (now in position middle) into position (low+1) */
		INT_SWAP(arr[middle], arr[low + 1]);

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do ll++; while (arr[low] > arr[ll]);
			do hh--; while (arr[hh] > arr[low]);

			if (hh < ll)
				break;

			INT_SWAP(arr[ll], arr[hh]);
		}

		/* Swap middle item (in position low) back into correct position */
		INT_SWAP(arr[low], arr[hh]);

		/* Re-set active partition */
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}
}
#undef ELEM_SWAP

void reNormalize(int** imageNG)
{
	int min = INT_MAX, max = INT_MIN, val;
	double ratio = 0.0;
	int i, j;
	for (j = 0; j < img_h; j++) {
		for (i = 0; i < img_w; i++) {
			val = imageNG[j][i];
			if (val < min) min = val;
			if (val > max) max = val;
		}
	}

	ratio = ((double)255) / (max - min);

#pragma omp parallel for
	for (j = 0; j < img_h; j++)
		for (i = 0; i < img_w; i++)
			imageNG[j][i] = (int)((imageNG[j][i]-min) * ratio);

}

void reBound(int** imageNG)
{
	int i, j;
	for (j = 0; j < img_h; j++)
	for (i = 0; i < img_w; i++)
		imageNG[j][i] = (imageNG[j][i] < 0) ? 0 : ((imageNG[j][i]) > 255 ? 255 : imageNG[j][i]);
}

u16** apply_filter(u16** src, filter_t* filter)
{
	int filter_size = filter->size;

	if (filter_size % 2 == 0) {
		error("Erreur : le filtre n'est pas de taille impaire !\n");
		return NULL;
	}
	if (filter_size > img_w) {
		error("Erreur : le filtre est plus grand que l'image source !\n");
		return NULL;
	}

	int i, j;
	int x, y;
	int tmp;

	int* median_array = NULL;
	if (filter->method == flt_m_Median) {
		median_array = (int*)calloc(filter_size * filter_size, sizeof(int));
		if (!median_array) return NULL;
	}

	uint median_array_idx;


	int** imageNG = (int**)malloc(img_w*img_h * sizeof(int*));
	if (!imageNG) return NULL;
	for (j = 0; j < img_h; j++) {
		imageNG[j] = (int*)calloc(img_w, sizeof(int));
		if (!imageNG[j]) return NULL;
	}

	int offset = filter_size >> 1;

	register int tmp_lbp;

	switch (filter->method) {
	case flt_m_LBP:
		for (j = offset; j < img_h - offset; j++)
		for (i = offset, tmp = 0; i < img_w - offset; i++, tmp = 0) {
			tmp_lbp = (int)(src[j][i]);
			for (y = -offset; y <= offset; y++)
			for (x = -offset; x <= offset; x++) 
				tmp += ((int)(src[j + y][x + i]) >= tmp_lbp) ? ((int)(filter->mask[y + offset][x + offset])) : 0;
			//error("tmp = %d\n", tmp);
			imageNG[j][i] = tmp;
		}
		break;
	case flt_m_Median:
		for (j = offset; j < img_h - offset; j++)
		for (i = offset, median_array_idx = 0; i < img_w - offset; i++, median_array_idx = 0) {
			for (y = -offset; y <= offset; y++)
			for (x = -offset; x <= offset; x++)
				median_array[median_array_idx++] = (int)(src[j + y][x + i]);
			imageNG[j][i] = getMedian(median_array, filter_size*filter_size);
		}
		secure_free(median_array);
		break;
	default:
#pragma omp parallel for
		for (j = offset; j < img_h - offset; j++)
		for (i = offset, tmp = 0; i < img_w - offset; i++)
			imageNG[j][i] = src[j][i];
		break;
	}

	if (filter->method == flt_m_LBP) {
		if (filter->reBoundOnly)
			reBound(imageNG);
		else
			reNormalize(imageNG);
	}

	u16** imageNG_16 = (u16**)malloc(img_h * sizeof(u16*));
	if (!imageNG_16) return NULL;
	for (j = 0; j < img_h; j++) {
		imageNG_16[j] = (u16*)calloc(img_w, sizeof(u16));
		if (!imageNG_16[j]) return NULL;
#pragma omp parallel for
		for (i = 0; i < img_w; i++)
			imageNG_16[j][i] = imageNG[j][i];
	}

#pragma omp parallel for
	for (j = 0; j < img_h; j++)
		secure_free(imageNG[j]);
	secure_free(imageNG);

	return imageNG_16;
}

void mask_copy_3(int** dest, int src[3][3])
{
	uint i, j;
	for (j = 0; j < 3; j++)
	for (i = 0; i < 3; i++)
		dest[j][i] = src[j][i];
}

void matrix_copy(int** dest, int** src, uint cols, uint rows)
{
	uint i, j;
	for (j = 0; j < rows; j++)
	for (i = 0; i < cols; i++)
		dest[j][i] = src[j][i];
}


filter_t** createFilters()
{
	uint i, j, size;

	filters = (filter_t**)malloc(NBR_FILTRES * sizeof(filter_t*));
	if (!filters) return NULL;

	int tmp_LBP[3][3] = { { 1, 2, 4 }, { 8, 0, 16 }, { 32, 64, 128 } };

	// default properties
	for (i = 0; i < NBR_FILTRES; i++)
	{
		filters[i] = (filter_t*)malloc(sizeof(filter_t));
		if (!filters[i]) return NULL;
		filters[i]->needNormalization = false;
		filters[i]->div = 0;
		filters[i]->reBoundOnly = false;
		size = 3;
		filters[i]->size = size;
		if (i != flt_Median) {
			filters[i]->mask = (int**)malloc(size * sizeof(int*));
			if (!filters[i]->mask) return NULL;
			for (j = 0; j < size; j++) {
				filters[i]->mask[j] = (int*)calloc(size, sizeof(int));
				if (!filters[i]->mask[j]) return NULL;
			}
		} else {
			filters[i]->mask = NULL;
		}
	}

	filters[flt_LBP]->method = flt_m_LBP;
	mask_copy_3(filters[flt_LBP]->mask, tmp_LBP);

	filters[flt_Median]->method = flt_m_Median;


	return filters;
}
