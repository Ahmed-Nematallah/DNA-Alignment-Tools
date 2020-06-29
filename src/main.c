
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>

#define UNUSED __attribute__ ((unused))

typedef short arr_type;

typedef enum {
	SECOND_GAP,
	FIRST_GAP,
	EQUAL,
	MISMATCH
} FIN_CASES;

arr_type max3(arr_type a, arr_type b, arr_type c) {
	if (a > b) {
		if (a > c) {
			return a;
		} else {
			return c;
		}
	} else {
		if (b > c) {
			return b;
		} else {
			return c;
		}
	}
}

void reverse_char_arr(char * arr, uint64_t len) {
	char t;
	for (uint64_t i = 0; i < len / 2; i++) {
		t = arr[i];
		arr[i] = arr[len - i - 1];
		arr[len - i - 1] = t;
	}
}

void print_similarities(const uint64_t start, const uint64_t len, const char * DNA1, const char * DNA2, char * fin, uint64_t fin_len, uint64_t * tmp1, uint64_t * tmp2, FILE * fp) {
	for (uint64_t i = start; (i < fin_len) && (i < start + len); i++) {
		switch (fin[i]) {
			case EQUAL:
				putc(DNA1[(*tmp1)++], fp);
				break;
			case MISMATCH:
				putc(DNA1[(*tmp1)++], fp);
				break;
			case FIRST_GAP:
				putc('-', fp);
				break;
			case SECOND_GAP:
				putc(DNA1[(*tmp1)++], fp);
				break;
			
			default:
				break;
		}
	}

	putc('\n', fp);

	for (uint64_t i = start; (i < fin_len) && (i < start + len); i++) {
		switch (fin[i]) {
			case EQUAL:
				putc('|', fp);
				break;
			
			default:
				putc(' ', fp);
				break;
		}
	}

	putc('\n', fp);

	for (uint64_t i = start; (i < fin_len) && (i < start + len); i++) {
		switch (fin[i]) {
			case EQUAL:
				putc(DNA2[(*tmp2)++], fp);
				break;
			case MISMATCH:
				putc(DNA2[(*tmp2)++], fp);
				break;
			case FIRST_GAP:
				putc(DNA2[(*tmp2)++], fp);
				break;
			case SECOND_GAP:
				putc('-', fp);
				break;
			
			default:
				break;
		}
	}

	putc('\n', fp);
}

void print_differences(const char * DNA1, const char * DNA2, char * fin, uint64_t fin_len, FILE * fp) {
	uint64_t x = 0;
	uint64_t y = 0;
	uint64_t j = 0;

	for (uint64_t i = 0; i < fin_len; i++) {
		if (fin[i] != EQUAL) {
			j = i + 1;
			// for (; fin[j] == fin[i]; j++);

			if (fin[i] == MISMATCH) {
				// fprintf(fp, "Difference found in offset %"PRId64" in first sequence and offset %"PRId64" in second sequence.\n", x, y);
				// fprintf(fp, "Length of difference is %"PRId64"\n", j - i);
				// fprintf(fp, "Data from first sequence is:\n%.*s\n", (int)(j - i), DNA1 + x);
				// fprintf(fp, "Data from second sequence is:\n%.*s\n", (int)(j - i), DNA2 + y);
				fprintf(fp, "%05"PRId64"\t%c>%c\n", x + 1, DNA1[x], DNA2[y]);

				x += j - i;
				y += j - i;
			} else if (fin[i] == FIRST_GAP) {
				// fprintf(fp, "Data not in first sequence found in offset %"PRId64" in second sequence.\n", y);
				// fprintf(fp, "Length of gap is %"PRId64"\n", j - i);
				// fprintf(fp, "Data from second sequence is:\n%.*s\n", (int)(j - i), DNA2 + y);
				fprintf(fp, "%05"PRId64"\t->%c\n", x + 1, DNA2[y]);

				y += j - i;
			} else {
				// fprintf(fp, "Data not in second sequence found in offset %"PRId64" in first sequence.\n", x);
				// fprintf(fp, "Length of gap is %"PRId64"\n", j - i);
				// fprintf(fp, "Data from first sequence is:\n%.*s\n", (int)(j - i), DNA1 + x);
				fprintf(fp, "%05"PRId64"\t%c>-\n", x + 1, DNA1[x]);

				x += j - i;
			}
			// printf("---------------------------------------------------------------\n");
			// printf("%d  %d\n", i, j);
			i = j - 1;
		} else {
			x++;
			y++;
		}
	}
}

void gen_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN) {
	for (uint64_t i = 0; i < DNA1_LEN + 1; i++) {
		arr[i * (DNA2_LEN + 1) + 0] = GAP_PEN * i;
	}

	for (uint64_t i = 0; i < DNA2_LEN + 1; i++) {
		arr[0 * (DNA2_LEN + 1) + i] = GAP_PEN * i;
	}

	for (uint64_t i = 1; i < DNA1_LEN + 1; i++) {
		for (uint64_t j = 1; j < DNA2_LEN + 1; j++) {
			arr[i * (DNA2_LEN + 1) + j] = max3(arr[(i - 1) * (DNA2_LEN + 1) + (j - 1)] + ((DNA1[i - 1] == DNA2[j - 1]) ? NOO_PEN : MIS_PEN), 
			                                   arr[i * (DNA2_LEN + 1) + (j - 1)] + GAP_PEN, 
										       arr[(i - 1) * (DNA2_LEN + 1) + j] + GAP_PEN);
		}
	}
}

void gen_fin_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, char * fin, uint64_t * fin_len, uint64_t * sim, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN) {
	uint64_t x = DNA1_LEN;
	uint64_t y = DNA2_LEN;

	while(1) {
		// printf("%d  %d\n", x, y);

		if ((x == 0) && (y == 0)) {
			break;
		}

		if (arr[x * (DNA2_LEN + 1) + y] == 
		    (arr[(x - 1) * (DNA2_LEN + 1) + (y - 1)] + ((DNA1[x - 1] == DNA2[y - 1]) ? NOO_PEN : MIS_PEN))) {
			fin[(*fin_len)++] = (DNA1[x - 1] == DNA2[y - 1]) ? EQUAL : MISMATCH;
			if (DNA1[x - 1] == DNA2[y - 1]) {
				(*sim)++;
			}
			x--; y--;
		} else if (arr[x * (DNA2_LEN + 1) + y] == (arr[(x - 1) * (DNA2_LEN + 1) + y] + GAP_PEN)) {
			fin[(*fin_len)++] = SECOND_GAP;
			x--;
		} else if (arr[x * (DNA2_LEN + 1) + y] == (arr[x * (DNA2_LEN + 1) + (y - 1)] + GAP_PEN)) {
			fin[(*fin_len)++] = FIRST_GAP;
			y--;
		} else {
			exit(-1);
		}
	}
	
	reverse_char_arr(fin, *fin_len);
}

void print_2d_arr(arr_type * arr, uint64_t dim1, uint64_t dim2, FILE * fp) {
	for (uint64_t i = 0; i < dim1 + 1; i++) {
		for (uint64_t j = 0; j < dim2 + 1; j++) {
			fprintf(fp, "%04d  ", arr[i * (dim2 + 1) + j]);
		}
		fputc('\n', fp);
	}
}

void shift_arr_left(char * arr, uint64_t start, uint64_t len, uint64_t shift) {
	for (uint64_t i = start; i < len - shift; i++) {
		arr[i] = arr[i + shift];
	}
}

uint64_t sanitize_input(char * input) {
	uint64_t len = strlen(input);
	uint64_t i;

	for (i = 0; i < len; i++)
		if (input[i] >= 'a')
			input[i] -= 0x20;

	for (i = 0; i < len; i++) {
		if (!((input[i] == 'A') || (input[i] == 'C') || (input[i] == 'T') || (input[i] == 'G'))) {
			shift_arr_left(input, i, len, 1);
			i--;
			len--;
		}
	}

	return len;
}

int main (UNUSED int argc, UNUSED char ** argv) {
	const char * FILE1_NAME = "DNA1_file.txt";
	const char * FILE2_NAME = "DNA2_file.txt";

	const char * FILE_OUT_NAME = "output.txt";

	FILE * fDNA1 = fopen64(FILE1_NAME, "r");
	FILE * fDNA2 = fopen64(FILE2_NAME, "r");

	fseek(fDNA1, 0L, SEEK_END);
	fseek(fDNA2, 0L, SEEK_END);

	uint64_t fDNA1_len = ftello64(fDNA1);
	uint64_t fDNA2_len = ftello64(fDNA2);

	rewind(fDNA1);
	rewind(fDNA2);

	char * DNA1 = calloc(fDNA1_len + 1, 1);
	char * DNA2 = calloc(fDNA2_len + 1, 1);

	fread(DNA1, 1, fDNA1_len, fDNA1);
	fread(DNA2, 1, fDNA2_len, fDNA2);

	fclose(fDNA1);
	fclose(fDNA2);

	FILE * fout;

	fout = fopen64(FILE_OUT_NAME, "w");
	// printf("%d\n", sizeof(uint64_t));

	uint64_t DNA1_LEN = sanitize_input(DNA1);
	uint64_t DNA2_LEN = sanitize_input(DNA2);

	// printf("%d   %d\n", DNA1_LEN, DNA2_LEN);

	arr_type GAP_PEN = -2;
	arr_type MIS_PEN = -1;
	arr_type NOO_PEN = 1;

	arr_type *arr = calloc((DNA1_LEN + 1), (DNA2_LEN + 1) * 4ULL);
	char *fin = calloc(sizeof(char) * (DNA1_LEN + DNA2_LEN + 2), 1);
	uint64_t fin_len = 0;

	if ((arr == NULL) || (fin == NULL)) {
		printf("Error: could not allocate memory\n");
		exit(-1);
	}

	clock_t c1;

	c1 = clock();

	gen_arr(DNA1, DNA2, DNA1_LEN, DNA2_LEN, arr, GAP_PEN, MIS_PEN, NOO_PEN);

	uint64_t sim = 0;

	gen_fin_arr(DNA1, DNA2, DNA1_LEN, DNA2_LEN, arr, fin, &fin_len, &sim, GAP_PEN, MIS_PEN, NOO_PEN);

	// print_2d_arr(arr, DNA1_LEN, DNA2_LEN, fout);

	free(arr);

	fprintf(fout, "Elapsed time: %lf\n", (double)(clock() - c1) / CLOCKS_PER_SEC);

	fprintf(fout, "Similarity percentage %lf%%\n", (double)sim / (double)DNA1_LEN * 100);

	// uint64_t t1 = 0;
	// uint64_t t2 = 0;

	// uint64_t print_per_screen = 80;
	// for (uint64_t i = 0; i < fin_len; i += print_per_screen) {
	// 	print_similarities(i, print_per_screen, DNA1, DNA2, fin, fin_len, &t1, &t2, fout);
	// 	putc('\n', fout);
	// }

	print_differences(DNA1, DNA2, fin, fin_len, fout);

	fclose(fout);
	free(fin);

	return 0;
}
