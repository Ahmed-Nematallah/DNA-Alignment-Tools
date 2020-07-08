
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>

#define UNUSED __attribute__ ((unused))

typedef short arr_type;
arr_type (*str_to_arr_type)(const char *) = atoll;

#define IS_BASIC_DNA_CHAR(x) (((x) == 'A') || ((x) == 'C') || ((x) == 'T') || ((x) == 'G'))
#define IS_EXTRA_DNA_CHAR(x) (((x) == 'R') || ((x) == 'Y') || ((x) == 'M') || ((x) == 'K') || \
                              ((x) == 'S') || ((x) == 'W') || ((x) == 'H') || ((x) == 'B') || \
							  ((x) == 'V') || ((x) == 'D') || ((x) == 'N'))

typedef enum {
	SECOND_GAP,
	FIRST_GAP,
	EQUAL,
	MISMATCH
} FIN_CASES;

typedef enum {
	NEEDLEMAN_WUNSCH,
	SMITH_WATERMAN,
	UNKNOWN
} ALIGNMENT_ALGORITHMS;

void reverse_char_arr(char * arr, uint64_t len);
void print_similarities(const uint64_t start, const uint64_t len, const char * DNA1, const char * DNA2, char * fin, uint64_t fin_len, uint64_t * tmp1, uint64_t * tmp2, FILE * fp);
void print_differences(const char * DNA1, const char * DNA2, char * fin, uint64_t fin_len, uint64_t x, uint64_t y, FILE * fp);
void nw_gen_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN);
void sw_gen_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN);
void nw_gen_fin_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, char * fin, uint64_t * fin_len, uint64_t * sim, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN);
void sw_gen_fin_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, char * fin, uint64_t * fin_len, uint64_t * sim, uint64_t * total, uint64_t start_x, uint64_t start_y, uint64_t * fin_x, uint64_t * fin_y, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN);
void print_2d_arr(arr_type * arr, uint64_t dim1, uint64_t dim2, FILE * fp);
void shift_arr_left(char * arr, uint64_t start, uint64_t len, uint64_t shift);
void find_max_val(arr_type * arr, uint64_t dim1, uint64_t dim2, arr_type * max_val);
uint64_t sanitize_input(char * input);
void print_usage(char * filename);

inline arr_type max2(arr_type a, arr_type b) {
	if (a > b) {
		return a;
	} else {
		return b;
	}
}

inline arr_type max3(arr_type a, arr_type b, arr_type c) {
	return max2(a, max2(b, c));
}

inline arr_type max4(arr_type a, arr_type b, arr_type c, arr_type d) {
	return max2(a, max3(b, c, d));
}

int main (int argc, char ** argv) {

	int output_file_idx = -1;
	int sequence_1_idx = -1;
	int sequence_2_idx = -1;

	int GAP_PEN_idx = -1;
	int MIS_PEN_idx = -1;
	int NOO_PEN_idx = -1;

	arr_type GAP_PEN = -2;
	arr_type MIS_PEN = -1;
	arr_type NOO_PEN = 1;

	ALIGNMENT_ALGORITHMS alg = NEEDLEMAN_WUNSCH;

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'S':
					if (argv[i][2] == '1') {
						sequence_1_idx = i + 1;
						i++;
					} else if (argv[i][2] == '2') {
						sequence_2_idx = i + 1;
						i++;
					} else {
						printf("Unknown argument %s\n", argv[i]);
						exit(-1);
					}

					break;
				case 'P':
					if (argv[i][2] == 'g') {
						GAP_PEN_idx = i + 1;
						i++;
					} else if (argv[i][2] == 'm') {
						MIS_PEN_idx = i + 1;
						i++;
					} else if (argv[i][2] == 'n') {
						NOO_PEN_idx = i + 1;
						i++;
					} else {
						printf("Unknown argument %s\n", argv[i]);
						exit(-1);
					}

					break;
				case 'A':
					if (argv[i][2] == 'n') {
						alg = NEEDLEMAN_WUNSCH;
					} else if (argv[i][2] == 's') {
						alg = SMITH_WATERMAN;
					} else {
						printf("Unknown argument %s\n", argv[i]);
						exit(-1);
					}

					break;
				case 'o':
					output_file_idx = i + 1;
					i++;
					break;
				case 'h':
					print_usage(argv[0]);
					exit(0);
				default:
					printf("Unknown argument %s\n", argv[i]);
					exit(-1);
					break;
			}
		} else {
			printf("Unknown argument %s\n", argv[i]);
			exit(-1);
		}
	}

	if (alg == UNKNOWN) {
		printf("Error, alignment algorithm not specified\n");
		print_usage(argv[0]);
	}

	if (GAP_PEN_idx >= 0) {
		GAP_PEN = str_to_arr_type(argv[GAP_PEN_idx]);
	}

	if (MIS_PEN_idx >= 0) {
		MIS_PEN = str_to_arr_type(argv[MIS_PEN_idx]);
	}

	if (NOO_PEN_idx >= 0) {
		NOO_PEN = str_to_arr_type(argv[NOO_PEN_idx]);
	}

	const char * FILE1_NAME = "DNA1_file.txt";
	const char * FILE2_NAME = "DNA2_file.txt";

	if (sequence_1_idx >= 0) {
		FILE1_NAME = argv[sequence_1_idx];
	}

	if (sequence_2_idx >= 0) {
		FILE2_NAME = argv[sequence_2_idx];
	}

	FILE * fDNA1 = fopen64(FILE1_NAME, "r");
	FILE * fDNA2 = fopen64(FILE2_NAME, "r");

	if (!fDNA1) {
		printf("Could not open first DNA file\n");
		exit(-1);
	}

	if (!fDNA2) {
		printf("Could not open second DNA file\n");
		exit(-1);
	}

	fseek(fDNA1, 0L, SEEK_END);
	fseek(fDNA2, 0L, SEEK_END);

	uint64_t fDNA1_len = ftello64(fDNA1);
	uint64_t fDNA2_len = ftello64(fDNA2);

	rewind(fDNA1);
	rewind(fDNA2);

	char * DNA1 = calloc(fDNA1_len + 1, 1);
	char * DNA2 = calloc(fDNA2_len + 1, 1);

	if ((DNA1 == NULL) || (DNA2 == NULL)) {
		printf("Error: could not allocate memory line %d file %s\n", __LINE__, __FILE__);
		exit(-1);
	}

	fread(DNA1, 1, fDNA1_len, fDNA1);
	fread(DNA2, 1, fDNA2_len, fDNA2);

	fclose(fDNA1);
	fclose(fDNA2);

	FILE * fout = stdout;

	if (output_file_idx >= 0) {
		fout = fopen64(argv[output_file_idx], "w");

		if (!fout) {
			printf("Could not open output file\n");
			exit(-1);
		}
	}

	uint64_t DNA1_LEN = sanitize_input(DNA1);
	uint64_t DNA2_LEN = sanitize_input(DNA2);

	// printf("%d   %d\n", DNA1_LEN, DNA2_LEN);

	arr_type *arr = calloc((DNA1_LEN + 1), (DNA2_LEN + 1) * 4ULL);
	char *fin = calloc(sizeof(char) * (DNA1_LEN + DNA2_LEN + 2), 1);

	if ((arr == NULL) || (fin == NULL)) {
		printf("Error: could not allocate memory line %d file %s\n", __LINE__, __FILE__);
		exit(-1);
	}

	clock_t c1;

	c1 = clock();

	if (alg == NEEDLEMAN_WUNSCH) {
		nw_gen_arr(DNA1, DNA2, DNA1_LEN, DNA2_LEN, arr, GAP_PEN, MIS_PEN, NOO_PEN);
	} else if (alg == SMITH_WATERMAN) {
		sw_gen_arr(DNA1, DNA2, DNA1_LEN, DNA2_LEN, arr, GAP_PEN, MIS_PEN, NOO_PEN);
	}

	// print_2d_arr(arr, DNA1_LEN + 1, DNA2_LEN + 1, stdout);

	if (alg == NEEDLEMAN_WUNSCH) {
		uint64_t sim = 0;
		uint64_t fin_len = 0;

		nw_gen_fin_arr(DNA1, DNA2, DNA1_LEN, DNA2_LEN, arr, fin, &fin_len, &sim, GAP_PEN, MIS_PEN, NOO_PEN);

		fprintf(fout, "Similarity percentage %lf%%\n", (double)sim / (double)DNA1_LEN * 100);

		// uint64_t t1 = 0;
		// uint64_t t2 = 0;

		// uint64_t print_per_screen = 80;
		// for (uint64_t i = 0; i < fin_len; i += print_per_screen) {
		// 	print_similarities(i, print_per_screen, DNA1, DNA2, fin, fin_len, &t1, &t2, fout);
		// 	putc('\n', fout);
		// }

		print_differences(DNA1, DNA2, fin, fin_len, 0, 0, fout);
	} else if (alg == SMITH_WATERMAN) {
		arr_type max_val = 0;
		find_max_val(arr, DNA1_LEN + 1, DNA2_LEN + 1, &max_val);

		if (max_val == 0) {
			fprintf(fout, "No local alignments found\n");
			exit(-1);
		}
		
		for (uint64_t i = 0; i < DNA1_LEN + 1; i++) {
			for (uint64_t j = 0; j < DNA2_LEN + 1; j++) {
				if (arr[(i * (DNA2_LEN + 1)) + j] == max_val) {
					uint64_t sim = 0;
					uint64_t total = 0;
					uint64_t fin_len = 0;
					uint64_t fin_x = 0;
					uint64_t fin_y = 0;

					sw_gen_fin_arr(DNA1, DNA2, DNA1_LEN, DNA2_LEN, arr, fin, &fin_len, &sim, &total, i, j, &fin_x, &fin_y, GAP_PEN, MIS_PEN, NOO_PEN);

					fprintf(fout, "Local alignment found in index %"PRId64" in reference and %"PRId64" in compared, length is %"PRId64"\n", fin_x, fin_y, i - fin_x);
					// printf("%d    %d    %d\n", i, j, arr[(i * (DNA2_LEN + 1)) + j]);
					fprintf(fout, "Similarity percentage %lf%%\n", (double)sim / (double)total * 100);

					// uint64_t t1 = 0;
					// uint64_t t2 = 0;

					// uint64_t print_per_screen = 80;
					// for (uint64_t i = 0; i < fin_len; i += print_per_screen) {
					// 	print_similarities(i, print_per_screen, DNA1, DNA2, fin, fin_len, &t1, &t2, fout);
					// 	putc('\n', fout);
					// }

					print_differences(DNA1, DNA2, fin, fin_len, fin_x, fin_y, fout);
				}
			}
		}
	}

	fprintf(fout, "Elapsed time: %lf\n", (double)(clock() - c1) / CLOCKS_PER_SEC);

	free(arr);
	fclose(fout);
	free(fin);

	return 0;
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

	printf("%06"PRId64"\t", *tmp1);

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

	printf("\t");
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

	printf("%06"PRId64"\t", *tmp2);
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

void print_differences(const char * DNA1, const char * DNA2, char * fin, uint64_t fin_len, uint64_t x, uint64_t y, FILE * fp) {
	// uint64_t x = 0;
	// uint64_t y = 0;
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

void nw_gen_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN) {
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

void sw_gen_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN) {
	for (uint64_t i = 0; i < DNA1_LEN + 1; i++) {
		arr[i * (DNA2_LEN + 1) + 0] = 0;
	}

	for (uint64_t i = 0; i < DNA2_LEN + 1; i++) {
		arr[0 * (DNA2_LEN + 1) + i] = 0;
	}

	for (uint64_t i = 1; i < DNA1_LEN + 1; i++) {
		for (uint64_t j = 1; j < DNA2_LEN + 1; j++) {
			arr[i * (DNA2_LEN + 1) + j] = max4(arr[(i - 1) * (DNA2_LEN + 1) + (j - 1)] + ((DNA1[i - 1] == DNA2[j - 1]) ? NOO_PEN : MIS_PEN), 
											   arr[i * (DNA2_LEN + 1) + (j - 1)] + GAP_PEN, 
											   arr[(i - 1) * (DNA2_LEN + 1) + j] + GAP_PEN, 
											   0);
		}
	}
}

void nw_gen_fin_arr(char * DNA1, char * DNA2, uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, char * fin, uint64_t * fin_len, uint64_t * sim, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN) {
	uint64_t x = DNA1_LEN;
	uint64_t y = DNA2_LEN;

	while(1) {
		if ((x == 0) && (y == 0)) {
			break;
		}

		if ((x > 0) && (y > 0) && (arr[x * (DNA2_LEN + 1) + y] == 
		    (arr[(x - 1) * (DNA2_LEN + 1) + (y - 1)] + ((DNA1[x - 1] == DNA2[y - 1]) ? NOO_PEN : MIS_PEN)))) {
			fin[(*fin_len)++] = (DNA1[x - 1] == DNA2[y - 1]) ? EQUAL : MISMATCH;
			if (DNA1[x - 1] == DNA2[y - 1]) {
				(*sim)++;
			}
			x--; y--;
		} else if ((x > 0) && (arr[x * (DNA2_LEN + 1) + y] == (arr[(x - 1) * (DNA2_LEN + 1) + y] + GAP_PEN))) {
			fin[(*fin_len)++] = SECOND_GAP;
			x--;
		} else if ((y > 0) && (arr[x * (DNA2_LEN + 1) + y] == (arr[x * (DNA2_LEN + 1) + (y - 1)] + GAP_PEN))) {
			fin[(*fin_len)++] = FIRST_GAP;
			y--;
		} else {
			printf("Error encountered at line %d file %s\n", __LINE__, __FILE__);
			exit(-1);
		}
	}
	
	reverse_char_arr(fin, *fin_len);
}

void sw_gen_fin_arr(char * DNA1, char * DNA2, UNUSED uint64_t DNA1_LEN, uint64_t DNA2_LEN, arr_type * arr, char * fin, uint64_t * fin_len, uint64_t * sim, uint64_t * total, uint64_t start_x, uint64_t start_y, uint64_t * fin_x, uint64_t * fin_y, arr_type GAP_PEN, arr_type MIS_PEN, arr_type NOO_PEN) {
	uint64_t x = start_x;
	uint64_t y = start_y;

	while(1) {
		if (arr[x * (DNA2_LEN + 1) + y] == 0) {
			break;
		}

		if ((x > 0) && (y > 0) && (arr[x * (DNA2_LEN + 1) + y] == 
		    (arr[(x - 1) * (DNA2_LEN + 1) + (y - 1)] + ((DNA1[x - 1] == DNA2[y - 1]) ? NOO_PEN : MIS_PEN)))) {
			fin[(*fin_len)++] = (DNA1[x - 1] == DNA2[y - 1]) ? EQUAL : MISMATCH;
			if (DNA1[x - 1] == DNA2[y - 1]) {
				(*sim)++;
			}
			(*total)++;
			x--; y--;
		} else if ((x > 0) && (arr[x * (DNA2_LEN + 1) + y] == (arr[(x - 1) * (DNA2_LEN + 1) + y] + GAP_PEN))) {
			fin[(*fin_len)++] = SECOND_GAP;
			x--;
		} else if ((y > 0) && (arr[x * (DNA2_LEN + 1) + y] == (arr[x * (DNA2_LEN + 1) + (y - 1)] + GAP_PEN))) {
			fin[(*fin_len)++] = FIRST_GAP;
			y--;
		} else {
			printf("Error encountered at line %d file %s\n", __LINE__, __FILE__);
			exit(-1);
		}
	}

	*fin_x = x;
	*fin_y = y;
	
	reverse_char_arr(fin, *fin_len);
}

void print_2d_arr(arr_type * arr, uint64_t dim1, uint64_t dim2, FILE * fp) {
	for (uint64_t i = 0; i < dim1; i++) {
		for (uint64_t j = 0; j < dim2; j++) {
			fprintf(fp, "%04d  ", arr[(i * dim2) + j]);
		}
		fputc('\n', fp);
	}
}

void shift_arr_left(char * arr, uint64_t start, uint64_t len, uint64_t shift) {
	for (uint64_t i = start; i < len - shift; i++) {
		arr[i] = arr[i + shift];
	}
}

void find_max_val(arr_type * arr, uint64_t dim1, uint64_t dim2, arr_type * max_val) {
	for (uint64_t i = 0; i < dim1; i++) {
		for (uint64_t j = 0; j < dim2; j++) {
			if (arr[(i * dim2) + j] > (*max_val)) {
				*max_val = arr[(i * dim2) + j];
			}
		}
	}
}

uint64_t sanitize_input(char * input) {
	uint64_t len = strlen(input);
	uint64_t i;

	for (i = 0; i < len; i++)
		if (input[i] >= 'a')
			input[i] -= 0x20;

	for (i = 0; i < len; i++) {
		if (!IS_BASIC_DNA_CHAR(input[i]) && !IS_EXTRA_DNA_CHAR(input[i])) {
			shift_arr_left(input, i, len, 1);
			i--;
			len--;
		}
	}

	return len;
}

void print_usage(char * filename) {
	printf("DNA Sequence alignment\n");
	printf("Usage: %s [OPTIONS]\n", filename);
	printf("OPTIONS:\n");
	printf("\t-h\t\t\tPrint this help message\n");
	printf("\t-o FILE\t\t\tSpecify output file. default: STDOUT\n");
	printf("\t-S1\t\t\tSpecify file as reference sequence. [default: DNA1_file.txt]\n");
	printf("\t-S2\t\t\tSpecify file as other sequence. [default: DNA2_file.txt]\n");
	printf("\t-Px N\t\t\tSpecify parameter for algorithm, N is the value:\n");
	printf("\t\tPossible values of x:\n");
	printf("\t\t\tg\tGap penalty. [default: -2]\n");
	printf("\t\t\tm\tMismatch penalty. [default: -1]\n");
	printf("\t\t\tn\tMatch penalty. [default: 1]\n");
	printf("\t-Ax\t\t\tSpecify desired algorithm\n");
	printf("\t\tPossible values of n:\n");
	printf("\t\t\tn\tNeedleman-Wunsch (Global alignment)\t\t[Default]\n");
	printf("\t\t\ts\tSmith-Waterman (Local alignment)\n");
}
