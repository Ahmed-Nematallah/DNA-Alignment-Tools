
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#define UNUSED __attribute__ ((unused))

typedef enum {
	SECOND_GAP,
	FIRST_GAP,
	EQUAL,
	MISMATCH
} FIN_CASES;

int max3(int a, int b, int c) {
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

int main (UNUSED int argc, UNUSED char ** argv) {
	const char * DNA1 = "STREETREEDTREE";
	const char * DNA2 = "TREEREEDSTREET";

	// printf("%d\n", sizeof(unsigned long long));

	unsigned long long DNA1_LEN = strlen(DNA1);
	unsigned long long DNA2_LEN = strlen(DNA2);

	int GAP_PEN = -2;
	int MIS_PEN = -1;
	int NOO_PEN = 1;

	int *arr = calloc((DNA1_LEN + 1), (DNA2_LEN + 1) * 4ULL);
	char *fin = malloc(sizeof(char) * (DNA1_LEN + DNA2_LEN + 2));
	unsigned long long fin_len = 0;

	if ((arr == NULL) || (fin == NULL)) {
		printf("Error: could not allocate memory\n");
		exit(-1);
	}

	clock_t c1;

	c1 = clock();


	for (unsigned long long i = 0; i < DNA1_LEN + 1; i++) {
		arr[i * (DNA2_LEN + 1) + 0] = GAP_PEN * i;
	}

	for (unsigned long long i = 0; i < DNA2_LEN + 1; i++) {
		arr[0 * (DNA2_LEN + 1) + i] = GAP_PEN * i;
	}

	for (unsigned long long i = 1; i < DNA1_LEN + 1; i++) {
		for (unsigned long long j = 1; j < DNA2_LEN + 1; j++) {
			arr[i * (DNA2_LEN + 1) + j] = max3(arr[(i - 1) * (DNA2_LEN + 1) + (j - 1)] + ((DNA1[i - 1] == DNA2[j - 1]) ? NOO_PEN : MIS_PEN), 
			                                   arr[i * (DNA2_LEN + 1) + (j - 1)] + GAP_PEN, 
											   arr[(i - 1) * (DNA2_LEN + 1) + j] + GAP_PEN);
		}
	}

	// for (int i = 0; i < DNA1_LEN + 1; i++) {
	// 	for (int j = 0; j < DNA2_LEN + 1; j++) {
	// 		printf("%06d  ", arr[i * (DNA2_LEN + 1) + j]);
	// 	}
	// 	putc('\n', stdout);
	// }

	unsigned long long x = DNA1_LEN;
	unsigned long long y = DNA2_LEN;

	while(1) {
		// printf("%d  %d\n", x, y);

		if ((x == 0) && (y == 0)) {
			break;
		}

		if (arr[x * (DNA2_LEN + 1) + y] == 
		    (arr[(x - 1) * (DNA2_LEN + 1) + (y - 1)] + ((DNA1[x - 1] == DNA2[y - 1]) ? NOO_PEN : MIS_PEN))) {
			fin[fin_len++] = (DNA1[x - 1] == DNA2[y - 1]) ? EQUAL : MISMATCH;
			x--; y--;
		} else if (arr[x * (DNA2_LEN + 1) + y] == (arr[(x - 1) * (DNA2_LEN + 1) + y] + GAP_PEN)) {
			fin[fin_len++] = FIRST_GAP;
			x--;
		} else if (arr[x * (DNA2_LEN + 1) + y] == (arr[x * (DNA2_LEN + 1) + (y - 1)] + GAP_PEN)) {
			fin[fin_len++] = SECOND_GAP;
			y--;
		} else {
			exit(-1);
		}
	}

	printf("Elapsed time: %lf\n", (double)(clock() - c1) / CLOCKS_PER_SEC);

	unsigned long long j = 0;
	for (unsigned long long i = 0; i < fin_len; i++) {
		switch (fin[i]) {
			case EQUAL:
				putc(DNA1[j++], stdout);
				break;
			case MISMATCH:
				putc(DNA1[j++], stdout);
				break;
			case FIRST_GAP:
				putc('-', stdout);
				break;
			case SECOND_GAP:
				putc(DNA1[j++], stdout);
				break;
			
			default:
				break;
		}
	}

	putc('\n', stdout);

	for (unsigned long long i = 0; i < fin_len; i++) {
		switch (fin[i]) {
			case EQUAL:
				putc('|', stdout);
				break;
			
			default:
				putc(' ', stdout);
				break;
		}
	}

	putc('\n', stdout);

	j = 0;
	for (unsigned long long i = 0; i < fin_len; i++) {
		switch (fin[i]) {
			case EQUAL:
				putc(DNA2[j++], stdout);
				break;
			case MISMATCH:
				putc(DNA2[j++], stdout);
				break;
			case FIRST_GAP:
				putc(DNA2[j++], stdout);
				break;
			case SECOND_GAP:
				putc('-', stdout);
				break;
			
			default:
				break;
		}
	}

	free(arr);
	free(fin);

	return 0;
}
