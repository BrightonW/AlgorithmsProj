/* DNA sequence mapper:
 *
 * Skeleton code written by Jianzhong Qi, April 2022
 * Edited by: Brighton Wankeaw 1263276
 * 
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#define STAGE_NUM_ONE 1						  /* stage numbers */ 
#define STAGE_NUM_TWO 2
#define STAGE_NUM_THREE 3
#define STAGE_NUM_FOUR 4
#define STAGE_NUM_FIVE 5
#define STAGE_HEADER "Stage %d\n==========\n" /* stage header format string */

#define MAX_READ_ID_LENGTH 100				  /* maximum read ID length */
#define MAX_READ_LENGTH 100					  /* maximum read length */
#define MAX_NUM_READS 100					  /* maximum number of reads */
#define MAX_REF_LENGTH 1000					  /* maximum reference DNA length */

typedef char read_id_t[MAX_READ_ID_LENGTH+1]; /* a read ID */
typedef char read_t[MAX_READ_LENGTH+1];		  /* a read */
typedef char score_t[MAX_READ_LENGTH+1];	  /* quality scores of a read */
typedef char ref_t[MAX_REF_LENGTH+1];		  /* a reference DNA sequence */

/****************************************************************/

/* function prototypes */
void take_one_read(read_t one_read, score_t scores_of_one_read);
void print_stage_header(int stage_num);
void take_reads(read_t reads[], score_t scores[], int *num_reads);
void ref_reads(ref_t ref_sequence);

void stage_one(read_t one_read, score_t scores_of_one_read);
void stage_two(read_t reads[], score_t scores[], int *num_reads);
void stage_three(read_t reads[], score_t scores[], int num_reads);
void stage_four(ref_t ref_sequence);
void stage_five(read_t reads[], score_t scores[], int num_reads, 
	ref_t ref_sequence);

/* add your own function prototypes here */
/****************************************************************/

/*algorithms are fun*/

/* main function controls all the action, do NOT modify this function */
int
main(int argc, char *argv[]) {
	/* to hold all input reads and quality scores */
	read_t reads[MAX_NUM_READS];	
	score_t scores[MAX_NUM_READS];
	/* to hold the number of input reads */
	int num_reads = 0;	
	/* to hold the input reference sequence */
	ref_t ref_sequence;
	
	/* stage 1: process one read */
	stage_one(reads[0], scores[0]); 
	num_reads++;
	
	/* stage 2: process all reads */
	stage_two(reads, scores, &num_reads);
	
	/* stage 3: mask bases with high error probability */ 
	stage_three(reads, scores, num_reads);
	
	/* stage 4: process reference sequence */
	stage_four(ref_sequence);
	
	/* stage 5: map reads to the reference sequence */
	stage_five(reads, scores, num_reads, ref_sequence);
	
	/* all done; take some rest */
	return 0;
}

/* print stage header given stage number */
void 
print_stage_header(int stage_num) {
	printf(STAGE_HEADER, stage_num);
}

/****************************************************************/
/* add your code below */

/* process a read record */
void 
take_one_read(read_t one_read, score_t scores_of_one_read) {
	/* modify this function for stage 2 */
	read_id_t id;
	scanf("%s", id);
	scanf("%s", one_read);
	getchar();
	getchar(); // skip '+' and '\n'
	scanf("%s", scores_of_one_read);
}

void
take_reads(read_t reads[], score_t scores[], int *num_reads) {
    int c;
    while ((c=getchar())!=EOF){
        read_id_t id;
        scanf("%s", id);
        if (!strcmp(id,"#####")){
            break;
        }
        scanf("%s", reads[*num_reads]);
        getchar();
        getchar(); // skip '+' and '\n'
        scanf("%s", scores[*num_reads]);
        (*num_reads)++;
    }
}

void
ref_reads(ref_t ref_sequence){
    scanf("%s", ref_sequence);
}
    
/* stage 1: process one read */
void 
stage_one(read_t one_read, score_t scores_of_one_read) {
	/* print stage header */
	print_stage_header(STAGE_NUM_ONE);
	
	/* add code for stage 1 */
	take_one_read(one_read, scores_of_one_read);
    int length, smallest_qual, stage1_index;
	char smallest_base;
	length = strlen(scores_of_one_read);
	smallest_qual = scores_of_one_read[0];
	for (int i = 1; i<length; i++){
		if (scores_of_one_read[i] < smallest_qual){
			smallest_qual = scores_of_one_read[i];
            smallest_base = one_read[i];
			stage1_index = i;
		}
	}
	/* process first read */
    printf("Base with the smallest quality score: %c\n", smallest_base);
    printf("Index: %d\n", stage1_index);
    printf("\n");
    }


/* stage 2: process all reads */
void 
stage_two(read_t reads[], score_t scores[], int *num_reads) {
    print_stage_header(STAGE_NUM_TWO);
	/* add code for stage 2 */
    take_reads(reads, scores, num_reads);
    int qual_index;
    double lowest_qual = 0;
    for (int i = 0; i<*num_reads; i++){
        int ascii_sum = 0;
        double char_count = 0.0;
        int score_len;
        score_len = strlen(scores[i]);
        for (int j = 0; j < score_len; j++){
            ascii_sum += scores[i][j];
            char_count++;
        }
        double avg = ascii_sum/char_count;
        if (lowest_qual == 0){
           lowest_qual = avg;
        }    
        else {
            if (lowest_qual > avg) {
                lowest_qual = avg;
                qual_index = i;
                }
  
        }
    }
    printf("Total number of reads: %d\n", *num_reads);
    printf("Smallest average quality score: %.2f\n", lowest_qual);
    printf("Read with the smallest average quality score: \n%s\n",
    reads[qual_index]);
    printf("\n");
}
    

/* stage 3: mask bases with high error probability */ 
void 
stage_three(read_t reads[], score_t scores[], int num_reads) {
	/* add code for stage 3 */
    print_stage_header(STAGE_NUM_THREE);
    double q_cutoff = 119/3.0;
	for (int i = 0; i < num_reads; i++){
            int score_len;
            score_len = strlen(scores[i]);
            for (int j = 0; j < score_len; j++){
                if (scores[i][j] < q_cutoff){
                    reads[i][j] = '*';
                }
		    }
		printf("%s\n", reads[i]);
	}
    printf("\n");
}

/* stage 4: process reference sequence */
void stage_four(ref_t ref_sequence) {
    print_stage_header(STAGE_NUM_FOUR);
	/* add code for stage 4 */
    ref_reads(ref_sequence);
    int A = 0, C = 0, G = 0, T = 0;
	int length;
	length = strlen(ref_sequence);
	for (int i = 0; i < length; i++){
		if (ref_sequence[i] == 'A'){
			A++;
		}
		else if (ref_sequence[i] == 'C'){
			C++;
		}
		else if (ref_sequence[i] == 'G'){
			G++;
		}
		else if (ref_sequence[i] == 'T'){
			T++;
		}
	}
	printf("Length of the reference sequence: %d\n", length);
	printf("Number of A bases: %d\n", A);
	printf("Number of G bases: %d\n", C);
	printf("Number of C bases: %d\n", G);
	printf("Number of T bases: %d\n", T);
    printf("\n");
	
}

/* stage 5: map reads to the reference sequence */
void 
stage_five(read_t reads[], score_t scores[], int num_reads, 
	ref_t ref_sequence) {
    print_stage_header(STAGE_NUM_FIVE);
	/* add code for stage 5  */
    for (int i = 0; i<num_reads; i++){
        double match_score = 0;
        int length = strlen(reads[i]);
        printf("length %d\n", length);
        char match[length];
        int end_of_ref = strlen(ref_sequence) - length;
        for (int start = 0; start < end_of_ref; start++){
            char str[length];
            int curr_score = 0;
            for (int j = start; j<(length + start); j++){
                if (reads[i][j] == ref_sequence[j]){
                    curr_score += log2(0.25);
                }
                else if (reads[i][j] == '*'){
                    curr_score += log2(0.25);
                }
                else  {
                    curr_score += log2(1);
                };
                str[j] = ref_sequence[j]; 
            }
            if (match_score < curr_score){
                match_score = curr_score;
                strcpy(match, str);
            }
        }
        printf("Reads: %s\n", reads[i]);
        printf("Match: %s\n", match);
    }
}
