/*
2021/4/4 debug end

Program for solve the traveling salesman problem by genetic algorithm
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <time.h>

#define		DLB_MAIN		"GA"		//data label
#define		NUM_CITY		12
#define		NUM_GENE		200

void main_GA();
void final_process();
void output_path_best_all();
void output_path_GA(long g);
void output_lng_GA(long g, double lng_best, double lng_worst, double lng_mid);
void mutate_gene();
void mutate_gene_each(long i, long j);
void kousa_gene();
void kousa_gene_pair(long p1[], long p2[]);
long search_vec(long v[], long n);
void sort_gene();
void koukan_gene(long i, long j);
void koukan_long(long *x, long *y);
double dist_gene(long i);
double dist_gene_each(long i, long j, long k);
void mk_city();
void init_gene();
void init_gene_each(long i, long j);
long check_same(long i, long j);
void my_fopen(FILE **fp, char fname[]);
void check_para_main();
void disp_main();
void c_lng_bwm(double *lng_best, double *lng_worst, double *lng_mid);
void c_gene_best_all(long g, double lng_best);
void select_gene();
void copy_gene(long i, long j);


//external variable
long Num_select, Num_kousa; 
double P_mut;
long G_fin;
long G_disp_GA, Nout_lng, Nout_path;

long Gene[NUM_GENE][NUM_CITY];

double City[NUM_CITY][2];		//x&y coordinates

double Gene_best_all[NUM_CITY], Lng_best_all, Lng_best_fin, Lng_worst_fin, Lng_mid_fin;

FILE *Fp_lng;
FILE *Fp_deb;

int main(void){
	if ( (Fp_deb = fopen("deb.dat", "w") ) == NULL ){ printf("err in fopen for Fp_deb"); exit(1); }
	
	Num_select  = NUM_GENE/2;
	Num_kousa	= NUM_GENE/4;
	P_mut  		= 5.e-2;			//probability of mutation				
	G_fin		= 500;
	G_disp_GA		= G_fin/10;
	Nout_lng		= 5;			//output frequency
	Nout_path		= G_fin/2;
	
	mk_city();
	main_GA();
	disp_main();
}

//screen output
void disp_main(){
	printf("Lng_best_all = %.7lf\n", Lng_best_all);
	printf("Lng_best_fin = %.7lf\n", Lng_best_fin);
	printf("Lng_worst_fin = %.7lf\n", Lng_worst_fin);
	printf("Lng_mid_fin %.7lf\n", Lng_mid_fin);
}

void main_GA(){
	double lng_best, lng_worst, lng_mid;
	long g;
	
	init_gene();
	
	for (g = 0; ; g++){
		if (g % G_disp_GA == 0) printf("g = %ld in main_GA\n", g);
		
		sort_gene();
		c_lng_bwm(&lng_best, &lng_worst, &lng_mid);
		c_gene_best_all(g, lng_best);
		
		if (g % Nout_lng == 0) output_lng_GA(g, lng_best, lng_worst, lng_mid);
		if (g % Nout_path == 0) output_path_GA(g);
		if (g == G_fin){ final_process(); return; }
		
		select_gene();
		sort_gene();
		kousa_gene();
		mutate_gene();
	}
}

//best gene
void c_gene_best_all(long g, double lng_best){
	long i;
	
	if (g == 0 || Lng_best_all > lng_best){
		Lng_best_all = lng_best;
		for (i = 0; i < NUM_CITY; i++) Gene_best_all[i] = Gene[0][i];
	}
}

void final_process(){
	c_lng_bwm(&Lng_best_fin, &Lng_worst_fin, &Lng_mid_fin);
	
	output_path_best_all();	
	
	fclose(Fp_lng);
}

//output the best path
void output_path_best_all(){
	FILE *fp;
	long id;
	long j;
	my_fopen(&fp, "path_best_of_all");
	
	fprintf(fp, "No_of_city, x_of_path(best_of_all)(%s), y_of_path(best_of_all)(%s)\n", DLB_MAIN, DLB_MAIN);
	
	for (j = 0; j <= NUM_CITY; j++){
		if (j != NUM_CITY) id = Gene_best_all[j];
		else id = Gene_best_all[0];
		
		fprintf(fp, "%ld, %.7le, %.7le\n", j, City[id][0], City[id][1]);
	}
	fclose(fp);
}

//output path
void output_path_GA(long g){
	char text1[200], text2[200];
	long id_b, id_w, id_m;
	long j;
	FILE *fp;
	
	sprintf(text1, "path_GA(g=%ld)", g);
	my_fopen(&fp, text1);
	
	sprintf(text2, "(g=%ld)(%s)", g, DLB_MAIN);
	fprintf(fp, "No_of_city, x_of_path(best)%s, y_of_path(best)%s, x_of_path(worst)%s, y_of_path(worst)%s, x_of_path(mid)%s, y_of_path(mid)%s\n",
		text2, text2, text2, text2, text2, text2);
	
	for (j = 0; j <= NUM_CITY; j++){
		
		if (j != NUM_CITY){
			id_b = Gene[0][j];
			id_w = Gene[NUM_GENE-1][j];
			id_m = Gene[NUM_GENE/2][j];
		}
		else{
			id_b = Gene[0][0];
			id_w = Gene[NUM_GENE-1][0];
			id_m = Gene[NUM_GENE/2][0];
		}
		
		fprintf(fp, "%ld, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le\n", 
				j, City[id_b][0], City[id_b][1], City[id_w][0], City[id_w][1], City[id_m][0], City[id_m][1]);
	}
	
	fclose(fp);
}

//output path length
void output_lng_GA(long g, double lng_best, double lng_worst, double lng_mid){
	if (g == 0){
		my_fopen(&Fp_lng, "lng_GA");
		fprintf(Fp_lng, "g, lng_best(%s), lng_worst(%s), lng_mid(%s)\n", DLB_MAIN, DLB_MAIN, DLB_MAIN);
	}
	
	fprintf(Fp_lng, "%ld, %.7f, %.7f, %.7f\n", g, lng_best, lng_worst, lng_mid);
}

void c_lng_bwm(double *lng_best, double *lng_worst, double *lng_mid){
	if (NUM_GENE%2 != 0){ printf("err in c-f-bwm\n"); exit(1); }
	
	*lng_best = dist_gene(0);
	*lng_worst = dist_gene(NUM_GENE-1);
	
	*lng_mid = (dist_gene(NUM_GENE/2-1) + dist_gene(NUM_GENE/2)) / 2.e0;
}

//mutation
void mutate_gene(){
	double u;
	long i, j;
	for (i = 0; i < NUM_GENE; i++){
		for (j = 0; j < NUM_CITY; j++){
			u = 1.e0 * rand()/RAND_MAX;
			if (u < P_mut) mutate_gene_each(i, j);
		}
	}
}

void mutate_gene_each(long i, long j){
	long n, m, k;
	n = Gene[i][j];
	for (;;){
		m = rand() % NUM_CITY;
		if (m != n) break;
	}
	
	Gene[i][j] = m;
	
	k = search_vec(Gene[i], m);
	Gene[i][k] = n;
}

//select
void select_gene(){
	long i;
	for (i = Num_select; i < NUM_GENE; i++){
		copy_gene(i, i%Num_select);	
	}
}

//copy
void copy_gene(long i, long j){
	long k;
	for (k = 0; k < NUM_CITY; k++){
		Gene[i][k] = Gene[j][k];
	}
}

//crossover
void kousa_gene(){
	long p1[NUM_GENE], p2[NUM_GENE];
	long i, j;
	if (Num_kousa%2 != 0){ printf("err: kousa-gene\n"); exit(1); }
	
	for (i = 0; i < Num_kousa/2; i++){
		
		for (j = 0; j < NUM_CITY; j++){
			p1[j] = Gene[2*i][j]; 
			p2[j] = Gene[2*i+1][j];
		}
		
		kousa_gene_pair(p1, p2);
		
		for (j = 0; j < NUM_CITY; j++){
			Gene[2*i][j] = p1[j];
			Gene[2*i+1][j] = p2[j];
		}
	}
}

void kousa_gene_pair(long p1[], long p2[]){
	long mask[NUM_CITY];
	long num_mask = 0;
	long k, n;
	
	n = 0;
	for (;;){
		mask[num_mask] = n; num_mask++;
		n = search_vec(p1, p2[n]);
		if (n == 0) break;
	}
	
	for (k = 0; k < num_mask; k++){
		n = mask[k];
		koukan_long(&(p1[n]), &(p2[n]));
	}
}

//search location
long search_vec(long v[], long n){
	long i;
	for (i = 0; i < NUM_CITY; i++){
		if (v[i] == n) return i;
	}
	printf("err in search-vec\n");
	exit(1);
}

//sort
void sort_gene(){
	long i, j;
	for (i = 0; i < NUM_GENE; i++){
		for (j = i+1; j < NUM_GENE; j++){
			if (dist_gene(i) > dist_gene(j)) koukan_gene(i, j);
		}
	}
}

void koukan_gene(long i, long j){
	long k;
	for (k = 0; k < NUM_CITY; k++){
		koukan_long(&(Gene[i][k]), &(Gene[j][k]));
	}
}

void koukan_long(long *x, long *y){
	long nwk1;
	nwk1 = *x;
	*x = *y;
	*y = nwk1;
}

//calculate the total distance of the path
double dist_gene(long i){
	long j;
	double sum;
	sum = 0.e0;
	for (j = 0; j < NUM_CITY-1; j++){
		sum += dist_gene_each(i, j, j+1);
	}
	sum += dist_gene_each(i, NUM_CITY-1, 0);
	return sum;
}

//calculate distance
double dist_gene_each(long i, long j, long k){
	long n1, n2;
	double dx, dy;
	n1 = Gene[i][j];
	n2 = Gene[i][k];
	
	dx = City[n1][0] - City[n2][0];
	dy = City[n1][1] - City[n2][1];
	
	return sqrt(dx*dx + dy*dy);
}

//city creation
void mk_city(){
	long i;
	FILE *fp;
	for (i = 0; i < NUM_CITY; i++){
		City[i][0] = 1.e0 * rand() / RAND_MAX;
		City[i][1] = 1.e0 * rand() / RAND_MAX;
	}
	
	my_fopen(&fp, "city");
	fprintf(fp, "i, City[i][0], City[i][1]\n");
	for (i = 0; i < NUM_CITY; i++){
		fprintf(fp, "%ld, %.7le, %.7le\n", i, City[i][0], City[i][1]);
	}
	fclose(fp);
}

void init_gene(){
	long i, j;
	for (i = 0; i < NUM_GENE; i++){
		for (j = 0; j < NUM_CITY; j++){
			init_gene_each(i, j);
		}
	}
}

//input the value of Gene[i][j]
void init_gene_each(long i, long j){
	long n_big = 10000;
	long k;
	for (k = 0; k < n_big; k++){
		Gene[i][j] = rand()%NUM_CITY;
		if (check_same(i, j) == 0) return;
	}
	printf("err in init-gene-each\n"); exit(1);
}

//check for city duplication in gene sequences
long check_same(long i, long j){
	long k;
	for (k = 0; k < j; k++){
		if (Gene[i][j] == Gene[i][k]) return 1;
	}
	return 0;
}

void my_fopen(FILE **fp, char fname[]){
	char text_wk1[200];
	
	sprintf(text_wk1, "out\\%s(%s).dat", fname, DLB_MAIN);
	
	if ( (*fp = fopen(text_wk1, "w") ) == NULL ){ 
		printf("err in fopcl for write mode open %s", text_wk1); 
		exit(1);
	}
}