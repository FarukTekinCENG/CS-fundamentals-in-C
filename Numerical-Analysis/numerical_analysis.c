#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define COL 10
#define MAX 1000
#define PI 3.14159265359
void menu(double matrix[][COL], float matrix2[][COL], float islemMatrisi[][COL], float birimMatris[][COL], char function[]);
int getFunct(double matrix[][COL], char function[]);
double calcFunct(double matrix[][COL], int h, double x);
int getMatrix(float matrix2[][COL]);
void bisection(double matrix[][COL], int h);
void regulaFalsi(double matrix[][COL], int h);
void newtonRaphson(double matrix[][COL], int h);
void inverseMatrix(float matrix2[][COL], float islemMatrisi[][COL], float birimMatris[][COL], int N);
void gaussElimination();
void gaussSeidel();
void rowExchanger(float matrix[][COL], float b[], int N, int row1, int row2);
double numericalDifferentiation(double matrix[][COL], int t, double xi);
void numericalDiff(double matrix[][COL], int t);
void simpsonsRule(double matrix[][COL], int t);
void trapezoidalRule(double matrix[][COL], int t);
void gregoryNewtonInterpolation();
//******************************************************************************************************************************
int main(){
	double matrix[6][COL]={{0}};
	float matrix2[COL][COL]={{0}};
	float islemMatrisi[COL][COL]={{0}};
	float birimMatris[COL][COL]={{0}};
	char function[MAX]="";
	menu(matrix, matrix2, islemMatrisi, birimMatris, function);
	return 0;
}
//******************************************************************************************************************************
void menu(double matrix[][COL], float matrix2[][COL], float islemMatrisi[][COL], float birimMatris[][COL], char function[]){
	int option;
	int h, N;
	do{
	printf("-1: Define a Function\n1:  Bisection\n2:  Regula-Falsi\n3:  Newton-Rapshon\n4:  Inverse Matrix\n5:  Gauss Elimination\n");
	printf("6:  Gauss-Seidal\n7:  Numerical Differentiation\n8:  Simpson's Rule\n9:  Trapezoidal Rule\n10: Gregory Newton\n0:  Quit\n");
	printf("Choice: ");
	scanf("%d",&option);
	printf("\n\n\n");
	switch(option){
		case -1:
			h=getFunct(matrix, function);
			break;
		case 0:
			printf("Quitting...");
			break;
		case 1:
			printf("Function:\n%s\n", function);
			bisection(matrix, h);
			break;
		case 2:
			printf("Function:\n%s\n", function);
			regulaFalsi(matrix, h);
			break;
		case 3:
			printf("Function:\n%s\n", function);
			newtonRaphson(matrix, h);
			break;
		case 4:
			N=getMatrix(matrix2);
			inverseMatrix(matrix2, islemMatrisi, birimMatris, N);
			break;
		case 5:
			gaussElimination();
			break;
		case 6:
			gaussSeidel();
			break;
		case 7:
			printf("Function:\n%s\n", function);
			numericalDiff(matrix, h);
			break;
		case 8:
			printf("Function:\n%s\n", function);
			simpsonsRule(matrix, h);
			break;
		case 9:
			printf("Function:\n%s\n", function);
			trapezoidalRule(matrix, h);
			break;
		case 10:
			gregoryNewtonInterpolation();
			break;
		default:
			printf("OLMADI\n\n");		
	}
	}while(option!=0);
}
//******************************************************************************************************************************
int getFunct(double matrix[][COL], char function[]){
	int countPoli, countExpo, countLoga, countTrigo, countInverseTrigo, i, count, index=0;
	double x_coef, x_exp, fn_coef, fn_exp, base;
	//*************************************
	printf(">>Polinomial count:");
	scanf("%d",&countPoli);
	i=0;
	while((countPoli!=0) && (i<countPoli)){
		printf("\n\nPolynomial : x_coef * x ^ x_exp\n");
		printf("\nx's cofactor (x_coef) : ");
		scanf("%lf",&x_coef);
		printf("x's exponent (x_exp) : ");
		scanf("%lf",&x_exp);
		printf("\nAdded : %lf * x ^ %lf",x_coef,x_exp);
		char sign = (x_coef < 0) ? '-' : '+';
		index += sprintf(function + index, "%c%.3lf x ^%.3lf ",sign, fabs(x_coef), x_exp);
		matrix[0][i]=1;
		matrix[1][i]=x_coef;
		matrix[2][i]=x_exp;
		i++;
	}
	//*************************************
	printf("\n\n>>Exponential count:");
	scanf("%d",&countExpo);
	count=0;
	while((countExpo!=0) && (count<countExpo)){
		printf("Exponential: fn_coef * (base ^ (x_coef * x ^ x_exp)) ^ fn_exp");
		printf("\nx's cofactor (x_coef) : ");
		scanf("%lf",&x_coef);
		printf("x's exponent (x_exp) : ");
		scanf("%lf",&x_exp);
		printf("Function cofactor (fn_coef) : ");
		scanf("%lf",&fn_coef);
		printf("Function exponent (fn_exp) : ");
		scanf("%lf",&fn_exp);
		printf("Base (base) : ");
		scanf("%lf",&base);
		printf("\nAdded : %lf * (%lf ^ (%lf * x ^ %lf)) ^ %lf",fn_coef, base, x_coef, x_exp, fn_exp);
		char sign = (fn_coef < 0) ? '-' : '+';
		index += sprintf(function + index, "%c%.3lf * [%.3lf^(%.3lfx^%.3lf)] ^ %.3lf ",sign,fabs(fn_coef), base, x_coef, x_exp, fn_exp);
		matrix[0][i]=2;
		matrix[1][i]=x_coef;
		matrix[2][i]=x_exp;
		matrix[3][i]=fn_coef;
		matrix[4][i]=fn_exp;
		matrix[5][i]=base;
		i++;
		count++;
	}
	//*************************************
	printf("\n\n>>Logarithmic count:");
	scanf("%d",&countLoga);
	count=0;
	while((countLoga!=0) && (count<countLoga)){
		printf("Logarithmic: fn_coef * (log _ base ^ (x_coef * x ^ x_exp)) ^ fn_exp");
		printf("\nx's cofactor (x_coef) : ");
		scanf("%lf",&x_coef);
		printf("x's exponent (x_exp) : ");
		scanf("%lf",&x_exp);
		printf("Function cofactor (fn_coef) : ");
		scanf("%lf",&fn_coef);
		printf("Function exponent (fn_exp) : ");
		scanf("%lf",&fn_exp);
		printf("Base (base) : ");
		scanf("%lf",&base);
		printf("\nAdded : %lf * [log _ %lf _(%lf * x ^ %lf)] ^ %lf",fn_coef,base, x_coef, x_exp, fn_exp);
		char sign = (fn_coef < 0) ? '-' : '+';
		index += sprintf(function + index, "%c%.3lf * [log_%.3lf_(%.3lfx^%.3lf)] ^ %.3lf ",sign,fabs(fn_coef), base, x_coef, x_exp, fn_exp);
		matrix[0][i]=3;
		matrix[1][i]=x_coef;
		matrix[2][i]=x_exp;
		matrix[3][i]=fn_coef;
		matrix[4][i]=fn_exp;
		matrix[5][i]=base;
		i++;
		count++;
	}
	//*************************************
	printf("\n\n>>Trigonometric count:");
	scanf("%d",&countTrigo);
	count=0;
	while((countTrigo!=0) && (count<countTrigo)){
		printf("\nfn_coef * <trig_fn>(x_coef * x ^ x_exp) ^ fn_exp");
		printf("\n<trig_fn>: sin: 0, cos: 1, tan: 2, cot: 3\n");
		scanf("%lf",&base);
		printf("\nx's cofactor (x_coef) : ");
		scanf("%lf",&x_coef);
		printf("x's exponent (x_exp) : ");
		scanf("%lf",&x_exp);
		printf("Function cofactor (fn_coef) : ");
		scanf("%lf",&fn_coef);
		printf("Function exponent (fn_exp) : ");
		scanf("%lf",&fn_exp);
		const char* trig = (base == 0.0) ? "sin" : (base == 1.0) ? "cos" : (base == 2.0) ? "tan" : (base == 3.0) ? "cot" : "";
		printf("\nAdded : %.3lf * %s(%.3lf * x ^ %.3lf) ^ %.3lf",fn_coef, trig, x_coef, x_exp, fn_exp);
		char sign = (fn_coef < 0) ? '-' : '+';
		index += sprintf(function + index, "%c%.3lf * %s(%.3lf x ^ %.3lf) ^ %.3lf ",sign,fabs(fn_coef),trig,x_coef,x_exp,fn_exp);
		matrix[0][i]=4;
		matrix[1][i]=x_coef;
		matrix[2][i]=x_exp;
		matrix[3][i]=fn_coef;
		matrix[4][i]=fn_exp;
		matrix[5][i]=base;
		i++;
		count++;
	}
	//*************************************
	printf("\n\n>>Inverse-Trigonometric count:");
	scanf("%d",&countInverseTrigo);
	count=0;
	while((countInverseTrigo!=0) && (count<countInverseTrigo)){
		printf("fn_coef * arc<trig_fn>(x_coef * x ^ x_exp) ^ fn_exp");
		printf("\narc<trig_fn>: arcsin: 0, arccos: 1, arctan: 2, arccot: 3 ");
		scanf("%lf",&base);
		printf("\nx's cofactor (x_coef) : ");
		scanf("%lf",&x_coef);
		printf("x's exponent (x_exp) : ");
		scanf("%lf",&x_exp);
		printf("Function cofactor (fn_coef) : ");
		scanf("%lf",&fn_coef);
		printf("Function exponent (fn_exp) : ");
		scanf("%lf",&fn_exp);
		const char* trig = (base == 0.0) ? "arcsin" : (base == 1.0) ? "arccos" : (base == 2.0) ? "arctan" : (base == 3.0) ? "arccot" : "";
		printf("\nAdded : %lf * %s(%lf * x ^ %lf) ^ %lf\n\n",fn_coef, trig, x_coef, x_exp, fn_exp);
		char sign = (fn_coef < 0) ? '-' : '+';
		index += sprintf(function + index, "%c%.3lf * %s(%.3lf x ^ %.3lf) ^ %.3lf ",sign,fabs(fn_coef),trig,x_coef,x_exp,fn_exp);
		matrix[0][i]=5;
		matrix[1][i]=x_coef;
		matrix[2][i]=x_exp;
		matrix[3][i]=fn_coef;
		matrix[4][i]=fn_exp;
		matrix[5][i]=base;
		i++;
		count++;
	}
	return i;
}

double calcFunct(double matrix[][COL], int h, double x){
	int i=0;
	double x_coef, x_exp, fn_coef, fn_exp, base, first, second, sum=0.0, temp;
	while(i<h){
		if(matrix[0][i]==1){//poli		
			x_coef=matrix[1][i];
			x_exp=matrix[2][i];
			sum += x_coef * pow(x, x_exp);
		}else if(matrix[0][i]==2){//expo
			x_coef=matrix[1][i];
			x_exp=matrix[2][i];
			fn_coef=matrix[3][i];
			fn_exp=matrix[4][i];
			base=matrix[5][i];
			first = x_coef * pow(x, x_exp);
			second = pow(base, first);
			sum += fn_coef * pow(second, fn_exp);
		}else if(matrix[0][i]==3){//loga
			x_coef=matrix[1][i];
			x_exp=matrix[2][i];
			fn_coef=matrix[3][i];
			fn_exp=matrix[4][i];
			base=matrix[5][i];
			first = log(x_coef * pow(x, x_exp))/log(base);
			sum += fn_coef * pow(first, fn_exp);
		}else if(matrix[0][i]==4){//trigo
			x_coef=matrix[1][i];
			x_exp=matrix[2][i];
			fn_coef=matrix[3][i];
			fn_exp=matrix[4][i];
			base=matrix[5][i];
			if(base==0.0){//sinx
				temp=x_coef*pow(x, x_exp);
				first= sin(temp);
			}else if(base==1.0){//cosx
				temp=x_coef*pow(x, x_exp);
				first=cos(temp);
			}else if(base==2.0){//tanx
				temp=x_coef*pow(x,x_exp);
				first=tan(temp);
			}else{//base==3 cotx
				temp=x_coef*pow(x, x_exp);
				first= 1.0/tan(temp);
			}
			sum += fn_coef * pow(first, fn_exp) ;
		}else if(matrix[0][i]==5){//inverse-trigo
			x_coef=matrix[1][i];
			x_exp=matrix[2][i];
			fn_coef=matrix[3][i];
			fn_exp=matrix[4][i];
			base=matrix[5][i];
			if(base==0.0){//arcsinx
				temp=x_coef*pow(x, x_exp);
				first=asin(temp);
			}else if(base==1.0){//arccosx
				temp=x_coef*pow(x, x_exp);
				first= acos(temp);
			}else if(base==2.0){//arctanx
				temp=x_coef*pow(x, x_exp);
				first= atan(temp);
			}else{//base==3 arccotx
				temp=x_coef*pow(x, x_exp);
				first= 1.0/atan(temp);
			}
			sum += fn_coef * pow(first, fn_exp);
		}	 
		i++;
	}
	return sum;
}

int getMatrix(float matrix2[][COL]){
	int i,j,N;
	printf("for NxN matrix enter N : ");
	scanf("%d",&N);
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			printf("[%d][%d]",i,j);
			scanf("%f",&matrix2[i][j]);
		}
	}
	printf("\nMatrix: [\n");
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			printf("%f ",matrix2[i][j]);
		}
		printf("\n");
	}	
	printf("]\n\n");
	return N;
}
//******************************************************************************************************************************
void bisection(double matrix[][COL], int h){
	double a, b, c, epsilon, f_a, f_b, f_c, yeni_aralik_ara_delta, x, f_x, control, f_turev, f_turev_2;
	int choice, maxIterations, iteration=1, exit=0, flag=0;
	while(exit==0){
		printf("\nstart: ");
		scanf("%lf",&a);
		printf("end: ");
		scanf("%lf",&b);
		if(b>a){
			exit=1;
		}else{
			printf("\nMust: start < end !!\n");
		}
	}
	printf("epsilon: ");
	scanf("%lf",&epsilon);
	do{
	printf("Stopping criterion:\n1: f(x) <= epsilon\n2: (end-start)/2^n <= epsilon\nChoice: ");
	scanf(" %d",&choice);
	}while(choice !=1 && choice != 2);
	c=(a+b)/2.0;
	f_c=calcFunct(matrix, h, c);
	if(choice==1){
		control=fabs(f_c);
	}else{
		control=(b-a)/pow(2, iteration);
	}
	printf("Max iterations: ");
	scanf("%d",&maxIterations);
	while((iteration<maxIterations) && (control>epsilon) && (flag==0)){
		f_a=calcFunct(matrix, h, a);
		f_b=calcFunct(matrix, h, b);
		printf("\nstart: %lf",a);
		printf("\nend: %lf",b);
		printf("\nmid: %lf",c);
		printf("\nf(start): %lf",f_a);
		printf("\nf(end): %lf",f_b);
		printf("\nf(mid): %lf",f_c);	
		printf("\niteration: %d\n\n",iteration);
		if((f_a>0.0) && (f_b<0.0)){
			c=(a+b)/2.0;
			if(f_c>0.0){
				a=c;
				c=(a+b)/2.0;
			}else if(f_c<0.0){
				b=c;
				c=(a+b)/2.0;
			}
		}else if((f_a<0.0) && (f_b>0.0)){
			if(f_c>0.0){
				b=c;
				c=(a+b)/2.0;
			}else if(f_c<0.0){
				a=c;
				c=(a+b)/2.0;
			}
		}else if((f_a==0.0) || (f_b==0.0)){
			if(f_a==0.0){
				f_c=f_a;
				c=a;
			}else{
				f_c=f_b;
				c=b;
			}
		}else if(f_a * f_b>0.0){
			printf("\nAra deger teoremine gore (f_a * f_b) > 0 oldugundan yeni araliklar aranacak!!!\n");
			yeni_aralik_ara_delta=(b-a)/100000;
			x=a+ yeni_aralik_ara_delta;
			f_x=calcFunct(matrix, h, x);
			f_turev=numericalDifferentiation(matrix, h, x);
			f_turev_2=f_turev;
			while((f_a * f_x > 0.0) && (x<b) && (f_turev*f_turev_2>0)){
				x=x+ yeni_aralik_ara_delta;
				f_x=calcFunct(matrix, h, x);
				f_turev_2=f_turev;
				f_turev=numericalDifferentiation(matrix, h, x);
			}
			if(f_turev*f_turev_2<=0 && fabs(f_x)<0.1){
				flag=1;
			}else if(x>=b){
				flag=2;
			}else if(f_a * f_x < 0.0 && x<b){
				b=x;
				f_b=f_x;
				c=(a+b)/2.0;
			}else{
				flag=3;
			}
		}
		f_c=calcFunct(matrix, h, c);
		iteration++;
		if(choice==1){
			control=fabs(f_c);
		}else{
			control=(b-a)/pow(2, iteration);
		}
	}
	
	if(flag==0){
		printf("\n\n########## RESULT: c : %lf ############\n\n", c);
		printf("\n\n########## RESULT: f_c : %lf ############\n\n", f_c);
	}else if(flag==1){
		printf("\n______CIFT KATLI KOK BULUNDU: %lf________\n",x);
	}else{
		printf("\n______VERILEN ARALIKTA KOK YOKTUR______\n");
	}
}

void regulaFalsi(double matrix[][COL], int h){
	double a, b, c, epsilon, f_a, f_b, f_c, yeni_aralik_ara_delta, x, f_x, control, f_turev, f_turev_2;
	int choice, maxIterations, iteration=1, exit=0, flag=0;
	while(exit==0){
		printf("\nstart: ");
		scanf("%lf",&a);
		printf("end: ");
		scanf("%lf",&b);
		if(b>a){
			exit=1;
		}else{
			printf("\nstart < end olacak sekilde dogru kullanimdir !!\n");
		}
	}
	printf("epsilon: ");
	scanf("%lf",&epsilon);
	do{
	printf("Stopping criterion:\n1: f(x) <= epsilon\n2: (end-start)/2^n <= epsilon\nChoice: ");
	scanf(" %d",&choice);
	}while(choice !=1 && choice != 2);
	f_a=calcFunct(matrix, h, a);
	f_b=calcFunct(matrix, h, b);
	c=(b*f_a-a*f_b)/(f_a-f_b);
	f_c=calcFunct(matrix, h, c);
	if(choice==1){
		control=fabs(f_c);
	}else{
		control=(b-a)/pow(2, iteration);
	}
	printf("Max iterations: ");
	scanf("%d",&maxIterations);
	while((iteration<maxIterations) && (control>epsilon) && (flag==0)){
		f_a=calcFunct(matrix, h, a);
		f_b=calcFunct(matrix, h, b);
		printf("\nstart: %lf",a);
		printf("\nend: %lf",b);
		printf("\nmid: %lf",c);
		printf("\nf(start): %lf",f_a);
		printf("\nf(end): %lf",f_b);
		printf("\nf(mid): %lf",f_c);	
		printf("\niteration: %d\n\n",iteration);
		if((f_a>0.0) && (f_b<0.0)){
			c=(b*f_a-a*f_b)/(f_a-f_b);
			f_c=calcFunct(matrix, h, c);
			if(f_c>0.0){
				a=c;
				c=(b*f_a-a*f_b)/(f_a-f_b);
			}else if(f_c<0.0){
				b=c;
				c=(b*f_a-a*f_b)/(f_a-f_b);
			}
		}else if((f_a<0.0) && (f_b>0.0)){
			if(f_c>0.0){
				b=c;
				c=(b*f_a-a*f_b)/(f_a-f_b);
			}else if(f_c<0.0){
				a=c;
				c=(b*f_a-a*f_b)/(f_a-f_b);
			}
		}else if((f_a==0.0) || (f_b==0.0)){
			if(f_a==0.0){
				f_c=f_a;
				c=a;
			}else{
				f_c=f_b;
				c=b;
			}
		}else if(f_a * f_b>0.0){
			printf("\nAra deger teoremine gore (f_a * f_b) > 0 oldugundan yeni araliklar aranacak!!!\n");
			yeni_aralik_ara_delta=(b-a)/100000;
			x=a+ yeni_aralik_ara_delta;
			f_x=calcFunct(matrix, h, x);
			f_turev=numericalDifferentiation(matrix, h, x);
			f_turev_2=f_turev;
			while((f_a * f_x > 0.0) && (x<b) && (f_turev*f_turev_2>0)){
				x=x+ yeni_aralik_ara_delta;
				f_x=calcFunct(matrix, h, x);
				f_turev_2=f_turev;
				f_turev=numericalDifferentiation(matrix, h, x);
			}
			if(f_turev*f_turev_2<=0 && fabs(f_x)<0.1){
				flag=1;
			}else if(x>=b){
				flag=2;
			}else if(f_a * f_x < 0.0 && x<b){
				b=x;
				f_b=f_x;
				c=(b*f_a-a*f_b)/(f_a-f_b);
			}else{
				flag=3;
			}
		}
		f_c=calcFunct(matrix, h, c);
		iteration++;
		if(choice==1){
			control=fabs(f_c);
		}else{
			control=(b-a)/pow(2, iteration);
		}	
	}
	
	if(flag==0){
		printf("\n\n########## RESULT: c : %lf ############\n\n", c);
		printf("\n\n########## RESULT: f_c : %lf ############\n\n", f_c);
	}else if(flag==1){
		printf("\n______CIFT KATLI KOK BULUNDU: f(%lf): %lf________\n",x,f_x);
	}else{
		printf("\n______VERILEN ARALIKTA KOK YOKTUR______\n");
	}
}

void newtonRaphson(double matrix[][COL], int h){
	double f_deriv_xn, xn, xn_1, f_xn, f_xn_1, epsilon;
	int iteration=0, maxIterations;
	printf("\nx0: ");
	scanf("%lf", &xn);
	printf("\nepsilon: ");
	scanf("%lf", &epsilon);
	printf("\nMax iterations: ");
	scanf("%d", &maxIterations);
	f_xn_1=calcFunct(matrix, h, xn);
	while((iteration<maxIterations) && (fabs(f_xn_1)>epsilon)){
		f_xn=calcFunct(matrix, h, xn);
		f_deriv_xn=numericalDifferentiation(matrix, h, xn);
		xn_1=xn-(f_xn/f_deriv_xn);
		f_xn_1=calcFunct(matrix, h, xn_1);
		printf("\n\nxn: %lf", xn);
		xn=xn_1;
		printf("\nxn+1: %lf", xn_1);
		printf("\nf(xn): %lf", f_xn);
		printf("\nf'(xn): %lf", f_deriv_xn);
		printf("\niteration: %d", iteration);
		iteration++;
	}
	printf("\n\nResult: f(%lf): %lf\n\n", xn_1, f_xn_1);
}

void inverseMatrix(float matrix2[][COL], float islemMatrisi[][COL], float birimMatris[][COL], int N){
	int i, j, k;
	float temp1, temp2;
	float factor, pivot;
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			islemMatrisi[i][j]=matrix2[i][j];
		}
	}
	for(i=0; i<N; i++){
		birimMatris[i][i]=1;
	}
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            if(i==j){
                birimMatris[i][j]=1.0;
            }else{
                birimMatris[i][j]=0.0;
            }
        }
    }
    for(i=0; i<N; i++){
        pivot=islemMatrisi[i][i];
        for(j=0; j<N; j++){
            islemMatrisi[i][j]/=pivot;
            birimMatris[i][j]/=pivot;
        }
        for(k=0; k<N; k++){
            if(k!=i){
                factor=islemMatrisi[k][i];
                for(j=0; j<N; j++) {
                    islemMatrisi[k][j] -=factor*islemMatrisi[i][j];
                    birimMatris[k][j] -=factor*birimMatris[i][j];
                }
            }
        }
    }	
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf(" %f", birimMatris[i][j]);
		}
		printf("\n");
	}	
}

void gaussElimination(){
	float matrix[COL][COL]={{0}};
	float birimMatris[COL]={0};
	float temp;
	float temp1, temp2;
	float factor, pivot;
	int N, i, j, k;
	printf("\nNumber of equations: \n");
	scanf("%d", &N);
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%d. equation's %d. unknown's cofactor: ", i+1, j+1);
			scanf("%f", &temp);
			matrix[i][j]=temp;
		}
		printf("%d. equation equals to: ", i);
		scanf("%f", &temp);
		birimMatris[i]=temp;
	}
	//
    for(i=0; i<N; i++){
        pivot=matrix[i][i];
        birimMatris[i]/=pivot;
        for(j=0; j<N; j++){
            matrix[i][j]/=pivot;
        }
        for(k=0; k<N; k++){
            if(k!=i){
                factor=matrix[k][i];
                birimMatris[k]-=factor*birimMatris[i];
                for(j=0; j<N; j++) {
                    matrix[k][j]-=factor*matrix[i][j];  
                }
            }
        }
    }	
	for(i=0; i<N; i++){
		printf("x%d : %f", i, birimMatris[i]);
		printf("\n");
	}
}
//*******************************************************************
void gaussSeidel(){
	int n;
    float epsilon;
	float A[COL][COL];
    float b[COL];
    float x[COL];
    int i, j, k;
    
    printf("Enter the Total Number of Equations: ");
    scanf("%d", &n);
    
    printf("Enter the error tolerance (epsilon): ");
    scanf("%f", &epsilon);

    printf("Enter the Coefficients\n");
    for (i = 0; i < n; i++) {
        printf("Equation %d:\n", i + 1);
        for (j = 0; j < n; j++) {
            printf("Matrix[%d][%d] = ", i, j);
            scanf("%f", &A[i][j]);
        }
        printf("Matrix[%d][%d] = ", i, n);
        scanf("%f", &b[i]);
    }
    char contiNu;
    int r1,r2;
    printf("\nDiagonal MAX?: y/n\n");
    scanf(" %c",&contiNu);
    while(contiNu == 'n'){
    	printf("\nWhich two rows to swap?: ");
    	printf("\nrow1(1 space " ")row2 ");
		scanf(" &d &d",&r1,&r2);
    	rowExchanger(A,b,n,r1,r2);
    	printf("\nDiagonal MAX yet?: y/n");
    	scanf(" %c",&contiNu);
	}
	printf("\n");
	for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("[%d][%d] = %f ", i, j, A[i][j]);
        }
        printf("[%d][%d]= %f\n",i,n,b[i]);
	}
	
    // Gauss-Seidel method
    float sum;
    float diff;
    int iterations = 0;
    do{
        diff = 0.0;
        for (i = 0; i < n; i++) {
            sum = 0.0;
            for (j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }
            float prevX = x[i];
            x[i] = (b[i] - sum) / A[i][i];
            diff = fabs(x[i] - prevX);
        }
        iterations++;
    }while(diff>epsilon && iterations < 100000);

    printf("\nSolution:\n");
    for (i = 0; i < n; i++) {
        printf("X[%d] = %f\n", i, x[i]);
    }
    printf("Number of iterations: %d\n", iterations);
}

void rowExchanger(float matrix[][COL], float b[], int N, int row1, int row2){
	float temp;
	int j;
	
	for(j=0;j<=N;j++){
		temp=matrix[row1][j];
		matrix[row1][j]=matrix[row2][j];
		matrix[row2][j]=temp;
	}
	temp=b[row1];
	b[row1]=b[row2];
	b[row2]=temp;
}
//*******************************************************************
double numericalDifferentiation(double matrix[][COL], int t, double xi){
	double f_deriv, f_xi_plus_h, f_xi_minus_h, h=0.000001;
	f_xi_plus_h=calcFunct(matrix, t, xi+h);
	f_xi_minus_h=calcFunct(matrix, t, xi-h);
	f_deriv=(f_xi_plus_h-f_xi_minus_h)/(2*h);
	return f_deriv;
}

void numericalDiff(double matrix[][COL], int t){
	double f_deriv, f_xi_plus_h, f_xi_minus_h, h, xi;
	printf("\nx0: ");
	scanf("%lf", &xi);
	printf("\nh:");
	scanf("%lf", &h);
	f_xi_plus_h=calcFunct(matrix, t, xi+h);
	f_xi_minus_h=calcFunct(matrix, t, xi-h);
	f_deriv=(f_xi_plus_h-f_xi_minus_h)/(2*h);
	printf("Merkezi farkla f'(%lf): %lf\n", xi,f_deriv);
}

void simpsonsRule(double matrix[][COL], int t){
	double sum=0.0, result=0.0, h, x0, xn, f_xi, temp, a, b;
	int i, n;
	printf("\nSimpson 1/3\n");
	printf("\nx0: ");
	scanf("%lf", &x0);
	printf("\nxn: ");
	scanf("%lf", &xn);
	do{
		printf("\nn(even): ");
		scanf("%d", &n);
	}while(n%2!=0);
	h=(xn-x0)/n;
	f_xi=calcFunct(matrix, t, x0);
	sum+=f_xi;
	f_xi=calcFunct(matrix, t, xn);
	sum+=f_xi;	
	for(i=1; i<n; i+=2){
		f_xi=calcFunct(matrix, t, x0+i*h);
		sum+=4*f_xi;
	}	
	for(i=2; i<n; i+=2){
		f_xi=calcFunct(matrix, t, x0+i*h);
		sum+=2*f_xi;
	}	
	result=sum*h/3;
	printf("\nResult(simpson 1/3): %lf\n", result);
	//
	sum=0.0;
	result=0.0;
	printf("\nSimpson 3/8\n");
	printf("\nx0: ");
	scanf("%lf", &x0);
	printf("\nxn: ");
	scanf("%lf", &xn);
	do{
		printf("\nn(3x): ");
		scanf("%d", &n);	
	}while(n%3!=0);
	h=(xn-x0)/n;
	for(i=0;i<n;i++){
		temp=x0+h*i;
		f_xi=calcFunct(matrix, t, temp);
		sum+=f_xi;
		
		temp=x0+h*i+(h/3.0);
		f_xi=calcFunct(matrix, t, temp);
		sum+=3*f_xi;	
		
		temp=x0+h*i+2*(h/3.0);
		f_xi=calcFunct(matrix, t, temp);
		sum+=3*f_xi;
		
		temp=x0+h*i+h;
		f_xi=calcFunct(matrix, t, temp);
		sum+=f_xi;
		
		result+=sum*(h/8.0);
		sum=0.0;
	}
	printf("\nResult(simpson 3/8): %lf\n", result);
}

void trapezoidalRule(double matrix[][COL], int t){
	double sum=0.0, result=0.0, h, x0, xn, f_xi, temp;
	int i, n;
	printf("\nTrapezoidal: \n");
	printf("\nx0: ");
	scanf("%lf", &x0);
	printf("\nxn: ");
	scanf("%lf", &xn);
	printf("\nn:");
	scanf("%d", &n);
	h=(xn-x0)/n;
	f_xi=calcFunct(matrix, t, x0);
	sum+=f_xi/2;
	f_xi=calcFunct(matrix, t, xn);
	sum+=f_xi/2;
	for(i=1; i<n; i++){
		temp=x0+i*h;
		f_xi=calcFunct(matrix, t, temp);
		sum+=f_xi;
	}
	result=h*sum;
	printf("\nResult: %lf\n", result);
}

void gregoryNewtonInterpolation(){
	float x0, xn, h, x, sum=0.0, k=1.0, hh=1.0;
	int i, j, n, col=2, factorial=1;
	float matrix[COL][COL]={{0}};
	printf("\nGregory-Newton Interpolation\n");
	printf("\nStart x0: ");
	scanf("%f",&x0);
	printf("\nEnd xn: ");
	scanf("%f",&xn);
	printf("\nNumber of points n: ");
	scanf("%d",&n);
	h=(xn-x0)/(n-1);
	for(i=0;i<n;i++){
		matrix[i][0]=x0+i*h;
		printf("f(%f): ",matrix[i][0]);
		scanf("%f",&matrix[i][1]);
	}
	do{
		printf("\nx (x0<x<xn): ");
		scanf("%f",&x);
	}while(x<x0 || x>xn);
	for(i=n-1;i>0;i--){
		for(j=0;j<i;j++){
			matrix[j][col]=matrix[j+1][col-1]-matrix[j][col-1];
		}
		col++;
	}
	sum+=matrix[0][1];
	for(i=1;i<n;i++){
		hh*=h;
		k*=(x-matrix[i-1][0]);
		factorial*=i;
		sum+=k*matrix[0][i+1]/(factorial*hh);
		printf("\nk: %f\ndelta_f: %f\nfactorial: %d\nsum: %f\n",k,matrix[0][i+1],factorial,sum);
	}
	printf("f(%f): %f\n",x,sum);
}
