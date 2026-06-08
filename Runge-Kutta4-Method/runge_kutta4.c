#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void enterG_x();
double calcG_x(double x);
double calcF_x(double xi, double yi);
struct {
	double a1, a0, x0, y0, xt, *ptrPolynomialG, *ptrTrigGx, *ptrLog, **ptrPay, **ptrPayda;
	int  NU, NA, MU, MA;
} ODE;
int main(){
	double k1, k2, k3, k4, xi, yi, h;
	int i;
	printf("Equation a1*y' + a0*y + g(x) = 0\n");
	printf("Enter a1:"); scanf("%lf",&ODE.a1);
	printf("Enter a0:"); scanf("%lf",&ODE.a0);
	ODE.ptrPolynomialG=(double*)calloc(50,sizeof(double));
	ODE.ptrTrigGx=(double*)calloc(75,sizeof(double));
	ODE.ptrLog=(double*)calloc(75,sizeof(double));
	ODE.ptrPay=(double**)calloc(75,sizeof(double*));
	ODE.ptrPayda=(double**)calloc(75,sizeof(double*));
	for(i=0;i<20;i++){
		ODE.ptrPay[i]=(double*)calloc(75,sizeof(double));
		ODE.ptrPayda[i]=(double*)calloc(75,sizeof(double));
	}
	enterG_x();
	printf("\nx0: "); scanf("%lf",&ODE.x0); xi=ODE.x0;
	printf("y0: "); scanf("%lf",&ODE.y0); yi=ODE.y0;
	printf("xt: "); scanf("%lf",&ODE.xt);
	printf("h: "); scanf("%lf",&h);
	while(xi<ODE.xt){
		printf("\ny(%.3lf): %lf",xi,yi);
		k1=calcF_x(xi, yi);
		k2=calcF_x(xi+0.5*h, yi+0.5*h*k1);
		k3=calcF_x(xi+0.5*h, yi+0.5*h*k2);
		k4=calcF_x(xi+h, yi+h*k3);
		xi=xi+h;
		yi=yi+(k1 + 2*k2 + 2*k3 + k4)*h/6.0;
	}
	printf("\ny(%.3lf): %lf",xi,yi);
	free(ODE.ptrPolynomialG);
	free(ODE.ptrTrigGx);
	free(ODE.ptrLog);
	for(i=0;i<20;i++){
		free(ODE.ptrPay[i]);
		free(ODE.ptrPayda[i]);
	}	
	free(ODE.ptrPay);
	free(ODE.ptrPayda);
	return 0;
}

void enterG_x(){ //kullanicidan g(x) polinomunu girdi olarak al
	int N, M, i, j;
	printf("\nSadece Polinom turundeki ifade sayisini giriniz:"); scanf("%d",&N);
	for(i=0;i<N;i++){
		printf("%d. Eleman katsayisini giriniz: ",i+1); scanf("%lf",&ODE.ptrPolynomialG[i*2]); //katsayilar cift 0 2 4 6 8 ... nolu indislerde
		printf("%d. Eleman derecesini giriniz: ",i+1); scanf("%lf",&ODE.ptrPolynomialG[i*2+1]); //Usler tek 1 3 5 7 9 ... nolu indislerde
	}
	printf("Polinom: ");
	for(i=0;i<N;i++){
		printf("+(%.2lf*x^%.2lf)",ODE.ptrPolynomialG[i*2],ODE.ptrPolynomialG[i*2+1]); //polinomu ekrana yazdir
	}
	printf("\nSadece Trigonometrik ifade sayisini giriniz:"); scanf("%d",&N);
	for(i=0;i<N;i++){
		printf("\n1:sin 2:cos 3:tan 4:cot"); scanf("%lf",&ODE.ptrTrigGx[i*3]);
		printf("%d. Eleman katsayisini giriniz: ",i+1); scanf("%lf",&ODE.ptrTrigGx[i*3+1]);
		printf("%d. Eleman derecesini giriniz: ",i+1); scanf("%lf",&ODE.ptrTrigGx[i*3+2]);
	}
	printf("\nSadece Logaritmik ifade sayisini giriniz:"); scanf("%d",&N);
	for(i=0;i<N;i++){
		printf("\nkatsayi_ log _ taban _ (x) ^ us\n %d.terim katsayisini giriniz:",i+1); scanf("%lf",&ODE.ptrLog[i*3]);
		printf("%d. terim taban giriniz: ",i+1); scanf("%lf",&ODE.ptrLog[i*3+1]);
		printf("%d. terim us giriniz: ",i+1); scanf("%lf",&ODE.ptrLog[i*3+2]);
	}
	printf("\nPay: toplam durumundaki ifade sayisini giriniz:"); scanf("%d",&N); ODE.NU=N;
	for(i=0;i<N;i++){//satir
		printf("\n%d.ifadede carpim durumunda kac terim var",i+1); scanf("%d",&M); ODE.MU=M;
		for(j=0;j<M;j++){//sutun
			printf("%d.ifadedeki %d. terim turu giriniz: (0:Polinom 1:sin 2:cos 3:tan 4:cot 5:log)   ",i+1,j+1);
			scanf("%lf",&ODE.ptrPay[i][3*j]);
			printf("%d. terim katsayisini giriniz: ",j+1); scanf("%lf",&ODE.ptrPay[i][3*j+1]);
			printf("%d. terim derecesini giriniz: ",j+1); scanf("%lf",&ODE.ptrPay[i][3*j+2]);		
		}
	}
	printf("\nPayda: toplam durumundaki ifade sayisini giriniz:"); scanf("%d",&N); ODE.NA=N;
	for(i=0;i<N;i++){//satir
		printf("\n%d.ifadede carpim durumunda kac terim var",i+1); scanf("%d",&M); ODE.MA=M;
		for(j=0;j<M;j++){//sutun
			printf("%d.ifadedeki %d. terim turu giriniz: (0:Polinom 1:sin 2:cos 3:tan 4:cot 5:log)   ",i+1,j+1);
			scanf("%lf",&ODE.ptrPayda[i][3*j]);
			printf("%d. terim katsayisini giriniz: ",j+1); scanf("%lf",&ODE.ptrPayda[i][3*j+1]);
			printf("%d. terim derecesini giriniz: ",j+1); scanf("%lf",&ODE.ptrPayda[i][3*j+2]);		
		}
	}		
}

double calcG_x(double x){
	double tmp1, tmp, us, katsayi, sum=0;
	int i, j;
	for(i=0;i<sizeof(ODE.ptrPolynomialG)/2;i++){
		us=ODE.ptrPolynomialG[2*i+1];
		katsayi=ODE.ptrPolynomialG[2*i];
		tmp=pow(x,us);
		sum+=katsayi*tmp;
	}
	for(i=0;i<sizeof(ODE.ptrTrigGx)/3;i++){
		if(ODE.ptrTrigGx[3*i]==1){
			tmp1=sin(x);
		}else if(ODE.ptrTrigGx[3*i]==2){
			tmp1=cos(x);
		}else if(ODE.ptrTrigGx[3*i]==3){
			tmp1=tan(x);
		}else{
			tmp1=1.0/tan(x);
		}
		katsayi=ODE.ptrTrigGx[3*i+1];
		us=ODE.ptrTrigGx[3*i+2];	
		tmp=pow(tmp1,us);
		sum+=katsayi*tmp;
	}
	for(i=0;i<sizeof(ODE.ptrLog)/3;i++){
		katsayi=ODE.ptrLog[3*i];
		us=ODE.ptrLog[3*i+2];
		tmp=log(x)/log(ODE.ptrLog[3*i+1]);//log_a_(x)=log(x)/log(a)
		tmp=pow(tmp,us);
		sum+=katsayi*tmp;
	}	
	double tmp1U, tmp1A, katsayiU, katsayiA, tmpU, tmpA, sumU=0, sumA=0, usU, usA, carpimxU, carpimxA;
	for(j=0;j<ODE.NU;j++){
		carpimxU=1;
		for(i=0;i<ODE.MU;i++){
			if(ODE.ptrPay[j][3*i]==0){
				katsayiU=ODE.ptrPay[j][3*i+1];
				usU=ODE.ptrPay[j][3*i+2];
				tmp1U=pow(x,usU);
				carpimxU*=katsayiU*tmp1U;
			}else{
				if(ODE.ptrPay[j][3*i]==1){
				tmp1U=sin(x);
				}else if(ODE.ptrPay[j][3*i]==2){
				tmp1U=cos(x);
				}else if(ODE.ptrPay[j][3*i]==3){
				tmp1U=tan(x);
				}else if(ODE.ptrPay[j][3*i]==4){
				tmp1U=1.0/tan(x);
				}else{
				tmp1U=log(x);
				}
				katsayiU=ODE.ptrPay[j][3*i+1];
				usU=ODE.ptrPay[j][3*i+2];	
				tmpU=pow(tmp1U,usU);
				carpimxU*=katsayiU*tmpU;
			}
		}
		sumU+=carpimxU;
	}
	for(j=0;j<ODE.NA;j++){
		carpimxA=1;
		for(i=0;i<ODE.MA;i++){
			if(ODE.ptrPayda[j][3*i]==0){
				katsayiA=ODE.ptrPayda[j][3*i+1];
				usA=ODE.ptrPayda[j][3*i+2];
				tmp1A=pow(x,usA);
				carpimxA*=katsayiA*tmp1A;
			}else{
				if(ODE.ptrPayda[j][3*i]==1){
				tmp1A=sin(x);
				}else if(ODE.ptrPayda[j][3*i]==2){
				tmp1A=cos(x);
				}else if(ODE.ptrPayda[j][3*i]==3){
				tmp1A=tan(x);
				}else if(ODE.ptrPayda[j][3*i]==4){
				tmp1A=1.0/tan(x);
				}else{
				tmp1A=log(x);	
				}
				katsayiA=ODE.ptrPayda[j][3*i+1];
				usA=ODE.ptrPayda[j][3*i+2];
				tmpA=pow(tmp1A,usA);
				carpimxA*=katsayiA*tmpA;
			}
		}
		sumA+=carpimxA;
	}
	if(ODE.NA!=0 && ODE.NU==0){
		sumU=1;
	}
	if(ODE.NA==0 && ODE.NU!=0){
		sumA=1;
	}
	if(ODE.NA==0 && ODE.NU==0){
		sumU=0; sumA=1;
	}
	sum+=sumU/sumA;
	return sum;
}

double calcF_x(double xi, double yi){
	double h=0.1, var=0, tmp=0;
	tmp=calcG_x(xi);
	var=-1*(ODE.a0*yi+tmp)/ODE.a1;
	return var;
}
