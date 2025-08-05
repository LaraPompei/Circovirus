#include <iostream>
#include <cmath>
#include <fstream>

#define eq 13 
#define t_store 1000        //Intervalo de pontos sendo salvos

//Parametros:

//#define pi_v 1.8e-1         //Taxa de replicacao viral
#define pi_v 0.8            //Taxa de replicacao viral
//#define c_v1 2.63           //Taxa de clareamento viral maximo pelo sistema inato
#define c_v1 32.63           //Taxa de clareamento viral maximo pelo sistema inato
//#define c_v2 6e-1           //Constante de meia saturacao
#define c_v2 8.1e-1           //Constante de meia saturacao
//#define k_v1 4.82e-5        //Taxa de neutralizacao do virus por unidade anticorpos neutralizantes
#define k_v1 1.6e-4        //Taxa de neutralizacao do virus por unidade anticorpos neutralizantes
//#define k_v2 7.48e-7        //Taxa de eliminacao do virus por unidade de celulas T CD8+
#define k_v2 3.2e-7        //Taxa de eliminacao do virus por unidade de celulas T CD8+
//#define k_v3 4.82e-5        //Taxa de neutralizacao do virus por unidade anticorpos IgG
#define k_v3 4.82e-7

//#define alpha_ap 2.5e-3     //Taxa de hosmeostase das APCs imaturas
#define alpha_ap 1.5e-3     //Taxa de hosmeostase das APCs imaturas
//#define beta_ap 5.5e-1      // Taxa de maturacao das APCs 
#define beta_ap 1e-1      // Taxa de maturacao das APCs

//#define c_ap1 8e-1          //Taxa de maturacao maxima das APCs
#define c_ap1 12e-1          //bem sensivel - desloca a curva dos anticorpos no eixo x
//#define c_ap2 4e1           //Constante de meia ativacao
#define c_ap2 2e2           //mexe no pico dos anticorpos

//#define delta_apm 5.38e-1   //Taxa de morte das APCs maduras
#define delta_apm 3.1e-1   //Taxa de morte das APCs maduras
//#define alpha_th 2.17e-4    //Taxa de gineistase das celulas T CD4+
#define alpha_th 2.17e-4    //Taxa de gineistase das celulas T CD4+
#define beta_th 1e-7        //Taxa de replicacao das celulas T CD4+ naive
#define pi_th 1e-8          //Taxa de replicacao das celulas T CD4+ efetoras
#define delta_th 2.2e-1     //Taxa de morte das celulas T CD4+ efetoras

//#define alpha_tk 2.17e-4    //Taxa de homeostase das celulas T CD8+
#define alpha_tk 2.17e-4    //Taxa de homeostase das celulas T CD8+
#define beta_tk 1e-5        //Taxa de ativacao das celulas T CD8+ naive
#define pi_tk 1e-8          //Taxa de replicacao das celulas T CD8+ efetoras
#define delta_tk 3e-4       //Taxa de morte das celulas T CD8+ efetoras

//#define alpha_b 6.0         //Taxa de homeostase das celulas B
#define alpha_b 5.5         //Taxa de homeostase das celulas B

//#define pi_b1 4.83e-6       //Taxa de ativacao das celulas B T-independente
#define pi_b1 4.83e-7       //Taxa de ativacao das celulas B T-independente
//#define pi_b2 1.27e-8       //Taxa de ativacao das celulas B T-dependentes
#define pi_b2 1.27e-6       //Taxa de ativacao das celulas B T-dependentes

//#define beta_ps 6.72e-4     //Taxa de diferenciacao das celulas B ativas em plasmocitos de vida curta
#define beta_ps 1.8522e-2     //Taxa de diferenciacao das celulas B ativas em plasmocitos de vida curta
//#define beta_pl 5.61e-6     //Taxa de diferenciacao das celulas B ativas em plasmocitos de vida longa
#define beta_pl 6.61e-3     //Taxa de diferenciacao das celulas B ativas em plasmocitos de vida longa

//#define beta_bm 1e-6        //Taxa de diferenciacao das celulas B ativas em celulas B de memoria
#define beta_bm 1e-6        //afeta a velocidade que o virus decai
//#define delta_ps 2.0        //Taxa de morte dos plasmocitos de vida curta
#define delta_ps 2.01        //Taxa de morte dos plasmocitos de vida curta
//#define delta_pl 2.4e-4     //Taxa de morte dos plasmocitos de vida longa
#define delta_pl 3.4e-3     //Taxa de morte dos plasmocitos de vida longa
//#define gama_bm 9.75e-4     //Taxa de diferenciacao das celulas B de memoria em plasmocitos de vida longa
#define gama_bm 4.75e-3     //Taxa de diferenciacao das celulas B de memoria em plasmocitos de vida longa

#define pi_bm1 1e-5         //Taxa de proliferacao das celulas B de memoria
#define pi_bm2 2.5e3        //Constante de crescimento maximo
//#define pi_ps 2e-3          //Taxa de secrecao de anticorpos por unidade de plasmocitos de vida curta

#define c_ps1 23.3e-3
#define c_ps2 23.5
#define c_ps3 50.0 
#define c_ps4 50.0
#define c_ps5 5.0

#define c_pl1 28.3e-5
#define c_pl2 44.5
#define c_pl3 60.0
#define c_pl4 80.0
#define c_pl5 21.0


#define pi_ps 19e-4      //Taxa de secrecao de anticorpos por unidade de plasmocitos de vida curta
//#define pi_pl 6.8e-4        //Taxa de secrecao de anticorpos por unidade de plasmocitos ded vida longa
#define pi_pl 1.0e-5        //Taxa de secrecao de anticorpos por unidade de plasmocitos ded vida longa
//#define delta_a 4e-2        //Taxa de morte de anticorpos
#define delta_IgM 4.42e-1        //Taxa de morte de IgM
#define delta_IgG 33.5e-2        //Taxa de morte de IgG

//Condicoes iniciais
//#define V0 724.0
//#define V0 275.28
//#define V0 79.39
//#define V0 15.47
//#define V0 4.16
#define V0 9570.81 
#define Ap0 0.83e6
//#define Ap0 10.6e6
#define Apm0 0.0
#define Thn0 1.56e6//1.56e6
//#define Thn0 12.8e6
#define The0 0.0
#define Tkn0 0.91e6//0.91e6
//#define Tkn0 10.4e6//28.8e6
#define Tke0 0.0
#define B0 1.39e6
//#define B0 12.4e6
#define Ps0 0.0
#define Pl0 0.0
#define Bm0 0.0
#define A0 0.0

using namespace std;

void SistemaTeste(double *y, double* dydt){
    double k = 9.0;
    double rmax = 0.0065;
    dydt[0] = rmax*((k-y[0])/k)*y[0];
}

void SistemaTeste2(double *y, double* dydt){
    double rmax = 0.5;
    double m = 0.50;
    double s = 3.0;
    double a1 = 0.3;
    double a2 = 0.7;
    ((rmax*y[1])-m)*y[0]*(1.0-y[0]);
    s*(1.0-y[1]) - (a1*y[0]*y[1]) - (a2*y[1]);
}

void SistemaTeste3(double *y, double* dydt){
    dydt[0] = 1.0*y[0] - 1.0*y[0]*y[1];
    dydt[1] = -1.0*y[1] + 2.0*y[0]*y[1];
}

void Sistema(double *y, double* dydt, double t){
    //V
    dydt[0] = pi_v*y[0] - (c_v1*y[0])/(c_v2+y[0]) - k_v1*y[0]*y[11] - k_v2*y[0]*y[6] - k_v3*y[0]*y[12];
    //Ap
    dydt[1] = alpha_ap*(Ap0 - y[1]) - beta_ap*y[1]*(c_ap1*(y[0])/(c_ap2 + y[0]));
    //Apm
    dydt[2] = beta_ap*y[1]*(c_ap1*(y[0])/(c_ap2 + y[0])) - delta_apm*y[2];
    //Thn
    dydt[3] = alpha_th*(Thn0 - y[3]) - beta_th*y[2]*y[3];
    //The
    dydt[4] = beta_th*y[2]*y[3] + pi_th*y[2]*y[4] - delta_th*y[4];
    //Tkn
    dydt[5] = alpha_tk*(Tkn0 - y[5]) - beta_tk*y[2]*y[5];
    //Tke
    dydt[6] = beta_tk*y[2]*y[5] + pi_tk*y[2]*y[6] - delta_tk*y[6];
    //B
    dydt[7] = alpha_b*(B0 - y[7]) + pi_b2*y[4]*y[7] - beta_ps*y[2]*y[7] - beta_pl*y[4]*y[7] - beta_bm*y[4]*y[7];
    //Ps
    dydt[8] = beta_ps*y[2]*y[7] - delta_ps*y[8];
    //Pl
    dydt[9] = beta_pl*y[4]*y[7] - delta_pl*y[9] + gama_bm*y[10];
    //Bm
    dydt[10] = beta_bm*y[4]*y[7] + pi_bm1*y[10]*(1 - (y[10]/pi_bm2)) - gama_bm*y[10];
    //IgM
    dydt[11] = c_ps1*y[8]*(1-exp(-pow((t/c_ps2),c_ps3))*exp(-pow((t/c_ps4),c_ps5))) - delta_IgM*y[11];
    //IgG
    //dydt[12] = pi_pl*y[9] - delta_IgG*y[12];
    dydt[12] = c_pl1*y[9]*(1-exp(-pow((t/c_pl2),c_pl3))*exp(-pow((t/c_pl4),c_pl5))) - delta_IgG*y[12];
}

void saveData (double** y, double* t, double h, int pont, int i_atual){
    fstream outputFile("output.csv",fstream::app);
    int x = (i_atual-1)/t_store;
    cout<<"x "<< x <<endl<<"pont "<<pont<<endl<<"i_atual "<<i_atual<<endl;
    for(int i=0; i<pont; i++){
        outputFile<<((x*t_store+i)*h);
        for(int j=0; j<eq;j++){
            outputFile<<","<<y[i][j];
        }
        outputFile<<endl;
    }
    outputFile.close();
}

void SaveK(double* k, int& lSave, double i_atual){

    fstream fileOut;
    char fileName[10];

    sprintf(fileName, "k%d.dat", lSave);
    fileOut.open(fileName, fstream::app);
    fileOut<<i_atual<<'\t';

    for(int i=0; i<eq; i++){
        fileOut<<k[i]<<'\t';
    }

    fileOut<<endl;
    fileOut.close();
    lSave++;
}

void RK5(double* t, double h, double** y, int inter){
    double* k1 = new double[eq];
    double* k2 = new double[eq];
    double* k3 = new double[eq];
    double* k4 = new double[eq];
    double* k5 = new double[eq];
    double* k6 = new double[eq];
    double* yk = new double[eq];
    int i,p;
    for(i=1; i<inter; i++){
        p = i%t_store; //posicao do vetor que estamos calculando
        int lSave = 0; //usado para salvar os ks]

        if(p==0){
            cout<<"Saving Data "<<i<<endl;
            saveData(y,t,h,t_store+p,i);
        }

        //Calculando K1=f(t[i],y[i])
        if(p==0){
            Sistema(y[t_store-1],k1,i*h); // passo de tempo anterior = ultima posicao do vetor
            cout<<"Tempo "<<i*h<<" dias"<<endl;
            cerr<<"Interacao "<<i<<endl;
        }else{
            Sistema(y[p-1],k1,i*h);
        }
        //SaveK(k1, lSave, i*h);

        //Calculando K2 = f(t[i]+(1/5)*h,y[i]+(1/5)*k1)
        for(int j=0;j<eq;j++){
            k1[j] *= h;
            yk[j] = (p == 0)?y[t_store-1][j] : y[p-1][j];
            yk[j] += (1.0/5.0)*k1[j];
        }
        Sistema(yk,k2,i*h);
        //SaveK(k2, lSave, i*h);

        //Calculando K3 = f(y[i]+(3/40)*k1+(9/40)*k2)
        for(int j=0;j<eq;j++){
            k2[j] *= h;
            yk[j] = (p == 0)?y[t_store-1][j] : y[p-1][j];
            yk[j] += (3.0/40.0)*k1[j]+(9.0/40.0)*k2[j];
        }
        Sistema(yk,k3,i*h);
        //SaveK(k3, lSave, i*h);

        //Calculando K4 = f(y[i]+(3/10)*k1-(9/10)*k2+(6/5)*k3)
        for(int j=0;j<eq;j++){
            k3[j] *= h;
            yk[j] = (p == 0)?y[t_store-1][j] : y[p-1][j];
            yk[j] += (3.0/10.0)*k1[j]-(9.0/10.0)*k2[j]+(6.0/5.0)*k3[j];
        }
        Sistema(yk,k4,i*h);
        //SaveK(k4, lSave, i*h);

        //Calculando K5 = f(y[i]-(11/54)*k1+(5/2)*k2-(70/27)*k3+(35/27)*k4)
        for(int j=0;j<eq;j++){
            k4[j] *= h;
            yk[j] = (p == 0)?y[t_store-1][j] : y[p-1][j];
            yk[j]+= -(11.0/54.0)*k1[j]+(5.0/2.0)*k2[j]-(70.0/27.0)*k3[j]+(35.0/27.0)*k4[j];
        }
        Sistema(yk,k5,i*h);
        //SaveK(k5, lSave, i*h);

        //Calculando K6 = f(y[i]+(1631/55296)*k1+(175/512)*k2+(575/13824)*k3+(44275/110592)*k4+(253/4096)*k5))
        for(int j=0;j<eq;j++){
            k5[j] *= h;
            yk[j] = (p == 0) ? y[t_store-1][j] : y[p-1][j];
            yk[j]+= (1631.0/55296.0)*k1[j]+(175.0/512.0)*k2[j]+(575.0/13824.0)*k3[j]+(44275.0/110592.0)*k4[j]+(253.0/4096.0)*k5[j];
        }
        Sistema(yk,k6,i*h);
        //SaveK(k6, lSave, i*h);

        //Calculando o y para cada eq h*(37*k1[j]/378 +0*k2[j] 250*k3[j]/621 + 125*k4[j]/594 +0*k5[j] + 512*k6[j]/1771);
        for(int j=0;j<eq;j++){
            k6[j] *= h;
            y[p][j] = (p==0)?y[t_store-1][j] : y[p-1][j];
            y[p][j] += ((37.0/378.0)*k1[j]+(250.0/621.0)*k3[j]+(125.0/594.0)*k4[j]+(512.0/1771.0)*k6[j]);
        }

        //passo de tempo
        t[p] = (p!=0) ? t[p-1] + h: t[t_store - 1] + h;

    }
    saveData(y,t,h,i%t_store,i);

}

int main(){

    double t0 = 0.0;
    
    //SistemaTeste 2:
    //double t_final = 50.0;
    //double h = 0.01;

    //SistemaTeste 3:
    //double t_final = 20.0;
    //double h = 0.00001;

    //Sistema:
    //double t_final = 45.0;
    double t_final = 600.0;
    //double t_final = 6000.0;
    double h = 0.001;
        
    int inter = int(abs(t_final-t0)/h);
    cout<<"Numero de interacoes: "<<inter<<endl;
    double t[t_store] = {t0};
    double** y = new double*[t_store];

    for(int i = 0; i<t_store; i++)
        y[i] = new double[eq];


//Sistema
    y[0][0] = V0;
    y[0][1] = Ap0;
    y[0][2] = Apm0;
    y[0][3] = Thn0;
    y[0][4] = The0;
    y[0][5] = Tkn0;
    y[0][6] = Tke0;
    y[0][7] = B0;
    y[0][8] = Ps0;
    y[0][9] = Pl0;
    y[0][10] = Bm0;
    y[0][11] = A0;

/*
//SistemaTeste 2:
    y[0][0] = 0.5;
    y[0][1] = 1.0;
*/

/*
    //SistemaTeste 3:
    y[0][0] = 1.0;
    y[0][1] = 1.0;
*/

    ofstream outputFile("output.csv");
    outputFile<<"t,V,Ap,Apm,Thn,The,Tkn,Tke,B,Ps,Pl,Bm,A"<<endl;
    outputFile.close();
    RK5(t,h,y,inter);

    for (int i = 0; i < t_store; ++i)
        delete [] y[i];
    delete [] y;

    return 0;
}
