#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>
#include<cstdio>
#include<iomanip>
double pi = 3.141592653589793;


//Définition de f :
double f(double x, double y) {
    return 0;
}

//Définition de V :
double V(double y, int b) {
    return 1 - cos((2 * pi * y) / b);
}

//fonction u:
double u(double x, double y) {
    return sin(2 * pi * x) * sin(2 * pi * y);
}

//fonction f :
double fxy(double x, double y) {
    return -8 * pow(pi, 2) * sin(2 * pi * x) * sin(2 * pi * y);
}

int main() {

    //Initialisation des paramétres :

    //int N = 12 ;
    //int M = 12;                                                          
    //int L = 100000;  
    int N = 12;
    int M = 12;
    int L = 100;
    double dx = 1. / (N - 1);
    double dy = 1. / (M - 1);
    double u0 = 0.5;
    double alpha = 0.5;
    double b = 1.;
    double a = 1.;



    //Allocation de la mémoire :
    double* sol = new double[N * M];
    double* solNew1 = new double[N * M];
    double* solNew2 = new double[N * M];
    double* U0 = new double[N];
    double* U1 = new double[N];
    double* U2 = new double[N];
    double* fx = new double[N * M];
    double* f_xy = new double[N * M];
    double* uex = new double[N * M];
    double err = 0;


    //Initialisation de fx :    pour f=0 (dans le problème)
    for (int i = 0; i < N * M; i++) {
        fx[i] = 0;
    }

    //Initialisation fxy :
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            f_xy[i * M + j] = fxy(i * dx, j * dy);
        }
    }

    ///Définition des conditions aux bord :


    //Initialisation U0 :
    for (int i = 0; i < N; i++) {
        U0[i] = u0;
    }

    //Initialisation U1 : (qui égale à U0 si N=M )
    for (int i = 0; i < M; i++) {
        U1[i] = u0 * (1 + alpha * V(i * dy, b));
    }

    //Initialisation de U2 : 
    for (int i = 0; i < M; i++) {
        U2[i] = u0;
    }

    //Initialisation de u(x,0)=U0  :
    for (int i = 0; i < N; i++) {

        solNew1[i * M] = U0[i];
        solNew2[i * M] = U0[i];
        sol[i * M] = U0[i];
        uex[i * M] = U0[i];
    }

    //Initialisation de u(x,b)=U0 : 
    for (int i = 1; i < N; i++) {
        solNew1[(M * i) - 1] = U0[i];
        solNew2[(M * i) - 1] = U0[i];
        sol[(M * i) - 1] = U0[i];
        uex[(M * i) - 1] = U0[i];
    }

    //Initialisation de u(a,y)=U0 :
    for (int i = 0; i < M; i++) {
        solNew1[(N - 1) * M + i] = U2[i];
        solNew2[(N - 1) * M + i] = U2[i];
        sol[(N - 1) * M + i] = U2[i];
        uex[(N - 1) * M + i] = U2[i];
    }

    //Initialisation de u(0,y)=U0(1+alpha*V(y)) :
    for (int i = 0; i < M; i++) {
        solNew1[i] = U1[i];
        solNew2[i] = U1[i];
        sol[i] = U1[i];
        uex[i] = U1[i];
    }

    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < M - 1; j++) {
            sol[i * (M)+j] = 0;
        }
    }

    //Méthode de Jacobi avec différence finis :


    double coef1, coef2, coef;
    coef1 = 1 / pow(dx, 2);
    coef2 = 1 / pow(dy, 2);
    coef = 1 / (2 * coef1 + 2 * coef2);
    for (int l = 0; l < L; l++) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < M - 1; j++) {
                solNew1[i * M + j] = coef * (coef1 * (sol[i * M + j - 1] + sol[i * M + j + 1]) + coef2 * (sol[(i + 1) * M + j] + sol[(i - 1) * M + j]) - f_xy[i * M + j]);
                uex[i * M + j] = u(i * dx, j * dy);
            }
        }

        double* tmp;
        tmp = sol;
        sol = solNew1;
        solNew1 = tmp;


    }
    //Méthode de Gauss avec Différences finis :
    for (int l = 0; l < L; l++) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < M - 1; j++) {
                solNew2[i * M + j] = coef * (coef1 * (solNew2[i * M + j - 1] + sol[i * M + j + 1]) + coef2 * (sol[(i + 1) * M + j] + solNew2[(i - 1) * M + j]) - f_xy[i * M + j]);
                uex[i * M + j] = u(i * dx, j * dy);
            }
        }
        double* tmp;
        tmp = sol;
        sol = solNew2;
        solNew2 = tmp;
    }



    //Affichage de la solution exacte :
    std::cout << " la solution exacte est :" << std::endl;
    std::cout << '[';
    for (int i = 0; i < N * M; i++) {
        std::cout << uex[i] << ',';
    }
    std::cout << ']' << std::endl;

    std::cout << "#########################################################################" << std::endl;
    std::cout << "La méthode de Jacobie : \n" << std::endl;
    //Calcul de l'erreur de Jacobi :
    for (int i = 0; i < N * M; i++) {
        err += pow(uex[i] - solNew1[i], 2);
    }
    err = sqrt(err);


    //Affichage de la solution approchée :
    std::cout << " la solution approche est :" << std::endl;
    std::cout << '[';
    for (int i = 0; i < N * M; i++) {
        std::cout << solNew1[i] << ',';
    }
    std::cout << ']' << std::endl;

    //Affichage de l'erreur :
    std::cout << "\n l\'erreur est : \n" << err << std::endl;


    // Print solution
    std::ofstream file;
    file.open("jacobi.dat");
    for (int n = 0; n < N * M; n++) {
        file << solNew1[n] << std::endl;
    }
    file.close();


    std::cout << "############################################################################" << std::endl;
    std::cout << "La methode de Gauss" << std::endl;
    //Calcul de l'erreur de Gauss :
    for (int i = 0; i < N * M; i++) {
        err += pow(uex[i] - solNew2[i], 2);
    }
    err = sqrt(err);


    //Affichage de la solution approchée :
    std::cout << "la solution approche est :" << std::endl;
    std::cout << '[';
    for (int i = 0; i < N * M; i++) {
        std::cout << solNew2[i] << ',';
    }
    std::cout << ']' << std::endl;

    //Affichage de l'erreur :
    std::cout << "\n l\'erreur est : " << err << std::endl;


    // Print solution
    std::ofstream file;
    file.open("Gauss.dat");
    for (int n = 0; n < N * M; n++) {
        file << solNew1[n] << std::endl;
    }
    file.close();

    // Memory deallocation
    delete[] sol;
    delete[] solNew1;
    delete[] solNew2;

    return 0;
}
