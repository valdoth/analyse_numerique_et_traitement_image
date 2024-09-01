/***********************************************
       NDIMBIARISOA VALDO TSIARO HASINA
         valdotsiarohasina@gmail.com

            DIMENSION FINIE 1D
***********************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cmath>

using namespace std;

// la classe de basse
class Choleski {
    private:
        float **matrice;
        float *vecteur;
        float *vecteurF;
        float *solution;
        int dim;
        float *x;
        float *y;
        float *x_final;

    public:
        // constructeur
        Choleski(string file);
        // destructeur
        ~Choleski();


        // getters
        float** getMatrice();
        float* getVecteur();
        float* getSolution();
        int getDim();

        // resolution de notre equation
        void triangularisation();
        void transposeMatrice();
        void resolution();
};

// afficher la matrice
void afficheMatrice(float **matrice, float* vecteur, int dim);
// afficher solution
void afficherSolution(float* solution, int dim);
float **initialiseMatrice(int rows, int cols);
float *initialiseVec(int n);

int main() {
    cout << "Méthode de Choleski" << endl;

    // initialisation de la matrice a partir du fichier data.txt
    Choleski matrices("data.txt");
    cout << "La matrice de départ :\n";
    afficheMatrice(matrices.getMatrice(), matrices.getVecteur(), matrices.getDim());

    // triangularisation de la matrice par la méthode de choleski
    cout << "\nLa matrice après triangularisation: " << endl;
    matrices.triangularisation();
    afficheMatrice(matrices.getMatrice(), matrices.getVecteur(), matrices.getDim());

    // résolution de léquation linéaire
    cout << "\nLa solution de notre équation linéaire est :\n";
    matrices.resolution();
    afficherSolution(matrices.getVecteur(), matrices.getDim());

    return 0;
}




void afficheMatrice(float** matrice, float* vecteur, int dim) {
    for (int i=0;i<dim;i++){
        for(int j=0; j<dim; j++)
            cout << setw(5) << matrice[i][j] << setw(5) << " ";
        cout << " = " << setw(5) << vecteur[i] << endl;
    }
}

void afficherSolution(float* solution, int dim){
    for(int i=0; i<dim; i++){
        cout << "x"<< i+1 << " = " << (double)solution[i] << endl;
    }
}


float **initialiseMatrice(int n) {
    float **matrice(NULL);
    matrice = new float*[n];

    for (int i=0;i<n;i++)
        matrice[i] = initialiseVec(n);
    return matrice;
}

// créer un tableau à n dimension
float *initialiseVec(int n) {
    return  new float [n];
}


Choleski::Choleski(string file) {
    double h = 0;
    double a,b;
    int N;
    cout << "Entrer borne inferieur de l'intervalle: " ;
    cin >> a;
    cout << "Entrer borne superieur de l'intervalle: ";
    cin >> b;
    cout << "Nombre de decomposition de l'intervalle: ";
    cin >> N;
    h = (double)(b - a) / (N + 1);

    dim = N;
    matrice = initialiseMatrice(dim);
    vecteur = initialiseVec(dim);
    vecteurF = initialiseVec(dim);
    solution = initialiseVec(dim);
    y = initialiseVec(dim);
    x = initialiseVec(dim);
    x_final = initialiseVec(dim);

    // generate the default matrice A
    double pas = a + h;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            double n = 0;
            if (i == j) {
                n = 2;
            } else if( i== j-1 || i == j+1){
                n = -1;
            }else{
                n = 0;
            }
            matrice[i][j] = n / (h * h);
        }
        x_final[i] = pas;
        vecteurF[i] = (2 - 4 * pas * pas) * (exp(-pas * pas));
        vecteur[i] = (2 - 4 * pas * pas) * (exp(-pas * pas));
        pas += h;
    }
    afficheMatrice(matrice, vecteur, dim);
}

Choleski::~Choleski() {
    for(int i = 0; i < dim; i++)    //delete colonne
        delete[] matrice[i];
    delete[] matrice;               //supprimé ligne
    delete[] vecteur;                     //supprimer le tableau contenant le second membre
    delete[] solution;              //supprimer le tableau contenant la solution
}

float** Choleski::getMatrice() {
    return matrice;
}
float* Choleski::getVecteur() {
    return vecteur;
}
float* Choleski::getSolution() {
    return solution;
}
int Choleski::getDim() {
    return dim;
}

void Choleski::triangularisation() {
    int n = dim;

    float ** lower;
    lower = initialiseMatrice(dim);

    float sum = 0;
    // utilisation de l'agorithme vu en cours
	for(int i(0); i<dim ; i++){

		for(int j(0); j<dim ; j++){
				sum = 0;
				if(j>i){
					matrice[i][j]=0;
				}
				else if(i==j){
					for(int k(0); k<j ; k++){
						sum += matrice[i][k]*matrice[i][k];
					}
					matrice[i][i] = sqrt(matrice[i][i] - sum);

				}else{
					for(int k(0);k<j ; k++){
						sum += matrice[i][k]*matrice[j][k];
					}
					matrice[i][j] = 1.0 / matrice[j][j] * (matrice[i][j] - sum);
			}

		}
	}

}

void Choleski::transposeMatrice() {
    float temp(0);
    // echange les lignes et les colonnes
    for(int i = 0; i<dim; i++) {
        for(int j = i; j<dim; j++) {
            temp = matrice[i][j];
            matrice[i][j] = matrice[j][i];
            matrice[j][i] = temp;
        }
    }
}


// resolution de la matrice triangulaire
void Choleski::resolution() {
     // Résolution de l'équation B * y = b, on obtient y
    for(int i=0;i<dim;i++){
		float s=0;
		for(int j=0;j<i;j++){
			s += matrice[i][j] * y[j];
		}
 		y[i]= ((vecteur[i]-s)/matrice[i][i]);
	}

    // On transpose la matrice B
    transposeMatrice();

	// On résout Bt * x = y, on obtient le résultat x
    for(int i=dim-1;i>=0;i--){
		float s=0;
		for(int j=i;j<dim;j++){
			s += matrice[i][j] * x[j];
		}
 		x[i]= (y[i]-s)/matrice[i][i];
 	}

    for(int i=0; i<dim; i++) {
        vecteur[i] = x[i];
    }

    // create file to save the result
    ofstream MyFile("curve.txt");
    MyFile << dim - 1 << endl;
    for (int i=0; i<dim; i++) {
        MyFile << x_final[i] << "," << vecteur[i] << endl;
    }
    ofstream MyFile1("curve1.txt");
    MyFile1 << dim - 1 << endl;
    for (int i=0; i<dim; i++) {
        MyFile1 << x_final[i] << "," << vecteurF[i] << endl;
    }
    MyFile.close();
    MyFile1.close();
}

double exp_funct (double x) {
    return exp( -1 * pow(x,2));
}
