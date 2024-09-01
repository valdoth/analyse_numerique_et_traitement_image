/*****************************************************
       NDIMBIARISOA VALDO TSIARO HASINA
         valdotsiarohasina@gmail.com

            Méthode de descente de gradient
******************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cmath>

using namespace std;

// la classe de basse
class DescenteGradient
{
private:
    float **matrice;
    float *vecteur;
    float *solution;
    int dim;
    float *x;
    float *y;

public:
    // constructeur
    DescenteGradient(string file);
    // destructeur
    ~DescenteGradient();

    // getters
    float **getMatrice();
    float *getVecteur();
    float *getSolution();
    int getDim();

    // resolution de notre equation
    void resolution();
};

// afficher la matrice
void afficheMatrice(float **matrice, float *vecteur, int dim);
// afficher solution
void afficherSolution(float *solution, int dim);
float **initialiseMatrice(int rows, int cols);
float *initialiseVec(int n);
// produit matrice
float *produitMatrice(float **matrice, float *vecteur, int dim);
// Calcul norme
float norm(float *vecteur, int dim);
// calcul scalaire
float scalaire(float *vec1, float *vec2, int dim);

int main()
{
    cout << "Méthode de DescenteGradient" << endl;

    // initialisation de la matrice a partir du fichier data.txt
    DescenteGradient matrices("data.txt");
    cout << "La matrice de départ :\n";
    afficheMatrice(matrices.getMatrice(), matrices.getVecteur(), matrices.getDim());


    // résolution de léquation linéaire
    cout << "\nLa solution de notre équation linéaire est :\n";
    matrices.resolution();
    // afficherSolution(matrices.getVecteur(), matrices.getDim());

    return 0;
}

void afficheMatrice(float **matrice, float *vecteur, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
            cout << setw(5) << matrice[i][j] << setw(5) << " ";
        cout << " = " << setw(5) << vecteur[i] << endl;
    }
}

void afficherSolution(float *solution, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        cout << "x" << i + 1 << " = " << solution[i] << endl;
    }
}

float **initialiseMatrice(int n)
{
    float **matrice(NULL);
    matrice = new float *[n];

    for (int i = 0; i < n; i++)
        matrice[i] = initialiseVec(n);
    return matrice;
}

// créer un tableau à n dimension
float *initialiseVec(int n)
{
    return new float[n];
}

DescenteGradient::DescenteGradient(string file)
{
    ifstream data(file);

    // recuperation des données dans notre fichier
    if (data)
    {
        data >> dim;
        matrice = initialiseMatrice(dim);
        vecteur = initialiseVec(dim);
        solution = initialiseVec(dim);
        y = initialiseVec(dim);
        x = initialiseVec(dim);

        // recuperer la valeur de la matrice
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                float temp = 0;
                data >> temp;
                matrice[i][j] = temp;
            }
        }

        // recuperation du vecteur
        for (int i = 0; i < dim; i++)
        {
            data >> vecteur[i];
        }
    }
    data.close();
}

DescenteGradient::~DescenteGradient()
{
    for (int i = 0; i < dim; i++) // delete colonne
        delete[] matrice[i];
    delete[] matrice;  // supprimé ligne
    delete[] vecteur;  // supprimer le tableau contenant le second membre
    delete[] solution; // supprimer le tableau contenant la solution
}

float **DescenteGradient::getMatrice()
{
    return matrice;
}
float *DescenteGradient::getVecteur()
{
    return vecteur;
}
float *DescenteGradient::getSolution()
{
    return solution;
}
int DescenteGradient::getDim()
{
    return dim;
}

// resolution de l'equation linéaire
void DescenteGradient::resolution()
{

    // algorithme
    // A = matrice
    // b = vecteur
    /*
    x quelconque donnée
    r = b - A.x
    while (|r| > eps)
        z = -A.r
        alpha = (||r||**2) / (z | r)
        x = x - alpha.r
        r = r - alpha.z
    */
    // les variables utiles
    float normeDeR;
    float epsilon = 10e-10;
    float *z = new float[dim];
    int iteration = 0;

    // initialiser vecteur x
    float *x = new float[dim];
    for(int i = 0; i < dim; i++) {
        x[i] = 1;
    }

    // calcul du résidu r = b - A.x
    float *ax = new float[dim];
    ax = produitMatrice(matrice, x, dim);
    float *r = new float[dim];
    for (int i = 0; i < dim; i++) {
        r[i] = vecteur[i] - ax[i];
    }

    // calcul du norme de r
    normeDeR = norm(r, dim);

    while (normeDeR > epsilon) {
        z = produitMatrice(matrice, r, dim);
        for (int i = 0; i < dim; i++) {
            z[i] = -z[i];
        }

        float alpha = (normeDeR*normeDeR) / scalaire(z, r, dim);
        for (int i = 0; i < dim; i++) {
            x[i] -=  alpha * r[i];
            r[i] -= alpha * z[i];
        }
        normeDeR = norm(r, dim);
        iteration++;
    }

    cout << "Après: " << iteration << " on obtient le résultat suivant: " << endl;

    for (int i = 0; i < dim; i++) {
        cout << "x" << i << " = " << x[i] << endl;
    }
}

float *produitMatrice(float **matrice, float *vecteur, int dim) {
    float *result = new float[dim];
    for(int i = 0; i < dim; i++) {
        float sum = 0;
        for (int j = 0; j < dim; j++) {
            sum += (matrice[i][j] * vecteur[j]);
        }
        result[i] = sum;
    }
    return result;
}

float norm(float *vecteur, int dim) {
    float norme = 0;
    for (int i = 0; i < dim ; i++ ) {
        norme += pow(vecteur[ i ], 2);
    }
    return sqrt(norme);
}

float scalaire(float *vec1, float *vec2, int dim) {
    float scal = 0;
    for(int i = 0; i < dim; i++) {
        scal += (vec1[i] * vec2[i]);
    }
    return scal;
}
