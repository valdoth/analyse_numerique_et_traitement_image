/*****************************************************
       NDIMBIARISOA VALDO TSIARO HASINA
         valdotsiarohasina@gmail.com

	Profillage des grands système matricielle
******************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

// la classe de basse
class Profil
{
private:
    vector<vector<float>> matrice;
    vector<float> vecteur;
    vector<float> solution;
    vector<float> ap;
    vector<float> ndiag;
    vector<float> colonnes;
    int dim;

public:
    // constructeur
    Profil(string file);
    // destructeur
    ~Profil();

    // getters
    vector<vector<float>> getMatrice();
    vector<float> getVecteur();
    int getDim();

    // créer notre profil
    void profillage();
    // sauver le profil dans un fichier
    void sauverProfil();

    // affiche matrice
    void afficheMatrice();
};

int main()
{
    cout << "Création de profil pour la resolution de grands systeme lineaire." << endl << endl;

    // initialisation de la matrice a partir du fichier data.txt
    Profil matrices("cuthill.txt");
    
    // affichage de notre matrice
    matrices.afficheMatrice();

    // profillage
    matrices.profillage();

    // Sauver le profil dans un fichier profil.txt
    matrices.sauverProfil();

    cout << "\nNotre profil est sauver dans le fichier \"profil.txt\"" << endl;

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


Profil::Profil(string file)
{
    ifstream data(file);

    // recuperation des données dans notre fichier
    if (data)
    {
        data >> dim;

        // recuperer la valeur de la matrice
        for (int i = 0; i < dim; i++)
        {
            vector<float> tmp;
            for (int j = 0; j < dim; j++)
            {
                float temp = 0;
                data >> temp;
                tmp.push_back(temp);
            }
            matrice.push_back(tmp);
        }

        // recuperation du vecteur
        for (int i = 0; i < dim; i++)
        {
            int temp;
            data >> temp;
            vecteur.push_back(temp);
        }
    }
    data.close();
}

Profil::~Profil() {}

vector<vector<float>> Profil::getMatrice()
{
    return matrice;
}
vector<float> Profil::getVecteur()
{
    return vecteur;
}
int Profil::getDim()
{
    return dim;
}

// resolution de l'equation linéaire
void Profil::profillage()
{
    int n = -1;
    for(int i = 0; i < dim; i++){
        int k = 0;
        for(int j = 0; j <= i; j++){
            if(matrice[i][j] != 0){
                k = j;
                goto here;
            }
        }
        here:
            for(int l = k; l <= i; l++){
                ap.push_back(matrice[i][l]);
                n++;
                if(l==i){
                    ndiag.push_back(n);
                }
            }
        
    }
    for (int i=0; i<dim; i++) {
        bool stop = true;
        for (int j=0; j<dim; j++) {
            if (matrice[i][j] != 0 && stop) {
                colonnes.push_back(j);
                stop = false;
            }
        }
    }
}

void Profil::sauverProfil() 
{
    ofstream fichier("profil.txt");
    if(fichier){
        fichier << ap.size() << endl;
        for(int i = 0; i < ap.size(); i++){
            fichier << ap[i] << "\t";
        }
        fichier << endl;
        fichier << ndiag.size() << endl;
        for(int i = 0; i < ndiag.size(); i++){
            fichier << ndiag[i] << "\t";
        }
        fichier << endl;
        for(int i = 0; i < ndiag.size(); i++){
            fichier << vecteur[i] << "\t";
        }
        fichier << endl;
        for(int i = 0; i < ndiag.size(); i++){
            fichier << colonnes[i] << "\t";
        }
    }
}

void Profil::afficheMatrice()
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j <= i; j++)
            cout << setw(5) << matrice[i][j] << setw(5) << " ";
        for (int j = i+1; j <= dim; j++)
            cout << setw(5) << " " << setw(5) << " ";
        cout << " = " << setw(5) << vecteur[i] << endl;
    }
}
