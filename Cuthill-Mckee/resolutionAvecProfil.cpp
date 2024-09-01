/************************************************************
       NDIMBIARISOA VALDO TSIARO HASINA
         valdotsiarohasina@gmail.com

  Resolution des grands système matricielle avec son profil
*************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

// la classe de basse
class Solve
{
private:
  vector<float> ap;
  vector<float> ndiag;
  vector<float> vecteur;
  vector<float> colonnes;
  vector<float> x;
  int dim;

public:
  // constructeur
  Solve(string file);
  // destructeur
  ~Solve();

  // getters
  vector<vector<float>> getMatrice();
  vector<float> getVecteur();
  int getDim();

  // factorisation matrice
  void factorisationMatrice();

  // affiche le vecteur
  void afficheVecteur();

  // resolution supérieur et inferieure
  void solveSuperieur();
  void solveInferieur();

  // affiche solution
  void afficheSolution();
};

int main()
{
  cout << "Resolution de grands systeme lineaire." << endl
       << endl;

  // initialisation de la matrice a partir du fichier profil.txt
  Solve matrices("profil.txt");

  // factorisation de notre matrice
  matrices.factorisationMatrice();

  // resolution de la matrice inférieur
  matrices.solveInferieur();
  cout << "\nLa valeur de x apres la resolution par rapport a la matrice inferieur" << endl;
  matrices.afficheVecteur();

  // resolution de la matrice supérieur
  matrices.solveSuperieur();
  cout << "\nLa valeur de x apres la resolution par rapport a la matrice superieur" << endl;
  matrices.afficheVecteur();

  // La solution de notre matrice
  cout << endl;
  matrices.afficheSolution();

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

Solve::Solve(string file)
{
  ifstream data(file);

  // recuperation des données dans notre fichier
  if (data)
  {
    float n;
    data >> n;
    float temp;
    for (int i = 0; i < n; i++)
    {
      data >> temp;
      ap.push_back(temp);
    }
    data >> dim;
    for (int i = 0; i < dim; i++)
    {
      x.push_back(0);
      data >> temp;
      ndiag.push_back(temp);
    }
    for (int i = 0; i < dim; i++)
    {
      data >> temp;
      vecteur.push_back(temp);
    }
    for (int i = 0; i < dim; i++)
    {
      data >> temp;
      colonnes.push_back(temp);
    }
  }
  data.close();
}

Solve::~Solve() {}

vector<float> Solve::getVecteur()
{
  return vecteur;
}
int Solve::getDim()
{
  return dim;
}

void Solve::factorisationMatrice()
{
  float somme = 0;
  for (int i = 0; i < dim; i++)
  {
    for (int j = colonnes[i]; j < i; j++)
    {

      for (int k = colonnes[j]; k < j; k++)
      {
        if (k >= colonnes[i])
        {
          somme += ap[ndiag[i] - i + k] * ap[ndiag[k]] * ap[ndiag[j] - j + k];
        }
      }
      cout << endl;
      ap[ndiag[i] - i + j] = 1.0 / ap[ndiag[j]] * (ap[ndiag[i] - i + j] - somme);
      somme = 0;
    }
    for (int k = colonnes[i]; k < i; k++)
    {
      somme += ap[ndiag[k]] * ap[ndiag[i] - i + k] * ap[ndiag[i] - i + k];
    }
    ap[ndiag[i]] -= somme;
    somme = 0;
  }
}

void Solve::solveInferieur()
{
  float s = 0;

  for (int i = 0; i < dim; i++)
  {
    s = 0;
    for (int j = colonnes[i]; j < i; j++)
    {
      s += ap[ndiag[i] - i + j] * x[j];
    }
    cout << endl;
    x[i] = vecteur[i] - s;
  }
}

void Solve::solveSuperieur()
{
  float s = 0;
  for (int i = 0; i < dim; i++)
  {
    x[i] /= ap[ndiag[i]];
  }
  for (int i = dim - 1; i >= 0; i--)
  {
    s = 0;
    for (int j = i + 1; j < dim; j++)
    {
      if (colonnes[j] <= i)
        s += ap[ndiag[j]-j+i] * x[j];
    }
    x[i] = x[i] - s;
  }
}

void Solve::afficheSolution()
{
  cout << "La solution de notre equation est: " << endl;

  for (int i(0); i < dim; i++)
  {
    cout << "x" << i << " = " << round(x[i]) << endl;
  }
}

void Solve::afficheVecteur()
{
  for (int i(0); i < dim; i++)
  {
    cout << "x" << i << " = " << x[i] << endl;
  }
}
