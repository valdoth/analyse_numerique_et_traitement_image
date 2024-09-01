#include <iostream>
#include <algorithm>
#include <queue>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

// la classe de basse
class Cuthill
{
private:
  vector<vector<int>> matrice;
  vector<int> vecteur;
  vector<int> order;
  int dim;

public:
  // constructeur
  Cuthill(string file);
  // destructeur
  ~Cuthill();

  // getters
  vector<vector<int>> getMatrice();
  vector<int> getVecteur();
  int getDim();

  // sauver le profil dans un fichier
  void sauverCuthill(vector<vector<int>> matrix);

  // affiche matrice
  void afficheMatrice();
};

// Structure représentant un sommet dans le graphe
struct Vertex
{
  int index;  // Indice du sommet
  int degree; // Degré du sommet

  bool operator<(const Vertex &other) const
  {
    return degree < other.degree;
  }
};

// Fonction pour réorganiser une matrice creuse avec la méthode Cuthill-McKee
vector<int> cuthillMckee(const std::vector<std::vector<int>> &adjacencyMatrix);

// Fonction utilitaire pour afficher une matrice
void printMatrix(const std::vector<std::vector<int>> &matrix);

int main()
{
  cout << "Création de profil pour la resolution de grands systeme lineaire." << endl
       << endl;

  // initialisation de la matrice a partir du fichier data.txt
  Cuthill matrices("data.txt");

  vector<vector<int>> matrix;
  matrix = matrices.getMatrice();

  // Réarrangement avec Cuthill-McKee
  vector<int> permutation = cuthillMckee(matrix);

  cout << "Matrice réarrangée : " << endl;
  vector<vector<int>> rearrangedMatrix(matrix.size(), std::vector<int>(matrix.size()));
  for (int i = 0; i < matrix.size(); ++i)
  {
    for (int j = 0; j < matrix.size(); ++j)
    {
      rearrangedMatrix[i][j] = matrix[permutation[i]][permutation[j]];
    }
  }
  printMatrix(rearrangedMatrix);

  matrices.sauverCuthill(matrix);

  cout << "\nNotre nouveau matrice est sauver dans le fichier \"cuthill.txt\"" << endl;

  return 0;
}

void afficheMatrice(int **matrice, int *vecteur, int dim)
{
  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
      cout << setw(5) << matrice[i][j] << setw(5) << " ";
    cout << " = " << setw(5) << vecteur[i] << endl;
  }
}

Cuthill::Cuthill(string file)
{
  ifstream data(file);

  // recuperation des données dans notre fichier
  if (data)
  {
    data >> dim;

    // recuperer la valeur de la matrice
    for (int i = 0; i < dim; i++)
    {
      vector<int> tmp;
      for (int j = 0; j < dim; j++)
      {
        int temp = 0;
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

Cuthill::~Cuthill() {}

vector<vector<int>> Cuthill::getMatrice()
{
  return matrice;
}
vector<int> Cuthill::getVecteur()
{
  return vecteur;
}
int Cuthill::getDim()
{
  return dim;
}

void Cuthill::sauverCuthill(vector<vector<int>> matrix)
{
  ofstream fichier("cuthill.txt");
  if (fichier)
  {
    fichier << dim << endl;
    for (int i = 0; i < matrix.size(); i++)
    {
      for (int j = 0; j < matrix[i].size(); j++)
      {
        fichier << matrix[i][j] << "\t";
      }
      fichier << endl;
    }
    for (int i = 0; i < vecteur.size(); i++)
    {
      fichier << vecteur[i] << endl;
    }
  }
}

void Cuthill::afficheMatrice()
{
  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j <= dim; j++)
      cout << setw(2) << matrice[i][j] << setw(2) << " ";
    for (int j = i + 1; j <= dim; j++)
      cout << setw(2) << " ";
    cout << " = " << setw(2) << vecteur[i] << endl;
  }
}

// Fonction utilitaire pour afficher une matrice
void printMatrix(const std::vector<std::vector<int>> &matrix)
{
  for (const auto &row : matrix)
  {
    for (int element : row)
    {
      std::cout << element << " ";
    }
    std::cout << std::endl;
  }
}

// Fonction pour réorganiser une matrice creuse avec la méthode Cuthill-McKee
std::vector<int> cuthillMckee(const std::vector<std::vector<int>> &adjacencyMatrix)
{
  int n = adjacencyMatrix.size();

  std::vector<int> permutation(n);     // Tableau de permutation
  std::vector<bool> visited(n, false); // Tableau pour marquer les sommets visités
  std::vector<Vertex> vertices(n);     // Tableau de sommets

  // Calculer le degré de chaque sommet
  for (int i = 0; i < n; ++i)
  {
    int degree = 0;
    for (int j = 0; j < n; ++j)
    {
      if (adjacencyMatrix[i][j] != 0)
      {
        ++degree;
      }
    }
    vertices[i] = {i, degree};
  }

  // Trier les sommets par degré croissant
  std::sort(vertices.begin(), vertices.end());

  // File d'attente pour le parcours en largeur
  std::queue<int> queue;

  int index = 0; // Indice pour la permutation

  // Parcours en largeur à partir de chaque sommet non visité
  for (int i = 0; i < n; ++i)
  {
    if (!visited[vertices[i].index])
    {
      queue.push(vertices[i].index);
      visited[vertices[i].index] = true;

      while (!queue.empty())
      {
        int vertex = queue.front();
        queue.pop();
        permutation[index++] = vertex;

        // Ajouter les voisins non visités à la file d'attente
        for (int j = 0; j < n; ++j)
        {
          if (adjacencyMatrix[vertex][j] != 0 && !visited[j])
          {
            queue.push(j);
            visited[j] = true;
          }
        }
      }
    }
  }

  return permutation;
}
