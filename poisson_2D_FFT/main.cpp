#include <iostream>
#include <cmath>

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr;

using namespace std;

// Fonction donnée
void sinft(double *y, int n);
void realft(double *data,int n, int isign);
void four1(double *data,int nn,int isign);

// fonction étudiée (fonction solution). 
double u(double x, double y);

// f est le Laplacien de u; Delta(U) = f
double f(double x, double y);

// Multiplcation d'une matrice avec un scalaire
double **scalarDotMatrix(double alpha, double **A, int N, int M);

// Calcul des valeurs propres de matrices B 
double *eigenvalueB(double dx, double dy, int M);

// Calcul de valeurs propres des deux matrices C 
double *eigenvalueC(int M);

//Calcul du second membre eta
double **matrixEta(double x0, double y0, double dx, double dy, int M, int N);

// La matrice H (transposé de eta)
double **matrixH(double **eta, int M, int N);

// Resolution des M systemes 
double **matrixV(int N, int M, double **h, double dx, double dy);

// La matrice U (transposé de v)
double **matrixU(double **v, int M, int N);

// La solution final du probleme
double **solution(double **U, int M, int N);

// La solution sur gnuplot
void gplot(double **sol, double xmin, double ymin, double dx, double dy, int M, int N);

// Affiche une matrice M
void displayMatrix(double **M, int rows, int col);

// Resoltion par tridiag
void tridiag(double *a, double *b, double *c, double *h, double *u, int n);

int main(){
    int N(64), M(64);
    double xmin(-3), xmax(3), ymin(-3), ymax(3),dx, dy;
    double *eigenB, *eigenC;
    double **S, **eta, **H, **U, **sol, **v;
    //Allocation dynamique
    eigenB = new double[M];
    eigenC = new double[M];

    S = new double*[N];
    for(int i = 0; i < N; i++) {
        S[i] = new double[M];
    }

    eta = new double*[N];
    for(int i = 0; i < N; i++) {
        eta[i] = new double[M];
    }

    H = new double*[M];
    for(int i = 0; i < M; i++) {
        H[i] = new double[N];
    }

    U = new double*[M];
    for(int i = 0; i < M; i++) {
        U[i] = new double[N];
    }

    v = new double*[N];
    for(int i = 0; i < N; i++) {
        v[i] = new double[M];
    }

    sol = new double*[N];
    for(int i = 0; i < N; i++) {
        sol[i] = new double[M];
    }

    dx = (xmax - xmin) / (double)M;
    dy = (ymax - ymin) / (double)N;

    cout << "Resolution du probleme de Poisson en 2D avec -Delta(u) = exp(-(x^2 + y^2)/2)" << endl;
    
    // etape de la fonction
    eigenB = eigenvalueB(dx, dy, M);
    eigenC = eigenvalueC(M);
    eta = matrixEta(xmin,ymin,dx,dy,M,N);
    // displayMatrix(eta, N, M);
    H = matrixH(eta, M, N);
    //  displayMatrix(H, N, M);
    v = matrixV(N, M, H, dx, dy);
    // displayMatrix(v,M,N);
    U = matrixU(v, M, N);
    // displayMatrix(U, N, M);
    sol = solution(U,M,N);
    // displayMatrix(sol, N, M);
    gplot(sol, xmin, ymin, dx, dy, M, N);
    return 0;
}



/*
    * La fonction solution
    * On prend sigma = 2.
    * u = exp(- (x*x + y*y) / 2).
*/
double u(double x, double y){
    return exp(- (x*x + y*y) / 2);
}

// f est le Laplacien de u; Delta(U) = f
double f(double x, double y){
    return (exp(- (x*x + y*y) / 2)) * (2 - (x*x + y*y));
}

// Multiplcation d'une matrice avec un scalaire
double **scalarDotMatrix(double alpha, double **A, int N, int M){
    double **res;
    res = new double*[N];
    for(int i = 0; i < N; i++) {
        res[i] = new double[M];
    }
    for (int i = 0; i < N; i++)
    {
        for(int j=0; j < M; j++){
            res[i][j] = A[i][j] * alpha;
        }
    }
    return res;
}

// Calcul des valeurs propres de matrices B 
double *eigenvalueB(double dx, double dy, int M){
    double *eigenB;
    eigenB = new double[M];
    for(int j = 0; j < M; j++)
        eigenB[j] = 2 + 2 * (dy / dx) * (dy / dx) * (1 - cos((M_PI * j) / (double)(M + 1)));
    return eigenB;
}

// Calcul de valeurs propres des deux matrices C 
double *eigenvalueC(int M){
    double *eigenC;
    for(int j = 0; j < M; j++ )
        eigenC[j] = -1;
    return eigenC;
}

// Affiche une matrice M
void displayMatrix(double **M, int rows, int col){
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < col; j++){
            cout << "   " << M[i][j];
        }
        cout << endl;
    }
}

// Calcul du second membre eta
double **matrixEta(double x0, double y0, double dx, double dy, int M, int N){
    double **S, **eta;
    double alpha = -1/(dx * dx);
    double gamma = -1/(dy * dy);
    S = new double*[N];
    for(int i = 0; i < N; i++) {
        S[i] = new double[M];
    }
    eta = new double*[N];
    for(int i = 0; i < N; i++) {
        eta[i] = new double[M];
    }
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            S[i][j] = f(x0 + ( dx * (j + 1)), y0 + (dy * (i + 1)));
        }
        S[i][0] -= alpha * u(x0, y0 + ( dy * (i + 1)));
        S[i][M - 1] -= alpha * u(x0 + (dx * (M)), y0 + ( dy * (i + 1)));
    }
    
    for(int j = 0; j < M; j++) {
        S[0][j] -= gamma * u(x0 + (dx * (j + 1)), y0);
        S[N-1][j] -= gamma * u(x0 + (dx * (j + 1)), y0 + (dy * M));
    }

    eta = scalarDotMatrix(( 2 / (double)M)*dy*dy , S , N , M );
    for(int j = 0; j < N; j++) {
        sinft(eta[j], M);
    }
    return eta;
}

// Calcul du H (transposé de eta)
double **matrixH(double **eta, int M, int N){
    double **H;
    H = new double*[M];
    for(int i = 0; i < M; i++) {
        H[i] = new double[N];
    }
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            H[j][i] = eta[i][j];
        }
    }
    return H;
};

// Resoltion par tridiag
void tridiag(double *a, double *b, double *c, double *h, double *u, int n)
{
    int j(0);
    double bet(b[0]), gam[n];

    u[0] = h[0] / bet;
    for (j = 1; j < n; j++)
    {
        gam[j] = c[j - 1] / bet;
        bet = b[j] - a[j] * gam[j];
        u[j] = (h[j] - a[j] * u[j - 1]) / bet;
    }
    for (j = n - 2; j >= 0; j--)
    {
        u[j] -= gam[j + 1] * u[j + 1];
    }
}

// Resolution des M systemes T[j].U[j]=h[j]
double **matrixV(int N, int M, double **h, double dx, double dy){
    double** v;
    v = new double*[N];
    for(int i = 0; i < N; i++) {
        v[i] = new double[M];
    }

    double a[N - 1];
    double b[N];
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            if (i < N - 1)
            {
                a[i] = -1;
                b[i] = eigenvalueB(dx,dy,M)[j];
            }
            else
            {
                b[i] = eigenvalueB(dx,dy,M)[j];
            }
        }
        tridiag(a, b, a, h[j], v[j], N);
    }
    return v;
}

// Calcul du U (transposé de v)
double **matrixU(double **v, int M, int N){
    double **U;
    U = new double*[M];
    for(int i = 0; i < M; i++) {
        U[i] = new double[N];
    }
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            U[j][i] = v[i][j];
        }
    }
    return U;
};

// La solution final du probleme
double **solution(double **U, int M, int N){
    double **sol;
    sol = new double*[N];
    for(int i = 0; i < N; i++) {
        sol[i] = new double[M];
    }
    for (int i = 0; i < N; i++)
    {
        sinft(U[i], M);
    }
    sol = U;
    return sol;
}
// La solution sur gnuplot
void gplot(double **sol, double xmin, double ymin, double dx, double dy, int M, int N) {
    FILE *GP = popen("gnuplot -persist", "w");

    if (GP) {
        fprintf(GP, "set term wxt title 'Probleme de Poisson 2D'\n");
        fprintf(GP, "set xlabel \"x(t)\"\n");
        fprintf(GP, "set ylabel \"y(t)\"\n");
        fprintf(GP, "set zlabel \"z(t)\"\n");
        fprintf(GP, "set surface\n");
        fprintf(GP, "set yrange [-3:3]\n");
        fprintf(GP, "set xrange [-3:3]\n");
        fprintf(GP, "set zrange [0:1.5]\n");
        fprintf(GP, "set style data points\n");
        fprintf(GP, "set decimalsign \".\"\n");
        
        fprintf(GP, "$data << EOF\n");
        int counter = 0;
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < M; j++) {
                fprintf(GP, "%f %f %f\n", xmin + (j + 1)*dx, ymin + (i + 1)*dy, sol[i][j]);
            }
        }
        fprintf(GP, "EOF\n");
        //plot u(x) = exp( - (x*x + y*y) / 2 )
        fprintf(GP, "splot  exp( - (x*x + y*y) / 2 ) with lines, $data using 1:2:3 with lines t 'Approximation'\n");
        
        /// commande GNUPLOT
        fflush(GP);
        pclose(GP);
    } 
    else {
        cout << "gnuplot not found..." << endl;
    }
}


// Fonction donnée
void sinft(double *y, int n){
    int  j,n2=n+2;
    double sum,y1,y2;
    double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;
    // cout << "Y[0] before : " << y[0];

    theta=3.14159265358979/(double) n;
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    y--;            // D�calage pour utiliser la num�rotation math
    y[1]=0.0;
    for (j=2;j<=(n>>1)+1;j++) {
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
        y1=wi*(y[j]+y[n2-j]);
        y2=0.5*(y[j]-y[n2-j]);
        y[j]=y1+y2;
        y[n2-j]=y1-y2;
    }
    realft(y,n,1);
    y[1]*=0.5;
    sum=y[2]=0.0;
    for (j=1;j<=n-1;j+=2) {
        sum += y[j];
        y[j]=y[j+1];
        y[j+1]=sum;
    }
}

void realft(double *data,int n, int isign){
    int i,i1,i2,i3,i4,np3;
    double c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta;

    theta=3.141592653589793/(double) (n>>1);
    if (isign == 1) {
        c2 = -0.5;
        four1(data,n>>1,1);
    }
    else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    np3=n+3;
    for (i=2;i<=(n>>2);i++) {
        i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    }
    else {
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        four1(data,n>>1,-1);
    }
}

void four1(double *data,int nn,int isign){
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;

    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=mmax << 1;
        theta=isign*(6.28318530717959/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}
