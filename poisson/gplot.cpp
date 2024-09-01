#include<vector>
#include<string>
#include<fstream>
#include<iostream>
using namespace std;

class Gplot {
    private:
        FILE* _pipe;
        unsigned int nbPts ;
        float xmin,xmax,ymin,ymax;
        float *xd,*yd;
        vector <string> script;
    public:
        Gplot();
        ~Gplot();
        bool isOpened() const;
        void setup(string filename , float xmin , float xmax , float ymin ,float ymax);
        void display(string xlab,string ylab , string title);
        void open();
        void flush();
        void close();
        void write(const char *line);
        void write(const string &line);
        void execute();
};

Gplot::Gplot(){
    _pipe= nullptr;
}

Gplot::~Gplot(){
    close();
}

bool Gplot::isOpened() const {
    return _pipe != nullptr;
}

void Gplot::setup(std::string filename, float xmin=0 , float xmax=1 , float ymin=-5, float ymax=5) {
	ifstream fichier(filename);
	std::string ligne;

	if(fichier.is_open()){
		fichier >> nbPts;
        // allocation et lecture des données
        xd = new float[nbPts]; yd = new float[nbPts];
        this->ymin = 1e8; this->ymax = -1e8;
        for(size_t i=0; i<nbPts ; i++){
            getline(fichier,ligne, ','); // jusqu'à la prochaine virgule
            xd[i] = atof(ligne.c_str()); // jusqu'à la fin de la ligne
            getline(fichier,ligne);
            yd[i] = atof(ligne.c_str());
            if(yd[i] < this->ymin) {
                this->ymin = yd[i];
            }
            if(yd[i] > this->ymax) {
                this->ymax = yd[i];
            }
		}
	} else {
        std::cout << "Impossible d'ouvrir le fichier " << endl;
    }
	this->xmin = xd[0]; this->xmax = xd[nbPts-1];
	if(this->ymin > ymin) this->ymin = ymin;
	if(this->ymax > ymax) this->ymax = ymax;

}

void Gplot::display(string xlab="x", string ylab="y", string title="Title") {
    open();
    script.push_back("reset");
    script.push_back("set term wxt size 640, 480");
    script.push_back("set title '" + title);
    script.push_back("set xlabel '" + xlab);
    script.push_back("set ylabel '" + ylab);
    script.push_back("set autoscale y");
    script.push_back("set nokey");
    script.push_back("set xzeroaxis");
    script.push_back("set xrange [" + to_string(xmin) + ":" + to_string(xmax) + "]");
    // script.push_back("set yrange [" + to_string(ymin) + ":" + to_string(ymax) + "]");

    // reusable internal data bloc
    script.push_back("$data << EOF");
    for (size_t i=0; i<nbPts; i++) {
        script.push_back(to_string(xd[i]) + " " + to_string(yd[i]));
    }
    script.push_back("EOF");
    script.push_back("plot $data using 1:2 w linespoints pt 5 , $data using 1:2 sm csplines");

    execute();
}

void Gplot::open() {
    close();
    _pipe = popen("/usr/bin/gnuplot -persist", "w");
}

void Gplot::flush() {
    if (isOpened()) {
        fflush(_pipe);
    }
}

void Gplot::close() {
    if (isOpened()) {
        pclose(_pipe);
        _pipe = nullptr;
    }
}

void Gplot::write(const char *line) {
    if (isOpened() && line != nullptr && line[0] != '\0') {
        fprintf(_pipe, "%s\n", line);
    }
}

void Gplot::write(const string &line) {
    if (!line.empty()) {
        write(line.c_str());
    }
}

void Gplot::execute() {
    if (isOpened()) {
        for (size_t i=0; i<script.size(); i++) {
            write(script[i]);
            flush();
        }
    }
}

int main() {
    Gplot plot;

    plot.setup("curve.txt", -4 , 4);
    plot.display("t", "f(t)", "approximation");
    plot.setup("curve1.txt", -4 , 4);
    plot.display("t", "f(t)", "initial");

    return 0;
}
