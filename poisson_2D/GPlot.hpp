//Gnuplot pipe
//pipe permet de faire un interaction entre deux programme
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

class GPlot {
    public:
        
        GPlot(){_pipe = nullptr; }
        virtual ~GPlot(){ close(); }
        bool isOpened() const { return _pipe != nullptr; }

        void setup(float xmin=0, float xmax=1, float ymin=0, float ymax=1, float zmin=5, float zmax=5);
        void display(std::string xlab="x", std::string ylab="y", std::string zlab="z", std::string title="Title");
        void open();
        void flush();
        void close();
        void write(const char *line);
        void write(const std::string &line);
        void execute();

    private:
        FILE* _pipe;
        float xmin, xmax, ymin, ymax, zmin, zmax;
        std::vector<std::string> script;

};