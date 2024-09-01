#ifndef DEF_GNUPLOT
#define DEF_GNUPLOT

#include "GPlot.hpp"

void GPlot::setup(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
    this->xmin = xmin; this->xmax = xmax;
    this->ymin = ymin; this->ymax = ymax;
    if(this->zmin > zmin) this->zmin = zmin;
    if(this->zmax < zmax) this->zmax = zmax;
}

void GPlot::display(std::string xlab, std::string ylab, std::string zlab, std::string title) {
    open();
    script.push_back("reset");
    script.push_back("set term wxt size 640,480");
    script.push_back("set title '" + title);
    script.push_back("set xlabel '" + xlab);
    script.push_back("set ylabel '" + ylab);
    script.push_back("set zlabel '" + zlab);
    script.push_back("set autoscale z");
    script.push_back("set nokey");
    script.push_back("set xzeroaxis");
    script.push_back("set yzeroaxis");
    script.push_back("set xrange [" + std::to_string(xmin) + ":" + std::to_string(xmax) + "]");
    script.push_back("set yrange [" + std::to_string(ymin) + ":" + std::to_string(ymax) + "]");
    script.push_back("splot 'surface.txt' u 1:2:3 w pm3d");

    execute();
}
void GPlot::open() {
    close();
    _pipe = popen("gnuplot -persist", "w");
}
void GPlot::flush() {
    if(isOpened()) {
        fflush(_pipe);
    }
}
void GPlot::close() {
    if(isOpened()) {
        pclose(_pipe);
        _pipe = nullptr;
    }
}
void GPlot::write(const char *line) {
    if(isOpened() && line != nullptr && line[0] != '\0') {
        fprintf(_pipe, "%s\n", line);
    }
}
void GPlot::write(const std::string &line) {
    if(!line.empty()) {
        write(line.c_str());
    }
}

void GPlot::execute() {
    if(isOpened()) {
        for(size_t i = 0; i < script.size(); i++) {
            write(script[i]);
            flush();
        }
    }
}

#endif //GNUPLOT