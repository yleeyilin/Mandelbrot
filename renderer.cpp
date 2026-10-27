#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

const int width = 1200;
const int height = 900;
int maxIterations = 10000;

/*
z_next = z ** 2 + c 
where z and c are complex numbers.

if abs(z) > 2, we have escaped the mandelbrot set! 
*/

struct Result {
    int iterations;
    double zr;
    double zi;
};

Result iterationsToEscapeMandelbrotSet(
    double cr, double ci, int maxIterations)
{
    double zr = 0.0;
    double zi = 0.0;
    int iteration = 0;

    while (zr * zr + zi * zi <= 4.0 &&
           iteration < maxIterations)
    {
        double nextR = zr * zr - zi * zi + cr;
        double nextI = 2.0 * zr * zi + ci;

        zr = nextR;
        zi = nextI;

        ++iteration;
    }

    return {iteration, zr, zi};
}


void colorPixel(int iterations, double zr, double zi, int maxIterations,
                int& red, int& green, int& blue)
{
    // Points inside the Mandelbrot set (did not escape) are pitch black.
    if (iterations == maxIterations) {
        red = 0;
        green = 0;
        blue = 0;
        return;
    }

    // Smooth potential calculation for the boundary
    double log_zn = std::log(zr * zr + zi * zi) / 2.0;
    double nu = std::log(log_zn / std::log(2.0)) / std::log(2.0);
    double smoothed = iterations + 1.0 - nu;

    double t = smoothed * 0.15;
    double edge = std::abs(std::sin(t));
    double thinLine = std::pow(edge, 4.0);

    // Dark indigo boundary lines on a pitch black background
    red   = static_cast<int>(40  * thinLine);
    green = static_cast<int>(25  * thinLine);
    blue  = static_cast<int>(130 * thinLine);
}



int main() {
    // Create an image file.
    std::ofstream image("mandelbrot_default.ppm");

    // PPM header: format, width, height, maximum color value.
    image << "P3\n";
    image << width << " " << height << "\n";
    image << "255\n"; // max value for each colour channel in a pixel vector 

    std::vector<int> pixels(width * height * 3);

    // Default full-set bounds scaled perfectly for 1200x900 (4:3 ratio)
    double xMin = -2.0;
    double xMax =  1.0;
    double yMin = -1.125;
    double yMax =  1.125;

    // Calculate the color of every pixel.
    #pragma omp parallel for schedule(dynamic)
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {

            double cr = xMin + (static_cast<double>(x) / (width - 1.0)) * (xMax - xMin);
            double ci = yMin + (static_cast<double>(y) / (height - 1.0)) * (yMax - yMin);

            Result res = iterationsToEscapeMandelbrotSet(cr, ci, maxIterations);

            // Next: use iterations and final z values to choose the pixel's color smoothly.
            int red, green, blue;
            colorPixel(res.iterations, res.zr, res.zi, maxIterations, red, green, blue);

            int index = (y * width + x) * 3;

            pixels[index] = red;
            pixels[index + 1] = green;
            pixels[index + 2] = blue;

        }
    }

    for (int i = 0; i < width * height * 3; ++i) {
        image << pixels[i] << ' ';

        if (i % 12 == 11) {
            image << '\n';
        }
    }

    image.close();

    std::cout << "Image saved to mandelbrot_default.ppm\n";
    return 0;
}