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

    // Smooth coloring algorithm to remove harsh color banding
    double log_zn = std::log(zr * zr + zi * zi) / 2.0;
    double nu = std::log(log_zn / std::log(2.0)) / std::log(2.0);
    double smoothed = iterations + 1.0 - nu;

    double t = smoothed * 0.05;

    // Much darker, moody bronze/gold outlines on a dark background
    red   = static_cast<int>(35  + 30  * std::sin(t * 1.2));
    green = static_cast<int>(18  + 15  * std::sin(t * 1.2 + 1.5));
    blue  = static_cast<int>(8   + 7   * std::sin(t * 0.8 + 3.0));
}



int main() {
    // Create an image file.
    std::ofstream image("mandelbrot_dark_gold.ppm");

    // PPM header: format, width, height, maximum color value.
    image << "P3\n";
    image << width << " " << height << "\n";
    image << "255\n"; // max value for each colour channel in a pixel vector 

    std::vector<int> pixels(width * height * 3);

    // Zoomed-in coordinate bounds for the spiral region (maintaining 4:3 ratio)
    double xMin = -0.744;
    double xMax = -0.743;
    double yMin =  0.1165;
    double yMax =  0.11725;

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

    std::cout << "Image saved to mandelbrot_dark_gold.ppm\n";
    return 0;
}