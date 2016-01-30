/******************************************************************
*
* PathTracing.cpp
*
* Description: This program demonstrates global illumination rendering
* based on the path tracing method. The intergral in the rendering
* equation is approximated via Monte-Carlo integration; explicit
* direct lighting is included to improve quality; the rendered image
* is saved in PPM format.
*
* The code is largely based on the software smallpt by Kevin Beason,
* released under the MIT License.
*
* Advanced Computer Graphics Proseminar WS 2015
*
* Interactive Graphics and Simulation Group
* Institute of Computer Science
* University of Innsbruck
*
*******************************************************************/

/* Standard includes */
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>

#include "Renderer.h"
#include "shape.h"
#include "Image.h"
#include "PathTracer.h"
#include "BiDirectionalPT.h"


/******************************************************************
* Main routine: Computation of path tracing image (2x2 subpixels)
* Key parameters
* - Image dimensions: width, height
* - Number of samples per subpixel (non-uniform filtering): samples
* Rendered result saved as PPM image file
*******************************************************************/
int main(int argc, char *argv[])
{
	std::string scene_filename = "../data/scene_new.obj";
	int samples = 1;
	time_t start, end;


	if (argc == 2)
		samples = atoi(argv[1]);
	else if (argc == 3)
	{
		samples = atoi(argv[1]);
		scene_filename = argv[2];
	}
	else
		std::cout << "Paramaters: <Nr. of Subsamples> <Scenefile>" << std::endl << "Usage Example: ./PathTracing 2 scene.obj" << std::endl;

	BiDirectionalPT p(scene_filename);
	//PathTracer p(scene_filename);
	Renderer r(1024, 768, samples);

	std::cout << "Starting to render..." << std::endl;
	time (&start);
	r.render_scene(&p).Save(std::string("image.ppm"));
	time(&end);

	std::cout << "Done Rendering, took: " << difftime(end, start) << " seconds." << std::endl;


}
