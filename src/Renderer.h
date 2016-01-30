#pragma once

#include "ray.h"
#include "PathTracer.h"
#include "Image.h"

class Renderer
{
	public: 
		/* Initialize Renderer used for Rendering the image.
		 * Width Heigt: pixels of the image
		 * Samples: # of Samples used per pixel
		 * enableDof: enable depth of field: Needs Samples > 1
		 * */
		Renderer(int _width, int _height, int _samples, bool _enableDoF = false) :
			width(_width), height(_height), samples(_samples), enableDoF(_enableDoF)
		{
			/* Set camera origin and viewing direction (negative z direction) */
			camera = Ray(Vector(0.0, 52.0, 295.6), Vector(0.0, -0.042612, -1.0).Normalized());

			/* Set depth of field variables aperture, focalLength and number of samples */
			aperture = 4;
			focalLength = 160;

			/* Image edge vectors for pixel sampling */
			cx = Vector(width * 0.5135 / height);
			cy = (cx.Cross(camera.dir)).Normalized() * 0.5135;
		}

		//Render the scene with the pathtracer p
		Image render_scene(PathTracer * p)
		{
			Image img(width, height);

			/* Loop over image rows */
#pragma omp parallel for
			for (int y = 0; y < height; y++)
			{
				srand(y * y * y);

				/* Loop over row pixels */
				for (int x = 0; x < width; x++)
				{
					img.setColor(x, y, Color());

					/* 2x2 subsampling per pixel */
					for (int sy = 0; sy < 2; sy++)
					{
						for (int sx = 0; sx < 2; sx++)
						{
							Color accumulated_radiance = Color();

							/* Compute radiance at subpixel using multiple samples */
							for (int s = 0; s < samples; s++)
							{
								const double r1 = 2.0 * drand48();
								const double r2 = 2.0 * drand48();

								/* Transform uniform into non-uniform filter samples */
								double dx;
								if (r1 < 1.0)
									dx = sqrt(r1) - 1.0;
								else
									dx = 1.0 - sqrt(2.0 - r1);

								double dy;
								if (r2 < 1.0)
									dy = sqrt(r2) - 1.0;
								else
									dy = 1.0 - sqrt(2.0 - r2);

								/* Ray direction into scene from camera through sample */
								Vector dir = cx * ((x + (sx + 0.5 + dx) / 2.0) / width - 0.5) +
											 cy * ((y + (sy + 0.5 + dy) / 2.0) / height - 0.5) +
											 camera.dir;

								/* Extend camera ray to start inside box */
								Vector start = camera.org + dir * 130.0;

								dir = dir.Normalized();

								Ray ray = Ray(start, dir);

								if(enableDoF)
								{
									//DoF
									double u1 = (drand48() * 2.0) - 1.0;
									double u2 = (drand48() * 2.0) - 1.0;

									double fac = (double) (2 * M_PI * u2);

									Vector offset = aperture * Vector(u1 * cos(fac), u1 * sin(fac), 0.0);
									Vector focalPlaneIntersection = ray.org + ray.dir * (focalLength / camera.dir.Dot(ray.dir));
									ray.org = ray.org + offset;
									ray.dir = (focalPlaneIntersection - ray.org).Normalized();

								}
									/* Accumulate radiance */
									accumulated_radiance = accumulated_radiance +
									p->calculatePixelColor(ray) / samples;// / dof_samples;
							}

							accumulated_radiance = accumulated_radiance.clamp() * 0.25;

							img.addColor(x, y, accumulated_radiance);
						}
					}
				}
			}

			return img;
		}

	private:
		int width, height;
		int samples;
		bool enableDoF;
		Ray camera;
		double aperture;
		double focalLength;
		Vector cx, cy;

};
