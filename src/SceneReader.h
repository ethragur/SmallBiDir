#pragma once

#include "scene.h"
#include "sphere.h"

#include <iostream>
#include <fstream>
#include <string>


class SceneReader
{
	public:
		void ReadScene(Scene & scene, const std::string & name)
		{
			std::ifstream f(name);
			std::cout << "Loading " << name << "..." << std::endl;
			if(f.is_open())
			{
				int i = 0;
				while ( !f.eof() )
				{
					i++;
					double radius;
					double posx;
					double posy;
					double posz;
					double emir;
					double emig;
					double emib;
					double colr;
					double colg;
					double colb;
					std::string reflection;

					f >> radius;
					f >> posx;
					f >> posy;
					f >> posz;
					f >> emir;
					f >> emig;
					f >> emib;
					f >> colr;
					f >> colg;
					f >> colb;
					f >> reflection;

					Refl_t refl;
					if(reflection.compare("DIFF") == 0)
							refl = DIFF;
					else if(reflection.compare("SPEC") == 0)
							refl = SPEC;
					else if(reflection.compare("REFR") == 0)
							refl = REFR;
					else if(reflection.compare("GLOSSY") == 0)
							refl = GLOSSY;
					else if(reflection.compare("TRANSL") == 0)
							refl = TRANSL;
					else
							refl = DIFF;

					Sphere *toAdd = new Sphere(radius, Vector(posx,posy,posz), Vector(emir, emig, emib), Vector(colr, colg, colb), refl);
					scene.add_sphere(toAdd);
					int tmp;
					//fix for reading
					f >> tmp;
				}
				f.close();
				std::cout << "# of spheres: " << i << std::endl;
			}

		}
};
