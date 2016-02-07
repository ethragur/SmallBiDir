#pragma once

#include "scene.h"
#include "sphere.h"

#include <iostream>
#include <fstream>
#include <string>


class SceneReader
{
	public:
		void ReadScene(const Scene & scene, const std::string & name)
		{
			std::ifstream f(name);
			if(f.isOpen())
			{
				while ( !f.eof() )
				{
					double radius;
					Vector position;
					Vector emission;
					Vector color;
					std::string reflection;

					f >> radius;
					f >> position.x;
					f >> position.y;
					f >> position.z;
					f >> emission.x;
					f >> emission.y;
					f >> emission.z;
					f >> color.x;
					f >> color.y;
					f >> color.z;
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

					scene.add_sphere(new Sphere(radius, position, emission, color, refl));
					
				}
				myfile.close();
			}

		}
};
