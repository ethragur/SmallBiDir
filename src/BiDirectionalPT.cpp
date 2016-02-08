#include "BiDirectionalPT.h"
/* constructor with scene name and the number of time the lightray should bounce through the scene */
BiDirectionalPT::BiDirectionalPT(const std::string & scene_name, unsigned int lightbounces) : PathTracer(scene_name, true)
{
	lightBounces = lightbounces;
	getLightEmitters();
}


/* constructor with the builded scene and the number of time the lightray should bounce through the scene */
BiDirectionalPT::BiDirectionalPT(const Scene & scene, unsigned int lightbounces) : PathTracer(scene, true)
{
	lightBounces = lightbounces;
	getLightEmitters();
}

BiDirectionalPT::BiDirectionalPT(const std::string & scene_name) : PathTracer(scene_name, true)
{
	lightBounces = 2;
	getLightEmitters();
}


BiDirectionalPT::BiDirectionalPT(const Scene & scene) : PathTracer(scene, true)
{
	lightBounces = 2;
	getLightEmitters();
}

/* save all the light emitters of the scene in an array */
void BiDirectionalPT::getLightEmitters()
{
	//save all light emitters and determine their size
	for (int i = 0; i < numObjects; i++)
	{
		Shape * cur_obj = our_scene.get_shapes()[i];
		if (cur_obj->emission.x <= 0 && cur_obj->emission.y <= 0 && cur_obj->emission.z <= 0)
			continue; /* Skip objects that are not light sources */
		else
		{
			lightEmitters[i] = cur_obj;
		}
	}
	noLights = lightEmitters.size();
}

/*  start tracing with path tracing */
Color BiDirectionalPT::calculatePixelColor(const Ray & ray) const
{
	std::vector<LightPath> allLightpaths = traceLightRays(lightBounces);
	return Radiance(ray, allLightpaths);
}


/* Shoot a ray from every light in the scene
 * and let it bounce  off a few objects in the scene
 * the accumulated color, the hitted object are stored in the LightPath array 
 *
 * This method looks similar to the Radiance() function.
 * */
std::vector<LightPath> BiDirectionalPT::traceLightRays(const int bounces) const
{
	std::vector<LightPath> lightPath;
	lightPath.reserve(bounces*noLights);

	//for(int i = 0; i < noLights; i++)
	std::map<int,Shape*>::const_iterator it;
	//c++11
	//for(const auto & light : lightEmitters) 
	for(it = lightEmitters.begin(); it != lightEmitters.end(); it++)
	{
		//get current light
		Shape* curLight = it->second;
		Vector origin = curLight->get_position();
		Vector nl = curLight->get_normal(origin);

		//create Random direction vector on a sphere
		Vector dir = uniformSampleSphere() ;
		//Shoot Ray from Lightsource with small eps
		Ray r = Ray(origin + dir * (curLight->radius + 0.1), dir);
		double t;
		int id = 0;
		Vector cf(1,1,1);
		Vector cl(0,0,0);

		
		for(int j = 0; j < bounces; j++)
		{
			if(!Intersect(r, t, id))
			{
				//that didn't hit anything, We should propably ignore it, but I'm not quite sure
				//lightPath[(noLights*i) + j] = LightPath(Vector(0,0,0), Color(0,0,0));
				continue;
			}
			//
			//so we hit something, lets find out what it is
			const Shape &hit_object = *our_scene.get_shapes()[id];	

			//where did we hit it?
			Vector hitpoint = r.org + r.dir *t;	

			//Normal of hitpoint
			Vector n = hit_object.get_normal(hitpoint);
			Vector nl = n;

			//flip it if we hit from within
			if (n.Dot(r.dir) >= 0)
				nl = nl * -1.0;

			//TODO Russian Roulette, it is biased right now
			Color col = hit_object.color;
			cf = cf.MultComponents(col);
			//cl = cl + cf.MultComponents(hit_object.emission);
			cl = cf.MultComponents(curLight->emission);

			//save the hitpoint and the color of the light ray
			//and let it bounce again :)
			Vector dir;
			if (hit_object.refl == DIFF)
			{
				lightPath.push_back(LightPath(hitpoint, cl, it->first, id));
				dir = MaterialBRDF::diffuseBRDF(nl);
			}
			else if (hit_object.refl == SPEC)
			{
				dir = MaterialBRDF::specularBRDF(r, nl);
			}
			else if(hit_object.refl == GLOSSY)
			{
				dir = MaterialBRDF::glossyBRDF(r, nl);
			}
			else if(hit_object.refl == TRANSL)
			{
				dir = MaterialBRDF::translBRDF(r, n, nl, cf);
			}
			else if(hit_object.refl == REFR)
			{
				dir = MaterialBRDF::refrBRDF(r, n, nl, cf);
			}
			else
			{
				//This case should not exist!!!
				dir = Vector(0,0,0);
			}

			r = Ray(hitpoint, dir);
		}
	}

	return lightPath;
}

/* iterative version of the Path Tracing algorithm, the recursive version with more documentation is in PathTracer.cpp */
Color BiDirectionalPT::Radiance(const Ray &ray, const std::vector<LightPath> & lp ) const
{
	int depth = 0;
	double t;	//distance to intersection
	int id = 0; //id of intersected object
	bool sampleLights = true; // if last hit was diffuse, don't sample lights
	Ray r = ray;
	Color cl(0,0,0); //accumulated color
	Color cf(1,1,1); //accumulated reflectance

	for(int i = 0;; i++)
	{
		if (!Intersect(r, t, id))   /* No intersection with scene */
			return cl;

		const Shape &obj = *our_scene.get_shapes()[id];

		Vector hitpoint = r.org + r.dir * t;	/* Intersection point */

		Vector n = obj.get_normal(hitpoint);

		Vector nl = n;
		/* Obtain flipped normal, if object hit from inside */
		if (n.Dot(r.dir) >= 0)
			nl = nl * -1.0;

		Color col = obj.color;

		/* Maximum RGB reflectivity for Russian Roulette */
		double p = col.Max();


		//when we didn't hit a diff surface save the color
		if(sampleLights || obj.refl != DIFF)
			cl = cl + cf.MultComponents(obj.emission);
		if (++depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
		{
			if (drand48() < p)			  /* Russian Roulette */
				col = col * (1 / p);		/* Scale estimator to remain unbiased */
			else
				return cl;  /* No further bounces, only return potential emission */
		}

		cf = cf.MultComponents(col);

		//direction Vector, calculated with the BRDFs
		Vector dir;

		if (obj.refl == DIFF)
		{
			dir = material.diffuseBRDF(nl);
			bool isNotLight = lightEmitters.find(id) == lightEmitters.end();
			if(!isNotLight && i == 0)
				cl = cl + cf.MultComponents(obj.emission);
			else if(isNotLight || sampleLights) 
				cl = cl + cf.MultComponents(explicitComputationOfDirectLight(hitpoint, obj, nl));
			sampleLights = false;
			cl = cl + cf.MultComponents(shootShadowRay(lp, obj, hitpoint, id));
		
		}
		else if (obj.refl == SPEC)
		{
			dir = MaterialBRDF::specularBRDF(r, nl);
			sampleLights = true;
		}
		else if(obj.refl == GLOSSY)
		{
			dir = MaterialBRDF::glossyBRDF(r, nl);
			sampleLights = true;
		}
		else if(obj.refl == TRANSL)
		{
			dir = MaterialBRDF::translBRDF(r, n, nl, cf);
			sampleLights = true;
		}
		else if(obj.refl == REFR)
		{
			dir = MaterialBRDF::refrBRDF(r, n, nl, cf);
			sampleLights = true;
		}
		else
		{
			//This case should not exist!!!
			dir = Vector(0,0,0);
			sampleLights = true;
		}

		r = Ray(hitpoint, dir);
	}
}

/* Direct Light Sample method */
Color BiDirectionalPT::explicitComputationOfDirectLight(const Vector & hitpoint, const Shape & object, const Vector & nl) const
{
	Vector e(0,0,0);	
	//for(int i = 0; i < noLights; i++)
	std::map<int,Shape*>::const_iterator it;
	//c++11
	for(it = lightEmitters.begin(); it != lightEmitters.end(); it++)
	{
		const Shape & curLight = *((it)->second);	
		/* Randomly sample spherical light source from surface intersection */
		/* Set up local orthogonal coordinate system su,sv,sw towards light source */
		Vector curCenter = curLight.get_position();
		Vector sw = (curCenter - hitpoint).Normalized(); 

		Vector su;

		if (fabs(sw.x) > 0.1)
			su = Vector(0.0, 1.0, 0.0);
		else
			su = Vector(1.0, 0.0, 0.0);

		su = (sw.Cross(su)).Normalized();
		Vector sv = sw.Cross(su);

		// next event estimation
		/* Create random sample direction l towards spherical light source */
		double cos_a_max = sqrt(1.0 - curLight.radius * curLight.radius /
									  (hitpoint - curCenter).Dot(hitpoint - curCenter));

		double eps1 = drand48();
		double eps2 = drand48();
		double cos_a = 1.0 - eps1 + eps1 * cos_a_max;
		double sin_a = sqrt(1.0 - cos_a * cos_a);
		double phi = 2.0 * M_PI * eps2;
		Vector l = su * cos(phi) * sin_a +
				   sv * sin(phi) * sin_a +
				   sw * cos_a;
		l = l.Normalized();
		/* Shoot shadow ray, check if intersection is with light source */
		double t;
		int id = 0;
		if (Intersect(Ray(hitpoint,l), t, id) && id == it->first)
		{
			double omega = 2.0 * M_PI * (1.0 - cos_a_max);
			/* Add diffusely reflected light from light source; note constant BRDF 1/PI */
			e = e + curLight.emission * l.Dot(nl) * omega * M_1_PI;
		}
	}

	return e;
}

/* this function returns a directional vector distributed on a sphere*/
Vector BiDirectionalPT::uniformSampleSphere() const
{
	//random angle between 0 and 2pi
	double phi = drand48() * 2 * M_PI;
	//random angle between -1 and 1
	double costheta = 2 * drand48() - 1; 

	double theta = std::acos(costheta);

	return Vector( std::sin(theta) * std::cos(phi),
				   std::sin(theta) * std::sin(phi),
				   std::cos(theta));
}


/* returns the probability that diffuse material reflects in a ceratin direction
 * nl : the normal of the Material
 * outRay : the outgoing direction
 */
double BiDirectionalPT::diffusePdf( const Vector & nl, const Vector & outRay )
{
	return (nl.Dot(outRay) / M_PI);
}



/* check if a point and an object see each other
 * no intersection between other things in the scene
 * p1: Point (used for ray origin)
 * p2: the object which is checked
 * id: the id of the object // bad design shouldn't be needed
 * */
bool BiDirectionalPT::checkVisibility(const Vector & p1, const Vector & p2, int id) const
{

	Ray r(p1, p2);
	int test_id = 0;
	//not needed...
	double d;
	if(Intersect(r, d, test_id))
	{
		return (test_id == id);
	}
	//this should never happen, couse we shoot a ray to an object in the scene
	return false;
}

//TODO pdf for all materials


/* Shoot shadow ray for each hitobject from camera (and bounces of camera) to each Path of the Light Ray
 * lightPath: The path which was calculated by the traceLightRays() function
 * object: the object that was hit by a camera ray
 * hitpoint: The point where the object was hit
 * id: the id where the object is stored in the scene
 * returns Color influence to Ray
 * */
Color BiDirectionalPT::shootShadowRay(const std::vector<LightPath> & lightPath, const Shape & object, const Vector & hitpoint, int id) const
{
	Vector e(0,0,0);

	if(object.refl == DIFF)
	{
		for(int i = 0; i < lightPath.size(); i++)
		{
			//create directional vector between light hitpoint and camera hitpoint
			Vector dir = hitpoint - lightPath[i].hitpoint ;
			Vector dirNorm = dir.Normalized();
			//is there a shadowray?
			if(checkVisibility(lightPath[i].hitpoint, dirNorm, id))
			{
				//get normal of enlightend Object
				const Vector & enlightenObjNormal = our_scene.get_shapes()[lightPath[i].objectid]->get_normal(lightPath[i].hitpoint);

				//the diffuse PDF for the enlightend object
				double diffEnlightenPdf = std::max(0.0, dirNorm.Dot(enlightenObjNormal)) * M_1_PI;
				//the diffuse PDF for the camera object
				double diffEyePdf = std::max(0.0,dirNorm.Dot(object.get_normal(hitpoint))) * M_1_PI;

				double distanceSqr = dir.Dot(dir);
				//First light constant which has to be added
				e = e + lightPath[i].absorbedColor * diffEnlightenPdf * diffEyePdf / distanceSqr;
				//e = e + (lightPath[i].absorbedColor * (dir.Dot(object.get_normal(hitpoint) * diffPdf))) / M_PI;

			}
			else
			{
				//Ignore if they don't see eachother
				continue;	
			}

		}
	}
	return e;
}

