#include "BiDirectionalPT.h"

BiDirectionalPT::BiDirectionalPT(const std::string & scene_name) : PathTracer(scene_name)
{
	getLightEmitters();
}


BiDirectionalPT::BiDirectionalPT(const Scene & scene) : PathTracer(scene)
{
	getLightEmitters();
}

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

Color BiDirectionalPT::calculatePixelColor(const Ray & ray) const
{
	std::vector<LightPath> allLightpaths = traceLightRays(2);
	return Radiance(ray, allLightpaths);
}


/* Shoot a ray from every light in the scene
 * and let it bounce  off a few objects in the scene
 *
 * This method looks similar to the Radiance() function... Maybe I should do something about that
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
		Vector cf(0,0,0);
		Vector cl(1,1,1);

		
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
			cl = cl + cf.MultComponents(hit_object.emission);

			//save the hitpoint and the color of the light ray
			//and let it bounce again :)
			Vector dir;
			//That is the accumulated reflectance needed for the Radiance(). But I have no idea what to do with it here
			if (hit_object.refl == DIFF)
			{
				lightPath.push_back(LightPath(hitpoint, cl, it->first, id));
				dir = diffuseBRDF(nl);
			}
			else if (hit_object.refl == SPEC)
			{
				dir = specularBRDF(r, nl);
			}
			else if(hit_object.refl == GLOSSY)
			{
				dir = glossyBRDF(r, nl);
			}
			else if(hit_object.refl == TRANSL)
			{
				dir = translBRDF(r, n, nl, cf);
			}
			else if(hit_object.refl == REFR)
			{
				dir = refrBRDF(r, n, nl, cf);
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

		//cl = cl + cf.MultComponents(obj.emission);

		if (++depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
		{
			if (drand48() < p)			  /* Russian Roulette */
				col = col * (1 / p);		/* Scale estimator to remain unbiased */
			else
				return cl;  /* No further bounces, only return potential emission */
		}

		cf = cf.MultComponents(col);

		Vector dir;
		if (obj.refl == DIFF)
		{
			dir = diffuseBRDF(nl);
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
			dir = specularBRDF(r, nl);
			sampleLights = true;
		}
		else if(obj.refl == GLOSSY)
		{
			dir = glossyBRDF(r, nl);
			sampleLights = true;
		}
		else if(obj.refl == TRANSL)
		{
			dir = translBRDF(r, n, nl, cf);
			sampleLights = true;
		}
		else if(obj.refl == REFR)
		{
			dir = refrBRDF(r, n, nl, cf);
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
//


/* Shoot shadow ray for each hitobject from camera (and bounces of camera) to each Path of the Light Ray
 * returns Color influence to Ray
 * */
Color BiDirectionalPT::shootShadowRay(const std::vector<LightPath> & lightPath, const Shape & object, const Vector & hitpoint, int id) const
{

	Vector e(0,0,0);

	if(object.refl == DIFF)
	{
		for(int i = 0; i < lightPath.size(); i++)
		{
			//TODO: copy id from object or check try direct comparison
			if(checkVisibility(lightPath[i].hitpoint, hitpoint, id))
			{
				const Vector dir = Vector(hitpoint - lightPath[i].hitpoint).Normalized();
				const double dot = object.get_normal(hitpoint).Dot(dir);
				double diffPdf;
				if(diffPdf <= 0)
					diffPdf = 0.0;
				else
					diffPdf = 1.0 / (2.0 * M_PI * dot);


				//First light constant which has to be added
				e = e + object.color.MultComponents(lightPath[i].absorbedColor * (dir.Dot(object.get_normal(hitpoint) * diffPdf))) / M_PI;
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

Vector BiDirectionalPT::diffuseBRDF(const Vector & nl) const
{
	/* Compute random reflection vector on hemisphere */
	double r1 = 2.0 * M_PI * drand48();
	double r2 = drand48();
	double r2s = sqrt(r2);

	/* Set up local orthogonal coordinate system u,v,w on surface */
	Vector w = nl;
	Vector u;

	if (fabs(w.x) > .1)
		u = Vector(0.0, 1.0, 0.0).Cross(w).Normalized();
	else
		u = (Vector(1.0, 0.0, 0.0).Cross(w)).Normalized();

	Vector v = w.Cross(u);

	/* Random reflection vector d */
	Vector d = (u * cos(r1) * r2s +
				v * sin(r1) * r2s +
				w * sqrt(1 - r2)).Normalized();

	/* Return potential light emission, direct lighting, and indirect lighting (via
	   recursive call for Monte-Carlo integration */
	//return obj.emission * E + e + col.MultComponents(Radiance(Ray(hitpoint, d), depth, 0));
	return d;

}


Vector BiDirectionalPT::glossyBRDF(const Ray & ray, const Vector & nl) const
{
	//cosine distribution??
	const double r1 = 2.0 * M_PI * drand48();
	const double r2 = pow(drand48(), 1.0 / (100.0 + 1.0));
	const double r2s = sqrt(1.0 - r2 * r2);

	//same as above
	Vector sw = (ray.dir - nl * 2.0 * nl.Dot(ray.dir)).Normalized();
	Vector su;
	if (fabs(sw.x) > 0.1)
		su = Vector(0.0, 1.0, 0.0);
	else
		su = Vector(1.0, 0.0, 0.0);

	su = su.Cross(sw).Normalized();

	Vector sv = sw.Cross(su);

	Vector d = (su * cos(r1) * r2s +
				sv * sin(r1) * r2s +
				sw * r2).Normalized();

	double orientation = nl.Dot(d);

	//if ray is beneath surface , reflect the resulting ray around the perfect mirror ray
	if(orientation < 0)
		d = d - sw * 2.0 * sw.Dot(d);

	return d;

}


Vector BiDirectionalPT::translBRDF(const Ray & ray, const Vector & n, const Vector & nl, Color & cf) const
{
	//cosine distribution??
	const double r1 = 2.0 * M_PI * drand48();
	//the exponent is the Transluency Factor, should be changable
	const double r2 = pow(drand48(), 1.0 / (1.0 + 100.0));
	const double r2s = sqrt(1.0 - r2 * r2);

	//same as above
	Vector sw_g = (ray.dir - nl * 2.0 * nl.Dot(ray.dir)).Normalized();
	Vector su_g;
	if (fabs(sw_g.x) > 0.1)
		su_g = Vector(0.0, 1.0, 0.0);
	else
		su_g = Vector(1.0, 0.0, 0.0);

	su_g = su_g.Cross(sw_g).Normalized();

	Vector sv_g = sw_g.Cross(su_g);
	Vector glos_refl_d = (su_g * cos(r1) * r2s +
				sv_g * sin(r1) * r2s +
				sw_g * r2).Normalized();
	
	/* Otherwise object transparent, i.e. assumed dielectric glass material */
	bool into = n.Dot(nl) > 0;		  /* Bool for checking if ray from outside going in */
	double nc = 1;						  /* Index of refraction of air (approximately) */
	double nt = 1.5;					  /* Index of refraction of glass (approximately) */
	double nnt;

	double Re, RP, TP, Tr;
	Vector tdir;
	
	if (into)	   /* Set ratio depending on hit from inside or outside */
		nnt = nc / nt;
	else
		nnt = nt / nc;
	
	double ddn = ray.dir.Dot(nl);
	double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
	
	/* Check for total internal reflection, if so only reflect */
	if (cos2t < 0)
	{
		return glos_refl_d;
	}
	
	/* Otherwise reflection and/or refraction occurs */

	//create Translucent refraction Vector
	const double r1p = 2.0 * M_PI * drand48();
	const double r2p= pow(drand48(), 1.0 / (1.0 + 100.0));
	const double r2sp = sqrt(1.0 - r2 * r2);

	Vector sw_p = (ray.dir);
	Vector su_p;
	if (fabs(sw_p.x) > 0.1)
		su_p = Vector(0.0, 1.0, 0.0);
	else
		su_p = Vector(1.0, 0.0, 0.0);

	su_p = su_p.Cross(sw_p).Normalized();

	Vector sv_p = sw_p.Cross(su_p);
	Vector transl_d = (su_p * cos(r1) * r2sp +
				sv_p * sin(r1p) * r2sp +
				sw_p * r2p).Normalized();


	/* Determine transmitted ray direction for refraction */
	if (into)
		tdir = (transl_d * nnt - n * (ddn * nnt + sqrt(cos2t))).Normalized();
	else
		tdir = (transl_d * nnt + n * (ddn * nnt + sqrt(cos2t))).Normalized();

	/* Determine R0 for Schlick's approximation */
	double a = nt - nc;
	double b = nt + nc;
	double R0 = a * a / (b * b);

	/* Cosine of correct angle depending on outside/inside */
	double c;
	if (into)
		c = 1 + ddn;
	else
		c = 1 - tdir.Dot(n);

	/* Compute Schlick's approximation of Fresnel equation */
	Re = R0 + (1 - R0) * c * c * c * c * c;	 /* Reflectance */
	Tr = 1 - Re;						/* Transmittance */

	/* Probability for selecting reflectance or transmittance */
	double P = .25 + .5 * Re;
	RP = Re / P;			/* Scaling factors for unbiased estimator */
	TP = Tr / (1 - P);

	if (drand48() < P)
	{
		cf = cf.MultComponents(RP);
		return glos_refl_d;
	}
	
	cf = cf.MultComponents(TP);
	return  tdir;

}


Vector BiDirectionalPT::refrBRDF(const Ray & ray, const Vector & n, const Vector & nl, Color & cf) const
{

	/* Object is transparent, i.e. assumed dielectric glass material */
	Vector reflectionV =  ray.dir - n * 2 * n.Dot(ray.dir);	/* Prefect reflection */
	bool into = n.Dot(nl) > 0;		  /* Bool for checking if ray from outside going in */
	double nc = 1;						  /* Index of refraction of air (approximately) */
	double nt = 1.5;					  /* Index of refraction of glass (approximately) */
	double nnt;

	double Re, RP, TP, Tr;
	Vector tdir;

	if (into)	   /* Set ratio depending on hit from inside or outside */
		nnt = nc / nt;
	else
		nnt = nt / nc;

	double ddn = ray.dir.Dot(nl);
	double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

	/* Check for total internal reflection, if so only reflect */
	if (cos2t < 0)
	{
		return reflectionV;
	}
	/* Otherwise reflection and/or refraction occurs */

	/* Determine transmitted ray direction for refraction */
	if (into)
		tdir = (ray.dir * nnt - n * (ddn * nnt + sqrt(cos2t))).Normalized();
	else
		tdir = (ray.dir * nnt + n * (ddn * nnt + sqrt(cos2t))).Normalized();

	/* Determine R0 for Schlick's approximation */
	double a = nt - nc;
	double b = nt + nc;
	double R0 = a * a / (b * b);

	/* Cosine of correct angle depending on outside/inside */
	double c;
	if (into)
		c = 1 + ddn;
	else
		c = 1 - tdir.Dot(n);

	/* Compute Schlick's approximation of Fresnel equation */
	Re = R0 + (1 - R0) * c * c * c * c * c;	 /* Reflectance */
	Tr = 1 - Re;						/* Transmittance */

	/* Probability for selecting reflectance or transmittance */
	double P = .25 + .5 * Re;
	RP = Re / P;			/* Scaling factors for unbiased estimator */
	TP = Tr / (1 - P);

	if (drand48() < P)
	{
		cf = cf.MultComponents(RP);
		return reflectionV;
	}
	
	cf = cf.MultComponents(TP);
	return tdir;

}

Vector BiDirectionalPT::specularBRDF(const Ray & ray, const Vector & nl) const
{
	return (ray.dir - nl * 2 * nl.Dot(ray.dir));
}

