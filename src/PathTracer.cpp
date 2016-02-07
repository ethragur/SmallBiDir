#include "PathTracer.h"
#include "SceneReader.h"

PathTracer::PathTracer(const std::string & scene_name, bool sampleLight_)
{
	our_scene.load_scene(std::string("../data/" + scene_name + ".obj"), std::string("../data/"));
//	Sphere *sp0 = new Sphere(16.5,	Vector(27, 16.5, -10), Vector(), Vector(1.0, 1.0, 1.0) * .999, REFR);
//	Sphere *sp1 = new Sphere(16.5,  Vector(-20, 16.5, 10), Vector(), Vector(1.0, 1.0, 1.0) * .999, REFR);
	//Sphere *sp2 = new Sphere(8.5,  Vector(-20, 8.5, 25), Vector(), Vector(1.0, 1.0, 0.0) * .999, DIFF);
//	Sphere *sp3 = new Sphere(2.0, Vector(20, 50, -15), Vector(80,80,80), Vector(1.0,1.0,1.0) * .999, DIFF);
//	our_scene.add_sphere(sp0);
//	our_scene.add_sphere(sp1);
//	our_scene.add_sphere(sp2);
//	our_scene.add_sphere(sp3);
	SceneReader s;
	s.ReadScene(our_scene, std::string("../data/" + scene_name + ".sce"));
	
	numObjects = our_scene.get_shapes().size();
	sampleLight = sampleLight_;
}


PathTracer::PathTracer(const Scene & scene, bool sampleLight_)
{
	our_scene = scene;
	numObjects = scene.get_shapes().size();
	sampleLight = sampleLight_;
	SceneReader s;
	s.ReadScene(our_scene, "../data/scene1.sce");
}


Color PathTracer::calculatePixelColor(const Ray & ray) const
{
	return Radiance(ray, 0, 1);
}


/******************************************************************
* Recursive path tracing for computing radiance via Monte-Carlo
* integration; only considers perfectly diffuse, specular or
* transparent materials;
* after 5 bounces Russian Roulette is used to possibly terminate rays;
* emitted light from light source only included on first direct hit
* (possibly via specular reflection, refraction), controlled by
* parameter E = 0/1;
* on diffuse surfaces light sources are explicitely sampled;
* for transparent objects, Schlick's approximation is employed;
* for first 3 bounces obtain reflected and refracted component,
* afterwards one of the two is chosen randomly
*******************************************************************/
Color PathTracer::Radiance(const Ray &ray, int depth, int E) const
{
	depth++;

	double t;
	int id = 0;

	if (!Intersect(ray, t, id))   /* No intersection with scene */
		return Color(0,0,0);

	const Shape &obj = *our_scene.get_shapes()[id];

	Vector hitpoint = ray.org + ray.dir * t;	/* Intersection point */

	Vector nl = obj.get_normal(hitpoint);

	/* Obtain flipped normal, if object hit from inside */
	if (nl.Dot(ray.dir) >= 0)
		nl = nl * -1.0;

	Color col = obj.color;

	/* Maximum RGB reflectivity for Russian Roulette */
	double p = col.Max();

	if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
	{
		if (drand48() < p)			  /* Russian Roulette */
			col = col * (1 / p);		/* Scale estimator to remain unbiased */
		else
			return obj.emission * E;  /* No further bounces, only return potential emission */
	}

	if (obj.refl == DIFF)
	{
		return diffuseBRDF(obj, ray, hitpoint, nl, depth, E);
	}
	else if (obj.refl == SPEC)
	{
		return specularBRDF(obj, ray, hitpoint, nl, depth);
	}
	else if(obj.refl == GLOSSY)
	{
		return glossyBRDF(obj, ray, hitpoint, nl, depth, E);
	}
	else if(obj.refl == TRANSL)
	{
		return translBRDF(obj, ray, hitpoint, nl, depth);
	}
	else if(obj.refl == REFR)
	{
		return refrBRDF(obj, ray, hitpoint, nl, depth);
	}
	
	return Color(0,0,0);

}



/******************************************************************
* Check for closest intersection of a ray with the scene;
* returns true if intersection is found, as well as ray parameter
* of intersection and id of intersected object
*******************************************************************/
bool PathTracer::Intersect(const Ray &ray, double &t, int &id) const
{
	t = 1e20;

	for (int i = 0; i < numObjects; i++)
	{
		double d = our_scene.get_shapes()[i]->intersect(ray);
		if (d > 0.0 && d < t)
		{
			t = d;
			id = i;
		}
	}
	return t < 1e20;
}

Color PathTracer::diffuseBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth, int E) const
{
	/* Compute random reflection vector on hemisphere */
	double r1 = 2.0 * M_PI * drand48();
	double r2 = drand48();
	double r2s = sqrt(r2);

	Color col = obj.color;
	Vector e(0,0,0);

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

	if(sampleLight)
	{
		/* Explicit computation of direct lighting */
		for (int i = 0; i < numObjects; i++)
		{
			const Shape &cur_obj = *our_scene.get_shapes()[i];
			if (cur_obj.emission.x <= 0 && cur_obj.emission.y <= 0 && cur_obj.emission.z <= 0)
				continue; /* Skip objects that are not light sources */

			/* Randomly sample spherical light source from surface intersection */

			/* Set up local orthogonal coordinate system su,sv,sw towards light source */
			Vector cur_cent = cur_obj.get_position();
			Vector sw = (cur_cent - hitpoint).Normalized(); 

			Vector su;

			if (fabs(sw.x) > 0.1)
				su = Vector(0.0, 1.0, 0.0);
			else
				su = Vector(1.0, 0.0, 0.0);

			su = (sw.Cross(su)).Normalized();
			Vector sv = sw.Cross(su);

			// next event estimation
			/* Create random sample direction l towards spherical light source */
			double cos_a_max = sqrt(1.0 - cur_obj.radius * cur_obj.radius /
										  (hitpoint - cur_cent).Dot(hitpoint - cur_cent));

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
			if (Intersect(Ray(hitpoint,l), t, id) && id==i)
			{
				double omega = 2 * M_PI * (1 - cos_a_max);

				/* Add diffusely reflected light from light source; note constant BRDF 1/PI */
				e = e + cur_obj.emission * l.Dot(nl) * omega * M_1_PI;
			}
		}
		return obj.emission * E + obj.color.MultComponents(e) + col.MultComponents(Radiance(Ray(hitpoint, d), depth, 0));
	}
	else
		return obj.emission + col.MultComponents(Radiance(Ray(hitpoint, d), depth, 0));

	
	/* Return potential light emission, direct lighting, and indirect lighting (via
	   recursive call for Monte-Carlo integration */

}

Color PathTracer::glossyBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth, int E) const
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

	Color col = obj.color;

	return obj.emission + col.MultComponents(Radiance(Ray(hitpoint, d), depth, E));

}

Color PathTracer::translBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth) const
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
	
	Vector obj_normal = obj.get_normal(hitpoint);
	
	/* Otherwise object transparent, i.e. assumed dielectric glass material */
	Ray reflRay = Ray(hitpoint, glos_refl_d);
	bool into = obj_normal.Dot(nl) > 0;		  /* Bool for checking if ray from outside going in */
	double nc = 1;						  /* Index of refraction of air (approximately) */
	double nt = 1.5;					  /* Index of refraction of glass (approximately) */
	double nnt;

	Color col = obj.color;
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
		return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1));
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
		tdir = (transl_d * nnt - obj_normal * (ddn * nnt + sqrt(cos2t))).Normalized();
	else
		tdir = (transl_d * nnt + obj_normal * (ddn * nnt + sqrt(cos2t))).Normalized();

	/* Determine R0 for Schlick's approximation */
	double a = nt - nc;
	double b = nt + nc;
	double R0 = a * a / (b * b);

	/* Cosine of correct angle depending on outside/inside */
	double c;
	if (into)
		c = 1 + ddn;
	else
		c = 1 - tdir.Dot(obj_normal);

	/* Compute Schlick's approximation of Fresnel equation */
	Re = R0 + (1 - R0) * c * c * c * c * c;	 /* Reflectance */
	Tr = 1 - Re;						/* Transmittance */

	/* Probability for selecting reflectance or transmittance */
	double P = .25 + .5 * Re;
	RP = Re / P;			/* Scaling factors for unbiased estimator */
	TP = Tr / (1 - P);

	if (depth < 3)
	{
		/* Initially both reflection and trasmission */
		return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * Re +
												 Radiance(Ray(hitpoint, tdir), depth, 1) * Tr);
	}
	else if (drand48() < P)
	{
		return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * RP);
	}
	
	return obj.emission + col.MultComponents(Radiance(Ray(hitpoint, tdir), depth, 1) * TP);

}


Color PathTracer::refrBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth) const
{
	Vector obj_normal = obj.get_normal(hitpoint);

	/* Otherwise object transparent, i.e. assumed dielectric glass material */
	Ray reflRay = Ray(hitpoint, ray.dir - obj_normal * 2 * obj_normal.Dot(ray.dir));	/* Prefect reflection */
	bool into = obj_normal.Dot(nl) > 0;		  /* Bool for checking if ray from outside going in */
	double nc = 1;						  /* Index of refraction of air (approximately) */
	double nt = 1.5;					  /* Index of refraction of glass (approximately) */
	double nnt;

	Color col = obj.color;
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
		return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1));
	}
	/* Otherwise reflection and/or refraction occurs */

	/* Determine transmitted ray direction for refraction */
	if (into)
		tdir = (ray.dir * nnt - obj_normal * (ddn * nnt + sqrt(cos2t))).Normalized();
	else
		tdir = (ray.dir * nnt + obj_normal * (ddn * nnt + sqrt(cos2t))).Normalized();

	/* Determine R0 for Schlick�s approximation */
	double a = nt - nc;
	double b = nt + nc;
	double R0 = a * a / (b * b);

	/* Cosine of correct angle depending on outside/inside */
	double c;
	if (into)
		c = 1 + ddn;
	else
		c = 1 - tdir.Dot(obj_normal);

	/* Compute Schlick�s approximation of Fresnel equation */
	Re = R0 + (1 - R0) * c * c * c * c * c;	 /* Reflectance */
	Tr = 1 - Re;						/* Transmittance */

	/* Probability for selecting reflectance or transmittance */
	double P = .25 + .5 * Re;
	RP = Re / P;			/* Scaling factors for unbiased estimator */
	TP = Tr / (1 - P);

	if (depth < 3)
	{
		/* Initially both reflection and trasmission */
		return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * Re +
												 Radiance(Ray(hitpoint, tdir), depth, 1) * Tr);
	}
	else if (drand48() < P)
	{
		return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * RP);
	}
	
	return obj.emission + col.MultComponents(Radiance(Ray(hitpoint, tdir), depth, 1) * TP);

}

Color PathTracer::specularBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth) const
{
	Color col = obj.color;
	/* Return light emission mirror reflection (via recursive call using perfect reflection vector) */
	return obj.emission + col.MultComponents(Radiance(Ray(hitpoint, ray.dir - nl * 2 * nl.Dot(ray.dir)), depth, 1));
}

