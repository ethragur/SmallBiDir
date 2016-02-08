#pragma once


#include "PathTracer.h"
#include "MaterialBRDF.h"
#include <map>

struct LightPath
{
	Vector hitpoint;
	Color absorbedColor;
	int lightid;
	int objectid;

	LightPath(Vector hitpoint_, Color absorbedColor_, int lightid_, int objectid_) :
		hitpoint(hitpoint_), absorbedColor(absorbedColor_), lightid(lightid_), objectid(objectid_)
	{}

	LightPath()
	{}
};

class BiDirectionalPT : public PathTracer 
{
	
	public:
		//functions
		BiDirectionalPT(const std::string & scene_name, unsigned int lightbounces);
		BiDirectionalPT(const Scene & scene, unsigned int lightbounces);
		BiDirectionalPT(const std::string & scene_name);
		BiDirectionalPT(const Scene & scene);
		Color calculatePixelColor(const Ray & ray) const;
	
	private:
		//variables:
		//std::vector<Shape*> lightEmitters;
		std::map<int,Shape*> lightEmitters;
		unsigned int noLights;
		unsigned int lightBounces;
		//functions:
		void getLightEmitters();
		std::vector<LightPath> traceLightRays(const int bounces) const;
		Color Radiance(const Ray &ray, const std::vector<LightPath> & lp) const;
		Vector uniformSampleSphere() const; 
		double diffusePdf( const Vector & inpDir, const Vector & outRay );
		/*
		Vector diffuseBRDF(const Vector & nl) const;
		Vector glossyBRDF(const Ray & ray, const Vector & nl) const;
		Vector translBRDF(const Ray & ray, const Vector & n, const Vector & nl, Color & cf) const;
		Vector refrBRDF(const Ray & ray, const Vector & n, const Vector & nl, Color & cf) const;
		Vector specularBRDF(const Ray & ray, const Vector & nl) const;
		*/
		bool checkVisibility(const Vector & p1, const Vector & p2, int id) const;
		Color shootShadowRay(const std::vector<LightPath> & lightPath, const Shape & object, const Vector & hitpoint, int id) const;
		Color explicitComputationOfDirectLight(const Vector & hitpoint, const Shape & object, const Vector & nl) const;
};
