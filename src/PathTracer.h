#pragma once


#include "shape.h"
#include "scene.h"

class PathTracer
{
	
	public:
		//functions
		PathTracer(const std::string & scene_name, bool sampleLight_);
		PathTracer(const Scene & scene, bool sampleLight_);
		virtual Color calculatePixelColor(const Ray & ray) const;
	
	private:
		//functions:
		Color Radiance(const Ray &ray, int depth, int E) const;
		Color refrBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth ) const;
		Color translBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth) const;
		Color glossyBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth, int E) const;
		Color diffuseBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth, int E) const;
		Color specularBRDF(const Shape &obj, const Ray & ray, const Vector & hitpoint, const Vector & nl, int depth) const;
		//
		bool sampleLight;

	protected:
		//variables:
		Scene our_scene;
		unsigned int numObjects;
		//functions:
		bool Intersect(const Ray &ray, double &t, int &id) const;

};
