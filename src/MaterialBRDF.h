class MaterialBRDF 
{
	public:
		Vector diffuseBRDF(const Vector & nl) const
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

		Vector glossyBRDF(const Ray & ray, const Vector & nl) const
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
		
		Vector translBRDF(const Ray & ray, const Vector & n, const Vector & nl, Color & cf) const
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
				cf = cf * (RP);
				return glos_refl_d;
			}
			
			cf = cf * (TP);
			return  tdir;
		}

		Vector refrBRDF(const Ray & ray, const Vector & n, const Vector & nl, Color & cf) const
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
				cf = cf * (RP);
				return reflectionV;
			}
			
			cf = cf * (TP);
			return tdir;
		}

		Vector specularBRDF(const Ray & ray, const Vector & nl) const
		{
			return (ray.dir - nl * 2 * nl.Dot(ray.dir));
		}

};
