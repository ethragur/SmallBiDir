#pragma once

/*------------------------------------------------------------------
| Struct holds pixels/colors of rendered image
------------------------------------------------------------------*/
class Image
{
public:
	int width, height;
	Color *pixels;

	Image(int _w, int _h) : width(_w), height(_h)
	{
		pixels = new Color[width * height];
	}

	Color getColor(int x, int y)
	{
		int image_index = (height - y - 1) * width + x;
		return pixels[image_index];
	}

	void setColor(int x, int y, const Color &c)
	{
		int image_index = (height - y - 1) * width + x;
		pixels[image_index] = c;
	}

	void addColor(int x, int y, const Color &c)
	{
		int image_index = (height - y - 1) * width + x;
		pixels[image_index] = pixels[image_index] + c;
	}

	int toInteger(double x)
	{
		/* Clamp to [0,1] */
		if (x < 0.0)
			x = 0.0;

		if (x > 1.0)
			x = 1.0;

		/* Apply gamma correction and convert to integer */
		return int(pow(x, 1 / 2.2) * 255 + .5);
	}

	void Save(const std::string &filename)
	{
		/* Save image in PPM format */
		FILE *f = fopen(filename.c_str(), "wb");
		fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
		for (int i = 0; i < width * height; i++)
			fprintf(f, "%d %d %d ", toInteger(pixels[i].x),
					toInteger(pixels[i].y),
					toInteger(pixels[i].z));
		fclose(f);
	}
};
