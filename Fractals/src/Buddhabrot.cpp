#include "Buddhabrot.h"
#include "stdlib.h"
#include "time.h"
using namespace std;

void Buddhabrot::gen_fractal()
{
  srand (time(NULL));

  double MinRe = -2.25;
  double MaxRe = 0.75;
  double MinIm = -1.50;
  double MaxIm = 1.50;
  int height = get_height();
  int width = get_width();
  double Re_factor = (MaxRe-MinRe)/(width-1);
  double Im_factor = (MaxIm-MinIm)/(height-1);

  int num_pixels = height*width;
  int rand_mult = ceil(num_pixels/double(RAND_MAX));

  int * mandel_set;
  mandel_set = new int[num_pixels];
  int * outer_array;
  //cout << num_pixels << endl;
  outer_array = new int[num_pixels];
  for (int bucket = 0; bucket < num_pixels; bucket++) {
	  mandel_set[bucket] = 0;
	  outer_array[bucket] = 0;
  }
  int * temp_array;
  temp_array = new int[num_pixels];

  #pragma omp parallel for
  for(int i = 0; i < height; i++) {
      double c_im = MaxIm - i*Im_factor;

      for(int j = 0; j < width; j++) {
		  double c_re = MinRe + j*Re_factor;
		  double Z_re = 0;
		  double Z_im = 0;
		  double old_re = 0;
		  double old_im = 0;
		  bool isInside = true;
		  int old_iter = 1;

		  for(int n = 1; n <= MAXITER; n++) {
              double Z_im2 = Z_im*Z_im;
			  Z_im = 2*Z_re*Z_im + c_im;
			  Z_re = Z_re*Z_re - Z_im2 + c_re;

			  if((Z_re*Z_re + Z_im*Z_im) > 4) {
				  isInside = false;
				  break;
			  }
			  
			  if((Z_re == old_re) && (Z_im == old_im))
				  break;

			  if(n == 2 * old_iter) {
				  old_re = Z_re;
				  old_im = Z_im;
				  old_iter = n;
			  }

		  if (isInside == true)
			  mandel_set[i*height + j] = 1;
		  }
	  }
  }

  #pragma omp parallel for
  for(int i = 0; i < (2*num_pixels); i++) {
	  for (int bucket = 0; bucket < num_pixels; bucket++)
		  temp_array[bucket] = 0;
	  // cout << i << endl;
	  int shift = -rand_mult + 1 + rand() % (2 * rand_mult - 1);
      int num = (rand() * rand_mult + shift) % num_pixels;
	  /*if (mandel_set[num] == 1)
		  continue;*/
      int x = num%height;
      int y = num/height;

      double c_re = MinRe+x*Re_factor;
      double c_im = MaxIm-y*Im_factor;
      double Z_re = 0;
      double Z_im = 0;
      bool isInside = true;

      for(int iter = 0; iter < MAXITER; iter++) {
		  if((Z_re*Z_re + Z_im*Z_im) > 4) {
			  isInside = false;
			  break;
		  }

		  double Z_im2 = Z_im*Z_im;
		  Z_im = 2*Z_re*Z_im+c_im;
		  Z_re= Z_re*Z_re-Z_im2+c_re;

		  int point = height * (MaxIm - Z_im) / Im_factor + (Z_re - MinRe) / Re_factor;
		  int reflect = height * (MaxIm - Z_im) / Im_factor + (width - 1) - (Z_re - MinRe) / Re_factor;
		  if (point >= 0 && point < num_pixels) {
			temp_array[point]++;
			// temp_array[reflect]++;
		  }
	  }

      if(isInside == false) {
		  for (int pos = 0; pos < num_pixels; pos++)
			  outer_array[pos] += temp_array[pos];
	  }
  }

  delete [] temp_array;
  delete [] mandel_set;
  double max = outer_array[0];
  #pragma omp parallel for
  for(int i = 1; i < num_pixels; i++) {
      if(outer_array[i] > max)
		max = outer_array[i];
  }

  #pragma omp parallel for
  for(int k = 0; k < num_pixels; k++) {
      m_bitmap[k*4] = (outer_array[k]/max)*255;
      m_bitmap[k*4 + 1] = (outer_array[k]/max)*255;
      m_bitmap[k*4 + 2] = (outer_array[k]/max)*255;
      m_bitmap[k*4 + 3] = 255;
  }

  delete [] outer_array;
}
