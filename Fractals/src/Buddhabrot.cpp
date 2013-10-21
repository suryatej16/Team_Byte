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
  int * outer_array;
  //cout << num_pixels << endl;
  outer_array = new int[num_pixels];
  for (int bucket = 0; bucket < num_pixels; bucket++)
	  outer_array[bucket] = 0;
  int * temp_array;
  temp_array = new int[num_pixels];
  #pragma omp parallel for
  for(int i = 0; i < (10*num_pixels); i++)
    {
      for (int bucket = 0; bucket < num_pixels; bucket++)
	{
	  temp_array[bucket] = 0;
	}
      // int temp_bucket[num_pixels];
      int num = rand() % num_pixels;
      //cout << num << endl;
      int x = num%height;
      int y = num/height;
      double c_re = MinRe+x*Re_factor;
      double c_im = MaxIm-y*Im_factor;
      double Z_re = 0;
      double Z_im = 0;
      bool isInside = true;
      for(int iter = 0; iter < MAXITER; iter++)
	{
	  double Z_im2 = Z_im*Z_im;
	  double Z_re2 = Z_re*Z_re;
	  if((Z_re2 + Z_im2) > 4)
	    {
	      isInside = false;
	      break;
	    }
	  Z_im = 2*Z_re*Z_im+c_im;
	  Z_re= Z_re2-Z_im2+c_re;
	  temp_array[(a+b)]++;
	  int array_pos = height * (MaxIm - Z_im) / Im_factor + (Z_re - MinRe) / Re_factor;
	  if (array_pos >= 0 && array_pos < num_pixels)
		temp_array[array_pos]++;
	}
      if(isInside == false)
	{
	  for (int pos = 0; pos < num_pixels; pos++)
	    {
	      outer_array[pos] += temp_array[pos];
	    }
	}
    }
  delete [] temp_array;
  double max;
  #pragma omp parallel for
  for(int i = 0; i < num_pixels; i++)
    {
      if(i == 0)
	{
	  max = outer_array[0];
	  continue;
	}
      if(outer_array[i] > max)
	{
	  max = outer_array[i];
	}
    }
  #pragma omp parallel for
  for(int i = 0; i < num_pixels; i++)
    {
      m_bitmap[k*4] = (outer_array[k]/max)*255;
      m_bitmap[k*4 + 1] = (outer_array[k]/max)*255;
      m_bitmap[k*4 + 2] = (outer_array[k]/max)*255;
      m_bitmap[k*4 + 3] = 255;
    }
  delete [] outer_array;
}
   
	// Real (-2.5, 1)
	// Imaginary (-1, 1)
    
    // Initialize a bucket array (one integer for each pixel) (this is the outer bucket array)

    // iterate over the following several thousand times (at least more times than # of pixels)
  
        // Create a temporary bucket array (one integer for each pixel
        //
        // Let C be a random point in the complex plane
        //
        // Trace the orbit of C, incrementing the temporary bucket that z falls in for each iteration
        // If Z is in the mandelbrot set, discard the temporary bucket
        // Else, merge the temporary bucket with the outer bucket array
        

     // Normalize the global bucket array by dividing each value by the maximum value
     // Color each pixel however you wish
     //
     // Parallelizing this function is tricky. It helps to have a list of temporary bucket arrays
     // Which are merged after the computation has finished.
     
     // Parallelizing is not required, but will save you a lot of time.
     


    

