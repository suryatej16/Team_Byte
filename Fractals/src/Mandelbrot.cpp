#include "Mandelbrot.h"
#include "math.h"

void Mandelbrot::gen_fractal()
{
  double MinRe = -2.5;
  double MaxRe = 1.0;
  double MinIm = -1.0;
  int height = get_height();
  int width = get_width();
  double MaxIm = MinIm+(MaxRe-MinRe)*height/width;
  double Re_factor = (MaxRe-MinRe)/(width-1);
  double Im_factor = (MaxIm-MinIm)/(height-1);
  #pragma omp parallel for
  for(int i = 0; i < height; i++)
    {
      double c_im = MaxIm - i*Im_factor;
      for(int j = 0; j < width; j++)
	{
	  double c_re = MinRe + j*Re_factor;
	  double Z_re = c_re;
	  double Z_im = c_im;
	  bool isInside = true;
	  int iter = 0;
	  for(int n = 0; n < MAXITER; n++)
	    {
	      if((Z_re*Z_re + Z_im*Z_im) > 4)
		{
		  isInside = false;
		  break;
		}
	      iter++;
	      double Z_im2 = Z_im*Z_im;
	      Z_im = 2*Z_re*Z_im + c_im;
	      Z_re = Z_re*Z_re - Z_im2 + c_re;
	    }
	  if(isInside == false)
	    {
	      m_bitmap[i*height*4 + j*4] = pow(double(iter/MAXITER),0.6)*255;
	      m_bitmap[i*height*4 + j*4 + 1] = pow(double(iter/MAXITER),0.5)*255;
	      m_bitmap[i*height*4 + j*4 + 2] = pow(double(iter/MAXITER), 0.4)*255;
	      m_bitmap[i*height*4 + j*4 + 3] = 255;
	    }
	}
    }
}
  
	// Real (-2.5, 1)
	// Imaginary (-1, 1)
    
    // Iterate over each pixel and determine the corresponding complex coordinate (or several complex coordinates if you're anti-aliasing)
    

        // For each pixel, compute the corresponding complex coordinate C (or multiple sub-coordinates, for AA)   

        // let z_r and z_i be 0
        // Begin iterating... while z is not infinity and not too many iterations have passed...
        
            // Z = Z^2 + C
            // increment an iteration counter

        // Color each pixel...


