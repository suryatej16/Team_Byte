#include "Mandelbrot.h"
#include "math.h"

void Mandelbrot::gen_fractal()
{
  double MinRe = -2.25;
  double MaxRe = 0.75;
  double MinIm = -1.50;
  double MaxIm = 1.50;
  int height = get_height();
  int width = get_width();
  double Re_factor = (MaxRe-MinRe)/(width-1);
  double Im_factor = (MaxIm-MinIm)/(height-1);
  #pragma omp parallel for
  for(int i = 0; i < height; i++)
    {
      double c_im = MaxIm - i*Im_factor;
      for(int j = 0; j < width; j++)
	{
	  double c_re = MinRe + j*Re_factor;
	  double Z_re = 0;
	  double Z_im = 0;
	  bool isInside = true;
	  int iter = 0;
	  for(int n = 1; n < MAXITER; n++)
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
	  double mu = iter - log(log(sqrt(pow(Z_re, 2) + pow(Z_im, 2)))) / log(2.0);
	  if(isInside == false)
	    {
	      m_bitmap[i*height*4 + j*4] = pow(mu/MAXITER,.60)*255;
	      m_bitmap[i*height*4 + j*4 + 1] = pow(mu/MAXITER,.40)*255;
	      m_bitmap[i*height*4 + j*4 + 2] = pow(mu/MAXITER,.30)*255;
	      m_bitmap[i*height*4 + j*4 + 3] = 255;
	    }
	  else
	    {
	      m_bitmap[i*height*4 + j*4] = 0;
	      m_bitmap[i*height*4 + j*4 + 1] = 0;
	      m_bitmap[i*height*4 + j*4 + 2] = 0;
	      m_bitmap[i*height*4 + j*4 + 3] = 255;
	    }
	}
    }
  int num_pixels = height*width;
  #pragma omp parallel for
  for(int i = 0; i < num_pixels; i++)
    {
      int x = i%height;
      int y = i/height;
      if((x-1) >= 0 && (y-1) >=0 && (x+1) < width && (y+1) < height)
	{
	  double temp1 = ( m_bitmap[(x+1)*height*4 + y*4]+ m_bitmap[(x-1)*height*4 + y*4] +  m_bitmap[x*height*4 + (y-1)*4] +  m_bitmap[x*height*4 + (y+1)*4])/4;
	  m_bitmap[x*height*4 + y*4] = temp1;

	  double temp2 = ( m_bitmap[(x+1)*height*4 + y*4 + 1]+ m_bitmap[(x-1)*height*4 + y*4 + 1] +  m_bitmap[x*height*4 + (y-1)*4 + 1] +  m_bitmap[x*height*4 + (y+1)*4 + 1])/4;
	  m_bitmap[x*height*4 + y*4 + 1] = temp2;
	  
	  double temp3 = ( m_bitmap[(x+1)*height*4 + y*4 + 2]+ m_bitmap[(x-1)*height*4 + y*4 + 2] +  m_bitmap[x*height*4 + (y-1)*4 + 2] +  m_bitmap[x*height*4 + (y+1)*4 + 2])/4;
	  m_bitmap[x*height*4 + y*4 + 2] = temp3;
	}
    }
}  