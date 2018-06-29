class complex
{
  public:
  double x,y;
  __host__ __device__ complex(double a,double b);
  __host__ __device__ complex();
};

__host__ __device__ complex::complex()
{
  x=0;
  y=0;
}

__host__ __device__ complex::complex(double a,double b)
{
  x=a;
  y=b;
}

__host__ __device__ complex make_complex(double a,double b)
{
  complex result;
  result.x = a;
  result.y = b;
  return result;
}

__host__ __device__ double real(const complex a)
{
  return a.x;
}

__host__ __device__ double imag(const complex a)
{
  return a.y;
}

__host__ __device__ complex conj(const complex a)
{
  return complex(a.x,
                      -a.y);
}

__host__ __device__ complex
operator + (const complex a,const complex b)
{
  return complex(a.x+b.x,
                      a.y+b.y);
}

__host__ __device__ complex
operator - (const double a,const complex b)
{
  return complex(a+b.x,
                      b.y);
}


__host__ __device__ complex
operator * (const complex a,const complex b)
{
  return complex(a.x*b.x-a.y*b.y,
                      a.x*b.y+a.y*b.y);
}


__host__ __device__ complex
operator * (const double a,const complex b)
{
  return complex(a*b.x,
                      a*b.y);
}

__host__ __device__ complex
operator * (const complex a,const double b)
{
  return complex(a.x*b,
                      a.y*b);
}

__host__ __device__ complex
operator / (const complex a,const double b)
{
  return complex(a.x/b,
                      a.y/b);
}

__host__ __device__ complex
operator / (const complex a,const complex b)
{
  double mod = b.x*b.x+b.y*b.y;
  return a*conj(b)/mod;
}

__host__ __device__ complex
operator + (const double a,const complex b)
{
  return complex(b.x+a,
                      b.y);
}

__host__ __device__ complex
operator / (const double a,const complex b)
{
  double mod = b.x*b.x+b.y*b.y;
  return a*conj(b)/mod;
}
