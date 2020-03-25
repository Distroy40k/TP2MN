
typedef struct {
  float real ;
  float imaginary ;
} complexe_float_t ;

typedef struct {
  double real ;
  double imaginary ;
} complexe_double_t ;


inline complexe_float_t init_complexe_float (const float real, const float img) {
  complexe_float_t r;
  r.real = real;
  r.imaginary = img;
  return r;
}

inline complexe_double_t init_complexe_double (const double real, const double img) {
  complexe_double_t r;
  r.real = real;
  r.imaginary = img;
  return r;
}

inline complexe_float_t add_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;

  return r ;
}

inline complexe_double_t add_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;

  return r ;
}

inline complexe_float_t mult_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = c1.real * c2.real - c1.imaginary * c2.imaginary ;
  r.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real  ;

  return r ;
}

inline complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;

  r.real = c1.real * c2.real - c1.imaginary * c2.imaginary ;
  r.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real  ;

  return r ;
}

inline complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
 complexe_float_t r ;
 float denominateur = c2.real * c2.real + c2.imaginary * c2.imaginary;

  r.real = (c1.real * c2.real + c1.imaginary * c2.imaginary) / denominateur ;
  r.imaginary = (c1.imaginary * c2.real - c1.real * c2.imaginary) / denominateur  ;

  return r ;
}

inline complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;
 double denominateur = c2.real * c2.real + c2.imaginary * c2.imaginary;

  r.real = (c1.real * c2.real + c1.imaginary * c2.imaginary) / denominateur ;
  r.imaginary = (c1.imaginary * c2.real - c1.real * c2.imaginary) / denominateur  ;

  return r ;
}
