#ifndef UTILS_H
#define UTILS_H

inline float absf (float x) {
  if (x < 0) return (-1 * x);
  return x;
}
inline double absd (double x){
  if (x < 0) return (-1 * x);
  return x;
}

inline float mn_squaref(float x){
  return x * x;
}
inline double mn_squared(double x){
  return x * x;
}

#endif
