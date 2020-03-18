float absf (float x) {
  if (x < 0) return (-1 * x);
  return x;
}
double absd (double x) {
  if (x < 0) return (-1 * x);
  return x;
}

float mn_squaref(float x) {
  return x * x;
}
inline double mn_squared(double x) {
  return x * x;
}
