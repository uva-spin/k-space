// Thin C ABI around the official LHAPDF interpolation, for ctypes clients.
#include <LHAPDF/LHAPDF.h>
#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
extern "C" {
void* kt_pdf_open(const char* data, const char* name, int member, char* err, int cap) {
  try { LHAPDF::setVerbosity(0); LHAPDF::setPaths(std::string(data));
        return LHAPDF::mkPDF(name,member); }
  catch (const std::exception& e) { std::snprintf(err,cap,"%s",e.what()); return nullptr; }
}
void kt_pdf_close(void* p) { delete static_cast<LHAPDF::PDF*>(p); }
int kt_pdf_eval(void* p, int pid, int n, const double* x, const double* q,
                double* out, char* err, int cap) {
  try { auto* pdf=static_cast<LHAPDF::PDF*>(p);
    if (!pdf || !pdf->hasFlavor(pid)) throw std::runtime_error("Unavailable PDF flavor");
    for (int i=0;i<n;++i) {
      if (!std::isfinite(x[i]) || !std::isfinite(q[i]) || x[i]<=0 || q[i]<=0 || !pdf->inRangeXQ(x[i],q[i]))
        throw std::runtime_error("PDF evaluation outside the tabulated domain");
      out[i]=pdf->xfxQ(pid,x[i],q[i])/x[i];
      if (!std::isfinite(out[i])) throw std::runtime_error("Nonfinite PDF value");
    } return 0;
  } catch (const std::exception& e) { std::snprintf(err,cap,"%s",e.what()); return 1; }
}
}
