"""Explicit local LHAPDF C-ABI adapter. No fallback values or extrapolation."""
import ctypes as c
import numpy as np

class LHAPDFProvider:
    def __init__(self,library,data_path,name='NNPDF40_nnlo_as_01180',member=0):
        self.lib=c.CDLL(str(library));self.handle=None
        ptr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')
        self.lib.kt_pdf_open.argtypes=[c.c_char_p,c.c_char_p,c.c_int,c.c_char_p,c.c_int]
        self.lib.kt_pdf_open.restype=c.c_void_p
        self.lib.kt_pdf_close.argtypes=[c.c_void_p];self.lib.kt_pdf_close.restype=None
        self.lib.kt_pdf_eval.argtypes=[c.c_void_p,c.c_int,c.c_int,ptr,ptr,ptr,c.c_char_p,c.c_int]
        self.lib.kt_pdf_eval.restype=c.c_int
        err=c.create_string_buffer(4096)
        self.handle=self.lib.kt_pdf_open(str(data_path).encode(),name.encode(),member,err,len(err))
        if not self.handle: raise RuntimeError(err.value.decode())
    def f(self,pid,x,mu):
        if not self.handle: raise RuntimeError('PDF provider is closed')
        xx,qq=np.broadcast_arrays(np.asarray(x,float),np.asarray(mu,float))
        xa=np.ascontiguousarray(xx.ravel());qa=np.ascontiguousarray(qq.ravel());out=np.empty_like(xa)
        err=c.create_string_buffer(4096)
        status=self.lib.kt_pdf_eval(self.handle,int(pid),xa.size,xa,qa,out,err,len(err))
        if status: raise ValueError(err.value.decode())
        result=out.reshape(xx.shape)
        return float(result) if result.ndim==0 else result
    def alpha(self,mu):
        if not self.handle: raise RuntimeError('PDF provider is closed')
        if not hasattr(self.lib,'kt_pdf_alpha'): raise RuntimeError('Bridge lacks explicit alpha_s anchor access')
        self.lib.kt_pdf_alpha.argtypes=[c.c_void_p,c.c_double,c.POINTER(c.c_double),c.c_char_p,c.c_int]
        self.lib.kt_pdf_alpha.restype=c.c_int
        out=c.c_double();err=c.create_string_buffer(4096)
        status=self.lib.kt_pdf_alpha(self.handle,float(mu),c.byref(out),err,len(err))
        if status: raise ValueError(err.value.decode())
        return out.value
    def close(self):
        if self.handle: self.lib.kt_pdf_close(self.handle);self.handle=None
    def __del__(self):
        if getattr(self,'handle',None): self.close()
