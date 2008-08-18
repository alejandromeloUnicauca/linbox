/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.27
 * 
 * This file is not intended to be easily readable and contains a number of 
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG 
 * interface file instead. 
 * ----------------------------------------------------------------------------- */


#ifdef __cplusplus
template<class T> class SwigValueWrapper {
    T *tt;
public:
    SwigValueWrapper() : tt(0) { }
    SwigValueWrapper(const SwigValueWrapper<T>& rhs) : tt(new T(*rhs.tt)) { }
    SwigValueWrapper(const T& t) : tt(new T(t)) { }
    ~SwigValueWrapper() { delete tt; } 
    SwigValueWrapper& operator=(const T& t) { delete tt; tt = new T(t); return *this; }
    operator T&() const { return *tt; }
    T *operator&() { return tt; }
private:
    SwigValueWrapper& operator=(const SwigValueWrapper<T>& rhs);
};
#endif

/***********************************************************************
 *
 *  This section contains generic SWIG labels for method/variable
 *  declarations/attributes, and other compiler dependent labels.
 *
 ************************************************************************/

/* template workaround for compilers that cannot correctly implement the C++ standard */
#ifndef SWIGTEMPLATEDISAMBIGUATOR
#  if defined(__SUNPRO_CC) && (__SUNPRO_CC <= 0x560)
#    define SWIGTEMPLATEDISAMBIGUATOR template
#  else
#    define SWIGTEMPLATEDISAMBIGUATOR 
#  endif
#endif

/* inline attribute */
#ifndef SWIGINLINE
# if defined(__cplusplus) || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#   define SWIGINLINE inline
# else
#   define SWIGINLINE
# endif
#endif

/* attribute recognised by some compilers to avoid 'unused' warnings */
#ifndef SWIGUNUSED
# if defined(__GNUC__) || defined(__ICC)
#   define SWIGUNUSED __attribute__ ((unused)) 
# else
#   define SWIGUNUSED 
# endif
#endif

/* internal SWIG method */
#ifndef SWIGINTERN
# define SWIGINTERN static SWIGUNUSED
#endif

/* internal inline SWIG method */
#ifndef SWIGINTERNINLINE
# define SWIGINTERNINLINE SWIGINTERN SWIGINLINE
#endif

/* exporting methods for Windows DLLs */
#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   define SWIGEXPORT
# endif
#endif

/* calling conventions for Windows */
#ifndef SWIGSTDCALL
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   define SWIGSTDCALL __stdcall
# else
#   define SWIGSTDCALL
# endif 
#endif



/* Fix for jlong on some versions of gcc on Windows */
#if defined(__GNUC__) && !defined(__INTELC__)
  typedef long long __int64;
#endif

/* Fix for jlong on 64-bit x86 Solaris */
#if defined(__x86_64)
# ifdef _LP64
#   undef _LP64
# endif
#endif

#include <jni.h>
#include <stdlib.h>
#include <string.h>


/* Support for throwing Java exceptions */
typedef enum {
  SWIG_JavaOutOfMemoryError = 1, 
  SWIG_JavaIOException, 
  SWIG_JavaRuntimeException, 
  SWIG_JavaIndexOutOfBoundsException,
  SWIG_JavaArithmeticException,
  SWIG_JavaIllegalArgumentException,
  SWIG_JavaNullPointerException,
  SWIG_JavaDirectorPureVirtual,
  SWIG_JavaUnknownError
} SWIG_JavaExceptionCodes;

typedef struct {
  SWIG_JavaExceptionCodes code;
  const char *java_exception;
} SWIG_JavaExceptions_t;


static void SWIG_JavaThrowException(JNIEnv *jenv, SWIG_JavaExceptionCodes code, const char *msg) {
  jclass excep;
  static const SWIG_JavaExceptions_t java_exceptions[] = {
    { SWIG_JavaOutOfMemoryError, "java/lang/OutOfMemoryError" },
    { SWIG_JavaIOException, "java/io/IOException" },
    { SWIG_JavaRuntimeException, "java/lang/RuntimeException" },
    { SWIG_JavaIndexOutOfBoundsException, "java/lang/IndexOutOfBoundsException" },
    { SWIG_JavaArithmeticException, "java/lang/ArithmeticException" },
    { SWIG_JavaIllegalArgumentException, "java/lang/IllegalArgumentException" },
    { SWIG_JavaNullPointerException, "java/lang/NullPointerException" },
    { SWIG_JavaDirectorPureVirtual, "java/lang/RuntimeException" },
    { SWIG_JavaUnknownError,  "java/lang/UnknownError" },
    { (SWIG_JavaExceptionCodes)0,  "java/lang/UnknownError" } };
  const SWIG_JavaExceptions_t *except_ptr = java_exceptions;

  while (except_ptr->code != code && except_ptr->code)
    except_ptr++;

  jenv->ExceptionClear();
  excep = jenv->FindClass(except_ptr->java_exception);
  if (excep)
    jenv->ThrowNew(excep, msg);
}


/* Contract support */

#define SWIG_contract_assert(nullreturn, expr, msg) if (!(expr)) {SWIG_JavaThrowException(jenv, SWIG_JavaIllegalArgumentException, msg); return nullreturn; } else


	#include <iostream>

extern bool det(std::istream& matrix_in, std::ostream& det_out);
extern char* detFiles(char *matfile);
extern bool rank(std::istream& matrix_in, std::ostream& rank_out);
extern char* rankFiles(char *matfile);
extern int estimateRankTime(char *matfile);
extern bool val(std::istream& matrix_in, std::ostream& val_out);
extern char* valFiles(char *matfile);
extern bool trace(std::istream& matrix_in, std::ostream& trace_out);
extern char* traceFiles(char *matfile);
extern bool smithNormalForm(std::istream& matrix_in, std::ostream& snf_out);
extern char* smithNormalFormFiles(char *matfile);


#ifdef __cplusplus
extern "C" {
#endif

JNIEXPORT jboolean JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_det(JNIEnv *jenv, jclass jcls, jlong jarg1, jlong jarg2) {
    jboolean jresult = 0 ;
    std::istream *arg1 = 0 ;
    std::ostream *arg2 = 0 ;
    bool result;
    
    (void)jenv;
    (void)jcls;
    arg1 = *(std::istream **)(void *)&jarg1;
    if(!arg1) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::istream & reference is null");
        return 0;
    } 
    arg2 = *(std::ostream **)(void *)&jarg2;
    if(!arg2) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::ostream & reference is null");
        return 0;
    } 
    result = (bool)det(*arg1,*arg2);
    
    jresult = (jboolean)result; 
    return jresult;
}


JNIEXPORT jstring JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_detFiles(JNIEnv *jenv, jclass jcls, jstring jarg1) {
    jstring jresult = 0 ;
    char *arg1 = (char *) 0 ;
    char *result;
    
    (void)jenv;
    (void)jcls;
    {
        arg1 = 0;
        if (jarg1) {
            arg1 = (char *)jenv->GetStringUTFChars(jarg1, 0);
            if (!arg1) return 0;
        }
    }
    result = (char *)detFiles(arg1);
    
    {
        if(result) jresult = jenv->NewStringUTF(result); 
    }
    {
        if (arg1) jenv->ReleaseStringUTFChars(jarg1, arg1); 
    }
    return jresult;
}


JNIEXPORT jboolean JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_rank(JNIEnv *jenv, jclass jcls, jlong jarg1, jlong jarg2) {
    jboolean jresult = 0 ;
    std::istream *arg1 = 0 ;
    std::ostream *arg2 = 0 ;
    bool result;
    
    (void)jenv;
    (void)jcls;
    arg1 = *(std::istream **)(void *)&jarg1;
    if(!arg1) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::istream & reference is null");
        return 0;
    } 
    arg2 = *(std::ostream **)(void *)&jarg2;
    if(!arg2) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::ostream & reference is null");
        return 0;
    } 
    result = (bool)rank(*arg1,*arg2);
    
    jresult = (jboolean)result; 
    return jresult;
}


JNIEXPORT jstring JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_rankFiles(JNIEnv *jenv, jclass jcls, jstring jarg1) {
    jstring jresult = 0 ;
    char *arg1 = (char *) 0 ;
    char *result;
    
    (void)jenv;
    (void)jcls;
    {
        arg1 = 0;
        if (jarg1) {
            arg1 = (char *)jenv->GetStringUTFChars(jarg1, 0);
            if (!arg1) return 0;
        }
    }
    result = (char *)rankFiles(arg1);
    
    {
        if(result) jresult = jenv->NewStringUTF(result); 
    }
    {
        if (arg1) jenv->ReleaseStringUTFChars(jarg1, arg1); 
    }
    return jresult;
}


JNIEXPORT jint JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_estimateRankTime(JNIEnv *jenv, jclass jcls, jstring jarg1) {
    jint jresult = 0 ;
    char *arg1 = (char *) 0 ;
    int result;
    
    (void)jenv;
    (void)jcls;
    {
        arg1 = 0;
        if (jarg1) {
            arg1 = (char *)jenv->GetStringUTFChars(jarg1, 0);
            if (!arg1) return 0;
        }
    }
    result = (int)estimateRankTime(arg1);
    
    jresult = (jint)result; 
    {
        if (arg1) jenv->ReleaseStringUTFChars(jarg1, arg1); 
    }
    return jresult;
}


JNIEXPORT jboolean JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_val(JNIEnv *jenv, jclass jcls, jlong jarg1, jlong jarg2) {
    jboolean jresult = 0 ;
    std::istream *arg1 = 0 ;
    std::ostream *arg2 = 0 ;
    bool result;
    
    (void)jenv;
    (void)jcls;
    arg1 = *(std::istream **)(void *)&jarg1;
    if(!arg1) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::istream & reference is null");
        return 0;
    } 
    arg2 = *(std::ostream **)(void *)&jarg2;
    if(!arg2) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::ostream & reference is null");
        return 0;
    } 
    result = (bool)val(*arg1,*arg2);
    
    jresult = (jboolean)result; 
    return jresult;
}


JNIEXPORT jstring JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_valFiles(JNIEnv *jenv, jclass jcls, jstring jarg1) {
    jstring jresult = 0 ;
    char *arg1 = (char *) 0 ;
    char *result;
    
    (void)jenv;
    (void)jcls;
    {
        arg1 = 0;
        if (jarg1) {
            arg1 = (char *)jenv->GetStringUTFChars(jarg1, 0);
            if (!arg1) return 0;
        }
    }
    result = (char *)valFiles(arg1);
    
    {
        if(result) jresult = jenv->NewStringUTF(result); 
    }
    {
        if (arg1) jenv->ReleaseStringUTFChars(jarg1, arg1); 
    }
    return jresult;
}


JNIEXPORT jboolean JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_trace(JNIEnv *jenv, jclass jcls, jlong jarg1, jlong jarg2) {
    jboolean jresult = 0 ;
    std::istream *arg1 = 0 ;
    std::ostream *arg2 = 0 ;
    bool result;
    
    (void)jenv;
    (void)jcls;
    arg1 = *(std::istream **)(void *)&jarg1;
    if(!arg1) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::istream & reference is null");
        return 0;
    } 
    arg2 = *(std::ostream **)(void *)&jarg2;
    if(!arg2) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::ostream & reference is null");
        return 0;
    } 
    result = (bool)trace(*arg1,*arg2);
    
    jresult = (jboolean)result; 
    return jresult;
}


JNIEXPORT jstring JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_traceFiles(JNIEnv *jenv, jclass jcls, jstring jarg1) {
    jstring jresult = 0 ;
    char *arg1 = (char *) 0 ;
    char *result;
    
    (void)jenv;
    (void)jcls;
    {
        arg1 = 0;
        if (jarg1) {
            arg1 = (char *)jenv->GetStringUTFChars(jarg1, 0);
            if (!arg1) return 0;
        }
    }
    result = (char *)traceFiles(arg1);
    
    {
        if(result) jresult = jenv->NewStringUTF(result); 
    }
    {
        if (arg1) jenv->ReleaseStringUTFChars(jarg1, arg1); 
    }
    return jresult;
}


JNIEXPORT jboolean JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_smithNormalForm(JNIEnv *jenv, jclass jcls, jlong jarg1, jlong jarg2) {
    jboolean jresult = 0 ;
    std::istream *arg1 = 0 ;
    std::ostream *arg2 = 0 ;
    bool result;
    
    (void)jenv;
    (void)jcls;
    arg1 = *(std::istream **)(void *)&jarg1;
    if(!arg1) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::istream & reference is null");
        return 0;
    } 
    arg2 = *(std::ostream **)(void *)&jarg2;
    if(!arg2) {
        SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "std::ostream & reference is null");
        return 0;
    } 
    result = (bool)smithNormalForm(*arg1,*arg2);
    
    jresult = (jboolean)result; 
    return jresult;
}


JNIEXPORT jstring JNICALL Java_samples_quickstart_service_adb_xsd_linboxfunctionsJNI_smithNormalFormFiles(JNIEnv *jenv, jclass jcls, jstring jarg1) {
    jstring jresult = 0 ;
    char *arg1 = (char *) 0 ;
    char *result;
    
    (void)jenv;
    (void)jcls;
    {
        arg1 = 0;
        if (jarg1) {
            arg1 = (char *)jenv->GetStringUTFChars(jarg1, 0);
            if (!arg1) return 0;
        }
    }
    result = (char *)smithNormalFormFiles(arg1);
    
    {
        if(result) jresult = jenv->NewStringUTF(result); 
    }
    {
        if (arg1) jenv->ReleaseStringUTFChars(jarg1, arg1); 
    }
    return jresult;
}


#ifdef __cplusplus
}
#endif

