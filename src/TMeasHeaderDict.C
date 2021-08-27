// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdITMeasHeaderDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/tracs/work/tracs-master/include/TMeasHeader.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TMeasHeader(void *p = 0);
   static void *newArray_TMeasHeader(Long_t size, void *p);
   static void delete_TMeasHeader(void *p);
   static void deleteArray_TMeasHeader(void *p);
   static void destruct_TMeasHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TMeasHeader*)
   {
      ::TMeasHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TMeasHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TMeasHeader", ::TMeasHeader::Class_Version(), "TMeasHeader.h", 24,
                  typeid(::TMeasHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TMeasHeader::Dictionary, isa_proxy, 4,
                  sizeof(::TMeasHeader) );
      instance.SetNew(&new_TMeasHeader);
      instance.SetNewArray(&newArray_TMeasHeader);
      instance.SetDelete(&delete_TMeasHeader);
      instance.SetDeleteArray(&deleteArray_TMeasHeader);
      instance.SetDestructor(&destruct_TMeasHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TMeasHeader*)
   {
      return GenerateInitInstanceLocal((::TMeasHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TMeasHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TMeasHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TMeasHeader::Class_Name()
{
   return "TMeasHeader";
}

//______________________________________________________________________________
const char *TMeasHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TMeasHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TMeasHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TMeasHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TMeasHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TMeasHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TMeasHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TMeasHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TMeasHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class TMeasHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TMeasHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(TMeasHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TMeasHeader(void *p) {
      return  p ? new(p) ::TMeasHeader : new ::TMeasHeader;
   }
   static void *newArray_TMeasHeader(Long_t nElements, void *p) {
      return p ? new(p) ::TMeasHeader[nElements] : new ::TMeasHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_TMeasHeader(void *p) {
      delete ((::TMeasHeader*)p);
   }
   static void deleteArray_TMeasHeader(void *p) {
      delete [] ((::TMeasHeader*)p);
   }
   static void destruct_TMeasHeader(void *p) {
      typedef ::TMeasHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TMeasHeader

namespace {
  void TriggerDictionaryInitialization_TMeasHeaderDict_Impl() {
    static const char* headers[] = {
"include/TMeasHeader.h",
0
    };
    static const char* includePaths[] = {
"/home/tracs/work/tracs-master/include/",
"/usr/local/root/include",
"/home/tracs/work/tracs-master/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TMeasHeaderDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Edge-TCT data header class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$include/TMeasHeader.h")))  TMeasHeader;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TMeasHeaderDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "include/TMeasHeader.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TMeasHeader", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TMeasHeaderDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TMeasHeaderDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TMeasHeaderDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TMeasHeaderDict() {
  TriggerDictionaryInitialization_TMeasHeaderDict_Impl();
}
