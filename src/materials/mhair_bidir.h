
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

/*

See the writeup "The Implementation of a Hair Scattering Model" at
http://pbrt.org/hair.pdf for a description of the implementation here.

*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_MHAIRBI_NEW_H
#define PBRT_MATERIALS_MHAIRBI_NEW_H

// materials/hair.h*
#include "material.h"
#include "pbrt.h"
#include "reflection.h"
#include <array>
#include<vector>
#include<mutex>
#include "hair.h"

namespace pbrt {
class hairSimBiBrdf{
  public:
    hairSimBiBrdf(const char* fileazi,const char* filelong);
    ~hairSimBiBrdf(){}
    Float getReflectA(Float it, Float ot);
    Float getReflectL(Float it, Float ot);
    std::vector<std::vector<RGB>> m_data_a;
    std::vector<std::vector<RGB>> m_data_l;
};

class SingBiBrdf{
  public:
  ~SingBiBrdf(){
    std::cout<<"sim Brdf destructor!"<<std::endl;
  }
  SingBiBrdf() = delete;
  SingBiBrdf& operator=(const SingBiBrdf&)= delete;
  static std::shared_ptr<hairSimBiBrdf> get_Instance(){
    if(m_instance_ptr==nullptr){
      std::lock_guard<std::mutex> lk(m_mutex);
      if(m_instance_ptr==nullptr){
        m_instance_ptr = std::shared_ptr<hairSimBiBrdf>(
          new hairSimBiBrdf("../table/HairMultilayerPerlinModel/incidenceLayer/(10nm,800cell)/TM/Ns/Reflection/","../table/HairMultilayerLongitudinalModel/incidenceLayer/(10nm,1200cell)/TM/Ns/Reflection/"));
      }
    }
    return m_instance_ptr;
  }
  private:
    static std::shared_ptr<hairSimBiBrdf> m_instance_ptr;
    static std::mutex m_mutex;
    
};



// MHairMaterial Declarations
class MHairNewBiMaterial : public HairMaterial {
  public:
    // MHairMaterial Public Methods
    MHairNewBiMaterial(const std::shared_ptr<Texture<Spectrum>> &sigma_a,
                 const std::shared_ptr<Texture<Spectrum>> &color,
                 const std::shared_ptr<Texture<Float>> &eumelanin,
                 const std::shared_ptr<Texture<Float>> &pheomelanin,
                 const std::shared_ptr<Texture<Float>> &eta,
                 const std::shared_ptr<Texture<Float>> &beta_m,
                 const std::shared_ptr<Texture<Float>> &beta_n,
                 const std::shared_ptr<Texture<Float>> &alpha)
        : HairMaterial(sigma_a, color, eumelanin, pheomelanin, eta, beta_m,beta_n,alpha) {}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;
  private:
    // MHairMaterial Private Data

};

MHairNewBiMaterial *CreateMHairNewBiMaterial(const TextureParams &mp);

// HairBSDF Declarations
class MHairNewBiBSDF : public HairBSDF {
  public:
    // HairBSDF Public Methods
    // just need the f function
    // pdf and sampled function are inherited
    MHairNewBiBSDF(Float h, Float eta, const Spectrum &sigma_a, Float beta_m,
             Float beta_n, Float alpha);
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;

  private:
    std::shared_ptr<hairSimBiBrdf> m_sim_brdf;
    // MHairBSDF Private data

};

}  // namespace pbrt

#endif  // PBRT_MATERIALS_HAIR_H
