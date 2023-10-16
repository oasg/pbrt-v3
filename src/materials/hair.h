
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

#ifndef PBRT_MATERIALS_HAIR_H
#define PBRT_MATERIALS_HAIR_H

// materials/hair.h*
#include "material.h"
#include "pbrt.h"
#include "reflection.h"
#include <array>

namespace pbrt {


// HairMaterial Declarations
class HairMaterial : public Material {
  public:
    // HairMaterial Public Methods
    HairMaterial(const std::shared_ptr<Texture<Spectrum>> &sigma_a,
                 const std::shared_ptr<Texture<Spectrum>> &color,
                 const std::shared_ptr<Texture<Float>> &eumelanin,
                 const std::shared_ptr<Texture<Float>> &pheomelanin,
                 const std::shared_ptr<Texture<Float>> &eta,
                 const std::shared_ptr<Texture<Float>> &beta_m,
                 const std::shared_ptr<Texture<Float>> &beta_n,
                 const std::shared_ptr<Texture<Float>> &alpha)
        : sigma_a(sigma_a),
          color(color),
          eumelanin(eumelanin),
          pheomelanin(pheomelanin),
          eta(eta),
          beta_m(beta_m),
          beta_n(beta_n),
          alpha(alpha) {}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // HairMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> sigma_a, color;
    std::shared_ptr<Texture<Float>> eumelanin, pheomelanin, eta;
    std::shared_ptr<Texture<Float>> beta_m, beta_n, alpha;
};

HairMaterial *CreateHairMaterial(const TextureParams &mp);

// HairBSDF Constants
static const int pMax = 3;
static const Float SqrtPiOver8 = 0.626657069f;

// HairBSDF Declarations
class HairBSDF : public BxDF {
  public:
    // HairBSDF Public Methods
    HairBSDF(Float h, Float eta, const Spectrum &sigma_a, Float beta_m,
             Float beta_n, Float alpha);
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;
    static Spectrum SigmaAFromConcentration(Float ce, Float cp);
    static Spectrum SigmaAFromReflectance(const Spectrum &c, Float beta_n);

  protected:  //share data with MHairMaterial friend class
    // HairBSDF Private Methods
    std::array<Float, pMax + 1> ComputeApPdf(Float cosThetaO) const;

    // HairBSDF Private Data
    const Float h, gammaO, eta;
    const Spectrum sigma_a;
    const Float beta_m, beta_n;
    Float v[pMax + 1];
    Float s;
    Float sin2kAlpha[3], cos2kAlpha[3];
};

// General Utility Functions
template <int n>
static Float Pow(Float v) {
    static_assert(n > 0, "Power can't be negative");
    Float n2 = Pow<n / 2>(v);
    return n2 * n2 * Pow<n & 1>(v);
}

template <>
inline Float Pow<1>(Float v) {
    return v;
}
template <>
inline Float Pow<0>(Float v) {
    return 1;
}

inline Float SafeASin(Float x) {
    CHECK(x >= -1.0001 && x <= 1.0001);
    return std::asin(Clamp(x, -1, 1));
}

inline Float SafeSqrt(Float x) {
    CHECK_GE(x, -1e-4);
    return std::sqrt(std::max(Float(0), x));
}

// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
static uint32_t Compact1By1(uint32_t x) {
    // TODO: as of Haswell, the PEXT instruction could do all this in a
    // single instruction.
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x &= 0x55555555;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    return x;
}

static Point2f DemuxFloat(Float f) {
    CHECK(f >= 0 && f < 1);
    uint64_t v = f * (1ull << 32);
    CHECK_LT(v, 0x100000000);
    uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
    return {bits[0] / Float(1 << 16), bits[1] / Float(1 << 16)};
}
//move some utility function to the header file
// to share with my own hair material
// and this utility function could only build once
inline Float Sqr(Float v) { return v * v; }
  inline Float I0(Float x) {
    Float val = 0;
    Float x2i = 1;
    int64_t ifact = 1;
    int i4 = 1;
    // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
    for (int i = 0; i < 10; ++i) {
        if (i > 1) ifact *= i;
        val += x2i / (i4 * Sqr(ifact));
        x2i *= x * x;
        i4 *= 4;
    }
    return val;
}
inline Float Phi(int p, Float gammaO, Float gammaT) {
    return 2 * p * gammaT - 2 * gammaO + p * Pi;
}

inline Float Logistic(Float x, Float s) {
    x = std::abs(x);
    return std::exp(-x / s) / (s * Sqr(1 + std::exp(-x / s)));
}

inline Float LogisticCDF(Float x, Float s) {
    return 1 / (1 + std::exp(-x / s));
}

inline Float TrimmedLogistic(Float x, Float s, Float a, Float b) {
    CHECK_LT(a, b);
    return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}
static std::array<Spectrum, pMax + 1> Ap(Float cosThetaO, Float eta, Float h,
                                         const Spectrum &T) {
    std::array<Spectrum, pMax + 1> ap;
    // Compute $p=0$ attenuation at initial cylinder intersection
    Float cosGammaO = SafeSqrt(1 - h * h);
    Float cosTheta = cosThetaO * cosGammaO;
    Float f = FrDielectric(cosTheta, 1.f, eta);
    ap[0] = f;

    // Compute $p=1$ attenuation term
    ap[1] = Sqr(1 - f) * T;

    // Compute attenuation terms up to $p=_pMax_$
    for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

    // Compute attenuation term accounting for remaining orders of scattering
    ap[pMax] = ap[pMax - 1] * f * T / (Spectrum(1.f) - T * f);
    return ap;
}
inline Float LogI0(Float x) {
    if (x > 12)
        return x + 0.5 * (-std::log(2 * Pi) + std::log(1 / x) + 1 / (8 * x));
    else
        return std::log(I0(x));
}
static Float Mp(Float cosThetaI, Float cosThetaO, Float sinThetaI,
                Float sinThetaO, Float v) {
    Float a = cosThetaI * cosThetaO / v;
    Float b = sinThetaI * sinThetaO / v;
    Float mp =
        (v <= .1)
            ? (std::exp(LogI0(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v))))
            : (std::exp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
    CHECK(!std::isinf(mp) && !std::isnan(mp));
    return mp;
}
inline Float Np(Float phi, int p, Float s, Float gammaO, Float gammaT) {
    Float dphi = phi - Phi(p, gammaO, gammaT);
    // Remap _dphi_ to $[-\pi,\pi]$
    while (dphi > Pi) dphi -= 2 * Pi;
    while (dphi < -Pi) dphi += 2 * Pi;
    return TrimmedLogistic(dphi, s, -Pi, Pi);
}





static Float SampleTrimmedLogistic(Float u, Float s, Float a, Float b) {
    CHECK_LT(a, b);
    Float k = LogisticCDF(b, s) - LogisticCDF(a, s);
    Float x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
    CHECK(!std::isnan(x));
    return Clamp(x, a, b);
}
// Hair Local Declarations
inline Float I0(Float x), LogI0(Float x);

// Hair Local Functions









}  // namespace pbrt

#endif  // PBRT_MATERIALS_HAIR_H
