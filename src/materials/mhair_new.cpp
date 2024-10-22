
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

// materials/hair.cpp*
#include "materials/mhair_new.h"

#include <array>
#include <numeric>
#include <fstream>
#include "interaction.h"
#include "mhair.h"
#include "paramset.h"
#include "reflection.h"
#include "sampling.h"
#include "spectrum.h"
#include "texture.h"
#include "textures/constant.h"
#include"../table/brdf_c5.h" //キューティクル多層
//#include"../table/brdf_nc5.h" //キューティクルなし
//#include"../tablebrdf_hk5.h" //キューティクル剥落
//#include"../table/brdf_ck5.h"//キューティクル欠け
//#include"../table/brdf_hkck5.h"//キューティクル欠け+剥落

namespace pbrt {

std::mutex SingBrdf::m_mutex;
std::shared_ptr<hairSimBrdf>  SingBrdf::m_instance_ptr = nullptr;
MHairNewMaterial *CreateMHairNewMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> sigma_a =
        mp.GetSpectrumTextureOrNull("sigma_a");
    std::shared_ptr<Texture<Spectrum>> color =
        mp.GetSpectrumTextureOrNull("color");
    std::shared_ptr<Texture<Float>> eumelanin =
        mp.GetFloatTextureOrNull("eumelanin");
    std::shared_ptr<Texture<Float>> pheomelanin =
        mp.GetFloatTextureOrNull("pheomelanin");
    if (sigma_a) {
        if (color)
            Warning(
                "Ignoring \"color\" parameter since \"sigma_a\" was provided.");
        if (eumelanin)
            Warning(
                "Ignoring \"eumelanin\" parameter since \"sigma_a\" was "
                "provided.");
        if (pheomelanin)
            Warning(
                "Ignoring \"pheomelanin\" parameter since \"sigma_a\" was "
                "provided.");
    } else if (color) {
        if (sigma_a)
            Warning(
                "Ignoring \"sigma_a\" parameter since \"color\" was provided.");
        if (eumelanin)
            Warning(
                "Ignoring \"eumelanin\" parameter since \"color\" was "
                "provided.");
        if (pheomelanin)
            Warning(
                "Ignoring \"pheomelanin\" parameter since \"color\" was "
                "provided.");
    } else if (eumelanin || pheomelanin) {
        if (sigma_a)
            Warning(
                "Ignoring \"sigma_a\" parameter since "
                "\"eumelanin\"/\"pheomelanin\" was provided.");
        if (color)
            Warning(
                "Ignoring \"color\" parameter since "
                "\"eumelanin\"/\"pheomelanin\" was provided.");
    } else {
        // Default: brown-ish hair.
        sigma_a = std::make_shared<ConstantTexture<Spectrum>>(
            HairBSDF::SigmaAFromConcentration(1.3, 0.));
    }

    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTexture("eta", 1.55f);
    std::shared_ptr<Texture<Float>> beta_m = mp.GetFloatTexture("beta_m", 0.3f);
    std::shared_ptr<Texture<Float>> beta_n = mp.GetFloatTexture("beta_n", 0.3f);
    std::shared_ptr<Texture<Float>> alpha = mp.GetFloatTexture("alpha", 2.f);

    return new MHairNewMaterial(sigma_a, color, eumelanin, pheomelanin, eta, beta_m,
                            beta_n, alpha);
}
MHairNewBSDF::MHairNewBSDF(Float h, Float eta, const Spectrum &sigma_a, Float beta_m,
        Float beta_n, Float alpha):HairBSDF(h,eta,sigma_a,beta_m,beta_n,alpha) {
        m_sim_brdf = SingBrdf::get_Instance();
    }
Spectrum MHairNewBSDF::f(const Vector3f &wo, const Vector3f &wi) const {
        // Compute hair coordinate system terms related to _wo_
    Float sinThetaO = wo.x;
    Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    Float phiO = std::atan2(wo.z, wo.y);

    // Compute hair coordinate system terms related to _wi_
    Float sinThetaI = wi.x;
    Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
    Float phiI = std::atan2(wi.z, wi.y);

    // Compute $\cos \thetat$ for refracted ray
    Float sinThetaT = sinThetaO / eta;
    Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    // Compute $\gammat$ for refracted ray
    Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
    Float sinGammaT = h / etap;
    Float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
    Float gammaT = SafeASin(sinGammaT);

    //角度を算出 θi, θo
    auto wo_theta = std::atan2(wo.x, std::sqrt(wo.y * wo.y + wo.z * wo.z));
    auto wo_ang = wo_theta * 180 / Pi;
    Float ot = std::abs(wo_ang);
    ot = ot + 90.0;

    auto wi_theta = std::atan2(wi.x, std::sqrt(wi.y * wi.y + wi.z * wi.z));
    auto wi_ang = wi_theta * 180 / Pi;
    Float it = std::abs(wi_ang);
    int x0 = floor(it);
    int x1 = ceil(it);

    int y0 = floor(ot);
    int y1 = ceil(ot);

    //remain
    Float alpha_x = it - x0;
    Float alpha_y = ot - y0;


        // 获取四个邻近点的颜色值
    auto c00 = m_sim_brdf->m_data[x0][y0];
    auto c10 = m_sim_brdf->m_data[x1][y0];
    auto c01 = m_sim_brdf->m_data[x0][y1];
    auto c11 = m_sim_brdf->m_data[x1][y1];


    // 对每个通道分别进行 x 轴插值
    Float r_x0 = (1 - alpha_x) * c00.r + alpha_x * c10.r;
    Float r_x1 = (1 - alpha_x) * c01.r + alpha_x * c11.r;
    Float g_x0 = (1 - alpha_x) * c00.g + alpha_x * c10.g;
    Float g_x1 = (1 - alpha_x) * c01.g + alpha_x * c11.g;
    Float b_x0 = (1 - alpha_x) * c00.b + alpha_x * c10.b;
    Float b_x1 = (1 - alpha_x) * c00.b + alpha_x * c11.b;

    // 对每个通道分别进行 y 轴插值，得到最终颜色值
    Float r_final = (1 - alpha_y) * r_x0 + alpha_y * r_x1;
    Float g_final = (1 - alpha_y) * g_x0 + alpha_y * g_x1;
    Float b_final = (1 - alpha_y) * b_x0 + alpha_y * b_x1;

    // 将结果除以10（根据你的原始代码）以得到最终结果
    const Float crgb[3] = {r_final / 10, g_final / 10, b_final / 10};

    // if(it > 90){
    //     it = it - 90;
    // }

    // if(ot < 0){
    //     ot = ot + 90;
    // }
    // it = it%90;
    // ot = (ot-90+180)%180;

    //std::cout<<"it:"<<it<<" ot:"<<ot<<std::endl;
    // // use angle in azimuthal
    // int phiI_ang = static_cast<int>(std::abs((std ::round( phiI * 180 / Pi))));
    // int phiO_ang = static_cast<int>(std::abs((std ::round(phiO * 180 / Pi))));
    // if (phiI_ang > 90) {
    //     phiI_ang = phiI_ang - 90;
    // }
    // if (phiO_ang < 0) {
    //     phiO_ang = phiO_ang + 90;
    // }

    // Compute the transmittance _T_ of a single path through the cylinder
    Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));

    // Evaluate hair BSDF
    Float phi = phiI - phiO;
    std::array<Spectrum, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
    Spectrum fsum(0.);
    //for p = 0
    // samplize the table
    // SampledSpectrum RN(0.);
    // for(int i =0;i<60;++i){
    //     RN[i] = BRDFTABLE5[it][ot][i] / 2.5;
    //     RN[i] *= SampledSpectrum::get_rgbIllum2SpectWhite()[i];
    // }
    Float sinThetaOp, cosThetaOp;
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
    cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
    //compute p =0 
    //fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[0]) * ap[0]*RN.ToRGBSpectrum();


    //std::cout<<rgb.r<< rgb.g<< rgb.b<<std::endl;
    RGBSpectrum reflect = RGBSpectrum::FromRGB(crgb);
    fsum += reflect;
    // // p >=1
    // for (int p = 1; p < pMax; ++p) {
    //     // Handle remainder of $p$ values for hair scale tilt
    //     if (p == 1) {
    //         sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
    //         cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
    //     } else if (p == 2) {
    //         sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
    //         cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
    //     } else {
    //         sinThetaOp = sinThetaO;
    //         cosThetaOp = cosThetaO;
    //     }

    //     // Handle out-of-range $\cos \thetao$ from scale adjustment
    //     cosThetaOp = std::abs(cosThetaOp);
    //     fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * ap[p] *
    //             Np(phi, p, s, gammaO, gammaT);
    // }

    // // Compute contribution of remaining terms after _pMax_
    // fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * ap[pMax] /
    //         (2.f * Pi);
    // if (AbsCosTheta(wi) > 0) fsum /= AbsCosTheta(wi);
    // CHECK(!std::isinf(fsum.y()) && !std::isnan(fsum.y()));
    return fsum;
}
void MHairNewMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                               MemoryArena &arena,
                                               TransportMode mode,
                                               bool allowMultipleLobes) const {
    Float bm = beta_m->Evaluate(*si);
    Float bn = beta_n->Evaluate(*si);
    Float a = alpha->Evaluate(*si);
    Float e = eta->Evaluate(*si);

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, e);

    Spectrum sig_a;
    if (sigma_a)
        sig_a = sigma_a->Evaluate(*si).Clamp();
    else if (color) {
        Spectrum c = color->Evaluate(*si).Clamp();
        sig_a = MHairNewBSDF::SigmaAFromReflectance(c, bn);
    } else {
        CHECK(eumelanin || pheomelanin);
        sig_a = MHairNewBSDF::SigmaAFromConcentration(
            std::max(Float(0), eumelanin ? eumelanin->Evaluate(*si) : 0),
            std::max(Float(0), pheomelanin ? pheomelanin->Evaluate(*si) : 0));
    }

    // Offset along width
    Float h = -1 + 2 * si->uv[1];
    si->bsdf->Add(ARENA_ALLOC(arena, MHairNewBSDF)(h, e, sig_a, bm, bn, a));

    
}
hairSimBrdf::hairSimBrdf(const char *file) {
    m_data = std::vector<std::vector<RGB>>(91,std::vector<RGB>(181));
    //read from file
    std::ifstream infile(file);
    if(!infile.is_open()){
        std::cout<<"can't open file:"<<file<<std::endl;
    }
    float a=0.1,b=0,c=0;
    int indeg = 0;
    int outdeg = 0;
    std::string line;
    while(std::getline(infile,line)){
        if(line == "save successed"){
                break;
            }
        if(std::abs(a)<1e-9){
            indeg++;
            outdeg = 0;
        }
        std::istringstream ss(line);
        std::string temp;
        if (std::getline(ss, temp, ',')) {
            
            std::stringstream tempStream(temp);
            tempStream >> a;
        }
        if (std::getline(ss, temp, ',')) {
            std::stringstream tempStream(temp);
            tempStream >> b;
        }
        if (std::getline(ss, temp, ',')) {
            std::stringstream tempStream(temp);
            tempStream >> c;
        }

        m_data[indeg][outdeg].r = a;
        m_data[indeg][outdeg].g = b;
        m_data[indeg][outdeg].b = c;
        outdeg++;
    }
    std::cout<<"gen hair brdf ok!!"<<file<<std::endl;
}
}  
