
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
#include "materials/mhair_bidir.h"

#include <array>
#include <numeric>
#include <fstream>
#include "interaction.h"
#include "hair.h"
#include "paramset.h"
#include "reflection.h"
#include "sampling.h"
#include "spectrum.h"
#include "texture.h"
#include "textures/constant.h"

namespace pbrt {

std::mutex SingBiBrdf::m_mutex;
std::shared_ptr<hairSimBiBrdf>  SingBiBrdf::m_instance_ptr = nullptr;

MHairNewBiMaterial *CreateMHairNewBiMaterial(const TextureParams &mp) {
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

    return new MHairNewBiMaterial(sigma_a, color, eumelanin, pheomelanin, eta, beta_m,
                            beta_n, alpha);
}

RGB GetRefData(Float it,Float ot,std::vector<std::vector<RGB>>& m_data){
    int x0 = floor(it);
    int x1 = ceil(it);

    int y0 = floor(ot);
    int y1 = ceil(ot);

    //remain
    Float alpha_x = it - x0;
    Float alpha_y = ot - y0;


        // 获取四个邻近点的颜色值
    auto c00 = m_data[x0][y0];
    auto c10 = m_data[x1][y0];
    auto c01 = m_data[x0][y1];
    auto c11 = m_data[x1][y1];


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
    return RGB{r_final,g_final,b_final};

}
MHairNewBiBSDF::MHairNewBiBSDF(Float h, Float eta, const Spectrum &sigma_a, Float beta_m,
        Float beta_n, Float alpha):HairBSDF(h,eta,sigma_a,beta_m,beta_n,alpha) {
        m_sim_brdf = SingBiBrdf::get_Instance();
    }
Spectrum MHairNewBiBSDF::f(const Vector3f &wo, const Vector3f &wi) const {
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

    auto phi_in = std::abs(phiI * 180 /Pi);
    if(phi_in > 90.0){
        phi_in-=90.0;
    }
    auto phi_out =std::abs(phiO*180/Pi);


    auto theta_in = std::abs(std::asin(sinThetaI));
    theta_in = theta_in*180/Pi;
    auto theta_out = std::abs(std::asin(sinThetaO));
    theta_out = theta_out*180/Pi;

    // 将结果除以10（根据你的原始代码）以得到最终结果

    auto resa = GetRefData(phi_in,phi_out,m_sim_brdf->m_data_a);
    auto resl = GetRefData(theta_in,theta_out,m_sim_brdf->m_data_l);
    const Float crgb[3] = {resa.r*0.5+resl.r*0.5 ,
        resa.g*0.5+resl.g*0.5, 
        resa.g*0.5+resl.g*0.5};

    // Compute the transmittance _T_ of a single path through the cylinder
    Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));

    // Evaluate hair BSDF
    Float phi = phiI - phiO;
    std::array<Spectrum, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
    Spectrum fsum(0.);
    //std::cout<<rgb.r<< rgb.g<< rgb.b<<std::endl;
    RGBSpectrum reflect = RGBSpectrum::FromRGB(crgb);
    fsum = reflect;
    if (AbsCosTheta(wi) > 0) fsum /= AbsCosTheta(wi);

    CHECK(!std::isinf(fsum.y()) && !std::isnan(fsum.y()));
    return fsum;
}
void MHairNewBiMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
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
        sig_a = MHairNewBiBSDF::SigmaAFromReflectance(c, bn);
    } else {
        CHECK(eumelanin || pheomelanin);
        sig_a = MHairNewBiBSDF::SigmaAFromConcentration(
            std::max(Float(0), eumelanin ? eumelanin->Evaluate(*si) : 0),
            std::max(Float(0), pheomelanin ? pheomelanin->Evaluate(*si) : 0));
    }

    // Offset along width
    Float h = -1 + 2 * si->uv[1];
    si->bsdf->Add(ARENA_ALLOC(arena, MHairNewBiBSDF)(h, e, sig_a, bm, bn, a));

    
}

void ReadFromFileAndTrans(std::vector<std::vector<RGB>>& m_data,const char* file){
    std::vector<std::vector<std::vector<double>>> data(91,std::vector<std::vector<double>>(181,std::vector<double>(60)));
    for (int a = 0; a < 91; a++) {
        // 内层循环遍历波长
        for (int w = 400; w <= 700; w+=5) {
            // 构造文件名
            std::string filename = std::string(file)+"(WL=" + std::to_string(w) + "nm,AOI=" + std::to_string(-a) + "deg).txt";
            std::ifstream file(filename);

            // 检查文件是否成功打开
            if (!file.is_open()) {
                std::cerr << "无法打开文件: " << filename << std::endl;
                continue; // 跳过当前文件，继续下一个
            }

            double number;
            // 逐行读取文件内容
            int i = 0;
            while (file >> number) {
                data[a][i][(w-400)/5]=number;
                i++;
            }
            // 关闭文件
            file.close();
        }
    }

    for(int it = 0;it<91;it++){
        for(int ot = 0;ot<181;ot++){
            SampledSpectrum RN(0.);
            for(int i =0;i<60;++i){
                RN[i] = data[it][ot][i];
                RN[i] *= SampledSpectrum::get_rgbIllum2SpectWhite()[i];
            }
            Float rgb[3];
            RN.ToRGB(rgb);
            m_data[it][ot] = {rgb[0],rgb[1],rgb[2]};
        }
    }
}
hairSimBiBrdf::hairSimBiBrdf(const char* fileazi,const char* filelong) {
    m_data_a = std::vector<std::vector<RGB>>(91,std::vector<RGB>(181));
    m_data_l = std::vector<std::vector<RGB>>(91,std::vector<RGB>(181));
    ReadFromFileAndTrans(m_data_a,fileazi);
    ReadFromFileAndTrans(m_data_l,filelong);
    std::cout<<"gen hair brdf ok!!"<<fileazi<<"and longitu "<<filelong<<std::endl;
}
}  
