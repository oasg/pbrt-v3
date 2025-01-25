#include "hairsimbrdf.h"
#include <fstream>
#include "spectrum.h"
std::mutex SingBrdf::m_mutex;
std::shared_ptr<hairSimBrdf>  SingBrdf::m_instance_ptr = nullptr;

hairSimBrdf::hairSimBrdf(const char *file) {
    m_data = std::vector<std::vector<RGB>>(91,std::vector<RGB>(181));
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
    pbrt::SampledSpectrum::Init();
    for(int it = 0;it<91;it++){
        for(int ot = 0;ot<181;ot++){
            pbrt::SampledSpectrum RN(0.);
            for(int i =0;i<60;++i){
                RN[i] = data[it][ot][i];
                RN[i] *= pbrt::SampledSpectrum::get_rgbIllum2SpectWhite()[i];
            }
            Float rgb[3];
            RN.ToRGB(rgb);
            m_data[it][ot] = {rgb[0],rgb[1],rgb[2]};
        }
    }

    std::cout<<"gen hair brdf ok!!"<<file<<std::endl;
}

RGB hairSimBrdf::getReflect(Float it, Float ot)
{
    int x0 = floor(it);
    int x1 = ceil(it);

    int y0 = floor(ot);
    int y1 = ceil(ot);
    Float alpha_x = it - x0;
    Float alpha_y = ot - y0;


        // 获取四个邻近点的颜色值
    auto c00 = this->m_data[x0][y0];
    auto c10 = this->m_data[x1][y0];
    auto c01 = this->m_data[x0][y1];
    auto c11 = this->m_data[x1][y1];


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
    return RGB{r_final , g_final , b_final};
}
