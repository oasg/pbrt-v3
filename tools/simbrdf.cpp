#include "hairsimbrdf.h"
#include <iostream>
#include<opencv4/opencv2/opencv.hpp>
cv::Mat gammaCorrection(const cv::Mat& image, double gamma) {
    cv::Mat lut(1, 256, CV_8UC1); // 生成查找表
    for (int i = 0; i < 256; i++) {
        lut.at<uchar>(i) = cv::saturate_cast<uchar>(pow(i / 255.0, gamma) * 255.0);
    }

    cv::Mat correctedImage;
    cv::LUT(image, lut, correctedImage); // 应用查找表
    return correctedImage;
}

int main(){
    auto brdf = SingBrdf::get_Instance();
    cv::Mat image(900, 900, CV_8UC3);
    for(int i = 0;i<900;i++){
        for(int j = 0;j<900;j++){
            Float it = (i/900.0)*90.0;
            Float ot = ((j+900)/1800.0)*180.0;
            auto ref = brdf->getReflect(it,ot);
            image.at<cv::Vec3b>(i, j)[0] =  int(ref.b*255)%255;
            image.at<cv::Vec3b>(i, j)[1] = int(ref.g*255)%255;
            image.at<cv::Vec3b>(i, j)[2] = int(ref.r*255)%255;
        }
    }
        // 在图像上绘制 X 和 Y 轴
    int axis_length = 200;  // 设置坐标轴的长度
    cv::Point origin(0, 0);  // 坐标轴的原点位置，位于图像中心


    double gamma = 0.4;
    cv::Mat brightenedImage = gammaCorrection(image, gamma);
    try {
    if (!cv::imwrite("output.png", brightenedImage)) {
        std::cerr << "Failed to write image!" << std::endl;
    } else {
        std::cout << "Image saved successfully!" << std::endl;
    }
    } catch (const cv::Exception& ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
    }
    
}