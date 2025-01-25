import cv2
import matplotlib.pyplot as plt
import os
def addaxis(filepath):
    # 读取现有图像
    image = cv2.imread(filepath)

    # 将图像从 BGR 转为 RGB（因为 OpenCV 默认使用 BGR）
    image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    height, width, _ = image.shape

    # 创建一个图形和坐标轴
    fig, ax = plt.subplots()
    ax.imshow(image_rgb)

    # 设置坐标轴范围，将最大值设置为图像的宽度和高度
    ax.set_xlim([0, width])
    ax.set_ylim([height, 0])  # Y轴从上到下反转

    # 设置坐标轴的刻度位置
    # 将坐标轴的最大值 90 映射到图像的边界位置
    xticks = [0, width // 4, width // 2, 3 * width // 4, width]
    yticks = [0, height // 4, height // 2, 3 * height // 4, height]
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.invert_yaxis()

    # 设置坐标轴的刻度标签
    ax.set_xticklabels([str(i * 90 // width) for i in xticks])  # 映射到 90 度
    ax.set_yticklabels([str(i * 90 // height) for i in yticks])  # 映射到 90 度

    # 添加 X 轴和 Y 轴标签
    ax.set_xlabel('Angle of incidence')
    ax.set_ylabel('Reflection angle')
    # 显示图像及坐标轴
    name = os.path.basename(filepath)
    plt.savefig(os.path.splitext(name)[0]+".png")
addaxis('../build/output_norm.png')
addaxis('../build/output_norm_rand.png')
addaxis('../build/output_damage.png')
addaxis('../build/output_damage_rand.png')
addaxis('../build/output_repaired.png')
addaxis('../build/output_repaired_normal.png')