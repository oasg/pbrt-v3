#!/bin/bash

# 定义要执行的模型数组
models=(
  #"HairMultilayerPerlinModel"
  #"HairPretteMultilayerPerlinModel"
  #"hairmodelDamagedLargeDisPerlin"
  #"HairMultilayerDamagedLackLayerPerlinModel"
  "HairMultilayerDamagedTiltPerlinModel"
  #"hairmodelRepairedPerlin"
)

# 定义一个函数来执行单个任务
run_model() {
  local model="$1"
  echo "Starting model: $model"
  ../build/pbrt  ../scenes/hair/hair/straight-hair-${model}.pbrt
  echo "Finished model: $model"
}

# 设置信号处理函数
trap cleanup CHLD

# 循环遍历模型数组
for model in "${models[@]}"; do
  # 执行任务
  run_model "$model" 
done

echo "complete"