# Motor Control Simulation: Command Smoothing Comparison

Python（control および numpy, matplotlib）を用いて、モータの多重制御ループ（速度ループ・電流ループ）における指令平滑化の有無による挙動の違いを検証・視覚化するためのシミュレーションスクリプトです。

## 概要 (Overview)
カスケード制御（外側：速度ループ、内側：電流ループ）を持つDCモータ制御系を想定し、指令平滑化あり/なしで同時シミュレーションして比較します。  
シミュレーションでは4次のルンゲ・クッタ法（RK4）を用いて連続時間系のモータプラント（状態方程式）を高精度に数値積分しています。

主な機能
---
制御周期の異なる速度ループ（例: $200\,\mu\text{s}$）と電流ループ（例: $50\,\mu\text{s}$）を、マイクロ秒オーダーの微小タイムステップ（ $5\,\mu\text{s}$）上で適切に分周・処理。  
実行すると、指令平滑化の有無による速度・電流・電圧の違いを比較したグラフウィンドウが立ち上がります。

## 実行方法 (Usage)

### 必要なライブラリ
```bash
pip install numpy matplotlib control
```
### 実行手順
```bash
python sample_command-smoothing.py
```
