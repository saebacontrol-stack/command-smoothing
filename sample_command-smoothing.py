import control as ct
import matplotlib.pyplot as plt
import numpy as np
import math

dt = 5e-6                   # タイムステップ 5μs
t_end = 0.1                 # シミュレーション時間 0.1s
steps = int(t_end/dt)       # ステップ数

Tsi = 50e-6         # 電流ループ制御周期[s]
Tsv = 200e-6        # 速度ループ制御周期[s]

Nsi = int(Tsi/dt)     # 電流ループ１ステップあたりの連続時間系ステップ数
Nsv = int(Tsv/dt)     # 速度ループ１ステップあたりの連続時間系ステップ数


t_hist = np.linspace(0, t_end, steps)   # 時間データ
x_hist1 = np.zeros((steps, 2))          # 状態変数データ１[steps x 2]
x_hist2 = np.zeros((steps, 2))          # 状態変数データ２[steps x 2]
V_hist1 = np.zeros((steps, 1))          # 電圧データ１
V_hist2 = np.zeros((steps, 1))          # 電圧データ２


Jm = 1e-4       # モータイナーシャ[kgm^2]
Rm = 0.2        # モータ巻線電気抵抗[Ω]
Lm = 2e-3       # モータ巻線インダクタンス[H]
Kt = 0.5        # モータトルク定数[Nm/A]
Ke = 0.5        # モータ逆起電力定数[V/(rad/s)]

InvKt = 1/Kt    # トルク定数の逆数[A/Nm]

wi = 2 * math.pi * 2000     # 電流制御応答周波数[rad/s]
wv = 2 * math.pi * 200      # 速度制御応答周波数[rad/s]

Kip = wi * Lm       # 電流比例ゲイン
Kii = wi * Rm       # 電流積分ゲイン

Kvp = Jm * wv               # 速度比例ゲイン
Kvi = 0.1 * Jm * wv ** 2    # 速度積分ゲイン

Ac = np.array(
        [
            [-Rm/Lm, -Ke/Lm],
            [Kt/Jm,  0.0,  ],
        ]
    )                               # 制御対象状態方程式A行列（連続時間系）
Bc = np.array([[1/Lm], [0.0]])      # 制御対象状態方程式B行列（連続時間系）


def sys_calc(x, u):
    """状態方程式 dx/dt = A@x + B@u 計算"""
    return Ac @ x + (Bc @ np.array([[u]])).flatten()


def solve_rk4(sys_calc, x, u, dt):
    """4次ルンゲ・クッタ法による1ステップ数値積分
    引数:
    -----------
    sys_calc : 状態方程式 dx/dt = A@x + B@u 計算用関数
    x : 現在の状態変数
    u : 現在の入力
    dt : 積分時間刻み幅[s]

    戻り値:
    --------
    x_next : dt 秒後の状態変数
    """
    k1 = sys_calc(x, u)
    k2 = sys_calc(x + 0.5 * dt * k1, u)
    k3 = sys_calc(x + 0.5 * dt * k2, u)
    k4 = sys_calc(x + dt * k3, u)

    x_next = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return x_next


def simulation():
  """シミュレーション実行関数"""
  global x_hist1, x_hist2
  
  x1 = np.array([0.0, 0.0])      # 初期状態の設定（状態変数データ１）
  x2 = np.array([0.0, 0.0])      # 初期状態の設定（状態変数データ２）

  Intg_i1 = 0       # 電流ループ内積分器（指令平滑化なし）
  Intg_i2 = 0       # 電流ループ内積分器（指令平滑化あり）

  Intg_v1 = 0       # 速度ループ内積分器（指令平滑化なし）
  Intg_v2 = 0       # 速度ループ内積分器（指令平滑化あり）

  cnt_i = Nsi       # カウント変数（電流ループ用）
  cnt_v = Nsv       # カウント変数（速度ループ用）

  vref = 0      # 速度指令（初期値０）

  Buf1r, Buf2r, Buf3r, Buf4r = 0, 0, 0, 0       # 指令平滑化（移動平均フィルタ）用バッファの初期化
  Buf1v, Buf2v, Buf3v, Buf4v = 0, 0, 0, 0       # 逆起電圧補償項平滑化（移動平均フィルタ）用バッファの初期化

  for k in range(steps):
    x_hist1[k] = x1
    x_hist2[k] = x2

    if cnt_v==Nsv:
       vref = 20*t_hist[k]      # 速度指令の生成
       if vref>1:
          vref = 1

       # 速度制御
       vfb1 = x1[1]     # 現在速度（指令平滑化なし）
       vfb2 = x2[1]     # 現在速度（指令平滑化あり）

       verr1 = vref - vfb1      # 速度偏差（指令平滑化なし）
       verr2 = vref - vfb2      # 速度偏差（指令平滑化あり）

       Intg_v1 += Tsv * verr1       # 速度偏差の積分計算（指令平滑化なし）
       Intg_v2 += Tsv * verr2       # 速度偏差の積分計算（指令平滑化あり）

       Tref1 = Kvp * verr1 + Kvi * Intg_v1  # PI制御によるトルク指令計算（指令平滑化なし）
       Tref2 = Kvp * verr2 + Kvi * Intg_v2  # PI制御によるトルク指令計算（指令平滑化あり）

       Iref1 = InvKt * Tref1    # トルク指令を電流指令に換算（指令平滑化なし）
       Iref2 = InvKt * Tref2    # トルク指令を電流指令に換算（指令平滑化あり）

       cnt_v = 0        # カウント変数の初期化

    if cnt_i==Nsi:
       # 電流制御
       Buf4r = Buf3r    # 指令平滑化用バッファ更新
       Buf3r = Buf2r    # 指令平滑化用バッファ更新
       Buf2r = Buf1r    # 指令平滑化用バッファ更新
       Buf1r = Iref2    # 指令平滑化用バッファ更新

       Buf4v = Buf3v    # 逆起電圧補償項平滑化用バッファ更新
       Buf3v = Buf2v    # 逆起電圧補償項平滑化用バッファ更新
       Buf2v = Buf1v    # 逆起電圧補償項平滑化用バッファ更新
       Buf1v = vfb2     # 逆起電圧補償項平滑化用バッファ更新

       Iref2_smth = 0.25 * (Buf1r + Buf2r + Buf3r + Buf4r)      # 移動平均フィルタの応用による指令平滑化
       vfb2_smth = 0.25 * (Buf1v + Buf2v + Buf3v + Buf4v)       # 移動平均フィルタの応用による逆起電圧補償用速度データ平滑化

       Ierr1 = Iref1 - x1[0]            # 電流偏差（指令平滑化なし）
       Ierr2 = Iref2_smth - x2[0]       # 電流偏差（指令平滑化あり）

       Intg_i1 += Tsi * Ierr1       # 電流偏差の積分計算（指令平滑化なし）
       Intg_i2 += Tsi * Ierr2       # 電流偏差の積分計算（指令平滑化あり）

       u1 = Kip * Ierr1 + Kii * Intg_i1 + Ke * vfb1             # PI制御＋逆起電圧補償による電圧計算（指令平滑化なし）
       u2 = Kip * Ierr2 + Kii * Intg_i2 + Ke * vfb2_smth        # PI制御＋逆起電圧補償による電圧計算（指令平滑化あり）

       cnt_i = 0        # カウント変数の初期化

    x1 = solve_rk4(sys_calc, x1, u1, dt)    # 1ステップ更新（指令平滑化なし）
    x2 = solve_rk4(sys_calc, x2, u2, dt)    # 1ステップ更新（指令平滑化あり）

    V_hist1[k] = u1
    V_hist2[k] = u2

    cnt_i += 1  # カウントアップ
    cnt_v += 1  # カウントアップ


# シミュレーション実行
simulation()

# グラフ描画（速度・電流・電圧をサブプロットで比較）
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

# 1. 速度プロット
axes[0].plot(t_hist, x_hist1[:, 1], label="without smoothing", color="tab:blue")
axes[0].plot(t_hist, x_hist2[:, 1], label="with smoothing", color="tab:orange")
axes[0].set_ylabel("Velocity [rad/s]")
axes[0].grid(True)
axes[0].legend()
axes[0].set_title("Comparison with/without Command Smoothing")

# 2. 電流プロット (x[0]が電流)
axes[1].plot(t_hist, x_hist1[:, 0], label="without smoothing", color="tab:blue")
axes[1].plot(t_hist, x_hist2[:, 0], label="with smoothing", color="tab:orange")
axes[1].set_ylabel("Current [A]")
axes[1].grid(True)
axes[1].legend()

# 3. 電圧プロット (u)
axes[2].plot(t_hist, V_hist1, label="without smoothing", color="tab:blue")
axes[2].plot(t_hist, V_hist2, label="with smoothing", color="tab:orange")
axes[2].set_ylabel("Voltage [V]")
axes[2].set_xlabel("Time [s]")
axes[2].grid(True)
axes[2].legend()

plt.tight_layout()
plt.show()