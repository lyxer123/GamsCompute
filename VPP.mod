# ================== 集合定义 ==================
set T ordered;          # 时间周期集合（24小时）
set S;                  # 场景集合（s1,s2,s3）
set W := {1,2,3};       # 产消者集合（对应3个VPP）
set C := {'HVAC', 'WM', 'DW'};  # 可控设备集合
set F;                  # 设备运行阶段集合（f1,f2）

# ================== 参数定义 ==================
# ---------- ESS参数 ----------
param eta_ESS_charge{W};        # ESS充电效率
param ESS_charge_rate{W};       # ESS充电速率（kW）
param ESS_discharge_rate{W};    # ESS放电速率（kW）
param SOC_ESS_initial{W};       # ESS初始SOC（kWh）
param SOC_ESS_max{W};           # ESS最大SOC（kWh）
param SOC_ESS_min{W};           # ESS最小SOC（kWh）

# ---------- EV参数 ----------
param eta_EV_charge{W};         # EV充电效率
param eta_EV_discharge{W};      # EV放电效率
param EV_charge_rate{W};        # EV充电速率（kW）
param EV_discharge_rate{W};     # EV放电速率（kW）
param SOC_EV_initial{W};        # EV初始SOC（kWh）
param SOC_EV_max{W};            # EV最大SOC（kWh）
param SOC_EV_min{W};            # EV最小SOC（kWh）
param EV_arrival{W};            # EV到达时间（时段t）
param EV_departure{W};          # EV离开时间（时段t）

# ---------- 光伏及负载参数 ----------
param P_PV{W,T};                # 光伏可用功率（kW）
param InfLoad{W,T};             # 不可调负载（kW）

# ---------- 可控设备参数 ----------
param N{W,S,C};                 # 设备运行时段数
param P_consumption{W,C,F};     # 设备阶段功率（kW）
param phase_duration{W,C,F};    # 阶段持续时间（时段数）

# ---------- 市场参数 ----------
param energy_price_buy{T};       # 购电价格（€/MWh）
param energy_price_sell{T};      # 售电价格（€/MWh）
param deltaT;                   # 时间间隔（小时）

# ================== 新增参数声明 ==================
param T_ref;                 # 参考温度（需在.dat文件中定义）
param T_actual{W,T};         # 实际温度（需在.dat文件中定义）

# ================== 变量定义 ==================
var P_grid{W,T} >= 0;           # 从电网采购功率（kW）
var P_community{W,T} >= 0;      # 本地社区采购功率（kW）
var ESS_charge{W,T} >= 0;       # ESS充电功率（kW）
var ESS_discharge{W,T} >= 0;    # ESS放电功率（kW）
var EV_charge{W,T} >= 0;        # EV充电功率（kW）
var EV_discharge{W,T} >= 0;     # EV放电功率（kW）
var SOC_ESS{W,T} >= 0;          # ESS的SOC（kWh）
var SOC_EV{W,T} >= 0;           # EV的SOC（kWh）

# ---------- 销售功率分配变量 ----------
var P_sold_PV{W,T} >= 0;        # 光伏出售功率（kW）
var P_sold_EV{W,T} >= 0;        # EV出售功率（kW）
var P_sold_ESS{W,T} >= 0;       # ESS出售功率（kW）

# ---------- 二进制变量 ----------
var y_community{T} binary;              # 社区从电网取电标志
var u_ESS_charge{W,T} binary;           # ESS充电标志
var u_ESS_discharge{W,T} binary;        # ESS放电标志（新增）
var u_EV_charge{W,T} binary;            # EV充电标志
var u_EV_discharge{W,T} binary;         # EV放电标志（新增）
var x_phase_start{W,T,C,F} binary;      # 阶段开始标志
var x_phase_active{W,T,C,F} binary;     # 阶段进行中标志
var x_phase_end{W,T,C,F} binary;        # 阶段结束标志

# ================== 目标函数 ==================
minimize TotalCost:
  sum{w in W, t in T} (P_grid[w,t] * energy_price_buy[t]/1000 * deltaT)  # 购电成本
  - sum{w in W, t in T} (P_sold_ESS[w,t] + P_sold_EV[w,t] + P_sold_PV[w,t]) * energy_price_sell[t]/1000 * deltaT;  # 售电收入

# ================== 约束条件 ==================
# ---------- 能源平衡（公式5） ----------
subject to EnergyBalance{w in W, t in T, s in S}:
  P_grid[w,t] + P_community[w,t] + ESS_discharge[w,t] + EV_discharge[w,t] + P_PV[w,t]
  >= InfLoad[w,t] + sum{c in C, f in F} (x_phase_active[w,t,c,f] * P_consumption[w,c,f]) 
  + ESS_charge[w,t] + EV_charge[w,t];

# ---------- ESS动态约束（公式26-28） ----------
subject to ESS_SOC_Update{w in W, t in T}:
  SOC_ESS[w,t] = if t == 1 then SOC_ESS_initial[w]
    else SOC_ESS[w,t-1] + (ESS_charge[w,t] * eta_ESS_charge[w] - ESS_discharge[w,t]/eta_ESS_charge[w]) * deltaT;

subject to ESS_SOC_Limits{w in W, t in T}:
  SOC_ESS_min[w] <= SOC_ESS[w,t] <= SOC_ESS_max[w];

# ---------- EV动态约束（公式19-22） ----------
subject to EV_SOC_Update{w in W, t in T: t > 1}:
  SOC_EV[w,t] = if t == EV_arrival[w] then SOC_EV_initial[w]
    else if t > EV_arrival[w] and t <= EV_departure[w] then
      SOC_EV[w,t-1] + (EV_charge[w,t] * eta_EV_charge[w] - EV_discharge[w,t]/eta_EV_discharge[w]) * deltaT
    else SOC_EV[w,t-1];

subject to EV_SOC_Init{w in W}:
  SOC_EV[w,1] = SOC_EV_initial[w];

subject to EV_SOC_Limits{w in W, t in T}:
  if t >= EV_arrival[w] and t <= EV_departure[w] then
    SOC_EV_min[w] <= SOC_EV[w,t] <= SOC_EV_max[w];

# ---------- 设备阶段连续性约束（公式10-15） ----------
subject to PhaseStart{w in W, t in T, c in C, f in F: t > 1}:
   x_phase_start[w,t,c,f] <= 1 - x_phase_active[w,t-1,c,f];

subject to PhaseActivation{w in W, t in T, c in C, f in F}:
  x_phase_active[w,t,c,f] >= x_phase_start[w,t,c,f] - x_phase_end[w,t,c,f];

subject to PhaseDuration{w in W, t in T, c in C, f in F: t + phase_duration[w,c,f] -1 <= card(T)}:
  sum{tau in t..t+phase_duration[w,c,f]-1} x_phase_active[w,tau,c,f] >= phase_duration[w,c,f] * x_phase_start[w,t,c,f];

# ---------- ESS充放电互斥约束（修正后） ----------
subject to ESS_ChargeDischarge_Mutex{w in W, t in T}:
  u_ESS_charge[w,t] + u_ESS_discharge[w,t] <= 1;

subject to ESS_Charge_Limit{w in W, t in T}:
  ESS_charge[w,t] <= ESS_charge_rate[w] * u_ESS_charge[w,t];

subject to ESS_Discharge_Limit{w in W, t in T}:
  ESS_discharge[w,t] <= ESS_discharge_rate[w] * u_ESS_discharge[w,t];

# ---------- EV充放电互斥约束（修正后） ----------
subject to EV_ChargeDischarge_Mutex{w in W, t in T}:
  u_EV_charge[w,t] + u_EV_discharge[w,t] <= 1;

subject to EV_Charge_Limit{w in W, t in T}:
  EV_charge[w,t] <= EV_charge_rate[w] * u_EV_charge[w,t];

subject to EV_Discharge_Limit{w in W, t in T}:
  EV_discharge[w,t] <= EV_discharge_rate[w] * u_EV_discharge[w,t];

# ---------- 其他约束 ----------
subject to LocalEnergyBalance{t in T, s in S}:
  sum{w in W} P_community[w,t] = sum{w in W} (P_sold_ESS[w,t] + P_sold_EV[w,t] + P_sold_PV[w,t]);

subject to PVPowerAllocation{w in W, t in T}:
  P_PV[w,t] = P_community[w,t] + P_sold_PV[w,t];

# ================== 数据加载 ==================
data VPP.dat;    # 显式加载数据文件

# ================== 求解器配置 ==================
option solver cplex;    # 显式指定使用CPLEX
solve;
