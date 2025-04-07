set T ordered;          
set S;                  
set W := {1,2,3};       
set C := {'HVAC', 'WM', 'DW'};  
set F;             

param eta_ESS_charge{W};        
param ESS_charge_rate{W};       
param ESS_discharge_rate{W};    
param SOC_ESS_initial{W};       
param SOC_ESS_max{W};           
param SOC_ESS_min{W};           

param eta_EV_charge{W};         
param eta_EV_discharge{W};      
param EV_charge_rate{W};        
param EV_discharge_rate{W};     
param SOC_EV_initial{W};        
param SOC_EV_max{W};            
param SOC_EV_min{W};            
param EV_arrival{W};            
param EV_departure{W};          

param P_PV{W,T};                
param InfLoad{W,T};             

param N{W,S,C};                 
param P_consumption{W,C,F};     
param phase_duration{W,C,F};    

param energy_price_buy{T};       
param energy_price_sell{T};      
param deltaT;                   

param T_ref;                 
param T_actual{W,T};         

var P_grid{W,T} >= 0;           
var P_community{W,T} >= 0;      
var ESS_charge{W,T} >= 0;       
var ESS_discharge{W,T} >= 0;    
var EV_charge{W,T} >= 0;        
var EV_discharge{W,T} >= 0;     
var SOC_ESS{W,T} >= 0;          
var SOC_EV{W,T} >= 0;           

var P_sold_PV{W,T} >= 0;        
var P_sold_EV{W,T} >= 0;        
var P_sold_ESS{W,T} >= 0;       

var y_community{T} binary;              
var u_ESS_charge{W,T} binary;           
var u_ESS_discharge{W,T} binary;        
var u_EV_charge{W,T} binary;            
var u_EV_discharge{W,T} binary;         
var x_phase_start{W,T,C,F} binary;      
var x_phase_active{W,T,C,F} binary;     
var x_phase_end{W,T,C,F} binary;     

# 输出变量
var Total_DW{W,T} >= 0;      
var Total_WM{W,T} >= 0;      
var Total_HVAC{W,T} >= 0;    
var Total_InfLoad{W,T} >= 0; 
var Total_ESS_charge{W,T} >= 0;
var Total_EV_charge{W,T} >= 0;
var Total_PV{W,T} >= 0;      
var Total_ESS_discharge{W,T} >= 0;
var Total_EV_discharge{W,T} >= 0;
var Total_Grid{W,T} >= 0;    

minimize TotalCost:
  sum{w in W, t in T} (P_grid[w,t] * energy_price_buy[t]/1000 * deltaT)
  - sum{w in W, t in T} (P_sold_ESS[w,t] + P_sold_EV[w,t] + P_sold_PV[w,t]) * energy_price_sell[t]/1000 * deltaT;

# 修改后的能量平衡约束（考虑所有产消者）
subject to EnergyBalance{w in W, t in T}:
  P_grid[w,t] + P_community[w,t] + ESS_discharge[w,t] + EV_discharge[w,t] + P_PV[w,t]
  = InfLoad[w,t] 
  + sum{c in C, f in F} (sum{s in S} x_phase_active[w,t,c,f] * P_consumption[w,c,f] * N[w,s,c])
  + ESS_charge[w,t] + EV_charge[w,t];

subject to ESS_SOC_Init{w in W}:
  SOC_ESS[w,1] = SOC_ESS_initial[w];

subject to ESS_SOC_Update{w in W, t in T: t > 1}:
  SOC_ESS[w,t] = SOC_ESS[w,t-1] + (ESS_charge[w,t] * eta_ESS_charge[w] - ESS_discharge[w,t]/eta_ESS_charge[w]) * deltaT;


# EV SOC更新（修正时间逻辑）
subject to EV_SOC_Init{w in W}:
  SOC_EV[w,EV_arrival[w]] = SOC_EV_initial[w];

subject to EV_SOC_Update{w in W, t in T: t > 1 and t > EV_arrival[w] and t <= EV_departure[w]}:
  SOC_EV[w,t] = SOC_EV[w,t-1] + (EV_charge[w,t] * eta_EV_charge[w] - EV_discharge[w,t]/eta_EV_discharge[w]) * deltaT;

subject to EV_SOC_OutsideRange{w in W, t in T: t < EV_arrival[w] or t > EV_departure[w]}:
  SOC_EV[w,t] = 0;


# 设备运行约束（修正相位激活逻辑）
subject to PhaseStart{w in W, t in T, c in C, f in F: t > 1}:
   x_phase_start[w,t,c,f] <= 1 - sum{tau in max(1,t-phase_duration[w,c,f])..t-1, ff in F} x_phase_active[w,tau,c,ff];

subject to PhaseActivation{w in W, t in T, c in C, f in F}:
  x_phase_active[w,t,c,f] >= x_phase_start[w,t,c,f] - x_phase_end[w,t,c,f];

subject to PhaseDuration{w in W, t in T, c in C, f in F: t + phase_duration[w,c,f] -1 <= card(T)}:
  sum{tau in t..t+phase_duration[w,c,f]-1} x_phase_active[w,tau,c,f] >= phase_duration[w,c,f] * x_phase_start[w,t,c,f];

# ESS充放电互斥约束
subject to ESS_ChargeDischarge_Mutex{w in W, t in T}:
  u_ESS_charge[w,t] + u_ESS_discharge[w,t] <= 1;

subject to ESS_Charge_Limit{w in W, t in T}:
  ESS_charge[w,t] <= ESS_charge_rate[w] * u_ESS_charge[w,t];

subject to ESS_Discharge_Limit{w in W, t in T}:
  ESS_discharge[w,t] <= ESS_discharge_rate[w] * u_ESS_discharge[w,t];

# EV充放电约束（修正可用时间）
subject to EV_ChargeDischarge_Mutex{w in W, t in T}:
  u_EV_charge[w,t] + u_EV_discharge[w,t] <= 1;

subject to EV_Charge_Limit{w in W, t in T}:
  EV_charge[w,t] <= if (t >= EV_arrival[w] and t <= EV_departure[w]) then EV_charge_rate[w] * u_EV_charge[w,t] else 0;

subject to EV_Discharge_Limit{w in W, t in T}:
  EV_discharge[w,t] <= if (t >= EV_arrival[w] and t <= EV_departure[w]) then EV_discharge_rate[w] * u_EV_discharge[w,t] else 0;

# 社区能量平衡
subject to LocalEnergyBalance{t in T}:
  sum{w in W} P_community[w,t] = sum{w in W} (P_sold_ESS[w,t] + P_sold_EV[w,t] + P_sold_PV[w,t]);

subject to PVPowerAllocation{w in W, t in T}:
  P_PV[w,t] = P_community[w,t] + P_sold_PV[w,t];
  
# 光伏时间验证（修正约束表达式）
subject to PV_Time_Validation{w in W, t in T: t < 8 or t > 17}:
  P_PV[w,t] = 0;
  
# 计算各设备总功率（修正计算方式）
subject to Calculate_DW{w in W, t in T}:
    Total_DW[w,t] = sum{s in S, f in F} (x_phase_active[w,t,'DW',f] * P_consumption[w,'DW',f] * N[w,s,'DW']);

subject to Calculate_WM{w in W, t in T}:
    Total_WM[w,t] = sum{s in S, f in F} (x_phase_active[w,t,'WM',f] * P_consumption[w,'WM',f] * N[w,s,'WM']);

subject to Calculate_HVAC{w in W, t in T}:
    Total_HVAC[w,t] = sum{s in S, f in F} (x_phase_active[w,t,'HVAC',f] * P_consumption[w,'HVAC',f] * N[w,s,'HVAC']);

subject to Calculate_InfLoad{w in W, t in T}:
    Total_InfLoad[w,t] = InfLoad[w,t];

subject to Calculate_ESS_charge{w in W, t in T}:
    Total_ESS_charge[w,t] = ESS_charge[w,t];

subject to Calculate_EV_charge{w in W, t in T}:
    Total_EV_charge[w,t] = EV_charge[w,t];

subject to Calculate_PV{w in W, t in T}:
    Total_PV[w,t] = P_PV[w,t];

subject to Calculate_ESS_discharge{w in W, t in T}:
    Total_ESS_discharge[w,t] = ESS_discharge[w,t];

subject to Calculate_EV_discharge{w in W, t in T}:
    Total_EV_discharge[w,t] = EV_discharge[w,t];

subject to Calculate_Grid{w in W, t in T}:
    Total_Grid[w,t] = P_grid[w,t];