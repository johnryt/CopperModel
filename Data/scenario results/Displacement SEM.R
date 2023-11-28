# Import Data and require packages
library(plm)
library(stargazer)
library(lmtest)
library(lavaan)
library(readxl)

# Set path, read data
wd <- ("C:/Users/ryter/Dropbox (MIT)/Group Research Folder_Olivetti/Displacement/00 Simulation/06 Module Integration/Data/scenario results")
setwd(wd)
df_all <- read_excel("SEM_df_all.xlsx")
df_high <- read_excel("SEM_df_high_opp.xlsx")
df_low <- read_excel("SEM_df_low_opp.xlsx")
panel_high <- pdata.frame(df_high,index='number', drop.index=TRUE, row.names=TRUE)

model <- '
# Between components - random intercepts
  RI_di =~ 1*di23 + 1*di22 + 1*di21 + 1*di20
  RI_mi =~ 1*mi23 + 1*mi22 + 1*mi21 + 1*mi20
  RI_de =~ 1*de23 + 1*de22 + 1*de21 + 1*de20
  RI_ss =~ 1*ss23 + 1*ss22 + 1*ss21 + 1*ss20
  RI_tc =~ 1*tc23 + 1*tc22 + 1*tc21 + 1*tc20
  RI_ca =~ 1*ca23 + 1*ca22 + 1*ca21 + 1*ca20
  RI_sh =~ 1*sh23 + 1*sh22 + 1*sh21 + 1*sh20

  RI_di ~~ RI_di + RI_mi + RI_de + RI_ss + RI_tc + RI_ca + RI_sh
  RI_mi ~~ RI_mi + RI_de + RI_ss + RI_tc + RI_ca + RI_sh
  RI_de ~~ RI_de + RI_ss + RI_tc + RI_ca + RI_sh
  RI_ss ~~ RI_ss + RI_tc + RI_ca + RI_sh
  RI_tc ~~ RI_tc + RI_ca + RI_sh
  RI_ca ~~ RI_ca + RI_sh
    
# Within-unit centered variables
  w_di23 =~ 1*di23
  w_di22 =~ 1*di22
  w_di21 =~ 1*di21
  w_di20 =~ 1*di20
  
  w_mi23 =~ 1*mi23
  w_mi22 =~ 1*mi22
  w_mi21 =~ 1*mi21
  w_mi20 =~ 1*mi20
  
  w_de23 =~ 1*de23
  w_de22 =~ 1*de22
  w_de21 =~ 1*de21
  w_de20 =~ 1*de20
  
  w_ss23 =~ 1*ss23
  w_ss22 =~ 1*ss22
  w_ss21 =~ 1*ss21
  w_ss20 =~ 1*ss20
  
  w_tc23 =~ 1*tc23
  w_tc22 =~ 1*tc22
  w_tc21 =~ 1*tc21
  w_tc20 =~ 1*tc20
  
  w_ca23 =~ 1*ca23
  w_ca22 =~ 1*ca22
  w_ca21 =~ 1*ca21
  w_ca20 =~ 1*ca20
  
  w_sh23 =~ 1*sh23
  w_sh22 =~ 1*sh22
  w_sh21 =~ 1*sh21
  w_sh20 =~ 1*sh20

# Lagged effects between within-unit centered variables
  w_di23 ~ a*w_di22 + b*w_mi22 + c*w_de22 + d*w_ss22 + e*w_tc22 + f*w_ca22 + g*w_sh22
  w_di22 ~ a*w_di21 + b*w_mi21 + c*w_de21 + d*w_ss21 + e*w_tc21 + f*w_ca21 + g*w_sh21
  w_di21 ~ a*w_di20 + b*w_mi20 + c*w_de20 + d*w_ss20 + e*w_tc20 + f*w_ca20 + g*w_sh20
  
  w_mi23 ~ h*w_mi22 + i*w_de22 + j*w_ss22 + k*w_tc22 + l*w_ca22 + m*w_sh22
  w_mi22 ~ h*w_mi21 + i*w_de21 + j*w_ss21 + k*w_tc21 + l*w_ca21 + m*w_sh21
  w_mi21 ~ h*w_mi20 + i*w_de20 + j*w_ss20 + k*w_tc20 + l*w_ca20 + m*w_sh20
  
  w_de23 ~ n*w_mi22 + o*w_de22 + p*w_ss22 + q*w_tc22 + r*w_ca22 + s*w_sh22
  w_de22 ~ n*w_mi21 + o*w_de21 + p*w_ss21 + q*w_tc21 + r*w_ca21 + s*w_sh21
  w_de21 ~ n*w_mi20 + o*w_de20 + p*w_ss20 + q*w_tc20 + r*w_ca20 + s*w_sh20
  
  w_ss23 ~ t*w_mi22 + u*w_de22 + v*w_ss22 + w*w_tc22 + x*w_ca22 + y*w_sh22
  w_ss22 ~ t*w_mi21 + u*w_de21 + v*w_ss21 + w*w_tc21 + x*w_ca21 + y*w_sh21
  w_ss21 ~ t*w_mi20 + u*w_de20 + v*w_ss20 + w*w_tc20 + x*w_ca20 + y*w_sh20
  
  w_tc23 ~ z*w_mi22 + aa*w_de22 + ab*w_ss22 + ac*w_tc22 + ad*w_ca22 + ae*w_sh22
  w_tc22 ~ z*w_mi21 + aa*w_de21 + ab*w_ss21 + ac*w_tc21 + ad*w_ca21 + ae*w_sh21
  w_tc21 ~ z*w_mi20 + aa*w_de20 + ab*w_ss20 + ac*w_tc20 + ad*w_ca20 + ae*w_sh20
  
  w_ca23 ~ af*w_mi22 + ag*w_de22 + ah*w_ss22 + ai*w_tc22 + aj*w_ca22 + al*w_sh22
  w_ca22 ~ af*w_mi21 + ag*w_de21 + ah*w_ss21 + ai*w_tc21 + aj*w_ca21 + al*w_sh21
  w_ca21 ~ af*w_mi20 + ag*w_de20 + ah*w_ss20 + ai*w_tc20 + aj*w_ca20 + al*w_sh20
  
  w_sh23 ~ am*w_mi22 + an*w_de22 + ao*w_ss22 + ap*w_tc22 + aq*w_ca22 + ar*w_sh22
  w_sh22 ~ am*w_mi21 + an*w_de21 + ao*w_ss21 + ap*w_tc21 + aq*w_ca21 + ar*w_sh21
  w_sh21 ~ am*w_mi20 + an*w_de20 + ao*w_ss20 + ap*w_tc20 + aq*w_ca20 + ar*w_sh20

# Covariance of within-unit centered variables - wave 1
  w_di20 ~~ w_di20 + w_mi20 + w_de20 + w_ss20 + w_tc20 + w_ca20 + w_sh20
  w_mi20 ~~ w_mi20 + w_de20 + w_ss20 + w_tc20 + w_ca20 + w_sh20
  w_de20 ~~ w_de20 + w_ss20 + w_tc20 + w_ca20 + w_sh20
  w_ss20 ~~ w_ss20 + w_tc20 + w_ca20 + w_sh20
  w_tc20 ~~ w_tc20 + w_ca20 + w_sh20
  w_ca20 ~~ w_ca20 + w_sh20

# Time-invariant covariance of within-unit centered variables - following waves
  w_di21 ~~ cov1*w_mi21 + cov2*w_de21 + cov3*w_ss21 + cov4*w_tc21 + cov5*w_ca21 + cov6*w_sh21
  w_mi21 ~~ cov7*w_de21 + cov8*w_ss21 + cov9*w_tc21 + cova*w_ca21 + covb*w_sh21
  w_de21 ~~ covc*w_ss21 + covd*w_tc21 + cove*w_ca21 + covf*w_sh21
  w_ss21 ~~ covg*w_tc21 + covh*w_ca21 + covi*w_sh21
  w_tc21 ~~ covj*w_ca21 + covk*w_sh21
  w_ca21 ~~ covl*w_sh21

  w_di22 ~~ cov1*w_mi22 + cov2*w_de22 + cov3*w_ss22 + cov4*w_tc22 + cov5*w_ca22 + cov6*w_sh22
  w_mi22 ~~ cov7*w_de22 + cov8*w_ss22 + cov9*w_tc22 + cova*w_ca22 + covb*w_sh22
  w_de22 ~~ covc*w_ss22 + covd*w_tc22 + cove*w_ca22 + covf*w_sh22
  w_ss22 ~~ covg*w_tc22 + covh*w_ca22 + covi*w_sh22
  w_tc22 ~~ covj*w_ca22 + covk*w_sh22
  w_ca22 ~~ covl*w_sh22
    
  w_di23 ~~ cov1*w_mi23 + cov2*w_de23 + cov3*w_ss23 + cov4*w_tc23 + cov5*w_ca23 + cov6*w_sh23
  w_mi23 ~~ cov7*w_de23 + cov8*w_ss23 + cov9*w_tc23 + cova*w_ca23 + covb*w_sh23
  w_de23 ~~ covc*w_ss23 + covd*w_tc23 + cove*w_ca23 + covf*w_sh23
  w_ss23 ~~ covg*w_tc23 + covh*w_ca23 + covi*w_sh23
  w_tc23 ~~ covj*w_ca23 + covk*w_sh23
  w_ca23 ~~ covl*w_sh23

# Time-invariant residual variance of within-unit centered variables
  w_di21 ~~ v_di*w_di21
  w_mi21 ~~ v_mi*w_mi21
  w_de21 ~~ v_de*w_de21
  w_ss21 ~~ v_ss*w_ss21
  w_tc21 ~~ v_tc*w_tc21
  w_ca21 ~~ v_ca*w_ca21
  w_sh21 ~~ v_sh*w_sh21

  w_di22 ~~ v_di*w_di22
  w_mi22 ~~ v_mi*w_mi22
  w_de22 ~~ v_de*w_de22
  w_ss22 ~~ v_ss*w_ss22
  w_tc22 ~~ v_tc*w_tc22
  w_ca22 ~~ v_ca*w_ca22
  w_sh22 ~~ v_sh*w_sh22

  w_di23 ~~ v_di*w_di23
  w_mi23 ~~ v_mi*w_mi23
  w_de23 ~~ v_de*w_de23
  w_ss23 ~~ v_ss*w_ss23
  w_tc23 ~~ v_tc*w_tc23
  w_ca23 ~~ v_ca*w_ca23
  w_sh23 ~~ v_sh*w_sh23

# Regression of observed variables on shock parameters
  di20 + di21 + di22 + di23 ~ ny_di*ny
  mi20 + mi21 + mi22 + mi23 ~ ny_mi*ny
  de20 + de21 + de22 + de23 ~ ny_de*ny
  ss20 + ss21 + ss22 + ss23 ~ ny_ss*ny
  tc20 + tc21 + tc22 + tc23 ~ ny_tc*ny
  ca20 + ca21 + ca22 + ca23 ~ ny_ca*ny
  sh20 + sh21 + sh22 + sh23 ~ ny_sh*ny
  
  # di20 + di21 + di22 + di23 ~ ns_di*ns
  # mi20 + mi21 + mi22 + mi23 ~ ns_mi*ns
  # de20 + de21 + de22 + de23 ~ ns_de*ns
  # ss20 + ss21 + ss22 + ss23 ~ ns_ss*ns
  # tc20 + tc21 + tc22 + tc23 ~ ns_tc*ns
  # ca20 + ca21 + ca22 + ca23 ~ ns_ca*ns
  # sh20 + sh21 + sh22 + sh23 ~ ns_sh*ns
  
  di20 + di21 + di22 + di23 ~ op_di*op
  mi20 + mi21 + mi22 + mi23 ~ op_mi*op
  de20 + de21 + de22 + de23 ~ op_de*op
  ss20 + ss21 + ss22 + ss23 ~ op_ss*op
  tc20 + tc21 + tc22 + tc23 ~ op_tc*op
  ca20 + ca21 + ca22 + ca23 ~ op_ca*op
  sh20 + sh21 + sh22 + sh23 ~ op_sh*op
  
  # di20 + di21 + di22 + di23 ~ os_di*os
  # mi20 + mi21 + mi22 + mi23 ~ os_mi*os
  # de20 + de21 + de22 + de23 ~ os_de*os
  # ss20 + ss21 + ss22 + ss23 ~ os_ss*os
  # tc20 + tc21 + tc22 + tc23 ~ os_tc*os
  # ca20 + ca21 + ca22 + ca23 ~ os_ca*os
  # sh20 + sh21 + sh22 + sh23 ~ os_sh*os  
'

model_small <- '
  RI_sh =~ 1*sh20 + 1*sh21 + 1*sh22 + 1*sh23
  RI_di =~ 1*di23 + 1*di22 + 1*di21 + 1*di20
  
  RI_sh ~~ RI_sh
  RI_sh ~~ RI_di
  RI_di ~~ RI_di

  w_di20 =~ 1*di20
  w_di21 =~ 1*di21
  w_di22 =~ 1*di22
  w_di23 =~ 1*di23
  w_sh20 =~ 1*sh20
  w_sh21 =~ 1*sh21
  w_sh22 =~ 1*sh22
  w_sh23 =~ 1*sh23

  w_di21 ~ w_di20 + w_sh20 # a* and b*
  w_di22 ~ w_di21 + w_sh21
  w_di23 ~ w_di22 + w_sh22

  w_sh21 ~ w_di20 + w_sh20 # c* and d*
  w_sh22 ~ w_di21 + w_sh21
  w_sh23 ~ w_di22 + w_sh22

  # w_sh20 ~~ w_di20
  # w_sh21 ~~ cov*w_di21
  # w_sh22 ~~ cov*w_di22
  # w_sh23 ~~ cov*w_di23

  # w_sh20 ~~ w_sh20
  # w_sh21 ~~ v_sh*w_sh21
  # w_sh22 ~~ v_sh*w_sh22
  # w_sh23 ~~ v_sh*w_sh23

  # w_di20 ~~ w_di20
  # w_di21 ~~ v_di*w_di21
  # w_di22 ~~ v_di*w_di22
  # w_di23 ~~ v_di*w_di23

  di20 + di21 + di22 + di23 ~ s1*ny
  #di20 + di21 + di22 + di23 ~ s2*op
  #sh20 + sh21 + sh22 + sh23 ~ s3*op
  sh20 + sh21 + sh22 + sh23 ~ s4*ny
  
'

another_try <- '
  sh2 ~ sh_op*op + sh_ny*ny
  sh3 ~ sh_op*op + sh_ny*ny + sh_sh*sh2
  sh4 ~ sh_op*op + sh_ny*ny + sh_sh*sh3
  sh5 ~ sh_op*op + sh_ny*ny + sh_sh*sh4
  sh6 ~ sh_op*op + sh_ny*ny + sh_sh*sh5
  sh7 ~ sh_op*op + sh_ny*ny + sh_sh*sh6
  sh8 ~ sh_op*op + sh_ny*ny + sh_sh*sh7
  sh9 ~ sh_op*op + sh_ny*ny + sh_sh*sh8
  sh10 ~ sh_op*op + sh_ny*ny + sh_sh*sh9

  mi5 ~ mi_sh*sh4 + mi_mi*mi4 
  mi6 ~ mi_sh*sh5 + mi_mi*mi5 
  mi7 ~ mi_sh*sh6 + mi_mi*mi6 + mi_de*de6
  mi8 ~ mi_sh*sh7 + mi_mi*mi7 + mi_de*de7
  mi9 ~ mi_sh*sh8 + mi_mi*mi8 + mi_de*de8
  mi10 ~ mi_sh*sh9 + mi_mi*mi9 + mi_de*de9

  ca

  de6 ~ de_mi*mi5 + de_sh*sh5
  de7 ~ de_mi*mi6 + de_sh*sh6 + de_de*de6
  de8 ~ de_mi*mi7 + de_sh*sh7 + de_de*de7
  de9 ~ de_mi*mi8 + de_sh*sh8 + de_de*de8
  de10 ~ de_mi*mi9 + de_sh*sh9 + de_de*de9

  eff_on_shock := sh_op
  eff_on_mining := sh_op*mi_sh + sh_op*mi_sh*mi_mi # unsure about whether to include the second term here
  eff_on_disp := -eff_on_mining/eff_on_shock

'

another_try <- '
 # sh2 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns
 # sh3 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh2
 # sh4 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh3
 # sh5 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh4
 # sh6 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh5
 # sh7 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh6
 # sh8 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh7
 # sh9 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh8
 # sh10 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh9
 # sh11 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh10
 # sh12 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh11
  sh13 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh12
  sh14 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh13
  sh15 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh14
  sh16 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh15
  sh17 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh16
  sh18 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh17
  sh19 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh18
  sh20 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh19
  sh21 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh20
  sh22 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh21
  sh23 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh22

  di23 ~ di_mi*mi23 + di_sh*sh23 + di_de*de23
  di22 ~ di_mi*mi22 + di_sh*sh22 + di_de*de22
  di21 ~ di_mi*mi21 + di_sh*sh21 + di_de*de21
  di20 ~ di_mi*mi20 + di_sh*sh20 + di_de*de20
  di19 ~ di_mi*mi19 + di_sh*sh19 + di_de*de19
  di18 ~ di_mi*mi18 + di_sh*sh18 + di_de*de18
  di17 ~ di_mi*mi17 + di_sh*sh17 + di_de*de17
  di16 ~ di_mi*mi16 + di_sh*sh16 + di_de*de16
  di15 ~ di_mi*mi15 + di_sh*sh15 + di_de*de15
  di14 ~ di_mi*mi14 + di_sh*sh14 + di_de*de14
 # di13 ~ di_mi*mi13 + di_sh*sh13 + di_de*de13
 # di12 ~ di_mi*mi12 + di_sh*sh12 + di_de*de12
 # di11 ~ di_mi*mi11 + di_sh*sh11 + di_de*de11
 # di10 ~ di_mi*mi10 + di_sh*sh10 + di_de*de10
 # di9 ~ di_mi*mi9 + di_sh*sh9 + di_de*de9
 # di8 ~ di_mi*mi8 + di_sh*sh8 + di_de*de8
 # di7 ~ di_mi*mi7 + di_sh*sh7 + di_de*de7
  
  
 # mi5 ~ mi_sh*sh4 + mi_mi*mi4 + mi_tc*tc4             + mi_op*op + mi_ny*ny
 # mi6 ~ mi_sh*sh5 + mi_mi*mi5 + mi_tc*tc5             + mi_op*op + mi_ny*ny
 # mi7 ~ mi_sh*sh6 + mi_mi*mi6 + mi_tc*tc6 + mi_de*de6 + mi_op*op + mi_ny*ny
 # mi8 ~ mi_sh*sh7 + mi_mi*mi7 + mi_tc*tc7 + mi_de*de7 + mi_op*op + mi_ny*ny
 # mi9 ~ mi_sh*sh8 + mi_mi*mi8 + mi_tc*tc8 + mi_de*de8 + mi_op*op + mi_ny*ny
 # mi10~ mi_sh*sh9 + mi_mi*mi9 + mi_tc*tc9 + mi_de*de9 + mi_op*op + mi_ny*ny
 # mi11~ mi_sh*sh10 + mi_mi*mi10 + mi_tc*tc10 + mi_de*de10 + mi_op*op + mi_ny*ny
 # mi12~ mi_sh*sh11 + mi_mi*mi11 + mi_tc*tc11 + mi_de*de11 + mi_op*op + mi_ny*ny
  mi13~ mi_sh*sh12 + mi_mi*mi12 + mi_tc*tc12 + mi_de*de12 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi14~ mi_sh*sh13 + mi_mi*mi13 + mi_tc*tc13 + mi_de*de13 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi15~ mi_sh*sh14 + mi_mi*mi14 + mi_tc*tc14 + mi_de*de14 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi16~ mi_sh*sh15 + mi_mi*mi15 + mi_tc*tc15 + mi_de*de15 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi17~ mi_sh*sh16 + mi_mi*mi16 + mi_tc*tc16 + mi_de*de16 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi18~ mi_sh*sh17 + mi_mi*mi17 + mi_tc*tc17 + mi_de*de17 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi19~ mi_sh*sh18 + mi_mi*mi18 + mi_tc*tc18 + mi_de*de18 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi20~ mi_sh*sh19 + mi_mi*mi19 + mi_tc*tc19 + mi_de*de19 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi21~ mi_sh*sh20 + mi_mi*mi20 + mi_tc*tc20 + mi_de*de20 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi22~ mi_sh*sh21 + mi_mi*mi21 + mi_tc*tc21 + mi_de*de21 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  mi23~ mi_sh*sh22 + mi_mi*mi22 + mi_tc*tc22 + mi_de*de22 + mi_op*op + mi_ny*ny + mi_os*os + mi_ns*ns
  
 # tc4 ~ tc_sh*sh3
 # tc5 ~ tc_sh*sh4 + tc_mi*mi4 + tc_tc*tc4
 # tc6 ~ tc_sh*sh5 + tc_mi*mi5 + tc_tc*tc5
 # tc7 ~ tc_sh*sh6 + tc_mi*mi6 + tc_tc*tc6
 # tc8 ~ tc_sh*sh7 + tc_mi*mi7 + tc_tc*tc7
 # tc9 ~ tc_sh*sh8 + tc_mi*mi8 + tc_tc*tc8
 # tc10 ~ tc_sh*sh9 + tc_mi*mi9 + tc_tc*tc9
 # tc11 ~ tc_sh*sh10 + tc_mi*mi10 + tc_tc*tc10
 # tc12 ~ tc_sh*sh11 + tc_mi*mi11 + tc_tc*tc11
  tc13 ~ tc_sh*sh12 + tc_mi*mi12 + tc_tc*tc12
  tc14 ~ tc_sh*sh13 + tc_mi*mi13 + tc_tc*tc13
  tc15 ~ tc_sh*sh14 + tc_mi*mi14 + tc_tc*tc14
  tc16 ~ tc_sh*sh15 + tc_mi*mi15 + tc_tc*tc15
  tc17 ~ tc_sh*sh16 + tc_mi*mi16 + tc_tc*tc16
  tc18 ~ tc_sh*sh17 + tc_mi*mi17 + tc_tc*tc17
  tc19 ~ tc_sh*sh18 + tc_mi*mi18 + tc_tc*tc18
  tc20 ~ tc_sh*sh19 + tc_mi*mi19 + tc_tc*tc19
  tc21 ~ tc_sh*sh20 + tc_mi*mi20 + tc_tc*tc20
  tc22 ~ tc_sh*sh21 + tc_mi*mi21 + tc_tc*tc21
  tc23 ~ tc_sh*sh22 + tc_mi*mi22 + tc_tc*tc22
  
  # de6 ~ de_mi*mi5 + de_sh*sh5             + de_op*op + de_os*os + de_ny*ny + de_ns*ns
 # de7 ~ de_mi*mi6 + de_sh*sh6 + de_de*de6 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
 # de8 ~ de_mi*mi7 + de_sh*sh7 + de_de*de7 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
 # de9 ~ de_mi*mi8 + de_sh*sh8 + de_de*de8 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
 # de10~ de_mi*mi9 + de_sh*sh9 + de_de*de9 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
 # de11~ de_mi*mi10 + de_sh*sh10 + de_de*de10 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
 # de12~ de_mi*mi11 + de_sh*sh11 + de_de*de11 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de13~ de_mi*mi12 + de_sh*sh12 + de_de*de12 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de14~ de_mi*mi13 + de_sh*sh13 + de_de*de13 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de15~ de_mi*mi14 + de_sh*sh14 + de_de*de14 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de16~ de_mi*mi15 + de_sh*sh15 + de_de*de15 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de17~ de_mi*mi16 + de_sh*sh16 + de_de*de16 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de18~ de_mi*mi17 + de_sh*sh17 + de_de*de17 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de19~ de_mi*mi18 + de_sh*sh18 + de_de*de18 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de20~ de_mi*mi19 + de_sh*sh19 + de_de*de19 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de21~ de_mi*mi20 + de_sh*sh20 + de_de*de20 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de22~ de_mi*mi21 + de_sh*sh21 + de_de*de21 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  de23~ de_mi*mi22 + de_sh*sh22 + de_de*de22 + de_op*op + de_os*os + de_ny*ny + de_ns*ns
  
  # opp_eff_on_shock := sh_op
  # opp_eff_on_mining := sh_op*mi_sh + sh_op*tc_sh*mi_tc + mi_op + de_op*mi_de
  # opp_eff_on_disp := -opp_eff_on_mining/opp_eff_on_shock
  # 
  # os_eff_on_shock := sh_os
  # os_eff_on_mining := sh_os*mi_sh + sh_os*tc_sh*mi_tc + de_os*mi_de
  # os_eff_on_disp := -os_eff_on_mining/os_eff_on_shock
  # 
  # ny_eff_on_shock := sh_ny
  # ny_eff_on_mining := sh_ny*mi_sh + sh_ny*mi_sh*mi_mi + sh_ny*tc_sh*mi_tc + sh_ny*tc_sh*mi_tc*mi_mi + mi_ny + de_ny*mi_de
  # ny_eff_on_disp := -ny_eff_on_mining/ny_eff_on_shock
  # 
  # ns_eff_on_shock := sh_ns
  # ns_eff_on_mining := sh_ns*mi_sh + sh_ns*tc_sh*mi_tc + de_ns*mi_de
  # ns_eff_on_disp := -ns_eff_on_mining/ns_eff_on_shock

  op_eff := sh_op*di_sh + sh_op*mi_sh*di_mi + sh_op*sh_sh*di_sh + sh_op*mi_sh*mi_mi*di_mi + sh_op*mi_sh*de_mi*di_de + sh_op*mi_sh*de_mi*de_de*di_de + sh_op*mi_sh*mi_mi*de_mi*di_de + sh_op*mi_sh*mi_mi*de_mi*de_de*di_de   + sh_op*mi_sh*tc_mi*mi_tc*di_mi + sh_op*tc_sh*mi_tc*di_mi + sh_op*sh_sh*tc_sh*mi_tc*di_mi + sh_op*sh_sh*tc_sh*tc_tc*mi_tc*di_mi + sh_op*tc_sh*tc_tc*mi_tc*di_mi + sh_op*tc_sh*mi_tc*mi_mi*di_mi + sh_op*sh_sh*tc_sh*mi_tc*mi_mi*di_mi + sh_op*mi_sh*tc_mi*mi_tc*mi_mi*di_mi + sh_op*mi_sh*tc_mi*mi_tc*de_mi*di_de + sh_op*mi_sh*tc_mi*mi_tc*de_mi*de_de*di_de + sh_op*mi_sh*mi_mi*tc_mi*mi_tc*de_mi*di_de + sh_op*mi_sh*mi_mi*tc_mi*mi_tc*de_mi*de_de*di_de       + de_op*di_de  + de_op*de_de*di_de + de_op*mi_de*di_mi + de_op*de_de*mi_de*di_mi + de_op*mi_de*mi_mi*di_mi + de_op*de_de*mi_de*mi_mi*di_mi + de_op*mi_de*tc_mi*mi_tc*di_mi + de_op*de_de*mi_de*tc_mi*mi_tc*di_mi + de_op*mi_de*mi_mi*tc_mi*mi_tc*di_mi + de_op*de_de*mi_de*mi_mi*tc_mi*mi_tc*di_mi +  de_op*de_de*mi_de*mi_mi*tc_mi*tc_tc*mi_tc*di_mi +  de_op*de_de*mi_de*mi_mi*tc_mi*mi_tc*mi_mi*di_mi +  de_op*de_de*mi_de*mi_mi*tc_mi*tc_tc*mi_tc*mi_mi*di_mi  + mi_op*mi_mi*di_mi + mi_op*di_mi + mi_op*tc_mi*mi_tc*di_mi + mi_op*mi_mi*tc_mi*mi_tc*di_mi + mi_op*tc_mi*tc_tc*mi_tc*di_mi + mi_op*mi_mi*tc_mi*tc_tc*mi_tc*mi_mi*di_mi + mi_op*de_mi*di_de + mi_op*de_mi*mi_de*di_mi + mi_op*mi_mi*de_mi*di_de + mi_op*de_mi*de_de*di_de
  os_eff := sh_os*di_sh + sh_os*mi_sh*di_mi + sh_os*sh_sh*di_sh + sh_os*mi_sh*mi_mi*di_mi + sh_os*mi_sh*de_mi*di_de + sh_os*mi_sh*de_mi*de_de*di_de + sh_os*mi_sh*mi_mi*de_mi*di_de + sh_os*mi_sh*mi_mi*de_mi*de_de*di_de   + sh_os*mi_sh*tc_mi*mi_tc*di_mi + sh_os*tc_sh*mi_tc*di_mi + sh_os*sh_sh*tc_sh*mi_tc*di_mi + sh_os*sh_sh*tc_sh*tc_tc*mi_tc*di_mi + sh_os*tc_sh*tc_tc*mi_tc*di_mi + sh_os*tc_sh*mi_tc*mi_mi*di_mi + sh_os*sh_sh*tc_sh*mi_tc*mi_mi*di_mi + sh_os*mi_sh*tc_mi*mi_tc*mi_mi*di_mi + sh_os*mi_sh*tc_mi*mi_tc*de_mi*di_de + sh_os*mi_sh*tc_mi*mi_tc*de_mi*de_de*di_de + sh_os*mi_sh*mi_mi*tc_mi*mi_tc*de_mi*di_de + sh_os*mi_sh*mi_mi*tc_mi*mi_tc*de_mi*de_de*di_de       + de_os*di_de  + de_os*de_de*di_de + de_os*mi_de*di_mi + de_os*de_de*mi_de*di_mi + de_os*mi_de*mi_mi*di_mi + de_os*de_de*mi_de*mi_mi*di_mi + de_os*mi_de*tc_mi*mi_tc*di_mi + de_os*de_de*mi_de*tc_mi*mi_tc*di_mi + de_os*mi_de*mi_mi*tc_mi*mi_tc*di_mi + de_os*de_de*mi_de*mi_mi*tc_mi*mi_tc*di_mi +  de_os*de_de*mi_de*mi_mi*tc_mi*tc_tc*mi_tc*di_mi +  de_os*de_de*mi_de*mi_mi*tc_mi*mi_tc*mi_mi*di_mi +  de_os*de_de*mi_de*mi_mi*tc_mi*tc_tc*mi_tc*mi_mi*di_mi  + mi_os*mi_mi*di_mi + mi_os*di_mi + mi_os*tc_mi*mi_tc*di_mi + mi_os*mi_mi*tc_mi*mi_tc*di_mi + mi_os*tc_mi*tc_tc*mi_tc*di_mi + mi_os*mi_mi*tc_mi*tc_tc*mi_tc*mi_mi*di_mi + mi_os*de_mi*di_de + mi_os*de_mi*mi_de*di_mi + mi_os*mi_mi*de_mi*di_de + mi_os*de_mi*de_de*di_de                                
  ny_eff := sh_ny*di_sh + sh_ny*mi_sh*di_mi + sh_ny*sh_sh*di_sh + sh_ny*mi_sh*mi_mi*di_mi + sh_ny*mi_sh*de_mi*di_de + sh_ny*mi_sh*de_mi*de_de*di_de + sh_ny*mi_sh*mi_mi*de_mi*di_de + sh_ny*mi_sh*mi_mi*de_mi*de_de*di_de   + sh_ny*mi_sh*tc_mi*mi_tc*di_mi + sh_ny*tc_sh*mi_tc*di_mi + sh_ny*sh_sh*tc_sh*mi_tc*di_mi + sh_ny*sh_sh*tc_sh*tc_tc*mi_tc*di_mi + sh_ny*tc_sh*tc_tc*mi_tc*di_mi + sh_ny*tc_sh*mi_tc*mi_mi*di_mi + sh_ny*sh_sh*tc_sh*mi_tc*mi_mi*di_mi + sh_ny*mi_sh*tc_mi*mi_tc*mi_mi*di_mi + sh_ny*mi_sh*tc_mi*mi_tc*de_mi*di_de + sh_ny*mi_sh*tc_mi*mi_tc*de_mi*de_de*di_de + sh_ny*mi_sh*mi_mi*tc_mi*mi_tc*de_mi*di_de + sh_ny*mi_sh*mi_mi*tc_mi*mi_tc*de_mi*de_de*di_de       + de_ny*di_de  + de_ny*de_de*di_de + de_ny*mi_de*di_mi + de_ny*de_de*mi_de*di_mi + de_ny*mi_de*mi_mi*di_mi + de_ny*de_de*mi_de*mi_mi*di_mi + de_ny*mi_de*tc_mi*mi_tc*di_mi + de_ny*de_de*mi_de*tc_mi*mi_tc*di_mi + de_ny*mi_de*mi_mi*tc_mi*mi_tc*di_mi + de_ny*de_de*mi_de*mi_mi*tc_mi*mi_tc*di_mi +  de_ny*de_de*mi_de*mi_mi*tc_mi*tc_tc*mi_tc*di_mi +  de_ny*de_de*mi_de*mi_mi*tc_mi*mi_tc*mi_mi*di_mi +  de_ny*de_de*mi_de*mi_mi*tc_mi*tc_tc*mi_tc*mi_mi*di_mi  + mi_ny*mi_mi*di_mi + mi_ny*di_mi + mi_ny*tc_mi*mi_tc*di_mi + mi_ny*mi_mi*tc_mi*mi_tc*di_mi + mi_ny*tc_mi*tc_tc*mi_tc*di_mi + mi_ny*mi_mi*tc_mi*tc_tc*mi_tc*mi_mi*di_mi + mi_ny*de_mi*di_de + mi_ny*de_mi*mi_de*di_mi + mi_ny*mi_mi*de_mi*di_de + mi_ny*de_mi*de_de*di_de
  ns_eff := sh_ns*di_sh + sh_ns*mi_sh*di_mi + sh_ns*sh_sh*di_sh + sh_ns*mi_sh*mi_mi*di_mi + sh_ns*mi_sh*de_mi*di_de + sh_ns*mi_sh*de_mi*de_de*di_de + sh_ns*mi_sh*mi_mi*de_mi*di_de + sh_ns*mi_sh*mi_mi*de_mi*de_de*di_de   + sh_ns*mi_sh*tc_mi*mi_tc*di_mi + sh_ns*tc_sh*mi_tc*di_mi + sh_ns*sh_sh*tc_sh*mi_tc*di_mi + sh_ns*sh_sh*tc_sh*tc_tc*mi_tc*di_mi + sh_ns*tc_sh*tc_tc*mi_tc*di_mi + sh_ns*tc_sh*mi_tc*mi_mi*di_mi + sh_ns*sh_sh*tc_sh*mi_tc*mi_mi*di_mi + sh_ns*mi_sh*tc_mi*mi_tc*mi_mi*di_mi + sh_ns*mi_sh*tc_mi*mi_tc*de_mi*di_de + sh_ns*mi_sh*tc_mi*mi_tc*de_mi*de_de*di_de + sh_ns*mi_sh*mi_mi*tc_mi*mi_tc*de_mi*di_de + sh_ns*mi_sh*mi_mi*tc_mi*mi_tc*de_mi*de_de*di_de       + de_ns*di_de  + de_ns*de_de*di_de + de_ns*mi_de*di_mi + de_ns*de_de*mi_de*di_mi + de_ns*mi_de*mi_mi*di_mi + de_ns*de_de*mi_de*mi_mi*di_mi + de_ns*mi_de*tc_mi*mi_tc*di_mi + de_ns*de_de*mi_de*tc_mi*mi_tc*di_mi + de_ns*mi_de*mi_mi*tc_mi*mi_tc*di_mi + de_ns*de_de*mi_de*mi_mi*tc_mi*mi_tc*di_mi +  de_ns*de_de*mi_de*mi_mi*tc_mi*tc_tc*mi_tc*di_mi +  de_ns*de_de*mi_de*mi_mi*tc_mi*mi_tc*mi_mi*di_mi +  de_ns*de_de*mi_de*mi_mi*tc_mi*tc_tc*mi_tc*mi_mi*di_mi  + mi_ns*mi_mi*di_mi + mi_ns*di_mi + mi_ns*tc_mi*mi_tc*di_mi + mi_ns*mi_mi*tc_mi*mi_tc*di_mi + mi_ns*tc_mi*tc_tc*mi_tc*di_mi + mi_ns*mi_mi*tc_mi*tc_tc*mi_tc*mi_mi*di_mi + mi_ns*de_mi*di_de + mi_ns*de_mi*mi_de*di_mi + mi_ns*mi_mi*de_mi*di_de + mi_ns*de_mi*de_de*di_de                 
  
  ns_o := mi_ns
  os_o := mi_os

  op ~~ 0*ny
  op ~~ 0*ns
  os ~~ 0*ny
  os ~~ 0*ns

'

just_shock <- '
sh2 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns
sh3 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh2
sh4 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh3
sh5 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh4
sh6 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh5
sh7 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh6
sh8 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh7
sh9 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh8
sh10 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh9
sh11 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh10
sh12 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh11
sh13 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh12
sh14 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh13
sh15 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh14
sh16 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh15
sh17 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh16
sh18 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh17
sh19 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh18
sh20 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh19
sh21 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh20
sh22 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh21
sh23 ~ sh_op*op + sh_os*os + sh_ny*ny + sh_ns*ns + sh_sh*sh22

di23 ~ di_sh*sh23 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di22 ~ di_sh*sh22 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di21 ~ di_sh*sh21 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di20 ~ di_sh*sh20 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di19 ~ di_sh*sh19 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di18 ~ di_sh*sh18 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di17 ~ di_sh*sh17 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di16 ~ di_sh*sh16 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di15 ~ di_sh*sh15 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di14 ~ di_sh*sh14 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di13 ~ di_sh*sh13 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di12 ~ di_sh*sh12 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di11 ~ di_sh*sh11 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di10 ~ di_sh*sh10 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di9 ~ di_sh*sh9 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di8 ~ di_sh*sh8 + di_op*op + di_os*os + di_ny*ny + di_ns*ns
di7 ~ di_sh*sh7 + di_op*op + di_os*os + di_ny*ny + di_ns*ns


op_eff := sh_op*di_sh + sh_op*sh_sh*di_sh + di_op
os_eff := sh_os*di_sh + sh_os*sh_sh*di_sh + di_os
ny_eff := sh_ny*di_sh + sh_ny*sh_sh*di_sh + di_ny
ns_eff := sh_ns*di_sh + sh_ns*sh_sh*di_sh + di_ns

sh_eff:= di_sh + sh_sh*di_sh

op ~~ 0*ny
op ~~ 0*ns
os ~~ 0*ny
os ~~ 0*ns

'

fit_high <- growth(model=another_try, data = df_high, check.gradient=FALSE)
fit_low <- growth(model=another_try, data = df_low, check.gradient=FALSE)
fit_all <- growth(model=another_try, data = df_all, check.gradient=FALSE)
fit_sh <- growth(model=just_shock, data = df_all, check.gradient=FALSE)
summary(fit_all)

varTable(fit_all)

x_all <- parameterEstimates(fit_all)
x_low <- parameterEstimates(fit_low)
x_high <- parameterEstimates(fit_high)
x_sh <- parameterEstimates(fit_sh)
