/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Mar 2020
*  LAST MODIFIED:    Mar 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      State dynamics for minimum-fuel bang-bang low thrust control using MEEs
*  REFERENCE:
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "const.hpp"
#include "c_functions.hpp"
#include "spacecraft_params.hpp"
#include "eclipse_model.hpp"

#include <boost/array.hpp>

typedef boost::array< double , 14 > state_type;

void states_twobody_thrust( const state_type &x , state_type &dxdt , const double t ){

  // // Temp
  double rho = 1.0;

  // MEE States
  double p = x[0];
  double f = x[1];
  double g = x[2];
  double h = x[3];
  double k = x[4];
  double L = x[5];
  double m = x[6];

  double Pa;
  double y[7] = {0.0};
  for (int i=0; i<=6; i++){
    y[i] = x[i];
  }

  eclipse_model(y,P,rho,true,&Pa);

  double Thr = Pa*A*eff/Isp/g0;  // Available Thrust (N)

  // MEE Costates
  double pco1 = x[7];
  double pco2 = x[8];
  double pco3 = x[9];
  double pco4 = x[10];
  double pco5 = x[11];
  double pco6 = x[12];
  double pco7 = x[13];

  // Computed Symbolically in MATLAB
  double t2 = cos(L);
  double t3 = sin(L);
  double t8 = f*t2;
  double t9 = g*t3;
  double t4 = t8+t9+1.0;
  double t5 = pow(t4,2.0);
  double t6 = 1.0/m;
  double t7 = h*t3;
  double t10 = sqrt(p);
  double t11 = 1.0/t4;
  double t12 = pow(h,2.0);
  double t13 = pow(k,2.0);
  double t14 = t12+t13+1.0;
  double t29 = k*t2;
  double t15 = t7-t29;
  double t30 = pco5*t3*t10*t11*t14*(1.0/2.0);
  double t31 = g*pco2*t10*t11*t15;
  double t32 = pco4*t2*t10*t11*t14*(1.0/2.0);
  double t16 = fabs(t30-t31+t32+pco6*t10*t11*t15+f*pco3*t10*t11*t15);
  double t17 = pow(p,1.5);
  double t18 = t8+t9+2.0;
  double t21 = pco1*t11*t17*2.0;
  double t22 = t2*t18;
  double t23 = f+t22;
  double t24 = pco2*t10*t11*t23;
  double t25 = t3*t18;
  double t26 = g+t25;
  double t27 = pco3*t10*t11*t26;
  double t28 = t21+t24+t27;
  double t19 = fabs(t28);
  double t35 = pco3*t2*t10;
  double t36 = pco2*t3*t10;
  double t37 = t35-t36;
  double t20 = fabs(t37);
  double t44 = pco6*t10*t11*(t7-t29);
  double t45 = f*pco3*t10*t11*(t7-t29);
  double t46 = t30-t31+t32+t44+t45;
  double t33 = fabs(t46);
  double t34 = pow(t19,2.0);
  double t38 = pow(t20,2.0);
  double t39 = 1.0/pow(p,4.0);
  double t40 = t7-t29;
  double t41 = 1.0/pow(t14,2.0);
  double t42 = pow(t5,2.0);
  double t43 = 1.0/rho;
  double t47 = pow(t33,2.0);
  double t48 = t34+t38+t47;
  double t49 = 1.0/sqrt(t48);
  double t50 = h*t2;
  double t51 = k*t3;
  double t52 = t50+t51;
  double t53 = J2*t15*t39*t41*t42*t52*12.0;
  double t54 = sqrt(t48);
  double t55 = c*eta*t6*t54;
  double t56 = pco7+t55-1.0;
  double t57 = t43*t56;
  double t58 = tanh(t57);
  double t59 = t58*(1.0/2.0);
  double t60 = t59+1.0/2.0;
  double t61 = t7-t29;
  double t62 = Thr*eta*t6*t37*t49*t60;
  double t63 = Thr*eta*t6*t28*t49*t60;
  double t64 = t53+t63;
  double t65 = t12+t13-1.0;
  double t66 = Thr*eta*t6*t46*t49*t60;
  double t68 = J2*t15*t39*t41*t42*t65*6.0;
  double t67 = t66-t68;
  double t69 = t7-t29;
  double t70 = 1.0/sqrt(p);
  double t71 = t7-t29;
  double t72 = 1.0/pow(p,4.5);
  double t73 = t7-t29;
  double t74 = t7-t29;
  double t75 = t7-t29;
  double t76 = t7-t29;
  double t77 = t7-t29;
  double t78 = 1.0/t14;
  double t79 = 1.0/pow(p,1.5);
  double t80 = 1.0/pow(p,3.5);
  double t81 = t7-t29;
  double t82 = 1.0/pow(t4,2.0);
  double t83 = t7-t29;
  double t84 = t7-t29;
  double t85 = pow(t2,2.0);
  double t86 = t7-t29;
  double t87 = t7-t29;
  double t88 = 1.0/pow(p,2.5);
  double t89 = t7-t29;
  double t90 = t2*t3*t10*t11*t64;
  double t91 = t7-t29;
  double t92 = pow(t3,2.0);
  double t93 = t7-t29;
  double t94 = t7-t29;
  double t95 = t7-t29;
  double t96 = -t7+t29;
  double t98 = pco6*t10*t11*t96;
  double t99 = f*pco3*t10*t11*t96;
  double t100 = g*pco2*t10*t11*t96;
  double t101 = t30+t32-t98-t99+t100;
  double t97 = fabs(t101);
  double t102 = pow(t97,2.0);
  double t103 = t34+t38+t102;
  double t104 = 1.0/pow(t14,3.0);
  double t105 = J2*t39*t41*t42*t65*t96*6.0;
  double t106 = sqrt(t103);
  double t107 = c*eta*t6*t106;
  double t108 = pco7+t107-1.0;
  double t109 = t43*t108;
  double t110 = tanh(t109);
  double t111 = t110*(1.0/2.0);
  double t112 = t111+1.0/2.0;
  double t113 = 1.0/sqrt(t103);
  double t114 = Thr*eta*t6*t101*t112*t113;
  double t115 = t105+t114;
  double t116 = J2*t3*t39*t41*t42*t65*6.0;
  double t117 = J2*h*t39*t42*t65*t96*t104*24.0;
  double t122 = J2*h*t39*t41*t42*t96*12.0;
  double t118 = t116+t117-t122;
  double t119 = J2*t3*t39*t41*t42*t52*12.0;
  double t120 = J2*h*t39*t42*t52*t96*t104*48.0;
  double t127 = J2*t2*t39*t41*t42*t96*12.0;
  double t121 = t119+t120-t127;
  double t123 = t3*t41*t96*24.0;
  double t124 = pow(t9,2.0);
  double t125 = h*t104*t124*48.0;
  double t126 = t123+t125;
  double t128 = J2*t2*t39*t41*t42*t65*6.0;
  double t129 = J2*k*t39*t41*t42*t96*12.0;
  double t134 = J2*k*t39*t42*t65*t96*t104*24.0;
  double t130 = t128+t129-t134;
  double t131 = J2*t2*t39*t41*t42*t52*12.0;
  double t132 = J2*t3*t39*t41*t42*t96*12.0;
  double t137 = J2*k*t39*t42*t52*t96*t104*48.0;
  double t133 = t131+t132-t137;
  double t135 = k*t104*t124*48.0;
  double t136 = t135-t2*t41*t96*24.0;
  double t138 = t41*t124*12.0;
  double t139 = t138-1.0;
  double t140 = f*t3;
  double t142 = g*t2;
  double t141 = t140-t142;
  double t143 = Thr*eta*t6*t28*t112*t113;
  double t151 = J2*t39*t41*t42*t52*t96*12.0;
  double t144 = t143-t151;
  double t145 = J2*t4*t5*t39*t139*t141*6.0;
  double t146 = J2*t39*t41*t42*t52*t96*36.0;
  double t147 = t145+t146;
  double t148 = J2*t39*t42*t139*(3.0/2.0);
  double t149 = Thr*eta*t6*t37*t112*t113;
  double t150 = t148+t149;
  double t152 = pow(t52,2.0);
  double t153 = J2*t39*t41*t42*t152*12.0;
  double t154 = J2*t4*t5*t39*t41*t52*t96*t141*48.0;
  double t155 = J2*t39*t41*t42*t52*t65*6.0;
  double t156 = J2*t4*t5*t39*t41*t65*t96*(t140-t142)*24.0;
  double t157 = t155+t156;
  double t158 = 1.0/pow(m,2.0);

  dxdt[0] = t11*t17*(t53+Thr*eta*t6*t28*t49*(tanh(t43*(pco7+c*eta*t6*sqrt(t34+t38+pow(t16,2.0))-1.0))*(1.0/2.0)+1.0/2.0))*-2.0;
  dxdt[1] = t3*t10*(t62+J2*t39*t42*(pow(t40,2.0)*t41*1.2e1-1.0)*(3.0/2.0))-t10*t11*t23*t64+g*t10*t11*t15*t67;
  dxdt[2] = -t2*t10*(t62+J2*t39*t42*(t41*pow(t61,2.0)*1.2e1-1.0)*(3.0/2.0))-t10*t11*t26*t64-f*t10*t11*t15*t67;
  dxdt[3] = t2*t10*t11*t14*t67*(-1.0/2.0);
  dxdt[4] = t3*t10*t11*t14*t67*(-1.0/2.0);
  dxdt[5] = t5*t79-t10*t11*t15*t67;
  dxdt[6] = -(Thr*t60)/c;

  dxdt[7] = pco6*(t5*t88*(3.0/2.0)+t11*t67*t70*(t7-t29)*(1.0/2.0)+J2*t4*t5*t41*t65*pow(t69,2.0)*t72*2.4e1)-pco2*(t3*t70*(t62+J2*t39*t42*(t41*pow(t71,2.0)*1.2e1-1.0)*(3.0/2.0))*(1.0/2.0)-t11*t23*t64*t70*(1.0/2.0)-J2*t3*t42*t72*(t41*pow(t73,2.0)*1.2e1-1.0)*6.0+g*t11*t15*t67*t70*(1.0/2.0)+J2*g*t4*t5*t41*t65*t72*pow(t74,2.0)*2.4e1+J2*t4*t5*t15*t23*t41*t52*t72*4.8e1)+pco3*(t2*t70*(t62+J2*t39*t42*(t41*pow(t75,2.0)*1.2e1-1.0)*(3.0/2.0))*(1.0/2.0)+t11*t26*t70*(t53+t63)*(1.0/2.0)-J2*t2*t42*t72*(t41*pow(t76,2.0)*1.2e1-1.0)*6.0+f*t11*t67*t70*(t7-t29)*(1.0/2.0)+J2*f*t4*t5*t41*t65*t72*pow(t77,2.0)*2.4e1-J2*t4*t5*t15*t26*t41*t52*t72*4.8e1)+pco1*t10*t11*(t53+t63)*3.0+pco4*t2*t11*t14*t67*t70*(1.0/4.0)+pco5*t3*t11*t14*t67*t70*(1.0/4.0)-J2*pco1*t4*t5*t15*t41*t52*t80*9.6e1+J2*pco4*t2*t4*t5*t65*t72*t78*(t7-t29)*1.2e1+J2*pco5*t3*t4*t5*t65*t72*t78*(t7-t29)*1.2e1;

  dxdt[8] = pco3*(t90+t10*t11*t15*t67-t2*t10*t26*t82*(t53+t63)-f*t2*t10*t67*t82*(t7-t29)+J2*t4*t5*t80*t85*(t41*pow(t86,2.0)*12.0-1.0)*6.0-J2*f*t2*t5*t41*t65*t80*pow(t87,2.0)*24.0+J2*t2*t5*t15*t26*t41*t52*t80*48.0)+pco2*(t10*t11*t64*(t85+1.0)-t2*t10*t23*t82*(t53+t63)+g*t2*t10*t15*t67*t82-J2*t2*t3*t4*t5*t80*(t41*pow(t83,2.0)*12.0-1.0)*6.0+J2*g*t2*t5*t41*t65*t80*pow(t84,2.0)*24.0+J2*t2*t5*t15*t23*t41*t52*t80*48.0)-pco6*(t2*t4*t79*2.0+t2*t10*t15*t67*t82+J2*t2*t5*t41*t65*t80*pow(t81,2.0)*24.0)-pco1*t2*t17*t64*t82*2.0-pco4*t10*t14*t67*t82*t85*(1.0/2.0)-pco5*t2*t3*t10*t14*t67*t82*(1.0/2.0)+J2*pco1*t2*t5*t41*t52*t88*(t7-t29)*96.0-J2*pco4*t5*t15*t65*t78*t80*t85*12.0-J2*pco5*t2*t3*t5*t15*t65*t78*t80*12.0;

  dxdt[9] = pco2*(t90-t10*t11*t67*(t7-t29)-t3*t10*t23*t82*(t53+t63)+g*t3*t10*t15*t67*t82-J2*t4*t5*t80*t92*(t41*pow(t89,2.0)*12.0-1.0)*6.0+J2*g*t3*t5*t41*t65*t80*pow(t91,2.0)*24.0+J2*t3*t5*t15*t23*t41*t52*t80*48.0)-pco6*(t3*t4*t79*2.0+t3*t10*t15*t67*t82+J2*t3*t5*t41*t65*t80*pow(t95,2.0)*24.0)+pco3*(t10*t11*t64*(t92+1.0)-t3*t10*t26*t82*(t53+t63)-f*t3*t10*t67*t82*(t7-t29)+J2*t2*t3*t4*t5*t80*(t41*pow(t93,2.0)*12.0-1.0)*6.0-J2*f*t3*t5*t41*t65*t80*pow(t94,2.0)*24.0+J2*t3*t5*t15*t26*t41*t52*t80*48.0)-pco1*t3*t17*t64*t82*2.0-pco5*t10*t14*t67*t82*t92*(1.0/2.0)-pco4*t2*t3*t10*t14*t67*t82*(1.0/2.0)+J2*pco1*t3*t5*t41*t52*t88*(t7-t29)*96.0-J2*pco5*t5*t15*t65*t78*t80*t92*12.0-J2*pco4*t2*t3*t5*t15*t65*t78*t80*12.0;

  dxdt[10] = pco3*(t10*t11*t26*t121-J2*t2*t42*t80*t126*(3.0/2.0)+f*t3*t10*t11*t115+f*t10*t11*t96*t118)+pco2*(t10*t11*t23*t121+J2*t3*t42*t80*t126*(3.0/2.0)-g*t3*t10*t11*t115-g*t10*t11*t96*t118)+pco6*(t3*t10*t11*t115+t10*t11*t96*t118)+pco1*t11*t17*t121*2.0+h*pco4*t2*t10*t11*t115+h*pco5*t3*t10*t11*t115-pco4*t2*t10*t11*t14*t118*(1.0/2.0)-pco5*t3*t10*t11*t14*t118*(1.0/2.0);

  dxdt[11] = -pco3*(t10*t11*t26*t133+J2*t2*t42*t80*t136*(3.0/2.0)+f*t2*t10*t11*t115+f*t10*t11*t96*t130)+pco2*(-t10*t11*t23*t133+J2*t3*t42*t80*t136*(3.0/2.0)+g*t2*t10*t11*t115+g*t10*t11*t96*t130)-pco6*(t2*t10*t11*t115+t10*t11*t96*t130)-pco1*t11*t17*t133*2.0+k*pco4*t2*t10*t11*t115+k*pco5*t3*t10*t11*t115+pco4*t2*t10*t11*t14*t130*(1.0/2.0)+pco5*t3*t10*t11*t14*t130*(1.0/2.0);

  dxdt[12] = pco6*(t4*t79*(t140-t142)*2.0+t10*t11*t52*t115+t10*t11*t96*t157-t10*t82*t96*t115*t141)+pco3*(-t2*t10*t147-t3*t10*t150+t10*t11*(t22-t3*t141)*(t143-t151)+t10*t11*t26*(t153+t154-J2*t39*t41*t42*t124*12.0)+t10*t26*t82*(t140-t142)*(t143-t151)+f*t10*t11*t52*t115+f*t10*t11*t96*t157-f*t10*t82*t96*t115*t141)+pco2*(t3*t10*t147-t2*t10*t150+t10*t11*t23*(t153+t154-J2*t39*t41*t42*t124*12.0)-t10*t11*t144*(t25+t2*t141)-g*t10*t11*t52*t115-g*t10*t11*t96*t157+t10*t23*t82*t141*t144+g*t10*t82*t96*t115*t141)+pco1*t11*t17*(t153+t154-J2*t39*t41*t42*t124*12.0)*2.0+pco1*t17*t82*(t140-t142)*(t143-t151)*2.0-pco4*t3*t10*t11*t14*t115*(1.0/2.0)+pco5*t2*t10*t11*t14*t115*(1.0/2.0)-pco4*t2*t10*t11*t14*t157*(1.0/2.0)-pco5*t3*t10*t11*t14*t157*(1.0/2.0)+pco4*t2*t10*t14*t82*t115*(t140-t142)*(1.0/2.0)+pco5*t3*t10*t14*t82*t115*(t140-t142)*(1.0/2.0);

  dxdt[13] = -pco3*(Thr*eta*t2*t10*t37*t112*t113*t158+Thr*eta*t10*t11*t26*t28*t112*t113*t158-Thr*eta*f*t10*t11*t96*t101*t112*t113*t158)-pco2*(-Thr*eta*t3*t10*t37*t112*t113*t158+Thr*eta*t10*t11*t23*t28*t112*t113*t158+Thr*eta*g*t10*t11*t96*t101*t112*t113*t158)-Thr*eta*pco1*t11*t17*t28*t112*t113*t158*2.0+Thr*eta*pco6*t10*t11*t96*t101*t112*t113*t158-Thr*eta*pco4*t2*t10*t11*t14*t101*t112*t113*t158*(1.0/2.0)-Thr*eta*pco5*t3*t10*t11*t14*t101*t112*t113*t158*(1.0/2.0);

}
