#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

double Tsup = 0.8; //周期
std::vector<double>sx = {0.3,0.3,0.3,0.3,0.0};//x歩幅
std::vector<double>sy = {0.2,0.2,0.2,0.2,0.2};//y歩幅   
double x = 0,y = 0,astpx = 0,astpy = 0,px = 0,py = 0,ppx,ppy ; //座標系
double T = 0,n = 0;
double xi = 0,yi = 0,dxi = 0,dyi = 0; //n歩目が始まる瞬間の重心速度と位置
double dx,dy;  //x,yの一回微分
double g = 9.8, zc = 0.8; //zc:高さ
double t = 0;
double dt = 4.0 / 1000.0;
double Tc,S,C;
double xr,yr; //歩行素片
double x_target,y_target,dx_target,dy_target; //目標状態
double a = 10, b = 1; //評価関数の重み
double D;



int main(){
    std::ofstream data("data.dat");  //data.dat用
//step3
    for(double step = 0;step < 4;step++){
        for (double period = 0 ; period < Tsup; period += dt, t += dt){
                Tc = std::sqrt(zc/g);
                S = std::sinh(period / Tc);
                C = std::cosh(period / Tc);
                x = (xi - astpx) * C + Tc * dxi * S + astpx;
                y = (yi - astpx) * C + Tc * dyi * S + astpy;
                dx = (xi - astpx) / Tc * S + dxi * C;
                dy = (yi - astpx) / Tc * S + dyi * C;
                data << t << " " << x << " " << y << std::endl;
                //data << dx << " " << dy << std::endl;
                //std::cout << dx << "," << dy << std::endl;
        }
//step4
        T = T + Tsup;
        n = n + 1;
//xi,dxi,yi,dyiの更新        
        xi = x;
        yi = y;
        dxi = dx;
        dyi = dy;
//step5
        px = astpx + sx.at(step);
        py = astpy + std::pow(-1, n) * sy.at(step);
//step6
        xr = sx.at(step) / 2;
        yr = std::pow(-1, n) * sy.at(step) / 2;
//step7
        x_target = px + xr;
        dx_target = (C + 1) / (Tc * S) * xr;
        y_target = py + yr;
        dy_target = (C - 1) / (Tc * S) * yr;
//step8
        D = a * (C -1) * (C - 1) + b * (S / Tc) * (S /Tc);
        astpx = -a * (C - 1) / D * (x_target - C * xi - Tc * S * dxi) - b * S / (Tc * D) * (dx_target - S * xi / Tc - C * dxi);
        astpy = -a * (C - 1) / D * (y_target - C  * yi - Tc * S * dyi) - b * S / (Tc * D) * (dy_target - S * yi / Tc - C * dyi);
    }
}
