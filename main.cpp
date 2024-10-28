#include <iostream>
#include <gsl/gsl_sf_bessel.h>
using namespace std;
int main(void) {
    double x = 5.0; // 定义一个双精度浮点数 x 并赋值为 5.0
    double y = gsl_sf_bessel_J0(x); // 计算贝塞尔函数 J0(5.0) 的值，并将结果存储在 y 中
    printf("J0(%g) = %.18e\n", x, y); // 打印结果，格式化输出 x 和 y 的值
    cout << "test " << endl;
    int i = 0; // 定义一个整型变量 i 并赋值为 
    return 0; // 返回 0，表示程序成功结束
}