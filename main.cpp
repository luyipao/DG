#include <iostream>
#include <gsl/gsl_sf_bessel.h>
using namespace std;
int main(void) {
    double x = 5.0; // ����һ��˫���ȸ����� x ����ֵΪ 5.0
    double y = gsl_sf_bessel_J0(x); // ���㱴�������� J0(5.0) ��ֵ����������洢�� y ��
    printf("J0(%g) = %.18e\n", x, y); // ��ӡ�������ʽ����� x �� y ��ֵ
    cout << "test " << endl;
    int i = 0; // ����һ�����ͱ��� i ����ֵΪ 
    return 0; // ���� 0����ʾ����ɹ�����
}