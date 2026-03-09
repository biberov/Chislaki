#include "NewTask.h"
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <locale.h>
#include <iostream>
using namespace std;

 double x0 = 0; // фиксированно 
 double xmax = 1; // фиксированно 
 double mu1 = 0.0; // фиксированно 
 double mu2 = 1.0; // фиксированно (для теста)
 double xi = 0.4; // фиксировано 

void outputMessage(double n, double max_diff, double max_x) {
	cout << "\nДля решения задачи использована равномерная сетка с числом разбиений n = " << n << endl;
	cout << "задача должна быть решена с точностью не более e = 0.5*10^–6" << endl;
	cout << "задача решена с точностью e2 = " << max_diff << endl;
	cout << "максимальная разность численных решений в общих узлах сетки наблюдается в точке x = " << max_x << endl;
}

void SolveTaskTestBalance() {
	cout << "== Test task ==" << endl;
	cout << "Введите число разбиений n: ";
	int n;
	cin >> n;

	double k1 = 1.4;
	double k2 = 0.4;
	double q1 = 0.4;
	double q2 = 0.16;
	double f1 = 0.4;
	double f2 = exp(-0.4);

	AnalyticSolver an(xi, k1, k2, q1, q2, f1, f2, mu1, mu2);
	std::vector<double> v = balance_method(n, 0, xi, k1, k2, q1, q2, f1, f2, mu1, mu2);

	double h = (xmax - x0) / n;
	double max_diff = 0.0;
	int max_index = 0;
	double max_x = 0.0;

	ofstream outFile("my_table1.csv");
	outFile.clear();
	if (!outFile.is_open()) cout << "Ошибка открытия файла!" << endl;
	outFile << "TEST" << endl;
	outFile << "i" << ";" << "x_i" << ";" << "U(x_i)" << ";" << "V(x_i)" << ";" << "|U-V|" << "\n";

	for (int i = 0; i <= n; ++i) {
		double x = x0 + i * h;
		double ua = an.u(x);
		double uv = v[i];
		double diff = fabs(ua - uv);

		// --- Фиксируем первую строку ---
		if (i == 0) {
			ua = 0.0;
			uv = 0.0;
			diff = 0.0;
		}

		// обновляем максимум
		if (diff > max_diff) {
			max_diff = diff;
			max_index = i;
			max_x = x;
		}

		outFile << i << ";" << x << ";" << ua << ";" << uv << ";" << diff << "\n";

	}

	outputMessage(n, max_diff, max_x);
	outFile.close();

}

void SolveTaskMainBalance() {
	cout << "======= MAIN TASK =======";
	double x0 = 0;
	double xmax = 1;
	double mu1 = 0.0;
	double mu2 = 1.0; // фиксированно для теста
	double xi = 0.4; // фиксировано
	int n = 10;   // можно брать из TextBox

	// Решение на сетке n
	std::vector<double> v_n = balance_method(n, 1, xi, 0, 0, 0, 0, 0, 0, mu1, mu2);
	// Решение на сетке 2n
	std::vector<double> v_2n = balance_method(2 * n, 1, xi, 0, 0, 0, 0, 0, 0, mu1, mu2);

	double h = 1.0 / n;
	double h2 = 1.0 / (2 * n);

	double max_diff = 0.0;
	int max_index = 0;
	double max_x = 0.0;

	ofstream outFile("my_table1.csv");
	outFile.clear();
	if(!outFile.is_open()) cout << "Ошибка открытия файла!" << endl;

	outFile << "MAIN" << endl;
	outFile << "i" << ";" << "x_i" << ";" << "V_n(x_i)" << ";" << "| V_2n(x_i) |" << ";" << "| V_n - V_2n |" <<"\n";

	for (int i = 0; i <= n; ++i) {
		double x = i * h;
		double v2_at_x = v_2n[2 * i]; // общий узел сетки 2n

		double diff = fabs(v_n[i] - v2_at_x);
		if (diff > max_diff) { max_diff = diff; max_index = i; }

		outFile << i << ";" << x << ";" << v_n[i] << ";" << v2_at_x << ";" << diff << "\n";

	}

	max_x = max_index * h;
	outputMessage(n, max_diff, max_x);
	outFile.close();

}

int main() {
	setlocale(LC_ALL, "RUS");
	bool choice;
	cout << "Test = 0 ; Main = 1\n";
	cin >> choice;
	choice == 1 ? SolveTaskMainBalance() : SolveTaskTestBalance();
}
