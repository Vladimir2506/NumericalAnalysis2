#include <iostream>
#include <iomanip>
#include <chrono>

#include "Integral.hpp"
#include "Ode.hpp"
#include "Utils.hpp"

using namespace std;
using namespace chrono;

int main()
{

#ifdef _DEBUG
	constexpr auto precision = 14;
#else
	constexpr auto precision = 6;
#endif 

	cout << fixed << setprecision(precision);

	auto beginPI1 = system_clock::now();
	auto MyPI1 = AtanSqrt3();
	auto endPI1 = system_clock::now();
	auto durationPI1 = duration_cast<nanoseconds>(endPI1 - beginPI1);
	cout << "PI (AtanSqrt3) = " << MyPI1 << endl;
	cout << "Duration for PI (AtanSqrt3):" << durationPI1.count() << "ns" << endl;

	cout << endl;

	auto beginPI2 = system_clock::now();
	auto MyPI2 = BBP();
	auto endPI2 = system_clock::now();
	auto durationPI2 = duration_cast<nanoseconds>(endPI2 - beginPI2);
	cout << "PI (BBP) = " << MyPI2 << endl;
	cout << "Duration for PI (BBP):" << durationPI2.count() << "ns" << endl;

	cout << endl;

	auto func1 = [](double t) {return 1.0 / t; };

	auto beginLnPI1 = chrono::system_clock::now();
	auto LnPI1 = SimpsonRule(func1)(1, MyPI2, 2944);
	auto endLnPI1 = system_clock::now();
	auto durationLnPI1 = duration_cast<chrono::nanoseconds>(endLnPI1 - beginLnPI1);
	cout << "ln PI (SimpsonRule) = " << LnPI1 << endl;
	cout << "Duration for ln PI (SimpsonRule):" << durationLnPI1.count() << "ns" << endl;

	cout << endl;

	auto beginLnPI2 = chrono::system_clock::now();
	auto LnPI2 = CotesRule(func1)(1, MyPI2, 165);
	auto endLnPI2 = system_clock::now();
	auto durationLnPI2 = duration_cast<chrono::nanoseconds>(endLnPI2 - beginLnPI2);
	cout << "ln PI (CotesRule) = " << LnPI2 << endl;
	cout << "Duration for ln PI (CotesRule):" << durationLnPI2.count() << "ns" << endl;

	cout << endl;

	auto x = 0.0;
	cout << "Please input x in [1, 10]:" << endl;
	cin >> x;

	cout << endl;
	
	auto func2 = [](double t, double y) {return y; };

	x = x > 10 ? 10 : x;
	x = x < 1 ? 1 : x;
	
	auto frac = x - int(x);

	auto result1 = 1.0, result11 = 1.0, result12 = 0.0;
	auto beginExp1 = system_clock::now();
	if (frac < 0.5)
	{
		for (int k = 1; k <= int(x); ++k) result11 *= MyPI2;
		auto fracLnPI = frac * LnPI2;
		result12 = TaylorExp(fracLnPI, 16);
		result1 = result11 * result12;
	}
	else
	{
		for (int k = 1; k <= int(x) + 1; ++k) result11 *= MyPI2;
		auto compLnPI = (1 - frac) * LnPI2;
		result12 = TaylorExp(compLnPI, 16);
		result1 = result11 / result12;
	}
	auto endExp1 = system_clock::now();
	auto durationExp1 = duration_cast<nanoseconds>(endExp1 - beginExp1);
	cout << "PI ^ x (TaylorExp) = " << result1 << endl;
	cout << "Duration for PI ^ x (TaylorExp):" << durationExp1.count() << "ns" << endl;

	cout << endl;

	auto result21 = 1.0, result22 = 1.0, result2 = 0.0;
	auto beginExp2 = system_clock::now();
	if (frac < 0.5)
	{
		for (int k = 1; k <= int(x); ++k) result21 *= MyPI2;
		auto fracLnPI = frac * LnPI2;
		result22 = RungeKutta4Method(func2)(0, fracLnPI, 333, 1);
		result2 = result21 * result22;
	}
	else
	{
		for (int k = 1; k <= int(x) + 1; ++k) result21 *= MyPI2;
		auto compLnPI = (1 - frac) * LnPI2;
		result22 = RungeKutta4Method(func2)(0, compLnPI, 333, 1);
		result2 = result21 / result22;
	}
	auto endExp2 = system_clock::now();
	auto durationExp2 = duration_cast<nanoseconds>(endExp2 - beginExp2);
	cout << "PI ^ x (RungeKutta4Method) = " << result2 << endl;
	cout << "Duration for PI ^ x (RungeKutta4Method):" << durationExp2.count() << "ns" << endl;

	system("pause");
	return 0;
}
