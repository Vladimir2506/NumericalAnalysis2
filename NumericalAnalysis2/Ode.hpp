#pragma once
#include <functional>

class OdeBase
{
protected:
	std::function<double(double, double)> fcnDerivative;
public:
	OdeBase(std::function<double(double, double)> fcn) : fcnDerivative(fcn) {}
	virtual double operator()(const double dBegin, const double dEnd, const int nInterval, const double dInit) = 0;
	virtual ~OdeBase() {}
};

class RungeKutta4Method : public OdeBase
{
public:
	RungeKutta4Method(std::function<double(double, double)> fcn) : OdeBase(fcn) {}
	virtual double operator()(const double dBegin, const double dEnd, const int nInterval, const double dInit)
	{
		double result = dInit;
		double h = (dEnd - dBegin) / nInterval;
		double xk = dBegin;
		for (int k = 0; k < nInterval; ++k)
		{
			double k1 = fcnDerivative(xk, result);
			double k2 = fcnDerivative(xk + 0.5 * h, result + 0.5 * h * k1);
			double k3 = fcnDerivative(xk + 0.5 * h, result + 0.5 * h * k2);
			double k4 = fcnDerivative(xk + h, result + h * k3);
			result += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);			
			xk += h;
		}
		return result;
	}
	virtual ~RungeKutta4Method() {}
};
 