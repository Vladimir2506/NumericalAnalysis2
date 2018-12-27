#pragma once
#include <functional>

class IntegralBase
{
protected:
	std::function<double(double)> fcnIntegral;
public:
	IntegralBase(std::function<double(double)> fcn) : fcnIntegral(fcn) {}
	virtual double operator()(const double dBegin, const double dEnd, const int nInterval) = 0;
	virtual ~IntegralBase() {}
};

class SimpsonRule : public IntegralBase
{
public:
	SimpsonRule(std::function<double(double)> fcn) : IntegralBase(fcn) {}
	virtual double operator()(const double dBegin, const double dEnd, const int nInterval)
	{
		double result = 0.0;
		double h = (dEnd - dBegin) / nInterval;
		double xk = dBegin;
		for (int k = 0; k < nInterval; ++k)
		{
			result += fcnIntegral(xk) + 4 * fcnIntegral(xk + h / 2) + fcnIntegral(xk + h);
			xk += h;
		}
		result *= h / 6;
		return result;
	}
	virtual ~SimpsonRule() {}
}; 

class CotesRule : public IntegralBase
{
public:
	CotesRule(std::function<double(double)> fcn) : IntegralBase(fcn) {}
	virtual double operator()(const double dBegin, const double dEnd, const int nInterval)
	{
		double result = 0.0;
		double h = (dEnd - dBegin) / nInterval;
		double xk = dBegin;
		for (int k = 0; k < nInterval; ++k)
		{
			result += 7 * fcnIntegral(xk)
				+ 32 * fcnIntegral(xk + 0.25 * h)
				+ 12 * fcnIntegral(xk + 0.5 * h)
				+ 32 * fcnIntegral(xk + 0.75 * h)
				+ 7 * fcnIntegral(xk + h);
			xk += h;
		}
		result *= h / 90;
		return result;
	}
	virtual ~CotesRule() {}
};