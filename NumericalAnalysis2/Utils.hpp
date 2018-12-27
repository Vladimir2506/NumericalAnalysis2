#pragma once

double MySqrt(double dSquare = 3.0, int nMaxStep = 5)
{
	double result = dSquare;
	for (int k = 0; k < nMaxStep; ++k)
	{
		result = (result + dSquare / result) * 0.5;
	}
	return result;
}

double Atan1(int nMaxStep)
{
	double result = 0.0;
	for (int k = 0; k < nMaxStep; ++k)
	{
		double p = 1.0 / (2 * k + 1);
		if (k % 2 == 0) result += p;
		else result -= p;
	}
	result *= 4;
	return result;
}

double BBP(int nMaxStep = 11)
{
	double result = 0.0;
	double power16 = 1.0;
	for (int k = 0; k < nMaxStep; ++k)
	{
		double p = 
			double(120 * k * k + 151 * k + 47) / 
			double(512 * k * k * k * k + 1024 * k * k * k + 712 * k * k + 194 * k + 15);
		p *= power16;
		power16 /= 16;
		result += p;
	}
	return result;
}

double GaussLegendre(int nMaxStep)
{
	double result = 0.0;
	double a = 1, b = MySqrt(0.5, 10), t = 0.25, p = 1;
	double aNext = a, bNext = b, tNext = t, pNext = p;
	for (int i = 0; i < nMaxStep; ++i)
	{
		aNext = (a + b) * 0.5;
		bNext = sqrt(a * b);
		tNext = t - p * (a - aNext) * (a - aNext);
		pNext = 2 * p;
		result = (aNext + bNext) * (aNext + bNext) / tNext * 0.25;
		a = aNext;
		b = bNext;
		t = tNext;
		p = pNext;
	}
	return result;
}

double AtanSqrt3(int nMaxStep = 31)
{
	double sqrt3 = MySqrt();
	double result = 0.0;
	double power3 = 1.0;
	for (int k = 0; k < nMaxStep; ++k)
	{
		double p = 1.0 / (2 * k + 1) * power3;
		if (k % 2 == 0) result += p;
		else result -= p;
		power3 /= 3;
	}
	result *= 2 * sqrt3;
	return result;
}

double TaylorExp(double x, int nIntervals)
{
	double result = 1.0;
	double dDenominator = 1.0;
	double dNumerator = x;
	for (int k = 1; k < nIntervals; ++k)
	{
		result += dNumerator * dDenominator;
		dNumerator *= x;
		dDenominator /= k + 1;
	}
	return result;
}
