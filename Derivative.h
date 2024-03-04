#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include <QFile>
#include <QIODevice>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

struct Data
{
    double t;
    double v;
    double e;
};

struct DerivativeData
{
    double t;
    double v;
};

class Derivative
{
public:
    Derivative(bool order, bool phase);
    ~Derivative();

    void WriteApproxFirstDerivative(QFile &FileSplineFirstDerivative);
    void WriteFirstDerivative(QFile &FileFirstDerivative);
    void WriteApprox(QFile &FileSpline);
    void WriteApproxTc(QFile &FileSplineTc);
    void GetData(QFile &FileInput);
    void WriteTc(QFile &FileTc);
    void GetFirstDerivative();
    void GetApproxFirstDerivative();
    void SetApproximation();

private:

    bool isco;
    bool isorder;
    double tc;
    double approx_tc;

    vector<Data> data;
    vector<Data> approx_data;
    vector<DerivativeData> approx_first_derivative;
    vector<DerivativeData> first_derivative;
};

#endif // DERIVATIVE_H
