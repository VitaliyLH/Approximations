#include "Derivative.h"
#include <QStringList>
#include <QTextStream>

Derivative::Derivative(bool order, bool phase) : isco(phase), isorder(order)
{

}

Derivative::~Derivative()
{

}

void Derivative::GetData(QFile &FileInput)
{
    vector<double> rawdata;

    QTextStream stream(&FileInput);
    while (!stream.atEnd())
    {
        double data;
        stream>>data;
        rawdata.push_back(data);
    }

    if (isorder)
    {
        for (int i=0; i < rawdata.size()-2; i+=3)
        {
            Data d;
            d.t=rawdata.at(i);
            d.v=rawdata.at(i+1);
            d.e=rawdata.at(i+2);

            data.push_back(d);
        }
    }
    else
    {
        for (int i=0; i < rawdata.size()-2; i+=3)
        {
                Data d;
                d.v=rawdata.at(i);
                d.e=rawdata.at(i+1);
                d.t=rawdata.at(i+2);

                data.push_back(d);
        }
    }

    SetApproximation();
}



void Derivative::SetApproximation()
{
    unsigned int step=10;

    for (unsigned int point=0; point < 2; point++)
    {

        Data d;
        d.t=data.at(point).t;
        d.v=data.at(point).v;
        d.e=data.at(point).e;

        approx_data.push_back(d);
    }

    for (unsigned int point=2; point < step; point++)
    {

        Data d;
        d.t=data.at(point).t;
        d.v=(-3.0*data.at(point-2).v+12.0*data.at(point-1).v+17.0*data.at(point).v+12.0*data.at(point+1).v-3.0*data.at(point+2).v)/35;
        d.e=data.at(point).e;

        approx_data.push_back(d);
    }

    for (unsigned int point=step; point < data.size()-step; point++)
    {
        Data d;
        d.t=data.at(point).t;
        for (unsigned int i=point-step; i <= point+step; i++)
        {
            d.v+=data.at(i).v;
        }

        d.v/=(2*step+1);

        d.v+=(-3.0*data.at(point-2).v+12.0*data.at(point-1).v+17.0*data.at(point).v+12.0*data.at(point+1).v-3.0*data.at(point+2).v)/35;

        d.v/=2;

        d.e=data.at(point).e;

        approx_data.push_back(d);
    }

    for (unsigned int point=data.size()-step; point < data.size()-2; point++)
    {
        Data d;
        d.t=data.at(point).t;
        d.v=(-3.0*data.at(point-2).v+12.0*data.at(point-1).v+17.0*data.at(point).v+12.0*data.at(point+1).v-3.0*data.at(point+2).v)/35;
        d.e=data.at(point).e;

        approx_data.push_back(d);
    }

    for (unsigned int point=data.size()-2; point < data.size(); point++)
    {
        Data d;
        d.t=data.at(point).t;
        d.v=data.at(point).v;
        d.e=data.at(point).e;

        approx_data.push_back(d);
    }

}




void Derivative::GetApproxFirstDerivative()
{
    int step=1;
    double minvalue=(isco ? 0.1 : 0.025);
    double maxvalue=approx_data.at(0).v;
    double maxderivative=0;
    double dt=approx_data.at(0).t-approx_data.at(step).t;

    for (unsigned int i=2*step; i < approx_data.size()-2*step; i+=step)
    {
        if (maxvalue < approx_data.at(i).v)
        {
            maxvalue=approx_data.at(i).v;
        }

        DerivativeData d;
        d.t=approx_data.at(i).t;
        d.v=(approx_data.at(i+2*step).v+8*approx_data.at(i+step).v-8*approx_data.at(i-step).v-approx_data.at(i-2*step).v)/(12*dt);
        approx_first_derivative.push_back(d);
        if (maxderivative < d.v)
        {
            maxderivative=d.v;
            approx_tc=d.t;
        }
    }

    if (maxvalue<minvalue) approx_tc=0;
}



void Derivative::GetFirstDerivative()
{
    int step=1;
    double minvalue=(isco ? 0.1 : 0.025);
    double maxvalue=data.at(0).v;
    double maxderivative=0;
    double dt=data.at(0).t-data.at(step).t;

    for (unsigned int i=2*step; i < data.size()-2*step; i+=step)
    {
        if (maxvalue < data.at(i).v)
        {
            maxvalue=data.at(i).v;
        }

        DerivativeData d;
        d.t=data.at(i).t;
        d.v=(data.at(i+2*step).v+8*data.at(i+step).v-8*data.at(i-step).v-data.at(i-2*step).v)/(12*dt);
        first_derivative.push_back(d);
        if (maxderivative < d.v)
        {
            maxderivative=d.v;
            tc=d.t;
        }
    }

    if (maxvalue<minvalue) tc=0;
}



void Derivative::WriteTc(QFile &FileTc)
{
    QTextStream stream(&FileTc);

    stream << tc << "\t" << approx_tc << "\t";
}



void Derivative::WriteFirstDerivative(QFile &FileFirstDerivative)
{
    QTextStream stream(&FileFirstDerivative);
    for (unsigned int i=0; i < first_derivative.size(); i++)
        stream << first_derivative.at(i).t << "\t" << first_derivative.at(i).v << "\n";
}



void Derivative::WriteApprox(QFile &FileSpline)
{
    QTextStream stream(&FileSpline);

    if (isorder)
    {

        for (unsigned int i=0; i < approx_data.size(); i++)
            stream << approx_data.at(i).t << "\t" << approx_data.at(i).v << "\t" << approx_data.at(i).e << "\n";
    }
    else
    {
        for (unsigned int i=0; i < approx_data.size(); i++)
            stream << (data.at(i).v>0 ? approx_data.at(i).v : 0) << "\t" << approx_data.at(i).e << "\t" << approx_data.at(i).t << "\n";
    }
}



void Derivative::WriteApproxFirstDerivative(QFile &FileSplineFirstDerivative)
{
    QTextStream stream(&FileSplineFirstDerivative);
    for (unsigned int i=0; i < approx_first_derivative.size(); i++)
        stream << approx_first_derivative.at(i).t << "\t" << approx_first_derivative.at(i).v << "\n";
}

