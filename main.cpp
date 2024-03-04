#include <QStringBuilder>
#include <iomanip>
#include "Derivative.h"
#include <string.h>


QString GetFileName(QString number, QString name, QString extension)
{
    int index = number.indexOf(extension);

    QString FileName(number);

    FileName.insert(index,name);

    return FileName;
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        puts("Input File (*.dat) not found!");
        return 1;
    }

    QFile FileInput(argv[3]);
    if (!FileInput.open(QIODevice::ReadOnly | QIODevice::Text))
        return EXIT_FAILURE;

    QFile FileApprox(GetFileName(argv[3],"Approx",".dat"));
    FileApprox.open(QIODevice::WriteOnly | QIODevice::Text);

    Derivative d(!strcmp(argv[1],"order"),!strcmp(argv[2],"co"));

    d.GetData(FileInput);

    d.WriteApprox(FileApprox);

    if (!strcmp(argv[1],"order"))
    {

        QFile FileTc(GetFileName(argv[3],"Tc",".dat"));
        FileTc.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileFirstDerivative(GetFileName(argv[3],"FirstDerivative",".dat"));
        FileFirstDerivative.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileApproxFirstDerivative(GetFileName(argv[3],"ApproxFirstDerivative",".dat"));
        FileApproxFirstDerivative.open(QIODevice::WriteOnly | QIODevice::Text);

        d.GetFirstDerivative();
        d.GetApproxFirstDerivative();

        d.WriteFirstDerivative(FileFirstDerivative);
        d.WriteApproxFirstDerivative(FileApproxFirstDerivative);

        d.WriteTc(FileTc);

        FileTc.close();
        FileFirstDerivative.close();
        FileApproxFirstDerivative.close();
    }

    FileInput.close();
    FileApprox.close();


    return 0;

}
