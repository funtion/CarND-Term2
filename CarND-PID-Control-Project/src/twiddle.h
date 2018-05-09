#ifndef TWIDDLE_H
#define TWIDDLE_H

#include "PID.h"

class Twiddle {

public:
    double delta[3];
    double* value[3];

    Twiddle(PID& pid);
    void Update(double error);

private:
    double bestError;
    int position;
    int cur;
    PID& pid;
};

#endif#endif