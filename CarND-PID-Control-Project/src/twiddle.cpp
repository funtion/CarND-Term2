#include "twiddle.h"
#include <limits>

Twiddle::Twiddle(PID& pid) : pid(pid) {

    this->delta[0] = 0.0001;
    this->delta[1] = 0.00001;
    this->delta[2] = 0.0001;
    
    this->value[0] = &this->pid.Kp;
    this->value[1] = &this->pid.Ki;
    this->value[2] = &this->pid.Kd;

    this->bestError = std::numeric_limits<double>::max();
    this->position = 1;
    this->cur = 0;

    *this->value[0] += this->delta[0];
}

void Twiddle::Update(double error) {
    if(this->delta[0] + this->delta[1] + this->delta[2] < 1e-6) {
        return;
    }

    if (this->position == 1) {
        if (error < this->bestError) {
            this->bestError = error;
            this->delta[cur] *= 1.05;
            this->cur = (this->cur + 1)%3;
            *this->value[cur] += this->delta[cur];
        } else {
            this->position = -1;
            *this->value[cur] += - 2 * this->delta[cur];
            if(*this->value[cur] < 0) {
                *this->value[cur] = 0;
            }
        }
    } else {
        if (error < this->bestError) {
            this->bestError = error;
            this->delta[cur] *= 1.1;
        } else {
            this->delta[cur] *= 0.9;
        }
        this->cur = (this->cur + 1)%3;
        *this->value[cur] += this->delta[cur];
        this->position = 1;
    }
}