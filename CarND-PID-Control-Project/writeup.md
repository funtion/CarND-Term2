1. Describe the effect each of the P, I, D components had in your implementation

The I, P, D are just three level of differential. Which means $$ di/dt = p, dp/dt = d $$. Thus, the effect of them is in decreasing order. To make the vehicle drive stable, we need to make sure that $$ i < p < d $$

2. Describe how the final hyperparameters were chosen

As discussed above, I choose some value that $$ i < p < d $$ and manually change the value. For example, if the car can not turn around successfully because it was turning too slow, I will increase the value of the p. Also, if the car is fluctuating, I may increase the d. Finally, I get $$ i = 0.000001, p = 0.15, d = 30$$

To even get a better result, I use the twiddle algorithm so fine tune the coefficient. It's implemented in `twiddle.cpp`.