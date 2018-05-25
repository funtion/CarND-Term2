# Write up

## The model

The model are using 6 states:

- x
- y
- psi
- v
- cte
- epsi

And two actuators:

- steering_angle
- throttle

The updates are implemented in `MPC.cpp`, from line 99 to 111.

## The hyper parameters

The `N` is 10 and `dt` is 0.1. The value is got by dozen of different experiments. Also I tune the weights for different parts of the loss function to make the vehicle drive more stable.

## Polynomial fitting

The waypoints are fitted with 3 order polynomial. They are transformed to car coordination fist to make the computation easier.

## Handling the latency

The latency is handled in `main.cpp`, from line 113 to 115. Basically, it assume that the vehicle will continue to drive under current state and predict where it will be after 100ms. The predicted values are sent to MPC for computation as the initial state.

