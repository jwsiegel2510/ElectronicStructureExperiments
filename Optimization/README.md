This folder contains optimization methods which are used in our electronic structure calculations.

In particular, we are exploring a direct optimization approach, so the self-consistent field iteration is not implemented.
Instead, we implement gradient descent, accelerated gradient descent, and preconditioned versions of both of these algorithms.
These methods can be found in the Methods folder.

Additionally, we have implemented tools for building preconditioners in the Preconditioners folder, and a variety of retractions
which are used by the gradient and accelerated gradient methods in the Retractions folder.

For calculations, we recommend preconditioned accelerated gradient descent.
