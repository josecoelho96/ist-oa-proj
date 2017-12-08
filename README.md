# Selecting reliable sensors via convex optimization

## Optimization and Alghoritms course project.

### Abstract
*Extraction of trusted and relevant information from sensor networks by fusing data from a multitude of heterogeneous, distinct, but possibly unreliable or irrelevant sensors. Recovering the desirable view of the environment from the maximum number of dependable sensors while specifying the unreliable ones. Performance of some methods is tested analytically, and through simulations.*

### Details
Fall Semester of 2017/2018.

Development using MATLAB and [CVX].

### MATLAB simulation notes
In order to generate the maximum number of results, MATLAB's [Parallel Computing Toolbox] was used. In order to make it work, a parallel pool must be running on the machine the scripts will be executed.

To start a new poll using the default local cluster:
```matlab
parpool
```

This command might be helpful if a parallel pool cannot be started.

```matlab
distcomp.feature( 'LocalUseMpiexec', false )
```

### Authors
José Coelho

Miguel de Moura

Gonçalo Pereira

   [CVX]: <http://cvxr.com/>
   [Parallel Computing Toolbox]: <https://www.mathworks.com/products/parallel-computing.html/>