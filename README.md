### Numerical Optimization

#### Task 1
Optimize
$$
\lambda = \arg\min\left \{ \int_0^1 e^{\lambda_0+\lambda_1x+\lambda_2x^2+\cdots +\lambda_n x^n}dx-(\lambda_0m_0+\lambda_1m_1+\cdots+\lambda_nx_n)\right\}
$$

##### Newton's Method
'''
[minf, lam_, errCode, itCount, fhist, xhist] = Newtons([0.5, 0.1], [2, 2], 0.0001, 100);  
trajectory
'''

Two evidence supporting the derivative rule
- Newtons method uses jacobian and works
- Algebra calculation in bivariable-case matches the partial derivative result

***

#### Task 2

$$
\frac{d}{d\lambda} \int_0^1 e^{\lambda x^2} dx
$$
