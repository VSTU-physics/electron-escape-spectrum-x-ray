HTML header: 	<script type="text/javascript"
				src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
				</script>
**Функция spe**

***

Будем численно решать уравнение диффузии неявным методом первого порядка.
Уравнение:
\\[
    \frac{\partial \Phi}{\partial t} =
    f(z, t) \frac{\partial^2 \Phi}{\partial z^2} + g(z, t)
\\]
Граничные условия:
\\[
    \begin{aligned}
        & \left(\alpha_1 \Phi + \alpha_2 \frac{\partial \Phi}{\partial z}\right)_{0, t} = \gamma_1(t) \\
        & \left(\beta_1 \Phi + \beta_2 \frac{\partial \Phi}{\partial z}\right)_{l, t} = \gamma_2(t) \\
    \end{aligned}
\\]