package rial

/** `math` package provides floating-point mathematical functions pipeline.
 *
 * It supports the following well-known math functions:
 *
 * <ul>
 *   <li> sqrt </li>
 *   <li> invsqrt </li>
 *   <li> reciprocal </li>
 *   <li> sin </li>
 *   <li> cos </li>
 *   <li> acos </li>
 *   <li> atan2 </li>
 *   <li> exp </li>
 *   <li> log </li>
 *   <li> sigmoid (`1/(1+exp(-x))`) </li>
 *   <li> softplus (`log(1+exp(x))`) </li>
 * </ul>
 *
 * All those functions are approximated by chebyshev polynomials after an
 * appropreate mathematical conversions.
 * The precision of the approximation functions provided here is typically
 * within 1 or 2 bits in the mantissa.
 * However, since the error depends on the function type, floating-point spec,
 * and other factors, it is necessary to verify that sufficient precision is
 * achieved in each specific use case.
 *
 * If necessary, the order of chebyshev polynomials used and the bit width of
 * temporaries while evaluating polynomiall can be configured via
 * [[rial.math.PolynomialSpec]].
 *
 * The [[rial.math.MathFunctions]] Module takes a list of math functions to
 * support and [[rial.math.MathFuncConfig]] provides you a value of signal to
 * select which function to be executed.
 *
 * The latency and register positions can be configured via
 * [[rial.math.MathFuncPipelineConfig]]. Since it is fully pipeplined, it takes
 * 1 argument at a cycle and never stalls.
 *
 * Since the complexity depends on functions, the optimal latency depends on
 * which functions are supported. For example, sin/cos requires `x/pi` before
 * calculating its main part and log requires `* ln2` after calculating `log2(x)`.
 * Those functions automatically adds a multiplier to [[rial.math.MathFunctions]]
 * and it requires a lot more latency compared to simpler functions like sqrt,
 * invsqrt, and reciprocal. If your list of functions contains those complicated
 * functions, the optimal latency of [[rial.math.MathFunctions]] increases and
 * simpler functions also takes longer cycles.
 *
 */
package object math {}
