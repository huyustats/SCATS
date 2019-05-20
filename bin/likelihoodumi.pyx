from libc.math cimport lgamma
from libc.math cimport exp
from libc.math cimport log
from libc.math cimport sqrt
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from math import pi
from scipy.integrate import quad
from scipy.integrate import dblquad

cdef double expit(double p):
    return 1.0 / (1 + exp(-p))

cdef double second_order_derivative(abkt_c, params_g, mu_cg, y_cg):
    cdef double alpha = abkt_c[0]
    cdef double beta = abkt_c[1]
    cdef double kappa = abkt_c[2]
    cdef double tau = abkt_c[3]
    cdef double theta_g = params_g[0]
    cdef double sigma_g = params_g[1]
    cdef double p_g = params_g[2]
    cdef double cg

    if sigma_g == 0:
        return float('inf')

    if y_cg == 0:
        cg = ((2 * p_g / (1 + exp(tau * mu_cg + kappa)) ** 3 * tau ** 2 * exp(tau * mu_cg + kappa) ** 2 - p_g / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau ** 2 * exp(tau * mu_cg + kappa) + 2 * p_g / (1 + exp(-tau * mu_cg - kappa)) ** 3 * exp(-exp(beta * mu_cg + alpha)) * tau ** 2 * exp(-tau * mu_cg - kappa) ** 2 - 2 * p_g / (1 + exp(-tau * mu_cg - kappa)) ** 2 * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - p_g / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau ** 2 * exp(-tau * mu_cg - kappa) - p_g / (1 + exp(-tau * mu_cg - kappa)) * beta ** 2 * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha)) + p_g / (1 + exp(-tau * mu_cg - kappa)) * beta ** 2 * exp(beta * mu_cg + alpha) ** 2 * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) / 2 - 2 * (-p_g / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau * exp(tau * mu_cg + kappa) + p_g / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - p_g / (1 + exp(-tau * mu_cg - kappa)) * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * (mu_cg - theta_g) / sigma_g ** 2 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) - (1 - p_g + p_g / (1 + exp(tau * mu_cg + kappa)) + p_g / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) / sigma_g ** 2 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) + 2 * (1 - p_g + p_g / (1 + exp(tau * mu_cg + kappa)) + p_g / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * (mu_cg - theta_g) ** 2 / sigma_g ** 4 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2)) / (1 - p_g + p_g / (1 + exp(tau * mu_cg + kappa)) + p_g / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * sqrt(pi * sigma_g ** 2) / exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) - ((-p_g / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau * exp(tau * mu_cg + kappa) + p_g / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - p_g / (1 + exp(-tau * mu_cg - kappa)) * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) / 2 - (1 - p_g + p_g / (1 + exp(tau * mu_cg + kappa)) + p_g / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * (mu_cg - theta_g) / sigma_g ** 2 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2)) / (1 - p_g + p_g / (1 + exp(tau * mu_cg + kappa)) + p_g / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) ** 2 * sqrt(2) * sqrt(pi * sigma_g ** 2) / exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) * (-p_g / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau * exp(tau * mu_cg + kappa) + p_g / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - p_g / (1 + exp(-tau * mu_cg - kappa)) * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha))) + 2 * ((-p_g / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau * exp(tau * mu_cg + kappa) + p_g / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - p_g / (1 + exp(-tau * mu_cg - kappa)) * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) / 2 - (1 - p_g + p_g / (1 + exp(tau * mu_cg + kappa)) + p_g / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * (mu_cg - theta_g) / sigma_g ** 2 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2)) / (1 - p_g + p_g / (1 + exp(tau * mu_cg + kappa)) + p_g / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * sqrt(pi * sigma_g ** 2) / exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) * (mu_cg - theta_g) / sigma_g ** 2
    else:
        cg = -(2 * beta ** 2 * exp(beta * mu_cg - tau * mu_cg + alpha - kappa) * sigma_g ** 2 + exp(beta * mu_cg - 2 * tau * mu_cg + alpha - 2 * kappa) * beta ** 2 * sigma_g ** 2 + tau ** 2 * exp(-tau * mu_cg - kappa) * sigma_g ** 2 + beta ** 2 * exp(beta * mu_cg + alpha) * sigma_g ** 2 + 2 * exp(-2 * tau * mu_cg - 2 * kappa) + 4 * exp(-tau * mu_cg - kappa) + 2) / sigma_g ** 2 / (1 + exp(-tau * mu_cg - kappa)) ** 2
    return cg


cdef double second_order_derivative_nob(abkt_c, params_g, mu_cg, y_cg):
    cdef double alpha = abkt_c[0]
    cdef double beta = abkt_c[1]
    cdef double kappa = abkt_c[2]
    cdef double tau = abkt_c[3]
    cdef double theta_g = params_g[0]
    cdef double sigma_g = params_g[1]
    cdef double cg

    if sigma_g == 0:
        return float('inf')

    if y_cg == 0:
        cg = ((2 / (1 + exp(tau * mu_cg + kappa)) ** 3 * tau ** 2 * exp(tau * mu_cg + kappa) ** 2 - 1 / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau ** 2 * exp(tau * mu_cg + kappa) + 2 / (1 + exp(-tau * mu_cg - kappa)) ** 3 * exp(-exp(beta * mu_cg + alpha)) * tau ** 2 * exp(-tau * mu_cg - kappa) ** 2 - 2 / (1 + exp(-tau * mu_cg - kappa)) ** 2 * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - 1 / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau ** 2 * exp(-tau * mu_cg - kappa) - 1 / (1 + exp(-tau * mu_cg - kappa)) * beta ** 2 * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha)) + 1 / (1 + exp(-tau * mu_cg - kappa)) * beta ** 2 * exp(beta * mu_cg + alpha) ** 2 * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) / 2 - 2 * (-1 / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau * exp(tau * mu_cg + kappa) + 1 / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - 1 / (1 + exp(-tau * mu_cg - kappa)) * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * (mu_cg - theta_g) / sigma_g ** 2 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) - (1 / (1 + exp(tau * mu_cg + kappa)) + 1 / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) / sigma_g ** 2 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) + 2 * (1 / (1 + exp(tau * mu_cg + kappa)) + 1 / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * (mu_cg - theta_g) ** 2 / sigma_g ** 4 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2)) / (1 / (1 + exp(tau * mu_cg + kappa)) + 1 / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * sqrt(pi * sigma_g ** 2) / exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) - ((-1 / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau * exp(tau * mu_cg + kappa) + 1 / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - 1 / (1 + exp(-tau * mu_cg - kappa)) * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) / 2 - (1 / (1 + exp(tau * mu_cg + kappa)) + 1 / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * (mu_cg - theta_g) / sigma_g ** 2 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2)) / (1 / (1 + exp(tau * mu_cg + kappa)) + 1 / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) ** 2 * sqrt(2) * sqrt(pi * sigma_g ** 2) / exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) * (-1 / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau * exp(tau * mu_cg + kappa) + 1 / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - 1 / (1 + exp(-tau * mu_cg - kappa)) * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha))) + 2 * ((-1 / (1 + exp(tau * mu_cg + kappa)) ** 2 * tau * exp(tau * mu_cg + kappa) + 1 / (1 + exp(-tau * mu_cg - kappa)) ** 2 * exp(-exp(beta * mu_cg + alpha)) * tau * exp(-tau * mu_cg - kappa) - 1 / (1 + exp(-tau * mu_cg - kappa)) * beta * exp(beta * mu_cg + alpha) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) / 2 - (1 / (1 + exp(tau * mu_cg + kappa)) + 1 / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * (pi * sigma_g ** 2) ** (-0.5) * (mu_cg - theta_g) / sigma_g ** 2 * exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2)) / (1 / (1 + exp(tau * mu_cg + kappa)) + 1 / (1 + exp(-tau * mu_cg - kappa)) * exp(-exp(beta * mu_cg + alpha))) * sqrt(2) * sqrt(pi * sigma_g ** 2) / exp(-(mu_cg - theta_g) ** 2 / sigma_g ** 2) * (mu_cg - theta_g) / sigma_g ** 2
    else:
        cg = -(2 * beta ** 2 * exp(beta * mu_cg - tau * mu_cg + alpha - kappa) * sigma_g ** 2 + exp(beta * mu_cg - 2 * tau * mu_cg + alpha - 2 * kappa) * beta ** 2 * sigma_g ** 2 + tau ** 2 * exp(-tau * mu_cg - kappa) * sigma_g ** 2 + beta ** 2 * exp(beta * mu_cg + alpha) * sigma_g ** 2 + 2 * exp(-2 * tau * mu_cg - 2 * kappa) + 4 * exp(-tau * mu_cg - kappa) + 2) / sigma_g ** 2 / (1 + exp(-tau * mu_cg - kappa)) ** 2
    return cg


cdef double log_dpois0(double log_mean):
    return -exp(log_mean)


cdef double log_dpois(double count, double log_mean):
    return count * log_mean - lgamma(long(count + 1.5)) - exp(log_mean)


cdef double log_expit(double x):
    return -log(1.0+exp(-x))


cdef double log_sum_exp(double a, double b, double c):
    cdef double max_el = max(a, b, c)
    return max_el + log(exp(a - max_el) + exp(b - max_el) + exp(c - max_el))


cdef double log_sum_exp2(double a, double b):
    cdef double max_el = max(a, b)
    return max_el + log(exp(a - max_el) + exp(b - max_el))

cdef double log_dnorm(double x, double mu, double sigma):
    if sigma == 0.0:
        if x == mu:
            return 0.0
        else:
            return -float('inf')
    else:
        return -0.918938533204672669540968854562379419803619384765625 - log(sigma) - (x-mu) * (x-mu) / sigma / sigma / 2


#### Beta distrubtion approach: NOT Completed #############################################################################
cdef double log_dbeta(double x, double alpha, double beta):
     return (alpha-1) * log(x) + (beta-1) * log(1-x) + lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta)

cdef double neg_log_single_complete_likelihood_nob_psi(double mu_cg, double psi_ce, params_e, params_g, abkt_c, y_ce1, y_ce0):
    cdef double theta_g = params_g[0]
    cdef double sigma_g = params_g[1]
    cdef double alpha_e = params_e[0]
    cdef double beta_e = params_e[1]
    cdef double a_c = abkt_c[0]
    cdef double b_c = abkt_c[1]
    
    return -(log_dpois(y_ce1, a_c + b_c * mu_cg * psi_ce) + log_dpois(y_ce0, a_c + b_c * mu_cg * (1-psi_ce) ) + log_dbeta(x=psi_ce, alpha=alpha_e, beta=beta_e))

cdef double neg_log_single_complete_likelihood_nob_psi_forminimize(double param_ce, double mu_cg, params_e, params_g, abkt_c, y_ce1, y_ce0):
    cdef double psi_ce = expit(param_ce)
    cdef double theta_g = params_g[0]
    cdef double sigma_g = params_g[1]
    cdef double alpha_e = params_e[0]
    cdef double beta_e = params_e[1]
    cdef double a_c = abkt_c[0]
    cdef double b_c = abkt_c[1]


    return -(log_dpois(y_ce1, a_c + b_c * mu_cg * psi_ce) + log_dpois(y_ce0, a_c + b_c * mu_cg * (1-psi_ce) ) + log_dbeta(x=psi_ce, alpha=alpha_e, beta=beta_e))

cdef double single_complete_likelihood_nob_psi(double psi_ce, double mu_cg, params_e, params_g, abkt_c, y_ce1, y_ce0, double scale_factor_ce):
    print "haha",-neg_log_single_complete_likelihood_nob_psi(mu_cg, psi_ce, params_e, params_g, abkt_c, y_ce1, y_ce0),scale_factor_ce,exp(-neg_log_single_complete_likelihood_nob_psi(mu_cg, psi_ce, params_e, params_g, abkt_c, y_ce1, y_ce0) + scale_factor_ce)
    return exp(-neg_log_single_complete_likelihood_nob_psi(mu_cg, psi_ce, params_e, params_g, abkt_c, y_ce1, y_ce0) + scale_factor_ce)


cdef double neg_log_single_marginal_likelihood_nob_psi(params_e, params_g, abkt_c, y_ce1, y_ce0, y_cg):
    # first get the min of the neg log-likelihood
    # use brent method
    cdef double min_val
    cdef double hessian
    cdef double lower_b
    cdef double upper_b
    min_neg_log = minimize_scalar(neg_log_single_complete_likelihood_nob_psi_forminimize, args=(9, params_e, params_g, abkt_c, y_ce1, y_ce0), method='brent')

    if min_neg_log.success:
        arg_min = min_neg_log.x
        min_val = min_neg_log.fun
        #print -neg_log_single_complete_likelihood_nob_psi(9, 0.5, params_e, params_g, abkt_c, y_ce1, y_ce0), min_val
        integral = quad(single_complete_likelihood_nob_psi, 0, 1, args = (9, params_e, params_g, abkt_c, y_ce1, y_ce0, min_val))
        print "hahaha",integral[0]
        return -(log(integral[0]) - min_val)
    else:
        return float('nan')
#############################################################################################################################

cdef double neg_log_single_complete_likelihood_nob(double mu_cg, params_g, abkt_c, y_cg):
    cdef double theta_g = params_g[0]
    cdef double sigma_g = params_g[1]
    cdef double a_c = abkt_c[0]
    cdef double b_c = abkt_c[1]
    cdef double k_c = abkt_c[2]
    cdef double t_c = abkt_c[3]
    if y_cg==0:
        return -(log_sum_exp2(log_expit(-(k_c + t_c * mu_cg)), log_expit(k_c + t_c * mu_cg) + log_dpois0(a_c + b_c * mu_cg)) + log_dnorm(x=mu_cg, mu=theta_g, sigma=sigma_g))
    else:
        return -(log_expit(k_c + t_c * mu_cg) + log_dpois(y_cg, a_c + b_c * mu_cg) + log_dnorm(x=mu_cg, mu=theta_g, sigma=sigma_g))


cdef double single_complete_likelihood_nob(double mu_cg, params_g, abkt_c, y_cg, double scale_factor_cg):

    return exp(-neg_log_single_complete_likelihood_nob(mu_cg, params_g, abkt_c, y_cg) + scale_factor_cg)




cdef double neg_log_single_marginal_likelihood_nob(params_g, abkt_c, y_cg):
    # first get the min of the neg log-likelihood
    # use brent method
    min_neg_log = minimize_scalar(neg_log_single_complete_likelihood_nob, args=(params_g, abkt_c, y_cg), method='brent')
    cdef double min_val
    cdef double hessian
    cdef double lower_b
    cdef double upper_b
    if min_neg_log.success:
        arg_min = min_neg_log.x
        min_val = min_neg_log.fun
        hessian = second_order_derivative_nob(abkt_c, params_g, arg_min, y_cg)
        lower_b = arg_min - 20 / sqrt(abs(hessian))
        upper_b = arg_min + 20 / sqrt(abs(hessian))
        integral = quad(single_complete_likelihood_nob, lower_b, upper_b, args = (params_g, abkt_c, y_cg, min_val))
        return -(log(integral[0]) - min_val)
    else:
        return float('nan')


def neg_log_sum_marginal_likelihood_nob(real_params_g, abkt, y_g):
    params_g = [real_params_g[0], exp(real_params_g[1])]
    cdef double sum_marginal_likelihood = 0
    for i in range(len(y_g)):
        sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g, abkt[i,:], y_g[i])
    return sum_marginal_likelihood


cdef double neg_log_single_complete_likelihood(double mu_cg, params_g, abkt_c, long y_cg):
    cdef double theta_g = params_g[0]
    cdef double sigma_g = params_g[1]
    cdef double p_g = params_g[2]
    cdef double a_c = abkt_c[0]
    cdef double b_c = abkt_c[1]
    cdef double k_c = abkt_c[2]
    cdef double t_c = abkt_c[3]
    if y_cg==0:
        return -(log_sum_exp(log(1-p_g), log(p_g) + log_expit(-(k_c + t_c * mu_cg)), log(p_g) + log_expit(k_c + t_c * mu_cg) + log_dpois0(a_c + b_c * mu_cg)) + log_dnorm(x=mu_cg, mu=theta_g, sigma=sigma_g))
    else:
        return -(log(p_g) + log_expit(k_c + t_c * mu_cg) + log_dpois(y_cg, a_c + b_c * mu_cg) + log_dnorm(x=mu_cg, mu=theta_g, sigma=sigma_g))

cdef double neg_log_single_complete_likelihood_umi(double mu_cg, params_g, abkt_c, long y_cg):
    cdef double theta_g = params_g[0]
    cdef double sigma_g = params_g[1]
    cdef double p_g = params_g[2]
    cdef double a_c = abkt_c[0]
    cdef double b_c = 1
    cdef double k_c = 1
    cdef double t_c = 1
    if y_cg==0:
        return -(log_sum_exp2(log(1-p_g), log(p_g) + log_dpois0(a_c + mu_cg)) + log_dnorm(x=mu_cg, mu=theta_g, sigma=sigma_g))
    else:
        return -(log(p_g) + log_dpois(y_cg, a_c + mu_cg) + log_dnorm(x=mu_cg, mu=theta_g, sigma=sigma_g))


cdef double single_complete_likelihood(double mu_cg, params_g, abkt_c, y_cg, scale_factor_cg):
    return exp(-neg_log_single_complete_likelihood(mu_cg, params_g, abkt_c, y_cg) + scale_factor_cg)

cdef double single_complete_likelihood_umi(double mu_cg, params_g, abkt_c, y_cg, scale_factor_cg):
    return exp(-neg_log_single_complete_likelihood_umi(mu_cg, params_g, abkt_c, y_cg) + scale_factor_cg)

cdef double neg_log_single_marginal_likelihood(params_g, abkt_c, y_cg):
    # first get the min of the neg log-likelihood
    # use brent method
    min_neg_log = minimize_scalar(neg_log_single_complete_likelihood, args = (params_g, abkt_c, y_cg), method='brent')
    cdef double min_val
    cdef double hessian
    cdef double lower_b
    cdef double upper_b
    if min_neg_log.success:
        arg_min = min_neg_log.x
        min_val = min_neg_log.fun
        hessian = second_order_derivative(abkt_c, params_g, arg_min, y_cg)
        lower_b = arg_min - 20 / sqrt(abs(hessian))
        upper_b = arg_min + 20 / sqrt(abs(hessian))
        integral = quad(single_complete_likelihood, lower_b, upper_b, args = (params_g, abkt_c, y_cg, min_val))
        return -(log(integral[0]) - min_val)
    else:
        return float('nan')

cdef double neg_log_single_marginal_likelihood_umi(params_g, abkt_c, y_cg):
    # first get the min of the neg log-likelihood
    # use brent method
    min_neg_log = minimize_scalar(neg_log_single_complete_likelihood_umi, args = (params_g, abkt_c, y_cg), method='brent')
    cdef double min_val
    cdef double hessian
    cdef double lower_b
    cdef double upper_b
    if min_neg_log.success:
        arg_min = min_neg_log.x
        min_val = min_neg_log.fun
        hessian = second_order_derivative(abkt_c, params_g, arg_min, y_cg)
        lower_b = arg_min - 20 / sqrt(abs(hessian))
        upper_b = arg_min + 20 / sqrt(abs(hessian))
        integral = quad(single_complete_likelihood_umi, lower_b, upper_b, args = (params_g, abkt_c, y_cg, min_val))
        return -(log(integral[0]) - min_val)
    else:
        return float('nan')


def neg_log_sum_marginal_likelihood(real_params_g, abkt, y_g):
    params_g = [real_params_g[0], exp(real_params_g[1]), expit(real_params_g[2])]
    cdef double sum_marginal_likelihood = 0
    for i in range(len(y_g)):
        sum_marginal_likelihood += neg_log_single_marginal_likelihood(params_g, abkt[i,:], y_g[i])
    return sum_marginal_likelihood

def neg_log_sum_marginal_likelihood_umi(real_params_g, abkt, y_g):
    params_g = [real_params_g[0], exp(real_params_g[1]), expit(real_params_g[2])]
    cdef double sum_marginal_likelihood = 0
    for i in range(len(y_g)):
        sum_marginal_likelihood += neg_log_single_marginal_likelihood_umi(params_g, abkt[i,:], y_g[i])
    return sum_marginal_likelihood

def neg_log_sum_marginal_likelihood_free_p(real_params_g, abkt, y_g, x_g):
    cdef double sum_marginal_likelihood = 0
    for i in range(len(y_g)):
        params_g = [real_params_g[0], exp(real_params_g[1]), expit(real_params_g[2]) * (1 - x_g[i]) + expit(real_params_g[3]) * x_g[i]]
        sum_marginal_likelihood += neg_log_single_marginal_likelihood(params_g, abkt[i,:], y_g[i])
    return sum_marginal_likelihood

def neg_log_sum_marginal_likelihood_free_theta(real_params_g, abkt, y_g, x_g):
    cdef double sum_marginal_likelihood = 0
    for i in range(len(y_g)):
        params_g = [real_params_g[0] * (1-x_g[i]) + real_params_g[1] * x_g[i], exp(real_params_g[2]), expit(real_params_g[3])]
        sum_marginal_likelihood += neg_log_single_marginal_likelihood(params_g, abkt[i,:], y_g[i])
    return sum_marginal_likelihood

def neg_log_sum_marginal_likelihood_free_theta_umi(real_params_g, abkt, y_g, x_g):
    cdef double sum_marginal_likelihood = 0
    for i in range(len(y_g)):
        params_g = [real_params_g[0] * (1-x_g[i]) + real_params_g[1] * x_g[i], exp(real_params_g[2]), expit(real_params_g[3])]
        sum_marginal_likelihood += neg_log_single_marginal_likelihood_umi(params_g, abkt[i,:], y_g[i])
    return sum_marginal_likelihood

def neg_log_sum_marginal_likelihood_free_both(real_params_g, abkt, y_g, x_g):
    cdef double sum_marginal_likelihood = 0
    for i in range(len(y_g)):
        params_g = [real_params_g[0] * (1-x_g[i]) + real_params_g[1] * x_g[i], exp(real_params_g[2]), expit(real_params_g[3]) * (1 - x_g[i])]
        sum_marginal_likelihood += neg_log_single_marginal_likelihood(params_g, abkt[i,:], y_g[i])
    return sum_marginal_likelihood

# testing BETE parameters NOT completed


def neg_log_sum_marginal_likelihood_psi_both(real_params_e, est_params_g, abkt, y_g, y_e1, y_e0, x_g):
    cdef double sum_marginal_likelihood = 0
    for i in range(len(y_g)):
        if y_e1[i] > 0 or y_e0[i] > 0:
            params_e = [exp(real_params_e[0]) * (1-x_g[i]) + exp(real_params_e[2]) * x_g[i], exp(real_params_e[1]) * (1-x_g[i]) + exp(real_params_e[3]) * x_g[i] ]
            params_g = [est_params_g[0] * (1-x_g[i]) + est_params_g[2] * x_g[i], est_params_g[1] * (1-x_g[i]) + est_params_g[3] * x_g[i]]
            #print "hahaha",neg_log_single_marginal_likelihood_nob_psi(params_e, params_g, abkt[i,:], y_e1[i], y_e0[i], y_g[i])
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob_psi(params_e, params_g, abkt[i,:], y_e1[i], y_e0[i], y_g[i])
    return sum_marginal_likelihood


# testing PSI

def neg_log_sum_marginal_likelihood_psi_equal_variance(real_params_g, abkt, y_ce1, y_ce0, x_g, theta_g1, theta_g2, sigma_g1, sigma_g2, group_status):
    cdef double sum_marginal_likelihood = 0
    cdef double psi_ce = expit(real_params_g[0])
    cdef double theta_e1_grp1 = theta_g1 * psi_ce
    cdef double theta_e0_grp1 = theta_g1 - theta_e1_grp1
    cdef double theta_e1_grp2 = theta_g2 * psi_ce
    cdef double theta_e0_grp2 = theta_g2 - theta_e1_grp2
    cdef double sigma_e1_grp1 = exp(real_params_g[1])
    cdef double sigma_e0_grp1 = exp(real_params_g[2])
    cdef double sigma_e1_grp2 = exp(real_params_g[1])
    cdef double sigma_e0_grp2 = exp(real_params_g[2])
    
    if group_status == "+\-":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g, abkt[i,:], y_ce1[i]) + neg_log_single_marginal_likelihood_nob(params_g1, abkt[i,:], y_ce0[i])
            
    if group_status == "+":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g, abkt[i,:], y_ce1[i])

    if group_status == "-":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g1, abkt[i,:], y_ce0[i])

    return sum_marginal_likelihood

def neg_log_sum_marginal_likelihood_psi_free_equal_variance(real_params_g, abkt, y_ce1, y_ce0, x_g, theta_g1, theta_g2, sigma_g1, sigma_g2, group_status):
    cdef double sum_marginal_likelihood = 0
    cdef double psi_ce_grp1 = expit(real_params_g[0])
    cdef double psi_ce_grp2 = expit(real_params_g[1])
    cdef double theta_e1_grp1 = theta_g1 * psi_ce_grp1
    cdef double theta_e0_grp1 = theta_g1 - theta_e1_grp1
    cdef double theta_e1_grp2 = theta_g2 * psi_ce_grp2
    cdef double theta_e0_grp2 = theta_g2 - theta_e1_grp2
    cdef double sigma_e1_grp1 = exp(real_params_g[2])
    cdef double sigma_e0_grp1 = exp(real_params_g[3])
    cdef double sigma_e1_grp2 = exp(real_params_g[2])
    cdef double sigma_e0_grp2 = exp(real_params_g[3])
    
    if group_status == "+\-":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g, abkt[i,:], y_ce1[i]) + neg_log_single_marginal_likelihood_nob(params_g1, abkt[i,:], y_ce0[i])
    
    if group_status == "+":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g, abkt[i,:], y_ce1[i])

    if group_status == "-":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g1, abkt[i,:], y_ce0[i])

    return sum_marginal_likelihood


### testing psi for UMI data

def neg_log_sum_marginal_likelihood_psi_equal_variance_umi(real_params_g, abkt, y_ce1, y_ce0, x_g, theta_g1, theta_g2, sigma_g1, sigma_g2, p_bursting, group_status):
    cdef double sum_marginal_likelihood = 0
    cdef double psi_ce = expit(real_params_g[0])
    cdef double theta_e1_grp1 = theta_g1 * psi_ce
    cdef double theta_e0_grp1 = theta_g1 - theta_e1_grp1
    cdef double theta_e1_grp2 = theta_g2 * psi_ce
    cdef double theta_e0_grp2 = theta_g2 - theta_e1_grp2
    cdef double sigma_e1_grp1 = exp(real_params_g[1])
    cdef double sigma_e0_grp1 = exp(real_params_g[2])
    cdef double sigma_e1_grp2 = exp(real_params_g[1])
    cdef double sigma_e0_grp2 = exp(real_params_g[2])

    if group_status == "+\-":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i], p_bursting]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i], p_bursting]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_umi(params_g, abkt[i,:], y_ce1[i]) + neg_log_single_marginal_likelihood_umi(params_g1, abkt[i,:], y_ce0[i])

    if group_status == "+":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g, abkt[i,:], y_ce1[i])

    if group_status == "-":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g1, abkt[i,:], y_ce0[i])

    return sum_marginal_likelihood



def neg_log_sum_marginal_likelihood_psi_free_equal_variance_umi(real_params_g, abkt, y_ce1, y_ce0, x_g, theta_g1, theta_g2, sigma_g1, sigma_g2, p_bursting, group_status):
    cdef double sum_marginal_likelihood = 0
    cdef double psi_ce_grp1 = expit(real_params_g[0])
    cdef double psi_ce_grp2 = expit(real_params_g[1])
    cdef double theta_e1_grp1 = theta_g1 * psi_ce_grp1
    cdef double theta_e0_grp1 = theta_g1 - theta_e1_grp1
    cdef double theta_e1_grp2 = theta_g2 * psi_ce_grp2
    cdef double theta_e0_grp2 = theta_g2 - theta_e1_grp2
    cdef double sigma_e1_grp1 = exp(real_params_g[2])
    cdef double sigma_e0_grp1 = exp(real_params_g[3])
    cdef double sigma_e1_grp2 = exp(real_params_g[2])
    cdef double sigma_e0_grp2 = exp(real_params_g[3])


    if group_status == "+\-":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i], p_bursting]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i], p_bursting]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_umi(params_g, abkt[i,:], y_ce1[i]) + neg_log_single_marginal_likelihood_umi(params_g1, abkt[i,:], y_ce0[i])

    if group_status == "+":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g, abkt[i,:], y_ce1[i])

    if group_status == "-":
        for i in range(len(y_ce1)):
            params_g = [theta_e1_grp1 * (1-x_g[i]) + (theta_e1_grp2 * x_g[i]), sigma_e1_grp1 * (1-x_g[i]) + sigma_e1_grp2 * x_g[i]]
            params_g1 = [theta_e0_grp1 * (1-x_g[i]) + (theta_e0_grp2 * x_g[i]), sigma_e0_grp1 * (1-x_g[i]) + sigma_e0_grp2 * x_g[i]]
            sum_marginal_likelihood += neg_log_single_marginal_likelihood_nob(params_g1, abkt[i,:], y_ce0[i])

    return sum_marginal_likelihood
