import celmech.disturbing_function as df
import numpy as np
import sympy as sy

def coef_E(j1,j2,j3,j4,j5,j6,nu1,nu2,nu3,nu4):
    m = j1+j3+j5
    p_i = int((j2+j4+1)/2)
    p_j = -int((j1+j3-1)/2)
    if p_i-(j2+j4+1)/2 != 0 or p_j-(j1+j3-1)/2 != 0:
        return 0
    if m != 0 and m != 1:
        return 0
    elif m == 0:
        k = 1
    else:
        k = 2
    coef_E = -k*np.math.factorial(1-m)/np.math.factorial(1+m)*df.KaulaF(1,m,p_i,nu1)*df.KaulaF(1,m,p_j,nu2)*df.HansenCoefficient_term(1,-j2-j4,-j2,nu3)*df.HansenCoefficient_term(-2,j1+j3,j1,nu4)
    return coef_E

def coef_I(j1,j2,j3,j4,j5,j6,nu1,nu2,nu3,nu4):
    m = j1+j3+j5
    pi = int((j2+j4+1)/2)
    pj = -int((j1+j3-1)/2)
    if pi-(j2+j4+1)/2 != 0 or pj-(j1+j3-1)/2 != 0:
        return 0
    if m != 0 and m != 1:
        return 0
    elif m == 0:
        k = 1
    else:
        k = 2
    coef_I = -k*np.math.factorial(1-m)/np.math.factorial(1+m)*df.KaulaF(1,m,pi,nu1)*df.KaulaF(1,m,pj,nu2)*df.HansenCoefficient_term(-2,-j2-j4,-j2,nu3)*df.HansenCoefficient_term(1,j1+j3,j1,nu4)
    return coef_I

def get_df_ctilde_coefficient_symbol_2(k1,k2,k3,k4,k5,k6,nu1,nu2,nu3,nu4,indexIn,indexOut):
    symbol_str  = r"C_{{({0}\,{1}\,{2}\,{3}\,{4}\,{5})}}".format(k1,k2,k3,k4,k5,k6)
    symbol_str += r"^{{({0}\,{1}\,{2}\,{3})}}".format(nu1,nu2,nu3,nu4)
    symbol_str += r"(\alpha_{{{0}\,{1}}})".format(indexIn,indexOut)
    return sy.symbols(symbol_str)

def funcion_perturbadora(A,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,mu_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_j,tipo=None,Alpha=None):
    R = 0
    if tipo == 'E' or tipo =='e':
        if (type(Alpha)==float):
            for k_vec,nu_vec in A:        
                CD = df.evaluate_df_coefficient_dict(df.df_coefficient_Ctilde(*k_vec, *nu_vec,  include_indirect=False), Alpha)
                R += directa(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_j,CD)
                CE = coef_E(*k_vec,*nu_vec)
                R += indirecta_externa(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_j,Alpha,CE)
        
        else:
            alpha = sy.Symbol('alpha')
            for k_vec,nu_vec in A:
                CD = get_df_ctilde_coefficient_symbol_2(*k_vec,*nu_vec,'i','j')
                R += directa(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_j,CD)
                CE = coef_E(*k_vec, *nu_vec)
                R += indirecta_externa(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_j,alpha,CE)

    elif tipo == 'I' or tipo == 'i':
        if (type(Alpha)==float):
            for k_vec,nu_vec in A:        
                CD = df.evaluate_df_coefficient_dict(df.df_coefficient_Ctilde(*k_vec, *nu_vec,  include_indirect=False), Alpha)
                R += directa(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_i,CD)
                CI = coef_I(*k_vec, *nu_vec)
                R += indirecta_interna(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_i,Alpha,CI)
        else:
            alpha = sy.Symbol('alpha')
            for k_vec,nu_vec in A:
                CD = get_df_ctilde_coefficient_symbol_2(*k_vec,*nu_vec,'i','j')
                R += directa(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_i,CD)
                CI = coef_I(*k_vec, *nu_vec)
                R += indirecta_interna(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu_i,alpha,CI)
    return R

def indirecta_externa(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu,alpha,CE):
    RE = alpha*mu/a_j*(sy.sin(i_j/2)**(abs(k_vec[4]) + 2*nu_vec[0]) * sy.sin(i_i/2)**(abs(k_vec[5]) + 2*nu_vec[1]) * e_j**(abs(k_vec[2]) + 2*nu_vec[2]) * e_i**(abs(k_vec[3]) + 2*nu_vec[3]))*CE*sy.cos(k_vec[1] * lambda_i + k_vec[0]*lambda_j+k_vec[2]*varpi_j+ k_vec[3]*varpi_i + k_vec[4]*Omega_j + k_vec[5]*Omega_i) 
    return RE

def indirecta_interna(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu,alpha,CI):
    RI = mu/(alpha**2*a_j)*(sy.sin(i_j/2)**(abs(k_vec[4]) + 2*nu_vec[0]) * sy.sin(i_i/2)**(abs(k_vec[5]) + 2*nu_vec[1]) * e_j**(abs(k_vec[2]) + 2*nu_vec[2]) * e_i**(abs(k_vec[3]) + 2*nu_vec[3]))*CI*sy.cos(k_vec[1] * lambda_i + k_vec[0]*lambda_j+k_vec[2]*varpi_j+ k_vec[3]*varpi_i + k_vec[4]*Omega_j + k_vec[5]*Omega_i) 
    return RI

def directa(k_vec,nu_vec,a_i,e_i,i_i,varpi_i,Omega_i,lambda_i,a_j,e_j,i_j,varpi_j,Omega_j,lambda_j,mu,CD):
    RD = mu/a_j*(sy.sin(i_i/2)**(abs(k_vec[4]) + 2*nu_vec[0]) * sy.sin(i_j/2)**(abs(k_vec[5]) + 2*nu_vec[1]) * e_i**(abs(k_vec[2]) + 2*nu_vec[2]) * e_j**(abs(k_vec[3]) + 2*nu_vec[3]))*CD*sy.cos(k_vec[0]*lambda_j + k_vec[1]*lambda_i + k_vec[2]*varpi_i+ k_vec[3]*varpi_j + k_vec[4]*Omega_i + k_vec[5]*Omega_j) 
    return RD