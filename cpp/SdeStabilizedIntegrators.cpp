#include "SdeStabilizedIntegrators.h"

SROCK2::SROCK2(Sde* sde, bool verb, bool dtadap, Real atol, 
               Real rtol, bool scalartol, int output_freq)
:OdeRungeKuttaIntegrator(sde,verb,dtadap,atol,rtol,scalartol,output_freq),
 SdeRungeKuttaIntegrator(sde,verb,dtadap,atol,rtol,scalartol,output_freq),
 DROCK2(sde,verb,dtadap,atol,rtol,scalartol,output_freq)
{
    integr5 = new Vector(sde->system_size());
    needDoubleInt = true;
    solver = "SROCK2";
}

SROCK2::~SROCK2()
{
    delete integr5;
}

void SROCK2::step(const Real t, const Real& h)
{
    Sde* sde = dynamic_cast<Sde*>(ode);
    if(sde->isDiagonal())
        this->step_diagonal_noise(t,h);
    else
        this->step_general_noise(t,h);
}

void SROCK2::step_diagonal_noise(const Real t, const Real& h)
{
    /**
     * Does an SROCK2 step.
     */
    
    //First convert ode to sde, so that we can call g and the brownian motion
    Sde* sde = dynamic_cast<Sde*>(ode);
    
    Real temp1,temp2,temp3;
    Real ci1,ci2,ci3;
    int mr,mz;
       
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
//             Initializations
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
    Vector*& yjm1= integr[0];
    Vector*& yjm2= integr[1];
    Vector*& yj= integr[2];
    Vector*& A = integr[3];
    Vector*& B = integr[4];
    Vector*& C = integr5;
    Vector* swap_ptr=0;
    
    Vector Ksm1;
    Vector fKsm1;
 
    mz=mp[0];//used to look in tables
    mr=mp[1];//as well
    const Real alpha = recalph[mz];
    const Real ha = alpha*h;
        
    temp1=ha*recf[mr];   //temp1=ha*mu_1
    ci1=t+temp1;        //ci1 = t+ha*mu_1 = t+ha*c_1
    ci2=ci1;            //ci2 = t+ha*mu_1
    ci3=t;              //ci3 = t
   
    //first step
    sde->f(t,*yn,*A);       // A = f(yn)
    *yjm2=*yn;              // yjm2 = K_0
    *yjm1=*yn+temp1*(*A);   // yjm1 = K_1
    
    //Stages for j=2...s (remember that for rock2 s is s-2)
    for(int j=2;j<=s+2;j++)
    {
        if(j<=s)
        {
            temp1= ha*recf[++mr];     //temp1 = ha*mu_j
            temp3= -recf[++mr];       //temp3 = -kappa_j
            temp2= 1.0-temp3;         //temp2 = -nu_j
            ci1=temp1+temp2*ci2+temp3*ci3; //ci1 = t+ha*c_j
        }
        else// j=s+1 or j=s+2
        {
            mr = mz*4 + 2*(j-s-1);
            temp1= ha*recf2[mr];      //temp1 = ha*mu_j
            temp3= -recf2[++mr];      //temp3 = -kappa_j
            temp2= 1.0-temp3;         //temp2 = -nu_j
            ci1=temp1+temp2*ci2+temp3*ci3; //ci1 = t+ha*c_j
        }

        //computing  yj = K_j
        sde->f(ci2,*yjm1,*yj);    // yj = f(K_{j-1)})
        *yj *= temp1;            
        *yj += temp2*(*yjm1);
        *yj += temp3*(*yjm2);     // yj = K_j

        //Shift the value "y" for the next stage
        if(j<s+2)
        {
            swap_ptr=yjm2;
            yjm2 = yjm1;    // yjm2 = K_{j-1}
            yjm1 = yj;      // yjm1 = K_j
            yj = swap_ptr;  // yj = free space previously occupied by yjm2
            
            ci3=ci2;        // ci3 = t+ha*c_{j-1}
            ci2=ci1;        // ci2 = t+ha*c_j
        }
    }
    // at this point we have yj   = K_s,     ci1 = t+ha*c_s
    //                       yjm1 = K_{s-1}, ci2 = t+ha*c_{s-1}
    //                       yjm2 = K_{s-2}, ci3 = t+ha*c_{s-2}
    
    
    // Begin of the two-stage finishing procedure
    // at the end yjm2 will contain the solution
    Real sigma, tau;
    sigma=fp1[mz]; 
    temp2=fp2[mz];  //temp2 = -sigma*(1-tau/sigma^2)
    tau = sigma*(temp2+sigma);
    
    Real sigma_a, tau_a;
    sigma_a = (1.-alpha)/2.+alpha*sigma;
    tau_a = (alpha-1.)*(alpha-1.)/2.+2.*alpha*(1.-alpha)*sigma+alpha*alpha*tau;
    
    // Stochastic integrals
    sde->sample(t,h);
    Vector& Ir = sde->getIr(); // Brownian increments
    Vector& Irr = sde->getIrr(); // double stochastic integrals
    Vector& Chir = sde->getChir();
            
    sde->g(ci1,*yj,Gv);     // G = ( g_1^1(K_s^1), g_2^2(K_s^2),...,g_r^r(K_s^r) )
    sde->f(ci3,*yjm2,*A);  // A = f(K_{s-2})
        
    *B = *yjm2;
    *B += 2.*tau_a*h*(*A);        // B = K_{s-1}^**
    sde->f(ci3+2.*tau_a*h,*B,*C); // C = f(K_{s-1}^**)
        
    // deterministic error estimation
    errD = 0.;
    for(int i=0;i<sde->system_size();i++)
    {
        temp1 = max(abs((*yn)(i)),abs((*yj)(i)));
        temp1 = temp1*r_tol+a_tol;
        errD += ((*A)(i)-(*C)(i))*((*A)(i)-(*C)(i))/temp1/temp1;
    }
    errD = abs(0.5*h*(1.-sigma_a*sigma_a/tau_a))*sqrt(errD)/sde->system_size();
    
    
    // yjm2 = K_{s-2}+(2.*sigma_a-1/2)*h*f(K_{s-2})
    *yjm2 += (2.*sigma_a-0.5)*h*(*A); 
        
    // yjm2 = K_{s-2}+(2.*sigma_a-1/2)*h*f(K_{s-2}) + first stochastic term
    *ynpu = Gv.cwiseProduct(Irr);
    *ynpu += *yj;
    sde->g(ci1,*ynpu,*A);
    *yjm2 += 0.5*(*A);
    *ynpu *= -1.;
    *ynpu += 2.*(*yj);
    sde->g(ci1,*ynpu,*A);
    *yjm2 -= 0.5*(*A);
    
    
    *yj = Gv.cwiseProduct(Ir);
    *B += *yj;                    // B = K_{s-1}^*   
    
    sde->f(ci3+2.*tau_a*h,*B,*A); // A = f(K_{s-1}^*)  
     
    *yjm2 += 0.5*h*(*A);          // yjm2 += 0.5*h*f(K_{s-1}^*) 
   
    // Computing C = 0.5*h*(f(K_{s-1}^**)-f(K_{s-1}^*)) + G*Ir
    // is part of stochastic error estimator
    *C -= *A;
    *C *= 0.5*h;
    *C += *yj;  
    
    *B = sqrt(h/2.)*Gv;//(Gv.cwiseProduct(Chir));
    *yj = *yjm1;
    *yj += *B;      // yj = K_{s-1} + \sqrt(h/2) \sum g^r(K_s)*chi_r

    sde->g(ci2,*yj,*A);
    *A *= 0.5;
    *ynpu = A->cwiseProduct(Ir);
    *yjm2 += *ynpu;
    *C -= *ynpu;
    
    *yj = *yjm1;
    *yj -= *B;      // yj = K_{s-1} - \sqrt(h)/2 \sum g^r(K_s)
    
    sde->g(ci2,*yj,*A);
    *A *= 0.5;
    *ynpu = A->cwiseProduct(Ir);
    *yjm2 += *ynpu;
    *C -= *ynpu;
    
    errS = 0.;
    for(int i=0;i<sde->system_size();i++)
    {
        temp1 = max(abs((*yn)(i)),abs((*yjm2)(i)));
        temp1 = temp1*r_tol+a_tol;
        errS += (*C)(i)*(*C)(i)/temp1/temp1;
    }
    errS = sqrt(errS)/sde->system_size();
        
    err = errD+pow(errS,4./3.); 
    err = 0.5*err;
    
    swap_ptr = ynpu;
    ynpu=yjm2;
    yjm2=swap_ptr;
    
    //update the number of right hand side evaluations
    n_f_eval=n_f_eval + s+5;
    n_g_eval= n_g_eval + 5*sde->brownian_size(); 
}  

void SROCK2::step_general_noise(const Real t, const Real& h)
{
    /**
     * Does an SROCK2 step.
     */
    
    //First convert ode to sde, so that we can call g and the brownian motion
    Sde* sde = dynamic_cast<Sde*>(ode);
    
    Real temp1,temp2,temp3;
    Real ci1,ci2,ci3;
    int mr,mz;
       
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** 
//             Initializations
// *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
    Vector*& yjm1= integr[0];
    Vector*& yjm2= integr[1];
    Vector*& yj= integr[2];
    Vector*& A = integr[3];
    Vector*& B = integr[4];
    Vector*& C = integr5;
    Vector* swap_ptr=0;
 
    mz=mp[0];//used to look in tables
    mr=mp[1];//as well
    const Real alpha = recalph[mz];
    const Real ha = alpha*h;
    
    
    temp1=ha*recf[mr];   //temp1=ha*mu_1
    ci1=t+temp1;        //ci1 = t+ha*mu_1 = t+ha*c_1
    ci2=ci1;            //ci2 = t+ha*mu_1
    ci3=t;              //ci3 = t
    
    
    //first step
    sde->f(t,*yn,*A);       // A = f(yn)
    *yjm2=*yn;              // yjm2 = K_0
    *yjm1=*yn+temp1*(*A);   // yjm1 = K_1

    //Stages for j=2...s (remember that for rock2 s is s-2)
    for(int j=2;j<=s+2;j++)
    {
        if(j<=s)
        {
            temp1= ha*recf[++mr];     //temp1 = ha*mu_j
            temp3= -recf[++mr];       //temp3 = -kappa_j
            temp2= 1.0-temp3;         //temp2 = -nu_j
            ci1=temp1+temp2*ci2+temp3*ci3; //ci1 = t+ha*c_j
        }
        else// j=s+1 or j=s+2
        {
            mr = mz*4 + 2*(j-s-1);
            temp1= ha*recf2[mr];      //temp1 = ha*mu_j
            temp3= -recf2[++mr];      //temp3 = -kappa_j
            temp2= 1.0-temp3;         //temp2 = -nu_j
            ci1=temp1+temp2*ci2+temp3*ci3; //ci1 = t+ha*c_j
        }

        //computing  yj = K_j
        sde->f(ci2,*yjm1,*yj);    // yj = f(K_{j-1)})
        *yj *= temp1;            
        *yj += temp2*(*yjm1);
        *yj += temp3*(*yjm2);     // yj = K_j

        //Shift the value "y" for the next stage
        if(j<s+2)
        {
            swap_ptr=yjm2;
            yjm2 = yjm1;    // yjm2 = K_{j-1}
            yjm1 = yj;      // yjm1 = K_j
            yj = swap_ptr;  // yj = free space previously occupied by yjm2
            
            ci3=ci2;        // ci3 = t+ha*c_{j-1}
            ci2=ci1;        // ci2 = t+ha*c_j
        }
    }
    // at this point we have yj   = K_s,     ci1 = t+ha*c_s
    //                       yjm1 = K_{s-1}, ci2 = t+ha*c_{s-1}
    //                       yjm2 = K_{s-2}, ci3 = t+ha*c_{s-2}
    
    // Begin of the two-stage finishing procedure
    // at the end yjm2 will contain the solution
    Real sigma, tau;
    sigma=fp1[mz]; 
    temp2=fp2[mz];  //temp2 = -sigma*(1-tau/sigma^2)
    tau = sigma*(temp2+sigma);
    
    Real sigma_a, tau_a;
    sigma_a = (1.-alpha)/2.+alpha*sigma;
    tau_a = (alpha-1.)*(alpha-1.)/2.+2.*alpha*(1.-alpha)*sigma+alpha*alpha*tau;
    
    // Stochastic integrals
    sde->sample(t,h);
    Vector& Ir = sde->getIr(); // Brownian increments
    vector<Vector>& Ipq = sde->getIpq(); // double stochastic integrals
    Vector& Chir = sde->getChir();
    Vector& ones = sde->getOnes();  // a vector of ones
    
    sde->g(ci1,*yj,G);     // G = ( g_1(K_s), g_2(K_s),...,g_r(K_s) )
    sde->f(ci3,*yjm2,*A);  // A = f(K_{s-2})
    
    *B = *yjm2;
    *B += 2.*tau_a*h*(*A);        // B = K_{s-1}^**
    sde->f(ci3+2.*tau_a*h,*B,*C); // C = f(K_{s-1}^**)
    
    // deterministic error estimation
    errD = 0.;
    for(int i=0;i<sde->system_size();i++)
    {
        temp1 = max(abs((*yn)(i)),abs((*yj)(i)));
        temp1 = temp1*r_tol+a_tol;
        errD += ((*A)(i)-(*C)(i))*((*A)(i)-(*C)(i))/temp1/temp1;
    }
    errD = abs(0.5*h*(1.-sigma_a*sigma_a/tau_a))*sqrt(errD)/sde->system_size();

    // yjm2 = K_{s-2}+(2.*sigma_a-1/2)*h*f(K_{s-2})
    *yjm2 += (2.*sigma_a-0.5)*h*(*A); 
        
    // yjm2 = K_{s-2}+(2.*sigma_a-1/2)*h*f(K_{s-2}) + first stochastic term
    for(int r=0;r<sde->brownian_size();r++)
    {
        *ynpu = *yj+G*Ipq[r];   // ynpu = K_s+\sum g^q(K_s)*Iqr
        sde->g(ci1,*ynpu,*A,r); // A = g^r(K_s+\sum g^q(K_s)*Iqr)
        *yjm2 += 0.5*(*A);
        *ynpu *= -1.;
        *ynpu += 2.*(*yj); // ynpu = yj-G*Ipq[r] (we want to avoid recomputation of G*Ipq[r])
        sde->g(ci1,*ynpu,*A,r);
        *yjm2 -= 0.5*(*A);
    }
        
    *yj = G*Ir; // yj = \sum g^r(K_s)Ir
    *B += *yj;                    // B = K_{s-1}^*
    sde->f(ci3+2.*tau_a*h,*B,*A); // A = f(K_{s-1}^*)
    
    *yjm2 += 0.5*h*(*A);          // yjm2 += 0.5*h*f(K_{s-1}^*) 
        
    // Computing C = 0.5*h*(f(K_{s-1}^**)-f(K_{s-1}^*)) + G*Ir
    // is part of stochastic error estimator
    *C -= *A;
    *C *= 0.5*h;
    *C += *yj;  
    
    *B = sqrt(h/2.)*G*Chir;
    *yj = *yjm1;
    *yj += *B;      // yj = K_{s-1} + \sqrt(h/2) \sum g^r(K_s)*chi_r
    
    for(int r=0;r<sde->brownian_size();r++)
    {
        sde->g(ci2,*yj,*A,r);
        *A *= 0.5*Ir(r);
        *yjm2 += *A;
        *C -= *A;
    }
        
    *yj = *yjm1;
    *yj -= *B;      // yj = K_{s-1} - \sqrt(h/2) \sum g^r(K_s)*chi_r
    
    for(int r=0;r<sde->brownian_size();r++)
    {
        sde->g(ci2,*yj,*A,r);
        *A *= 0.5*Ir(r);
        *yjm2 += *A;
        *C -= *A;
    }
    
    errS = 0.;
    for(int i=0;i<sde->system_size();i++)
    {
        temp1 = max(abs((*yn)(i)),abs((*yjm2)(i)));
        temp1 = temp1*r_tol+a_tol;
        errS += (*C)(i)*(*C)(i)/temp1/temp1;
    }
    errS = sqrt(errS)/sde->system_size();
        
    err = errD + pow(errS,4./3.);
    err = 0.5*err;
    
    swap_ptr = ynpu;
    ynpu=yjm2;
    yjm2=swap_ptr;
    
    //update the number of right hand side evaluations
    n_f_eval=n_f_eval + s+5;
    n_g_eval= n_g_eval + 5*sde->brownian_size(); 
}  