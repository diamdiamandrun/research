# Used in CM3192

def debye(mu,dmu,sigma,dsigma,em,K=None,error=None,scale=None,NP_name='None'):
    # first defining lambda
    if em == 'CuKa':
        Lambda = 1.5418e-10 # Cu K\alpha

    if K == 'None':
        K = 0.89
    elif K == 'sphere':
        K = 0.98
    
    if error == None:
        FWHM = 0.06
    elif error != None:
        FWHM = error
    
    # calculating the crystallite size
    if isinstance(mu and sigma and dmu and dsigma, (int,float)):
        dtheta,dbeta = (dmu/2)*np.pi/180, dsigma*(sigma/np.sqrt(sigma**2-FWHM**2))*np.pi/180
        beta = (np.pi/180) * np.sqrt(sigma**2-FWHM**2)
        theta = (np.pi/180) * 1/2 * mu
        D = (K * Lambda) / (beta * np.cos(theta))
        dD = np.sqrt(D**2 * ((dbeta/beta)**2 + (dtheta*np.tan(theta))**2))
        #return [D,dD]
        print(f'----------Calculated Crystallite Size of {NP_name}----------')
        print(f'D = {D:} \u00B1 {dD}')

    elif isinstance(mu and sigma and dmu and dsigma, list):
        dtheta_intm,dbeta_intm = [dmu[i]/2 * np.pi/180 for i in range(len(dmu))],[dsigma[i]*(sigma[i]/np.sqrt(sigma[i]**2-FWHM**2))*np.pi/180 for i in range(len(dsigma))]
        beta_intm = [(np.pi/180) * np.sqrt(sigma[i]**2-FWHM**2) for i in range(len(sigma))]
        theta_intm = [(np.pi/180) * mu[i] / 2 for i in range(len(mu))]
        D_intm = [(K * Lambda) / (beta_intm[i] * np.cos(theta_intm[i])) for i in range(len(beta_intm))] # also = len(theta_intm)
        #dD = [np.sqrt(D_intm**2 * (dbeta**2 + (dtheta*np.sin(i))**2)) for i in theta]
        dD = [np.sqrt(D_intm[i]**2 * ((dbeta_intm[i]/beta_intm[i])**2 + (dtheta_intm[i]*np.tan(theta_intm[i]))**2)) for i in range(len(D_intm))]

        D,beta,theta=[],[],[]
        for i in range(len(D_intm) and len(dD)):
            D.append([])
            D[i].append(D_intm[i])
            D[i].append(dD[i])

            beta.append([])
            beta[i].append(beta_intm[i])
            beta[i].append(dbeta_intm[i])

            theta.append([])
            theta[i].append(theta_intm[i])
            theta[i].append(dtheta_intm[i])

        if scale in (None,'nm','nano'):
            print(f'----------Calculated Crystallite Size of {NP_name}----------')
            for i in range(len(D)):
                print(f'Signal {i+1}: D = {D[i][0]*(1e9):.3f} \u00B1 {D[i][1]*(1e9):.3f} nm')
                print(f'\u03B8 = {theta[i][0]:.4f} \u00B1 {theta[i][1]:.4f} rad')
                print(f'\u03B2 = {beta[i][0]:.4f} \u00B1 {beta[i][1]:.4f} rad')
        
            Dbar,dDbar=0,0
            for i in range(len(D)):
                Dbar += D[i][0]
                dDbar += D[i][1]  
            Dbar,dDbar=Dbar/len(D),dDbar/len(D)
            print(f'Average Crystallite Size = {Dbar*1e9:.3f} \u00B1 {dDbar*1e9:.3f} nm')
    
        elif scale in ('AA','Ang','ang','Angstrom','Angstrom'):
            print(f'----------Calculated Crystallite Size of {NP_name}----------')
            for i in range(len(D)):
                print(f'Signal {i+1}: D = {D[i][0]*(1e10):.3f} \u00B1 {D[i][1]*(1e10):.3f} Å')
                print(f'\u03B8 = {theta[i][0]:.4f} \u00B1 {theta[i][1]:.4f} rad')
                print(f'\u03B2 = {beta[i][0]:.4f} \u00B1 {beta[i][1]:.4f} rad')
        
            Dbar,dDbar=0,0
            for i in range(len(D)):
                Dbar += D[i][0]
                dDbar += D[i][1]  
            Dbar,dDbar=Dbar/len(D),dDbar/len(D)
            print('-----------------------------------------------------------------')
            print(f'---------Average Crystallite Size = {Dbar*1e10:.3f} \u00B1 {dDbar*1e10:.3f} Å-----------')
            print('-----------------------------------------------------------------')
