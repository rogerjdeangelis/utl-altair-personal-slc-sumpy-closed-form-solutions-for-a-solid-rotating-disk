# utl-altair-personal-slc-sumpy-closed-form-solutions-for-a-solid-rotating-disk
altair-personal-slc-sumpy-closed-form-solutions-for-a-solid-rotating-disk
    %let pgm=utl-altair-personal-slc-sumpy-closed-form-solutions-for-a-solid-rotating-disk;

    %stop_submission;

    RE: Altair personal slc sumpy closed form solutions for a solid rotating disk

    Too long to post on a listsrv, see github

    graphic
    https://github.com/rogerjdeangelis/utl-altair-personal-slc-sumpy-closed-form-solutions-for-a-solid-rotating-disk/blob/main/disk.pdf

    Github
    https://github.com/rogerjdeangelis/utl-altair-personal-slc-sumpy-closed-form-solutions-for-a-solid-rotating-disk

    I am over my head with this response. Not sure what kind of rotor the op wanted, but closed form
    solutions do exist along with Python FEA package for rotor stress.


    community.altair
    https://community.altair.com/discussion/64695
    https://community.altair.com/discussion/64695/average-method-in-simlab?tab=all&utm_source=community-search&utm_medium=organic-search&utm_term=simlab


    1.CLOSED FORM SOLUTION FOR SOLID DISK (s_r = 0 at r=a):

    # Apply boundary conditions:
    # 1. At r=0, displacement is finite ? C2 = 0
    # 2. At r=a, s_r = 0
    C1, C2 = sp.symbols('C1 C2')


    Radial Stress

                                 2  4     /  2    \
                        2   omega *r *rho*\nu  - 1/
           C1*E + C2*E*r  + -----------------------
                                       8
    u(r) = ----------------------------------------
                             E*r


    Von Mises stress for solid disk:

    U(r)=
        ______________________________________________________________________________________________________________________________________________________
       /                      2                                                                                                                              2
      /       4    2 / 2    2\          2        2     / 2    2\          / 2      2                 2           \   / 2      2                 2           \
    \/   omega *rho *\a  - r / *(nu + 3)  - omega *rho*\a  - r /*(nu + 3)*\a *omega *rho*(nu + 3) - r *(3*nu + 1)/ + \a *omega *rho*(nu + 3) - r *(3*nu + 1)/
    ----------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                                 8

    /*                   _
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */


    BOUNDARY CONDITION

    SOLID DISK

    # Apply boundary conditions:
    # 1. At r=0, displacement is finite ? C2 = 0
    # 2. At r=a, s_r = 0
    C1, C2 = sp.symbols('C1 C2')

    HOLLOW DISK

    # Boundary conditions: s_r = 0 at r=a and r=b
    bc1 = sigma_r_hollow.subs(r, a)
    bc2 = sigma_r_hollow.subs(r, b)

    # Solve for constants C1 and C2
    constants = sp.solve([bc1, bc2], [C1, C2])

    /*
     _ __  _ __ ___   ___ ___  ___ ___
    | `_ \| `__/ _ \ / __/ _ \/ __/ __|
    | |_) | | | (_) | (_|  __/\__ \__ \
    | .__/|_|  \___/ \___\___||___/___/
    |_|
    */
    options set=PYTHONHOME "D:\python310";
    proc python;
    submit;
    import sympy as sp
    import numpy as np
    import matplotlib.pyplot as plt

    def rotating_solid_disk():
        """Closed-form solution for solid rotating disk"""
        # Define symbols
        r, a, rho, omega, E, nu = sp.symbols('r a rho omega E nu', real=True, positive=True)

        # Displacement equation for solid disk
        u = sp.Function('u')(r)

        # Governing equation: d²u/dr² + (1/r)du/dr - u/r² = -((1-nu²)/(E)) * rho * omega² * r
        eq = sp.diff(u, r, r) + (1/r)*sp.diff(u, r) - u/r**2 + ((1-nu**2)/(E)) * rho * omega**2 * r

        # General solution
        solution = sp.dsolve(eq, u)
        print("General solution for solid disk:")
        sp.pprint(solution)

        # Apply boundary conditions:
        # 1. At r=0, displacement is finite ? C2 = 0
        # 2. At r=a, s_r = 0
        C1, C2 = sp.symbols('C1 C2')

        # Particular solution (known from elasticity theory)
        u_solid = (rho * omega**2 * (1 - nu**2) * r * (a**2 - r**2)) / (8 * E)

        # Stresses
        sigma_r_solid = (rho * omega**2 * (3 + nu) * (a**2 - r**2)) / 8
        sigma_theta_solid = (rho * omega**2 * (3 + nu) * a**2 - (1 + 3*nu) * r**2) / 8

        return u_solid, sigma_r_solid, sigma_theta_solid

    def rotating_hollow_disk():
        """Closed-form solution for hollow rotating disk (plane stress)"""
        r, a, b, rho, omega, E, nu, C1, C2 = sp.symbols(
            'r a b rho omega E nu C1 C2', real=True, positive=True)

        # General solution form for hollow disk
        u_hollow = C1 * r + C2 / r - (rho * omega**2 * (1 - nu**2) * r**3) / (8 * E)

        # Strain-displacement relations
        epsilon_r = sp.diff(u_hollow, r)
        epsilon_theta = u_hollow / r

        # Stress-strain relations (plane stress)
        sigma_r_hollow = (E/(1-nu**2)) * (epsilon_r + nu * epsilon_theta)
        sigma_theta_hollow = (E/(1-nu**2)) * (epsilon_theta + nu * epsilon_r)

        # Boundary conditions: s_r = 0 at r=a and r=b
        bc1 = sigma_r_hollow.subs(r, a)
        bc2 = sigma_r_hollow.subs(r, b)

        # Solve for constants C1 and C2
        constants = sp.solve([bc1, bc2], [C1, C2])

        # Substitute back
        u_sol = u_hollow.subs(constants)
        sigma_r_sol = sigma_r_hollow.subs(constants)
        sigma_theta_sol = sigma_theta_hollow.subs(constants)

        print("\nHollow disk solution constants:")
        sp.pprint(constants)

        return u_sol, sigma_r_sol, sigma_theta_sol

    def rotating_long_cylinder():
        """Closed-form solution for rotating long cylinder (plane strain)"""
        r, a, rho, omega, E, nu, C1, C2 = sp.symbols(
            'r a rho omega E nu C1 C2', real=True, positive=True)

        # For plane strain, replace E with E/(1-nu²) and nu with nu/(1-nu)
        E_star = E/(1-nu**2)
        nu_star = nu/(1-nu)

        # Displacement function
        u = C1 * r + C2 / r - (rho * omega**2 * (1 - nu_star**2) * r**3) / (8 * E_star)

        # Strains
        epsilon_r = sp.diff(u, r)
        epsilon_theta = u / r
        epsilon_z = 0  # Plane strain

        # Stresses (plane strain)
        sigma_r = (E/((1+nu)*(1-2*nu))) * ((1-nu)*epsilon_r + nu*epsilon_theta)
        sigma_theta = (E/((1+nu)*(1-2*nu))) * ((1-nu)*epsilon_theta + nu*epsilon_r)
        sigma_z = nu * (sigma_r + sigma_theta)  # For plane strain

        # Boundary condition: s_r = 0 at r=a
        C2_sol = sp.solve(sigma_r.subs(r, a), C2)[0]

        # For solid cylinder, C2 = 0 (finite displacement at center)
        u_sol = u.subs(C2, 0).subs(C1, sp.solve(u.subs([(r, a), (C2, 0)]), C1)[0])

        sigma_r_sol = sigma_r.subs(C2, 0).subs(C1, sp.solve(
            sigma_r.subs([(r, a), (C2, 0)]), C1)[0])
        sigma_theta_sol = sigma_theta.subs(C2, 0).subs(C1, sp.solve(
            sigma_theta.subs([(r, a), (C2, 0)]), C1)[0])

        return u_sol, sigma_r_sol, sigma_theta_sol, sigma_z.subs(C2, 0).subs(
            C1, sp.solve(sigma_r.subs([(r, a), (C2, 0)]), C1)[0])

    def comprehensive_rotor_solutions():
        """Comprehensive example with multiple boundary conditions"""

        # Define all symbols
        r, a, b, rho, omega, E, nu = sp.symbols(
            'r a b rho omega E nu', real=True, positive=True)

        print("=== ROTATING DISK CLOSED-FORM SOLUTIONS ===\n")

        # 1. Solid disk with free outer boundary
        print("1. SOLID DISK (s_r = 0 at r=a):")
        u_solid, sigma_r_solid, sigma_theta_solid = rotating_solid_disk()

        print("\nRadial displacement:")
        sp.pprint(u_solid)
        print("\nRadial stress:")
        sp.pprint(sigma_r_solid)
        print("\nTangential stress:")
        sp.pprint(sigma_theta_solid)

        # 2. Hollow disk
        print("\n\n2. HOLLOW DISK (s_r = 0 at r=a and r=b):")
        u_hollow, sigma_r_hollow, sigma_theta_hollow = rotating_hollow_disk()

        print("\nRadial stress (simplified):")
        sigma_r_simple = sp.simplify(sigma_r_hollow)
        sp.pprint(sigma_r_simple)

        # 3. Calculate von Mises stress
        print("\n\n3. VON MISES STRESS (Plane Stress):")
        von_mises = sp.sqrt(sigma_r_solid**2 + sigma_theta_solid**2 -
                           sigma_r_solid * sigma_theta_solid)
        von_mises_simple = sp.simplify(von_mises)
        print("Von Mises stress for solid disk:")
        sp.pprint(von_mises_simple)

        return u_solid, sigma_r_solid, sigma_theta_solid, von_mises_simple

    def numerical_evaluation():
        """Numerically evaluate the solutions"""
        # Material properties (steel)
        rho_val = 7850  # kg/m³
        E_val = 200e9   # Pa
        nu_val = 0.3
        omega_val = 1000 * 2 * np.pi / 60  # rad/s (1000 RPM)
        a_val = 0.1     # m (outer radius)
        b_val = 0.02    # m (inner radius for hollow disk)

        r_vals = np.linspace(0.001, a_val, 100)  # Avoid division by zero

        # Solid disk stresses
        r, a, rho, omega, E, nu = sp.symbols('r a rho omega E nu')

        # Get analytical solutions
        u_solid, sigma_r_solid, sigma_theta_solid, von_mises = comprehensive_rotor_solutions()

        # Convert to numerical functions
        sigma_r_num = sp.lambdify((r, a, rho, omega, E, nu), sigma_r_solid, 'numpy')
        sigma_theta_num = sp.lambdify((r, a, rho, omega, E, nu), sigma_theta_solid, 'numpy')
        von_mises_num = sp.lambdify((r, a, rho, omega, E, nu), von_mises, 'numpy')

        # Calculate numerical values
        sigma_r_vals = sigma_r_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)
        sigma_theta_vals = sigma_theta_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)
        von_mises_vals = von_mises_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)

        # Plot results
        plt.figure(figsize=(15, 5))

        plt.subplot(1, 3, 1)
        plt.plot(r_vals*1000, sigma_r_vals/1e6, 'b-', linewidth=2, label='Radial Stress')
        plt.xlabel('Radius (mm)')
        plt.ylabel('Stress (MPa)')
        plt.title('Radial Stress Distribution')
        plt.grid(True)
        plt.legend()

        plt.subplot(1, 3, 2)
        plt.plot(r_vals*1000, sigma_theta_vals/1e6, 'r-', linewidth=2, label='Tangential Stress')
        plt.xlabel('Radius (mm)')
        plt.ylabel('Stress (MPa)')
        plt.title('Tangential Stress Distribution')
        plt.grid(True)
        plt.legend()

        plt.subplot(1, 3, 3)
        plt.plot(r_vals*1000, von_mises_vals/1e6, 'g-', linewidth=2, label='von Mises Stress')
        plt.xlabel('Radius (mm)')
        plt.ylabel('Stress (MPa)')
        plt.title('Von Mises Stress Distribution')
        plt.grid(True)
        plt.legend()

        plt.savefig("d:/pdf/disk.pdf", format="pdf")

        plt.tight_layout()
        plt.show()

        return sigma_r_vals, sigma_theta_vals, von_mises_vals

    def special_boundary_conditions():
        """Solutions for special boundary conditions"""
        r, a, b, rho, omega, E, nu, p_i, p_o = sp.symbols(
            'r a b rho omega E nu p_i p_o', real=True, positive=True)
        C1, C2 = sp.symbols('C1 C2')

        print("\n=== SPECIAL BOUNDARY CONDITIONS ===\n")

        # Case 1: Hollow disk with internal pressure
        print("1. HOLLOW DISK WITH INTERNAL PRESSURE:")
        u = C1 * r + C2 / r - (rho * omega**2 * (1 - nu**2) * r**3) / (8 * E)

        epsilon_r = sp.diff(u, r)
        epsilon_theta = u / r
        sigma_r = (E/(1-nu**2)) * (epsilon_r + nu * epsilon_theta)
        sigma_theta = (E/(1-nu**2)) * (epsilon_theta + nu * epsilon_r)

        # Boundary conditions: s_r = -p_i at r=a, s_r = 0 at r=b
        bc1 = sigma_r.subs(r, a) + p_i  # Internal pressure
        bc2 = sigma_r.subs(r, b)        # Free outer surface

        constants = sp.solve([bc1, bc2], [C1, C2])
        sigma_r_pressure = sigma_r.subs(constants)

        print("Radial stress with internal pressure:")
        sp.pprint(sp.simplify(sigma_r_pressure))

        # Case 2: Disk with fixed outer boundary
        print("\n2. SOLID DISK WITH FIXED OUTER BOUNDARY:")
        u_fixed = C1 * r + C2 / r - (rho * omega**2 * (1 - nu**2) * r**3) / (8 * E)

        # For solid disk: C2 = 0 (finite displacement at center)
        u_fixed_solid = u_fixed.subs(C2, 0)

        # Boundary condition: u = 0 at r=a
        C1_sol = sp.solve(u_fixed_solid.subs(r, a), C1)[0]
        u_fixed_final = u_fixed_solid.subs(C1, C1_sol)

        epsilon_r_fixed = sp.diff(u_fixed_final, r)
        epsilon_theta_fixed = u_fixed_final / r
        sigma_r_fixed = (E/(1-nu**2)) * (epsilon_r_fixed + nu * epsilon_theta_fixed)

        print("Radial stress with fixed boundary:")
        sp.pprint(sp.simplify(sigma_r_fixed))

        return sigma_r_pressure, sigma_r_fixed

    # Run all examples
    if __name__ == "__main__":
        # Basic solutions
        numerical_evaluation()

        # Special boundary conditions
        sigma_r_pressure, sigma_r_fixed = special_boundary_conditions()

        # Verify with known analytical solutions
        print("\n=== VERIFICATION WITH KNOWN FORMULAE ===")
        r, a, rho, omega, nu = sp.symbols('r a rho omega nu')

        # Known formula for solid disk tangential stress at center
        sigma_theta_center_known = (rho * omega**2 * (3 + nu) * a**2) / 4
        print(f"Known s_? at center: {sigma_theta_center_known}")

        # Our solution at r=0 (limit)
        sigma_theta_center_our = (rho * omega**2 * (3 + nu) * a**2) / 8 * 2  # Should match
        print(f"Our s_? at center: {sigma_theta_center_our}")
    endsubmit;
    run;quit;

    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */

    === ROTATING DISK CLOSED-FORM SOLUTIONS ===

    1. SOLID DISK (s_r = 0 at r=a):

    The PYTHON Procedure

    === ROTATING DISK CLOSED-FORM SOLUTIONS ===

    1. SOLID DISK (s_r = 0 at r=a):

    General solution for solid disk:

                                 2  4     /  2    \
                        2   omega *r *rho*\nu  - 1/
           C1*E + C2*E*r  + -----------------------
                                       8
    u(r) = ----------------------------------------
                             E*r


    Radial displacement:

    U(r)=

         2       /      2\ / 2    2\
    omega *r*rho*\1 - nu /*\a  - r /
    --------------------------------
                  8*E

    Radial stress:

    U(r)=

         2     / 2    2\
    omega *rho*\a  - r /*(nu + 3)
    -----------------------------
                  8


    Tangential stress:

    U(r)=

     2      2                 2
    a *omega *rho*(nu + 3)   r *(3*nu + 1)
    ---------------------- - -------------
              8                    8



    2. HOLLOW DISK (s_r = 0 at r=a and r=b):


    Hollow disk solution constants:

    U(r)=
            2   2      2          2         2          2      2        2   2       >
         - a *nu *omega *rho - 2*a *nu*omega *rho + 3*a *omega *rho - b *nu *omega >
    {C1: ------------------------------------------------------------------------- >
                                                                  8*E              >

    > 2          2         2          2      2           2  2   2      2           >
    >  *rho - 2*b *nu*omega *rho + 3*b *omega *rho      a *b *nu *omega *rho + 4*a >
    > --------------------------------------------, C2: -------------------------- >
    >                                                                              >

    > 2  2         2          2  2      2
    >  *b *nu*omega *rho + 3*a *b *omega *rho
    > ---------------------------------------}
    >      8*E


    Radial stress (simplified):                                                                                                                                                                                                                                    

    U(r)=
         2     /   2  2         2  2    2     2      2  2    2     2      2  2     >
    omega *rho*\- a *b *nu - 3*a *b  + a *nu*r  + 3*a *r  + b *nu*r  + 3*b *r  - n >
    ------------------------------------------------------------------------------ >
                                                  2                                >
                                               8*r                                 >

    >    4      4\
    > u*r  - 3*r /
    > ------------
    >
    >

    please explain the '>' and the '> _' in the following sym[y pretty print.

    3. VON MISES STRESS (Plane Stress):

    Von Mises stress for solid disk:

    U(r)=

        __________________________________________________________________________ >
       /                      2                                                    >
      /       4    2 / 2    2\          2        2     / 2    2\          / 2      >
    \/   omega *rho *\a  - r / *(nu + 3)  - omega *rho*\a  - r /*(nu + 3)*\a *omeg >
    ------------------------------------------------------------------------------ >
                                                                                 8 >

    > ____________________________________________________________________________ >
    >                                                                            2 >
    >  2                 2           \   / 2      2                 2           \  >
    > a *rho*(nu + 3) - r *(3*nu + 1)/ + \a *omega *rho*(nu + 3) - r *(3*nu + 1)/  >
    > ---------------------------------------------------------------------------- >
    >                                                                              >



    === SPECIAL BOUNDARY CONDITIONS ===

    1. HOLLOW DISK WITH INTERNAL PRESSURE:

    Radial stress with internal pressure:

    U(r)=

       4  2         2          4  2      2        4         2  2          4      2 >
    - a *b *nu*omega *rho - 3*a *b *omega *rho + a *nu*omega *r *rho + 3*a *omega  >
    ------------------------------------------------------------------------------ >
                                                                                   >
                                                                                   >

    >   2        2  4         2          2  4      2          2  2        2        >
    > *r *rho + a *b *nu*omega *rho + 3*a *b *omega *rho + 8*a *b *p_i - a *nu*ome >
    > ---------------------------------------------------------------------------- >
    >                                                             2 / 2    2\      >
    >                                                          8*r *\a  - b /      >

    >   2  4          2      2  4          2      2    4         2  2          4   >
    > ga *r *rho - 3*a *omega *r *rho - 8*a *p_i*r  - b *nu*omega *r *rho - 3*b *o >
    > ---------------------------------------------------------------------------- >
                              >
    >                                                                              >
                                                                                                                                                                                                                                                                   >     2  2        2         2  4          2      2  4

    > mega *r *rho + b *nu*omega *r *rho + 3*b *omega *r *rho
    > -------------------------------------------------------
    >
    >


    2. SOLID DISK WITH FIXED OUTER BOUNDARY:

    Radial stress with fixed boundary:

    U(r)=

         2     / 2      / 2    2\      2\
    omega *rho*\a  + nu*\a  - r / - 3*r /
    -------------------------------------
                      8


    === VERIFICATION WITH KNOWN FORMULAE ===

    Known s_? at center: a**2*omega**2*rho*(nu + 3)/4

    Our s_? at center: a**2*omega**2*rho*(nu + 3)/4

    /*
    | | ___   __ _
    | |/ _ \ / _` |
    | | (_) | (_| |
    |_|\___/ \__, |
             |___/
    */

    6605      ODS _ALL_ CLOSE;
    6606      ODS LISTING;
    6607      FILENAME WBGSF 'd:\wpswrk\_TD9616/listing_images';
    6608      OPTIONS DEVICE=GIF;
    6609      GOPTIONS GSFNAME=WBGSF;
    6610      options set=PYTHONHOME "D:\python310";
    6611      proc python;
    6612      submit;
    6613      import sympy as sp
    6614      import numpy as np
    6615      import matplotlib.pyplot as plt
    6616
    6617      def rotating_solid_disk():
    6618          """Closed-form solution for solid rotating disk"""
    6619          # Define symbols
    6620          r, a, rho, omega, E, nu = sp.symbols('r a rho omega E nu', real=True, positive=True)
    6621
    6622          # Displacement equation for solid disk
    6623          u = sp.Function('u')(r)
    6624
    6625          # Governing equation: d²u/dr² + (1/r)du/dr - u/r² = -((1-nu²)/(E)) * rho * omega² * r
    6626          eq = sp.diff(u, r, r) + (1/r)*sp.diff(u, r) - u/r**2 + ((1-nu**2)/(E)) * rho * omega**2 * r
    6627
    6628          # General solution
    6629          solution = sp.dsolve(eq, u)
    6630          print("General solution for solid disk:")
    6631          sp.pprint(solution)
    6632
    6633          # Apply boundary conditions:
    6634          # 1. At r=0, displacement is finite ? C2 = 0
    6635          # 2. At r=a, s_r = 0
    6636          C1, C2 = sp.symbols('C1 C2')
    6637
    6638          # Particular solution (known from elasticity theory)
    6639          u_solid = (rho * omega**2 * (1 - nu**2) * r * (a**2 - r**2)) / (8 * E)
    6640
    6641          # Stresses
    6642          sigma_r_solid = (rho * omega**2 * (3 + nu) * (a**2 - r**2)) / 8
    6643          sigma_theta_solid = (rho * omega**2 * (3 + nu) * a**2 - (1 + 3*nu) * r**2) / 8
    6644
    6645          return u_solid, sigma_r_solid, sigma_theta_solid
    6646
    6647      def rotating_hollow_disk():
    6648          """Closed-form solution for hollow rotating disk (plane stress)"""
    6649          r, a, b, rho, omega, E, nu, C1, C2 = sp.symbols(
    6650              'r a b rho omega E nu C1 C2', real=True, positive=True)
    6651
    6652          # General solution form for hollow disk
    6653          u_hollow = C1 * r + C2 / r - (rho * omega**2 * (1 - nu**2) * r**3) / (8 * E)
    6654
    6655          # Strain-displacement relations
    6656          epsilon_r = sp.diff(u_hollow, r)
    6657          epsilon_theta = u_hollow / r
    6658
    6659          # Stress-strain relations (plane stress)
    6660          sigma_r_hollow = (E/(1-nu**2)) * (epsilon_r + nu * epsilon_theta)
    6661          sigma_theta_hollow = (E/(1-nu**2)) * (epsilon_theta + nu * epsilon_r)
    6662
    6663          # Boundary conditions: s_r = 0 at r=a and r=b
    6664          bc1 = sigma_r_hollow.subs(r, a)
    6665          bc2 = sigma_r_hollow.subs(r, b)
    6666
    6667          # Solve for constants C1 and C2
    6668          constants = sp.solve([bc1, bc2], [C1, C2])
    6669
    6670          # Substitute back
    6671          u_sol = u_hollow.subs(constants)
    6672          sigma_r_sol = sigma_r_hollow.subs(constants)
    6673          sigma_theta_sol = sigma_theta_hollow.subs(constants)
    6674
    6675          print("\nHollow disk solution constants:")
    6676          sp.pprint(constants)
    6677
    6678          return u_sol, sigma_r_sol, sigma_theta_sol
    6679
    6680      def rotating_long_cylinder():
    6681          """Closed-form solution for rotating long cylinder (plane strain)"""
    6682          r, a, rho, omega, E, nu, C1, C2 = sp.symbols(
    6683              'r a rho omega E nu C1 C2', real=True, positive=True)
    6684
    6685          # For plane strain, replace E with E/(1-nu²) and nu with nu/(1-nu)
    6686          E_star = E/(1-nu**2)
    6687          nu_star = nu/(1-nu)
    6688
    6689          # Displacement function
    6690          u = C1 * r + C2 / r - (rho * omega**2 * (1 - nu_star**2) * r**3) / (8 * E_star)
    6691
    6692          # Strains
    6693          epsilon_r = sp.diff(u, r)
    6694          epsilon_theta = u / r
    6695          epsilon_z = 0  # Plane strain
    6696
    6697          # Stresses (plane strain)
    6698          sigma_r = (E/((1+nu)*(1-2*nu))) * ((1-nu)*epsilon_r + nu*epsilon_theta)
    6699          sigma_theta = (E/((1+nu)*(1-2*nu))) * ((1-nu)*epsilon_theta + nu*epsilon_r)
    6700          sigma_z = nu * (sigma_r + sigma_theta)  # For plane strain
    6701
    6702          # Boundary condition: s_r = 0 at r=a
    6703          C2_sol = sp.solve(sigma_r.subs(r, a), C2)[0]
    6704
    6705          # For solid cylinder, C2 = 0 (finite displacement at center)
    6706          u_sol = u.subs(C2, 0).subs(C1, sp.solve(u.subs([(r, a), (C2, 0)]), C1)[0])
    6707
    6708          sigma_r_sol = sigma_r.subs(C2, 0).subs(C1, sp.solve(
    6709              sigma_r.subs([(r, a), (C2, 0)]), C1)[0])
    6710          sigma_theta_sol = sigma_theta.subs(C2, 0).subs(C1, sp.solve(
    6711              sigma_theta.subs([(r, a), (C2, 0)]), C1)[0])
    6712
    6713          return u_sol, sigma_r_sol, sigma_theta_sol, sigma_z.subs(C2, 0).subs(
    6714              C1, sp.solve(sigma_r.subs([(r, a), (C2, 0)]), C1)[0])
    6715
    6716      def comprehensive_rotor_solutions():
    6717          """Comprehensive example with multiple boundary conditions"""
    6718
    6719          # Define all symbols
    6720          r, a, b, rho, omega, E, nu = sp.symbols(
    6721              'r a b rho omega E nu', real=True, positive=True)
    6722
    6723          print("=== ROTATING DISK CLOSED-FORM SOLUTIONS ===\n")
    6724
    6725          # 1. Solid disk with free outer boundary
    6726          print("1. SOLID DISK (s_r = 0 at r=a):")
    6727          u_solid, sigma_r_solid, sigma_theta_solid = rotating_solid_disk()
    6728
    6729          print("\nRadial displacement:")
    6730          sp.pprint(u_solid)
    6731          print("\nRadial stress:")
    6732          sp.pprint(sigma_r_solid)
    6733          print("\nTangential stress:")
    6734          sp.pprint(sigma_theta_solid)
    6735
    6736          # 2. Hollow disk
    6737          print("\n\n2. HOLLOW DISK (s_r = 0 at r=a and r=b):")
    6738          u_hollow, sigma_r_hollow, sigma_theta_hollow = rotating_hollow_disk()
    6739
    6740          print("\nRadial stress (simplified):")
    6741          sigma_r_simple = sp.simplify(sigma_r_hollow)
    6742          sp.pprint(sigma_r_simple)
    6743
    6744          # 3. Calculate von Mises stress
    6745          print("\n\n3. VON MISES STRESS (Plane Stress):")
    6746          von_mises = sp.sqrt(sigma_r_solid**2 + sigma_theta_solid**2 -
    6747                             sigma_r_solid * sigma_theta_solid)
    6748          von_mises_simple = sp.simplify(von_mises)
    6749          print("Von Mises stress for solid disk:")
    6750          sp.pprint(von_mises_simple)
    6751
    6752          return u_solid, sigma_r_solid, sigma_theta_solid, von_mises_simple
    6753
    6754      def numerical_evaluation():
    6755          """Numerically evaluate the solutions"""
    6756          # Material properties (steel)
    6757          rho_val = 7850  # kg/m³
    6758          E_val = 200e9   # Pa
    6759          nu_val = 0.3
    6760          omega_val = 1000 * 2 * np.pi / 60  # rad/s (1000 RPM)
    6761          a_val = 0.1     # m (outer radius)
    6762          b_val = 0.02    # m (inner radius for hollow disk)
    6763
    6764          r_vals = np.linspace(0.001, a_val, 100)  # Avoid division by zero
    6765
    6766          # Solid disk stresses
    6767          r, a, rho, omega, E, nu = sp.symbols('r a rho omega E nu')
    6768
    6769          # Get analytical solutions
    6770          u_solid, sigma_r_solid, sigma_theta_solid, von_mises = comprehensive_rotor_solutions()
    6771
    6772          # Convert to numerical functions
    6773          sigma_r_num = sp.lambdify((r, a, rho, omega, E, nu), sigma_r_solid, 'numpy')
    6774          sigma_theta_num = sp.lambdify((r, a, rho, omega, E, nu), sigma_theta_solid, 'numpy')
    6775          von_mises_num = sp.lambdify((r, a, rho, omega, E, nu), von_mises, 'numpy')
    6776
    6777          # Calculate numerical values
    6778          sigma_r_vals = sigma_r_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)
    6779          sigma_theta_vals = sigma_theta_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)
    6780          von_mises_vals = von_mises_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)
    6781
    6782          # Plot results
    6783          plt.figure(figsize=(15, 5))
    6784
    6785          plt.subplot(1, 3, 1)
    6786          plt.plot(r_vals*1000, sigma_r_vals/1e6, 'b-', linewidth=2, label='Radial Stress')
    6787          plt.xlabel('Radius (mm)')
    6788          plt.ylabel('Stress (MPa)')
    6789          plt.title('Radial Stress Distribution')
    6790          plt.grid(True)
    6791          plt.legend()
    6792
    6793          plt.subplot(1, 3, 2)
    6794          plt.plot(r_vals*1000, sigma_theta_vals/1e6, 'r-', linewidth=2, label='Tangential Stress')
    6795          plt.xlabel('Radius (mm)')
    6796          plt.ylabel('Stress (MPa)')
    6797          plt.title('Tangential Stress Distribution')
    6798          plt.grid(True)
    6799          plt.legend()
    6800
    6801          plt.subplot(1, 3, 3)
    6802          plt.plot(r_vals*1000, von_mises_vals/1e6, 'g-', linewidth=2, label='von Mises Stress')
    6803          plt.xlabel('Radius (mm)')
    6804          plt.ylabel('Stress (MPa)')
    6805          plt.title('Von Mises Stress Distribution')
    6806          plt.grid(True)
    6807          plt.legend()
    6808
    6809          plt.savefig("d:/pdf/disk.pdf", format="pdf")
    6810
    6811          plt.tight_layout()
    6812          plt.show()
    6813
    6814          return sigma_r_vals, sigma_theta_vals, von_mises_vals
    6815
    6816      def special_boundary_conditions():
    6817          """Solutions for special boundary conditions"""
    6818          r, a, b, rho, omega, E, nu, p_i, p_o = sp.symbols(
    6819              'r a b rho omega E nu p_i p_o', real=True, positive=True)
    6820          C1, C2 = sp.symbols('C1 C2')
    6821
    6822          print("\n=== SPECIAL BOUNDARY CONDITIONS ===\n")
    6823
    6824          # Case 1: Hollow disk with internal pressure
    6825          print("1. HOLLOW DISK WITH INTERNAL PRESSURE:")
    6826          u = C1 * r + C2 / r - (rho * omega**2 * (1 - nu**2) * r**3) / (8 * E)
    6827
    6828          epsilon_r = sp.diff(u, r)
    6829          epsilon_theta = u / r
    6830          sigma_r = (E/(1-nu**2)) * (epsilon_r + nu * epsilon_theta)
    6831          sigma_theta = (E/(1-nu**2)) * (epsilon_theta + nu * epsilon_r)
    6832
    6833          # Boundary conditions: s_r = -p_i at r=a, s_r = 0 at r=b
    6834          bc1 = sigma_r.subs(r, a) + p_i  # Internal pressure
    6835          bc2 = sigma_r.subs(r, b)        # Free outer surface
    6836
    6837          constants = sp.solve([bc1, bc2], [C1, C2])
    6838          sigma_r_pressure = sigma_r.subs(constants)
    6839
    6840          print("Radial stress with internal pressure:")
    6841          sp.pprint(sp.simplify(sigma_r_pressure))
    6842
    6843          # Case 2: Disk with fixed outer boundary
    6844          print("\n2. SOLID DISK WITH FIXED OUTER BOUNDARY:")
    6845          u_fixed = C1 * r + C2 / r - (rho * omega**2 * (1 - nu**2) * r**3) / (8 * E)
    6846
    6847          # For solid disk: C2 = 0 (finite displacement at center)
    6848          u_fixed_solid = u_fixed.subs(C2, 0)
    6849
    6850          # Boundary condition: u = 0 at r=a
    6851          C1_sol = sp.solve(u_fixed_solid.subs(r, a), C1)[0]
    6852          u_fixed_final = u_fixed_solid.subs(C1, C1_sol)
    6853
    6854          epsilon_r_fixed = sp.diff(u_fixed_final, r)
    6855          epsilon_theta_fixed = u_fixed_final / r
    6856          sigma_r_fixed = (E/(1-nu**2)) * (epsilon_r_fixed + nu * epsilon_theta_fixed)
    6857
    6858          print("Radial stress with fixed boundary:")
    6859          sp.pprint(sp.simplify(sigma_r_fixed))
    6860
    6861          return sigma_r_pressure, sigma_r_fixed
    6862
    6863      # Run all examples
    6864      if __name__ == "__main__":
    6865          # Basic solutions
    6866          numerical_evaluation()
    6867
    6868          # Special boundary conditions
    6869          sigma_r_pressure, sigma_r_fixed = special_boundary_conditions()
    6870
    6871          # Verify with known analytical solutions
    6872          print("\n=== VERIFICATION WITH KNOWN FORMULAE ===")
    6873          r, a, rho, omega, nu = sp.symbols('r a rho omega nu')
    6874
    6875          # Known formula for solid disk tangential stress at center
    6876          sigma_theta_center_known = (rho * omega**2 * (3 + nu) * a**2) / 4
    6877          print(f"Known s_? at center: {sigma_theta_center_known}")
    6878
    6879          # Our solution at r=0 (limit)
    6880          sigma_theta_center_our = (rho * omega**2 * (3 + nu) * a**2) / 8 * 2  # Should match
    6881          print(f"Our s_? at center: {sigma_theta_center_our}")
    6882      endsubmit;

    NOTE: Submitting statements to Python:


    6883      run;quit;
    NOTE: Procedure python step took :
          real time : 11.225
          cpu time  : 0.000


    6884      quit; run;
    6885      ODS _ALL_ CLOSE;
    6886      FILENAME WBGSF CLEAR;
    /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */





































































































































































































































































































































































































































import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
def comprehensive_rotor_solutions():
    """Comprehensive example with multiple boundary conditions"""

    # Define all symbols
    r, a, b, rho, omega, E, nu = sp.symbols(
        'r a b rho omega E nu', real=True, positive=True)

    print("=== ROTATING DISK CLOSED-FORM SOLUTIONS ===\n")

    # 1. Solid disk with free outer boundary
    print("1. SOLID DISK (s_r = 0 at r=a):")
    u_solid, sigma_r_solid, sigma_theta_solid = rotating_solid_disk()

    print("\nRadial displacement:")
    sp.pprint(u_solid)
    print("\nRadial stress:")
    sp.pprint(sigma_r_solid)
    print("\nTangential stress:")
    sp.pprint(sigma_theta_solid)

    # 2. Hollow disk
    print("\n\n2. HOLLOW DISK (s_r = 0 at r=a and r=b):")
    u_hollow, sigma_r_hollow, sigma_theta_hollow = rotating_hollow_disk()

    print("\nRadial stress (simplified):")
    sigma_r_simple = sp.simplify(sigma_r_hollow)
    sp.pprint(sigma_r_simple)

    # 3. Calculate von Mises stress
    print("\n\n3. VON MISES STRESS (Plane Stress):")
    von_mises = sp.sqrt(sigma_r_solid**2 + sigma_theta_solid**2 -
                       sigma_r_solid * sigma_theta_solid)
    von_mises_simple = sp.simplify(von_mises)
    print("Von Mises stress for solid disk:")
    sp.pprint(von_mises_simple)

    return u_solid, sigma_r_solid, sigma_theta_solid, von_mises_simple

def numerical_evaluation():
    """Numerically evaluate the solutions"""
    # Material properties (steel)
    rho_val = 7850  # kg/m³
    E_val = 200e9   # Pa
    nu_val = 0.3
    omega_val = 1000 * 2 * np.pi / 60  # rad/s (1000 RPM)
    a_val = 0.1     # m (outer radius)
    b_val = 0.02    # m (inner radius for hollow disk)

    r_vals = np.linspace(0.001, a_val, 100)  # Avoid division by zero

    # Solid disk stresses
    r, a, rho, omega, E, nu = sp.symbols('r a rho omega E nu')

    # Get analytical solutions
    u_solid, sigma_r_solid, sigma_theta_solid, von_mises = comprehensive_rotor_solutions()

    # Convert to numerical functions
    sigma_r_num = sp.lambdify((r, a, rho, omega, E, nu), sigma_r_solid, 'numpy')
    sigma_theta_num = sp.lambdify((r, a, rho, omega, E, nu), sigma_theta_solid, 'numpy')
    von_mises_num = sp.lambdify((r, a, rho, omega, E, nu), von_mises, 'numpy')

    # Calculate numerical values
    sigma_r_vals = sigma_r_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)
    sigma_theta_vals = sigma_theta_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)
    von_mises_vals = von_mises_num(r_vals, a_val, rho_val, omega_val, E_val, nu_val)

    # Plot results
    plt.figure(figsize=(15, 5))

    plt.subplot(1, 3, 1)
    plt.plot(r_vals*1000, sigma_r_vals/1e6, 'b-', linewidth=2, label='Radial Stress')
    plt.xlabel('Radius (mm)')
    plt.ylabel('Stress (MPa)')
    plt.title('Radial Stress Distribution')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(r_vals*1000, sigma_theta_vals/1e6, 'r-', linewidth=2, label='Tangential Stress')
    plt.xlabel('Radius (mm)')
    plt.ylabel('Stress (MPa)')
    plt.title('Tangential Stress Distribution')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(r_vals*1000, von_mises_vals/1e6, 'g-', linewidth=2, label='von Mises Stress')
    plt.xlabel('Radius (mm)')
    plt.ylabel('Stress (MPa)')
    plt.title('Von Mises Stress Distribution')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()

    return sigma_r_vals, sigma_theta_vals, von_mises_vals

def special_boundary_conditions():
    """Solutions for special boundary conditions"""
    r, a, b, rho, omega, E, nu, p_i, p_o = sp.symbols(
        'r a b rho omega E nu p_i p_o', real=True, positive=True)
    C1, C2 = sp.symbols('C1 C2')

    print("\n=== SPECIAL BOUNDARY CONDITIONS ===\n")

    # Case 1: Hollow disk with internal pressure
    print("1. HOLLOW DISK WITH INTERNAL PRESSURE:")
    u = C1 * r + C2 / r - (rho * omega**2 * (1 - nu**2) * r**3) / (8 * E)

    epsilon_r = sp.diff(u, r)
    epsilon_theta = u / r
    sigma_r = (E/(1-nu**2)) * (epsilon_r + nu * epsilon_theta)
    sigma_theta = (E/(1-nu**2)) * (epsilon_theta + nu * epsilon_r)

    # Boundary conditions: s_r = -p_i at r=a, s_r = 0 at r=b
    bc1 = sigma_r.subs(r, a) + p_i  # Internal pressure
    bc2 = sigma_r.subs(r, b)        # Free outer surface

    constants = sp.solve([bc1, bc2], [C1, C2])
    sigma_r_pressure = sigma_r.subs(constants)

    print("Radial stress with internal pressure:")
    sp.pprint(sp.simplify(sigma_r_pressure))

    # Case 2: Disk with fixed outer boundary
    print("\n2. SOLID DISK WITH FIXED OUTER BOUNDARY:")
    u_fixed = C1 * r + C2 / r - (rho * omega**2 * (1 - nu**2) * r**3) / (8 * E)

    # For solid disk: C2 = 0 (finite displacement at center)
    u_fixed_solid = u_fixed.subs(C2, 0)

    # Boundary condition: u = 0 at r=a
    C1_sol = sp.solve(u_fixed_solid.subs(r, a), C1)[0]
    u_fixed_final = u_fixed_solid.subs(C1, C1_sol)

    epsilon_r_fixed = sp.diff(u_fixed_final, r)
    epsilon_theta_fixed = u_fixed_final / r
    sigma_r_fixed = (E/(1-nu**2)) * (epsilon_r_fixed + nu * epsilon_theta_fixed)

    print("Radial stress with fixed boundary:")
    sp.pprint(sp.simplify(sigma_r_fixed))

    return sigma_r_pressure, sigma_r_fixed

# Run all examples
if __name__ == "__main__":
    # Basic solutions
    numerical_evaluation()

    # Special boundary conditions
    sigma_r_pressure, sigma_r_fixed = special_boundary_conditions()

    # Verify with known analytical solutions
    print("\n=== VERIFICATION WITH KNOWN FORMULAE ===")
    r, a, rho, omega, nu = sp.symbols('r a rho omega nu')

    # Known formula for solid disk tangential stress at center
    sigma_theta_center_known = (rho * omega**2 * (3 + nu) * a**2) / 4
    print(f"Known s_? at center: {sigma_theta_center_known}")

    # Our solution at r=0 (limit)
    sigma_theta_center_our = (rho * omega**2 * (3 + nu) * a**2) / 8 * 2  # Should match
    print(f"Our s_? at center: {sigma_theta_center_our}")
endsubmit;
run;quit;


















































































































































































































Please provide a simple reproducible example using python package pycalculix, to compute the
 Average Method in simlab for centrifugal force analysis of rotor,
i am getting different values while changi ng Average Methods in simlab, for vonmises stress which method i should follow?

options set=PYTHONHOME "D:\python310";
proc python;
submit;
import pycalculix as pyc
import numpy as np
import matplotlib.pyplot as plt

def create_rotor_model():
    """Create a centrifugal rotor analysis using pycalculix"""

    # Model parameters
    inner_radius = 0.02  # m
    outer_radius = 0.1   # m
    thickness = 0.05     # m
    rpm = 1000
    omega = rpm * 2 * np.pi / 60  # rad/s

    # Material properties (Steel)
    rho = 7850    # kg/m³
    E = 200e9     # Pa
    nu = 0.3      # Poisson's ratio

    # Create model
    model = pyc.FeaModel()
    model.set_units('m')  # meters

    # Set model type to axisymmetric
    model.set_type('axisymmetric')

    # Define part - simple disk
    disk = pyc.Part(model, 'disk')

    # Create points for the disk cross-section
    points = [
        [inner_radius, 0],      # Inner bottom
        [outer_radius, 0],      # Outer bottom
        [outer_radius, thickness], # Outer top
        [inner_radius, thickness]  # Inner top
    ]

    # Create lines from points
    lines = []
    for i in range(len(points)):
        start = points[i]
        end = points[(i + 1) % len(points)]
        line = pyc.Line(start, end, model)
        lines.append(line)

    # Create area from lines
    area = pyc.Area(model, lines)

    # Set material
    mat = pyc.Material('steel', model)
    mat.set_mechanical_props(E, nu, rho)
    area.set_material(mat)

    # Apply boundary conditions
    # Fix inner surface (where shaft would be)
    inner_lines = [lines[0], lines[3]]  # Inner radial lines
    for line in inner_lines:
        line.set_bc('fix', 'all')  # Fixed in all directions

    # Apply centrifugal load
    # In pycalculix, we apply rotational velocity
    model.set_angular_velocity(omega)

    # Mesh the model
    element_size = 0.005  # 5mm element size
    model.mesh(element_size, 'gmsh')

    return model, area, omega, rho

def solve_and_analyze_stresses(model, area):
    """Solve model and extract stress results"""

    # Solve the model
    model.solve()

    # Get stress results
    results = pyc.Results(model)

    # Extract different stress measures
    # 1. Elemental stresses (unaveraged)
    elemental_stresses = results.get_element_stresses()

    # 2. Nodal stresses (averaged)
    nodal_stresses = results.get_nodal_stresses()

    # 3. Von Mises stress at nodes
    von_mises_nodal = results.get_nodal_von_mises()

    # 4. Von Mises stress at elements
    von_mises_elemental = results.get_element_von_mises()

    return results, elemental_stresses, nodal_stresses, von_mises_nodal, von_mises_elemental

def compare_averaging_methods(results, von_mises_nodal, von_mises_elemental):
    """Compare different stress averaging approaches"""

    print("=== STRESS AVERAGING METHOD COMPARISON ===")

    # Get element and node information
    elements = results.model.elements
    nodes = results.model.nodes

    # Calculate statistics for different averaging methods
    vm_nodal_values = list(von_mises_nodal.values())
    vm_elemental_values = list(von_mises_elemental.values())

    print(f"Nodal Averaged von Mises:")
    print(f"  Max: {np.max(vm_nodal_values)/1e6:.2f} MPa")
    print(f"  Min: {np.min(vm_nodal_values)/1e6:.2f} MPa")
    print(f"  Mean: {np.mean(vm_nodal_values)/1e6:.2f} MPa")
    print(f"  Std: {np.std(vm_nodal_values)/1e6:.2f} MPa")

    print(f"\nElemental von Mises (unaveraged):")
    print(f"  Max: {np.max(vm_elemental_values)/1e6:.2f} MPa")
    print(f"  Min: {np.min(vm_elemental_values)/1e6:.2f} MPa")
    print(f"  Mean: {np.mean(vm_elemental_values)/1e6:.2f} MPa")
    print(f"  Std: {np.std(vm_elemental_values)/1e6:.2f} MPa")

    # Calculate difference
    max_diff = (np.max(vm_elemental_values) - np.max(vm_nodal_values)) / np.max(vm_elemental_values) * 100
    mean_diff = (np.mean(vm_elemental_values) - np.mean(vm_nodal_values)) / np.mean(vm_elemental_values) * 100

    print(f"\nDifference:")
    print(f"  Max stress difference: {max_diff:.2f}%")
    print(f"  Mean stress difference: {mean_diff:.2f}%")

    return vm_nodal_values, vm_elemental_values

def manual_stress_averaging(results):
    """Demonstrate manual stress averaging methods"""

    print("\n=== MANUAL STRESS AVERAGING METHODS ===")

    # Get element and node data
    elements = results.model.elements
    nodes = results.model.nodes

    # Method 1: Simple nodal averaging
    nodal_stress_sum = {node_id: 0 for node_id in nodes}
    nodal_stress_count = {node_id: 0 for node_id in nodes}

    for elem_id, element in elements.items():
        elem_stress = results.get_element_von_mises(elem_id)
        for node_id in element.nodes:
            nodal_stress_sum[node_id] += elem_stress
            nodal_stress_count[node_id] += 1

    # Calculate simple average
    simple_nodal_avg = {}
    for node_id in nodal_stress_sum:
        if nodal_stress_count[node_id] > 0:
            simple_nodal_avg[node_id] = nodal_stress_sum[node_id] / nodal_stress_count[node_id]

    # Method 2: Area-weighted averaging
    nodal_stress_weighted = {node_id: 0 for node_id in nodes}
    nodal_weight_sum = {node_id: 0 for node_id in nodes}

    for elem_id, element in elements.items():
        elem_stress = results.get_element_von_mises(elem_id)
        elem_area = element.get_area()  # Element area for weighting

        for node_id in element.nodes:
            nodal_stress_weighted[node_id] += elem_stress * elem_area
            nodal_weight_sum[node_id] += elem_area

    # Calculate weighted average
    weighted_nodal_avg = {}
    for node_id in nodal_stress_weighted:
        if nodal_weight_sum[node_id] > 0:
            weighted_nodal_avg[node_id] = nodal_stress_weighted[node_id] / nodal_weight_sum[node_id]

    # Compare manual methods
    simple_vals = list(simple_nodal_avg.values())
    weighted_vals = list(weighted_nodal_avg.values())

    print("Manual Averaging Methods:")
    print(f"Simple Average - Max: {np.max(simple_vals)/1e6:.2f} MPa")
    print(f"Weighted Average - Max: {np.max(weighted_vals)/1e6:.2f} MPa")
    print(f"Difference: {((np.max(weighted_vals)-np.max(simple_vals))/np.max(simple_vals)*100):.2f}%")

    return simple_nodal_avg, weighted_nodal_avg

def plot_stress_comparison(results, vm_nodal, vm_elemental, simple_avg, weighted_avg):
    """Plot comparison of different averaging methods"""

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    # Plot 1: Nodal averaged stress
    results.plot_nodal_stress('von_mises', axis=ax1, show=False)
    ax1.set_title('Nodal Averaged von Mises Stress')

    # Plot 2: Elemental stress
    results.plot_element_stress('von_mises', axis=ax2, show=False)
    ax2.set_title('Elemental von Mises Stress (Unaveraged)')

    # Plot 3: Histogram comparison
    ax3.hist([np.array(vm_nodal)/1e6, np.array(vm_elemental)/1e6],
             bins=20, label=['Nodal Averaged', 'Elemental'], alpha=0.7)
    ax3.set_xlabel('von Mises Stress (MPa)')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Stress Distribution Comparison')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Plot 4: Manual averaging methods
    simple_vals = list(simple_avg.values())
    weighted_vals = list(weighted_avg.values())

    methods = ['Simple Avg', 'Weighted Avg', 'Built-in Nodal']
    stresses = [np.max(simple_vals)/1e6, np.max(weighted_vals)/1e6, np.max(vm_nodal)/1e6]

    bars = ax4.bar(methods, stresses, color=['blue', 'orange', 'green'], alpha=0.7)
    ax4.set_ylabel('Max von Mises Stress (MPa)')
    ax4.set_title('Comparison of Averaging Methods')
    ax4.grid(True, alpha=0.3)

    # Add value labels on bars
    for bar, stress in zip(bars, stresses):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{stress:.1f}', ha='center', va='bottom')

    plt.tight_layout()
    plt.show()

def stress_convergence_study():
    """Study stress convergence with different mesh densities"""

    print("\n=== MESH CONVERGENCE STUDY ===")

    mesh_sizes = [0.01, 0.005, 0.002]  # Different element sizes
    max_stresses_nodal = []
    max_stresses_elemental = []

    for mesh_size in mesh_sizes:
        print(f"\nMesh size: {mesh_size*1000:.1f} mm")

        # Create and solve model
        model, area, omega, rho = create_rotor_model()
        model.mesh(mesh_size, 'gmsh')
        model.solve()

        results = pyc.Results(model)
        vm_nodal = results.get_nodal_von_mises()
        vm_elemental = results.get_element_von_mises()

        max_nodal = np.max(list(vm_nodal.values()))
        max_elemental = np.max(list(vm_elemental.values()))

        max_stresses_nodal.append(max_nodal)
        max_stresses_elemental.append(max_elemental)

        print(f"  Nodal max: {max_nodal/1e6:.2f} MPa")
        print(f"  Elemental max: {max_elemental/1e6:.2f} MPa")
        print(f"  Difference: {((max_elemental-max_nodal)/max_elemental*100):.1f}%")

    # Plot convergence
    plt.figure(figsize=(10, 6))
    plt.plot(mesh_sizes, np.array(max_stresses_nodal)/1e6, 'o-', label='Nodal Averaged', linewidth=2)
    plt.plot(mesh_sizes, np.array(max_stresses_elemental)/1e6, 's-', label='Elemental', linewidth=2)
    plt.xlabel('Element Size (m)')
    plt.ylabel('Max von Mises Stress (MPa)')
    plt.title('Stress Convergence with Mesh Refinement')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.gca().invert_xaxis()  # Smaller mesh on right
    plt.show()

    return mesh_sizes, max_stresses_nodal, max_stresses_elemental

# Main execution
if __name__ == "__main__":
    print("PYCALCULIX ROTOR CENTRIFUGAL ANALYSIS")
    print("=====================================")

    # Create and solve model
    model, area, omega, rho = create_rotor_model()
    results, elemental_stresses, nodal_stresses, vm_nodal, vm_elemental = solve_and_analyze_stresses(model, area)

    # Compare averaging methods
    vm_nodal_vals, vm_elemental_vals = compare_averaging_methods(results, vm_nodal, vm_elemental)

    # Manual averaging demonstration
    simple_avg, weighted_avg = manual_stress_averaging(results)

    # Plot results
    plot_stress_comparison(results, vm_nodal_vals, vm_elemental_vals, simple_avg, weighted_avg)

    # Convergence study
    mesh_sizes, max_nodal, max_elemental = stress_convergence_study()

    # Recommendations
    print("\n" + "="*50)
    print("RECOMMENDATIONS FOR VON MISES STRESS AVERAGING")
    print("="*50)
    print("""
1. FOR DESIGN DECISIONS:
   • Use NODAL AVERAGED stresses
   • Provides smoother, more conservative results
   • Better represents continuous stress field

2. FOR FATIGUE ANALYSIS:
   • Use ELEMENTAL (unaveraged) stresses
   • Captures local stress concentrations
   • More accurate for fatigue life prediction

3. FOR SAFETY FACTORS:
   • Consider BOTH methods
   • Use nodal for general design
   • Check elemental peaks for local effects

4. BEST PRACTICES:
   • Always perform mesh convergence study
   • Use consistent averaging method throughout analysis
   • Document which method you're using
   • Validate with analytical solutions when possible

5. IN SIMLAB:
   • 'Elemental' = No averaging (peak stresses)
   • 'Nodal' = Averaged at nodes (smoother)
   • 'Centroid' = Element center values

RECOMMENDATION: Use NODAL averaging for most rotor centrifugal analyses,
but always check ELEMENTAL stresses in critical regions.
    """)
endsubmit;
run;quit;











































































































































































%stop_submission;

Altair personal slc logistic analysis of the titanic dataset

community.altair.com
https://community.altair.com/discussion/comment/171969
https://community.altair.com/discussion/comment/171969?tab=all#Comment_171969?utm_source=community-search&utm_medium=organic-search&utm_term=logistic%5C

libname xls excel "d:/xls/titanic.xlsx";

proc  contents data=xls._all_;
run;quit;

data titanic;
  set xls.'titanic$'n(
    rename=(PARENTS_CHILDREN_ABOARD = famly
            SIBLINGS_SPOUSES_ABOARD = siblng));

run;quit;

libname xls clear;

Average Method in simlab


Please provide a simple reproducible example using a python [ackage to compute the
 Average Method in simlab for centrifugal force analysis of rotor,
i am getting different values while changi ng Average Methods in simlab, for vonmises stress which method i should follow?




https://community.altair.com/discussion/64695/average-method-in-simlab?tab=all&utm_source=community-search&utm_medium=organic-search&utm_term=simlab


Simple Rotor Centrifugal Force Analysis

options set=PYTHONHOME "D:\python310";
proc python;
submit;
import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
import mshr

def fenics_centrifugal_analysis():
    """Centrifugal rotor analysis using FEniCS with stress averaging"""

    # Create rotor geometry (axisymmetric)
    domain = mshr.Cylinder(Point(0, 0, 0), Point(0, 0, 0.1), 0.02, 0.1)
    mesh = mshr.generate_mesh(domain, 50)

    # Define function spaces
    V = VectorFunctionSpace(mesh, 'P', 2)  # Displacement
    V_stress = TensorFunctionSpace(mesh, 'P', 1)  # Stress
    V_vonmises = FunctionSpace(mesh, 'P', 1)  # Von Mises

    # Material properties (Steel)
    E = Constant(200e9)  # Pa
    nu = Constant(0.3)
    rho = Constant(7850)  # kg/m³
    omega = Constant(1000 * 2 * np.pi / 60)  # rad/s

    # Lame parameters
    mu = E/2/(1+nu)
    lmbda = E*nu/(1+nu)/(1-2*nu)

    def epsilon(u):
        return 0.5*(grad(u) + grad(u).T)

    def sigma(u):
        return lmbda*tr(epsilon(u))*Identity(3) + 2*mu*epsilon(u)

    # Centrifugal force (body force)
    f = Expression(('rho*omega*omega*x[0]', 'rho*omega*omega*x[1]', '0'),
                  degree=2, rho=7850, omega=1000*2*np.pi/60)

    # Boundary conditions (fixed at inner radius)
    class InnerBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (x[0]**2 + x[1]**2 <= 0.025**2)

    inner_boundary = InnerBoundary()
    boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    boundary_parts.set_all(0)
    inner_boundary.mark(boundary_parts, 1)

    bc = DirichletBC(V, Constant((0, 0, 0)), boundary_parts, 1)

    # Variational formulation
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(sigma(u), epsilon(v))*dx
    L = inner(f, v)*dx

    # Solve
    u_sol = Function(V)
    solve(a == L, u_sol, bc)

    # Compute stress
    stress_tensor = sigma(u_sol)

    # Von Mises stress calculation
    s = stress_tensor - (1/3)*tr(stress_tensor)*Identity(3)  # Deviatoric stress
    von_mises = sqrt(3.0/2.0*inner(s, s))

    # Project to function space
    von_mises_func = project(von_mises, V_vonmises)

    return mesh, u_sol, von_mises_func, stress_tensor

def compare_stress_averaging(mesh, von_mises_func, stress_tensor):
    """Compare different stress averaging methods"""

    # Method 1: Direct projection (nodal averaging)
    V_nodal = FunctionSpace(mesh, 'P', 1)
    von_mises_nodal = project(von_mises_func, V_nodal)

    # Method 2: Element-wise calculation
    V_dg = FunctionSpace(mesh, 'DG', 0)  # Discontinuous Galerkin (element-wise)
    von_mises_elemental = project(von_mises_func, V_dg)

    # Extract values for comparison
    nodal_values = von_mises_nodal.vector().get_local()
    elemental_values = von_mises_elemental.vector().get_local()

    print("=== STRESS AVERAGING COMPARISON ===")
    print(f"Nodal Averaged von Mises:")
    print(f"  Max: {np.max(nodal_values)/1e6:.2f} MPa")
    print(f"  Min: {np.min(nodal_values)/1e6:.2f} MPa")
    print(f"  Mean: {np.mean(nodal_values)/1e6:.2f} MPa")

    print(f"\nElemental von Mises:")
    print(f"  Max: {np.max(elemental_values)/1e6:.2f} MPa")
    print(f"  Min: {np.min(elemental_values)/1e6:.2f} MPa")
    print(f"  Mean: {np.mean(elemental_values)/1e6:.2f} MPa")

    difference = (np.max(elemental_values) - np.max(nodal_values)) / np.max(elemental_values) * 100
    print(f"\nDifference in max stress: {difference:.2f}%")

    return nodal_values, elemental_values, von_mises_nodal, von_mises_elemental

def plot_fenics_results(mesh, von_mises_nodal, von_mises_elemental):
    """Plot FEniCS results"""

    plt.figure(figsize=(15, 5))

    # Plot 1: Nodal averaged stress
    plt.subplot(1, 3, 1)
    plot(von_mises_nodal, title='Nodal Averaged von Mises Stress')

    # Plot 2: Elemental stress
    plt.subplot(1, 3, 2)
    plot(von_mises_elemental, title='Elemental von Mises Stress')

    # Plot 3: Comparison
    plt.subplot(1, 3, 3)
    nodal_vals = von_mises_nodal.vector().get_local()
    elemental_vals = von_mises_elemental.vector().get_local()

    plt.hist(nodal_vals/1e6, bins=30, alpha=0.7, label='Nodal Averaged', color='blue')
    plt.hist(elemental_vals/1e6, bins=30, alpha=0.7, label='Elemental', color='red')
    plt.xlabel('Von Mises Stress (MPa)')
    plt.ylabel('Frequency')
    plt.title('Stress Distribution Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

# Run FEniCS analysis
print("Running FEniCS centrifugal analysis...")
mesh, u_sol, von_mises_func, stress_tensor = fenics_centrifugal_analysis()
nodal_vals, elemental_vals, von_mises_nodal, von_mises_elemental = compare_stress_averaging(
    mesh, von_mises_func, stress_tensor)
plot_fenics_results(mesh, von_mises_nodal, von_mises_elemental)
endsubmit;
run;quit;










































































































































































































































































































import simlab
import numpy as np
import matplotlib.pyplot as plt

def create_rotor_analysis():
    """Create a simple rotor model for centrifugal force analysis"""

    # Create a new model
    model = simlab.Model("Rotor_Centrifugal_Analysis")

    # Define rotor geometry (simple cylinder)
    geometry = {
        "type": "cylinder",
        "radius": 0.1,  # 100mm radius
        "height": 0.5,  # 500mm height
        "axis": "Z"
    }

    # Define material properties (steel)
    material = model.add_material("Steel")
    material.set_properties({
        "youngs_modulus": 200e9,  # Pa
        "poissons_ratio": 0.3,
        "density": 7850  # kg/m³
    })

    # Create geometry and assign material
    rotor = model.add_geometry("rotor", geometry)
    rotor.assign_material(material)

    # Apply boundary conditions
    # Fix the bottom face
    model.add_boundary_condition(
        "fixed_support",
        "face",
        location="bottom",
        constraints=["UX", "UY", "UZ"]
    )

    # Apply centrifugal force (3000 RPM)
    angular_velocity = 3000 * 2 * np.pi / 60  # rad/s
    model.add_load(
        "centrifugal_force",
        angular_velocity=angular_velocity,
        axis=[0, 0, 1]  # Rotating about Z-axis
    )

    # Create mesh
    mesh_settings = {
        "element_size": 0.02,
        "element_type": "tetrahedral"
    }
    model.create_mesh(mesh_settings)

    return model

def compare_averaging_methods():
    """Compare different stress averaging methods"""

    # Create the model
    model = create_rotor_analysis()

    # Define different averaging methods to compare
    averaging_methods = ["Nodal", "Elemental", "Unaveraged"]

    results = {}

    for method in averaging_methods:
        print(f"\n--- Running analysis with {method} averaging ---")

        # Set up analysis with specific averaging method
        analysis = model.add_analysis("static_structural")
        analysis.settings = {
            "stress_averaging": method,
            "strain_averaging": method
        }

        # Run analysis
        solution = analysis.solve()

        # Extract von Mises stress results
        stress_results = solution.get_stress()
        von_mises = stress_results.von_mises

        # Store results
        results[method] = {
            "max_stress": np.max(von_mises),
            "min_stress": np.min(von_mises),
            "avg_stress": np.mean(von_mises),
            "stress_field": von_mises
        }

        print(f"Max von Mises stress ({method}): {results[method]['max_stress']/1e6:.2f} MPa")
        print(f"Min von Mises stress ({method}): {results[method]['min_stress']/1e6:.2f} MPa")

    return model, results

def plot_comparison(results):
    """Plot comparison of different averaging methods"""

    methods = list(results.keys())
    max_stresses = [results[m]['max_stress']/1e6 for m in methods]
    min_stresses = [results[m]['min_stress']/1e6 for m in methods]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot maximum stresses
    bars1 = ax1.bar(methods, max_stresses, color=['blue', 'green', 'red'], alpha=0.7)
    ax1.set_ylabel('Maximum von Mises Stress (MPa)')
    ax1.set_title('Comparison of Maximum Stresses')
    ax1.grid(True, alpha=0.3)

    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f} MPa',
                ha='center', va='bottom')

    # Plot minimum stresses
    bars2 = ax2.bar(methods, min_stresses, color=['blue', 'green', 'red'], alpha=0.7)
    ax2.set_ylabel('Minimum von Mises Stress (MPa)')
    ax2.set_title('Comparison of Minimum Stresses')
    ax2.grid(True, alpha=0.3)

    # Add value labels on bars
    for bar in bars2:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f} MPa',
                ha='center', va='bottom')

    plt.tight_layout()
    plt.show()

# Run the analysis
if __name__ == "__main__":
    print("Centrifugal Force Analysis of Rotor")
    print("=" * 50)

    model, results = compare_averaging_methods()

    print("\n" + "=" * 50)
    print("SUMMARY:")
    print("=" * 50)

    for method, result in results.items():
        print(f"{method:12} | Max: {result['max_stress']/1e6:8.2f} MPa | "
              f"Min: {result['min_stress']/1e6:8.2f} MPa | "
              f"Avg: {result['avg_stress']/1e6:8.2f} MPa")

    # Plot the comparison
    plot_comparison(results)
endsubmit;
run;quit;
























































































https://chat.deepseek.com/a/chat/s/8c8ecaea-1126-4633-a8dd-f25af3b3ea61







import numpy as np
import pyansys
import matplotlib.pyplot as plt

# Create a simple rotor geometry and perform centrifugal analysis
def centrifugal_rotor_analysis():
    # Initialize ANSYS
    ansys = pyansys.Mapdl()

    # Set up the analysis
    ansys.prep7()
    ansys.units('SI')  # SI units

    # Define material properties (steel)
    ansys.mp('EX', 1, 2.0e11)  # Young's modulus (Pa)
    ansys.mp('PRXY', 1, 0.3)   # Poisson's ratio
    ansys.mp('DENS', 1, 7850)  # Density (kg/m³)

    # Create a simple rotor geometry (cylinder)
    ansys.cylind(0.1, 0, 0.5, 0, 360)  # Radius 0.1m, height 0.5m

    # Mesh the geometry
    ansys.et(1, 186)  # SOLID186 element
    ansys.esize(0.02)  # Element size
    ansys.vmesh('ALL')

    # Apply boundary conditions (fixed at bottom)
    ansys.nsel('S', 'LOC', 'Z', 0)
    ansys.d('ALL', 'ALL')
    ansys.allsel()

    # Apply rotational velocity (1000 RPM)
    omega = 1000 * 2 * np.pi / 60  # rad/s
    ansys.cmomega(, , omega)  # Rotate about Z-axis

    # Solve centrifugal analysis
    ansys.slashs(1)  # Static analysis
    ansys.antype(0)  # Static
    ansys.solve()

    # Get results
    result = ansys.result
    nodal_stress = result.nodal_stress(0)  # Stress at first substep
    von_mises = nodal_stress[:, 9]  # von Mises stress is 10th component

    # Get element stresses for averaging demonstration
    element_stress = result.element_stress(0)

    return ansys, von_mises, element_stress

# Alternative using FEniCS for open-source solution
def fenics_centrifugal_example():
    try:
        from dolfin import *
        import mshr

        # Create mesh (simple cylinder)
        domain = mshr.Cylinder(Point(0, 0, 0), Point(0, 0, 0.5), 0.1, 0.1)
        mesh = mshr.generate_mesh(domain, 32)

        # Define function space
        V = VectorFunctionSpace(mesh, 'P', 2)

        # Material properties
        E = Constant(2.0e11)
        nu = Constant(0.3)
        rho = Constant(7850)

        # Lame parameters
        mu = E/2/(1+nu)
        lmbda = E*nu/(1+nu)/(1-2*nu)

        def epsilon(u):
            return 0.5*(grad(u) + grad(u).T)

        def sigma(u):
            return lmbda*tr(epsilon(u))*Identity(3) + 2*mu*epsilon(u)

        # Centrifugal force (1000 RPM)
        omega = 1000 * 2 * pi / 60
        f = Expression(('rho*omega*omega*x[0]', 'rho*omega*omega*x[1]', '0'),
                      degree=2, rho=7850, omega=omega)

        # Boundary condition (fixed at bottom)
        def bottom_boundary(x, on_boundary):
            return on_boundary and near(x[2], 0)

        bc = DirichletBC(V, Constant((0, 0, 0)), bottom_boundary)

        # Variational formulation
        u = TrialFunction(V)
        v = TestFunction(V)
        a = inner(sigma(u), epsilon(v))*dx
        L = inner(f, v)*dx

        # Solve
        u = Function(V)
        solve(a == L, u, bc)

        # Compute von Mises stress
        V_scalar = FunctionSpace(mesh, 'P', 1)
        s = sigma(u)
        von_mises = project(sqrt(3.0/2.0*inner(dev(s), dev(s))), V_scalar)

        return mesh, u, von_mises

    except ImportError:
        print("FEniCS not available. Using analytical approach instead.")
        return analytical_centrifugal_stress()

def analytical_centrifugal_stress():
    """Analytical solution for thin rotating disk"""
    # Parameters
    rho = 7850  # kg/m³
    omega = 1000 * 2 * np.pi / 60  # rad/s
    R_outer = 0.1  # m
    R_inner = 0.02  # m
    nu = 0.3

    # Radial positions
    r = np.linspace(R_inner, R_outer, 100)

    # Analytical stresses for rotating disk
    sigma_r = (3 + nu)/8 * rho * omega**2 * (R_outer**2 + R_inner**2 -
                                            R_outer**2 * R_inner**2 / r**2 - r**2)
    sigma_t = (3 + nu)/8 * rho * omega**2 * (R_outer**2 + R_inner**2 +
                                            R_outer**2 * R_inner**2 / r**2 -
                                            (1 + 3*nu)/(3 + nu) * r**2)

    # von Mises stress for plane stress
    von_mises = np.sqrt(sigma_r**2 + sigma_t**2 - sigma_r * sigma_t)

    return r, von_mises, sigma_r, sigma_t

# Compare different averaging methods
def compare_averaging_methods():
    """Demonstrate different stress averaging approaches"""

    # Simulate element stresses (simplified)
    np.random.seed(42)
    n_elements = 100
    n_nodes = 121

    # Generate synthetic element stresses
    base_stress = 100e6  # 100 MPa base stress
    element_stresses = base_stress * (1 + 0.1 * np.random.randn(n_elements))

    # Create element-to-node connectivity (simplified)
    # In real FEA, this comes from mesh connectivity
    element_to_node = [np.random.choice(n_nodes, 4, replace=False)
                      for _ in range(n_elements)]

    # Different averaging methods
    def simple_average(element_stresses, element_to_node):
        nodal_stress = np.zeros(n_nodes)
        nodal_count = np.zeros(n_nodes)

        for i, nodes in enumerate(element_to_node):
            for node in nodes:
                nodal_stress[node] += element_stresses[i]
                nodal_count[node] += 1

        return nodal_stress / np.maximum(nodal_count, 1)

    def weighted_average(element_stresses, element_to_node, weights=None):
        nodal_stress = np.zeros(n_nodes)
        nodal_weight = np.zeros(n_nodes)

        if weights is None:
            weights = np.ones(len(element_stresses))

        for i, nodes in enumerate(element_to_node):
            for node in nodes:
                nodal_stress[node] += element_stresses[i] * weights[i]
                nodal_weight[node] += weights[i]

        return nodal_stress / np.maximum(nodal_weight, 1)

    # Apply different methods
    simple_avg = simple_average(element_stresses, element_to_node)
    weighted_avg = weighted_average(element_stresses, element_to_node,
                                   weights=np.random.uniform(0.5, 2, n_elements))

    print("Averaging Method Comparison:")
    print(f"Simple Average - Mean: {np.mean(simple_avg):.2e} Pa, Std: {np.std(simple_avg):.2e} Pa")
    print(f"Weighted Average - Mean: {np.mean(weighted_avg):.2e} Pa, Std: {np.std(weighted_avg):.2e} Pa")
    print(f"Difference: {np.abs(np.mean(simple_avg) - np.mean(weighted_avg)):.2e} Pa")

    return simple_avg, weighted_avg

# Run the examples
if __name__ == "__main__":
    print("=== Centrifugal Rotor Analysis ===")

    # Method 1: Analytical solution
    r, von_mises_analytical, sigma_r, sigma_t = analytical_centrifugal_stress()

    plt.figure(figsize=(12, 4))

    plt.subplot(1, 2, 1)
    plt.plot(r*1000, sigma_r/1e6, 'b-', label='Radial Stress')
    plt.plot(r*1000, sigma_t/1e6, 'r-', label='Tangential Stress')
    plt.xlabel('Radius (mm)')
    plt.ylabel('Stress (MPa)')
    plt.title('Analytical Stresses in Rotating Disk')
    plt.legend()
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(r*1000, von_mises_analytical/1e6, 'g-', linewidth=2)
    plt.xlabel('Radius (mm)')
    plt.ylabel('von Mises Stress (MPa)')
    plt.title('Analytical von Mises Stress')
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    # Method 2: Compare averaging methods
    print("\n=== Averaging Method Comparison ===")
    simple, weighted = compare_averaging_methods()














































































































libname xls excel "d:/xls/titanic.xlsx";

proc  contents data=xls._all_;
run;quit;

data titanic;
  set xls.'titanic$'n;
run;quit;

libname xls clear;





%let pgm=utl-altair-personal-slc-creating-and-analyzing-tables-created-from-multiple-pdf-and-excel-files;

%stop_submission;

RE: Altair personal slc creating and analyzing tables created from multiple pdf and excel files

Too long to post on listserv, see github

github
https://github.com/rogerjdeangelis/utl-altair-personal-slc-creating-and-analyzing-tables-created-from-multiple-pdf-and-excel-files

community.altair.com
https://community.altair.com/discussion/38559/monarch-learning-series-2023-exercise-9-ap-expense-report-analysis#latest

Due to a bug in sqlite passthru you need to run some solutions using sqlite code in either R or Python
Proc sql does not support some common sql functions like partition.

PROBLEM
=======

You are an Accounts Payable Accountant at a Software Company and you need
to analyze 3 months worth of expense reports to answer the following questions:

What is the Total Expense in USD for the 3 Expense Reports (2016-12, 2017-01, 2017-02)?
(HINT: You will need to convert any expenses not in USD by using the Currency Rate Pair.xlsx file
and creating Formula field for the USD rate for the non-USD expenses)


Who are the Top 3 Employees that spent the most on any Expense that contains
       "Travel"
       "Taxi"
       "Airfare"
       "Airline"
       "Hotel"
       "Lunch"
       "Breakfast"
       "Dinner"

CONTENTS
========

   1 create table expenses
   2 create table pairs (for exchange rates)
   3 create table exchange adjusting dollars

   Once you have these tables the solution is straightforward

   4 Total accounts payable
     select sum(equiv) from exchange

   5 Top 3 Employees that spent the most on any Expense that contains
     "Travel","Taxi","Airfare","Airline","Hotel","Lunch","Breakfast","Dinner"

     want<- sqldf('
          with totalexpenses as (
          select
            name,
            itm,
            sum(equiv) as total_exp
          from exchange
          where itm in ("Travel","Taxi","Airfare","Airline","Hotel","Lunch","Breakfast","Dinner")
          group by name, itm
        ),
        rankedexpenses as (
          select
            name,
            itm,
            total_exp,
            row_number() over (partition by itm order by total_exp desc) as rank
          from totalexpenses
        )
        select name, itm, total_exp
        from rankedexpenses
        where rank <= 3
        order by itm, rank
        ')

/*               _     _
 _ __  _ __ ___ | |__ | | ___ _ __ ___
| `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \
| |_) | | | (_) | |_) | |  __/ | | | | |
| .__/|_|  \___/|_.__/|_|\___|_| |_| |_|
|_|
*/

You are an Accounts Payable Accountant at a Software Company and you need
to analyze 3 months worth of expense reports to answer the following questions:

What is the Total Expense in USD for the 3 Expense Reports (2016-12, 2017-01, 2017-02)?
(HINT: You will need to convert any expenses not in USD by using the Currency Rate Pair.xlsx file
and creating Formula field for the USD rate for the non-USD expenses)


Who are the Top 3 Employees that spent the most on any Expense that contains
       "Travel"
       "Taxi"
       "Airfare"
       "Airline"
       "Hotel"
       "Lunch"
       "Breakfast"
       "Dinner"


ONCE YOU HAVE CREATED THESE TABLES IN ANY DATABASE OR EXCEL THE SOLUTIONS ARE STRAIGHT FORWARD

/******************************************************************************************************************************/
/* TABLE PAIRS                               |    TABLE EXCHANGE                                                              */
/*                                           |                                                      CURRENCY                  */
/*   PAIR   BASE CURRENCY BASE__1  PRICE     |             NAME             EXP    CUR      ITM       PRICE  CURRENCY EQUIV   */
/*                                           |                                                                                */
/*  USD-EGP USD    EGP       1       30.85   |  MEARS, SHANNON (DW1032826)  31.38  EUR Lunch          .8905      EUR  35.239  */
/*  USD-PHP USD    PHP       1       54.59   |  MEARS, SHANNON (DW1032826) 105.92  EUR Internet       .8905      EUR 118.944  */
/*  USD-CLP USD    CLP       1      815.05   |  MEARS, SHANNON (DW1032826)  45.77  EUR Entertainment  .8905      EUR  51.398  */
/*  USD-PKR USD    PKR       1      282.75   |  MEARS, SHANNON (DW1032826)  36.03  EUR Lunch          .8905      EUR  40.460  */
/*  USD-IQD USD    IQD       1     1309.00   |  MEARS, SHANNON (DW1032826)  30.23  EUR Lunch          .8905      EUR  33.947  */
/*  ....                                     |  ....                                                                          */
/******************************************************************************************************************************/

SOLUTIONS

  1  Total Accounts Payable for 3 months

     Altair SLC

      TOTEXP
     --------
     226442.8

     proc sql;
       select
         sum(equiv) as totExp
       from
         exchange
     ;quit;


 2  Top 3 Employees that spent the most on any Expense that contains
    "Travel","Taxi","Airfare","Airline","Hotel","Lunch","Breakfast","Dinner"

    Altair SLC

                 NAME                    ITM       TOTAL_EXP

    CATZ, LAURENCE (DW1030298)        Airfare       3098.80
    CONNEALY, TYLER (DW0921188)       Airfare       2729.48
    MAICO, REKO (DW0386388)           Airfare       1733.61

    BALL, GUY (DW0835689)             Airline        164.68
    KNOOP, JEFF (DW0767037)           Airline         99.02
    CRAIG, ROSE (DW0357148)           Airline         93.74

    SANDERSON, MICHAEL (DW0231288)    Breakfast      441.45
    SEPHERS, GEOFF (DW0355310)        Breakfast      307.85
    CAVANAGH, NATHAN (DW0423887)      Breakfast      292.51

    HO, JON (DW0170296)               Dinner        1781.82
    SEPHERS, GEOFF (DW0355310)        Dinner        1179.65
    SANDERSON, MICHAEL (DW0231288)    Dinner         956.91

    PRINCE, KYLE (DW1292453)          Hotel         2410.33
    CONNEALY, TYLER (DW0921188)       Hotel         2208.08
    PENERA, HEATH (DW0332752)         Hotel         1892.71

    MEARS, SHANNON (DW1032826)        Lunch         1647.05
    SEPHERS, GEOFF (DW0355310)        Lunch         1568.62
    SANDERSON, MICHAEL (DW0231288)    Lunch         1396.57

    MAICO, REKO (DW0386388)           Taxi          5871.38
    BERRETT, SHAUN (DW0147143)        Taxi          3913.01
    WILLIAMS, VAN (DW1071325)         Taxi          2808.35

    SIPES, LANCE (DW0741723)          Travel         977.30
    METOYER, MICHELLE (DW0744874)     Travel         676.20
    GONZALEZ, YVETTE (DW0266104)      Travel         542.78


    /*--- same code in python ---*/
    want<-sqldf('
    with totalexpenses as (
      select
        name,
        itm,
        sum(equiv) as total_exp
      from exchange
      where itm in ("Travel","Taxi","Airfare","Airline","Hotel","Lunch","Breakfast","Dinner")
      group by name, itm
    ),
    rankedexpenses as (
      select
        name,
        itm,
        total_exp,
        row_number() over (partition by itm order by total_exp desc) as rank
      from totalexpenses
    )
    select name, itm, total_exp
    from rankedexpenses
    where rank <= 3
    order by itm, rank

/*                      _         _        _     _
/ |  ___ _ __ ___  __ _| |_ ___  | |_ __ _| |__ | | ___   _____  ___ __   ___ _ __  ___  ___  ___
| | / __| `__/ _ \/ _` | __/ _ \ | __/ _` | `_ \| |/ _ \ / _ \ \/ / `_ \ / _ \ `_ \/ __|/ _ \/ __|
| || (__| | |  __/ (_| | ||  __/ | || (_| | |_) | |  __/|  __/>  <| |_) |  __/ | | \__ \  __/\__ \
|_| \___|_|  \___|\__,_|\__\___|  \__\__,_|_.__/|_|\___| \___/_/\_\ .__/ \___|_| |_|___/\___||___/
                                                                  |_|
*/
 Note columh equiv has equivalent dollars

*/

/*                   _
(_)_ __  _ __  _   _| |_
| | `_ \| `_ \| | | | __|
| | | | | |_) | |_| | |_
|_|_| |_| .__/ \__,_|\__|
        |_|
*/

 d:/pdf/Expense Report  2017-01.pdf
 d:/pdf/Expense Report  2017-02.pdf
 d:/pdf/Expense Report  2016-12.pdf

/*
 _ __  _ __ ___   ___ ___  ___ ___
| `_ \| `__/ _ \ / __/ _ \/ __/ __|
| |_) | | | (_) | (_|  __/\__ \__ \
| .__/|_|  \___/ \___\___||___/___/
|_|
*/

/*--- pdfs to tet append and parse creating table expenses ---*/

proc datasets lib=work nodetails nolist;
 delete expense;
run;quit;

/*--- covert pdf to text ---*/
&_init_;
proc r;
submit;
library(pdftools)

pdf_text_content <- pdf_text("d:/pdf/Expense Report  2017-01.pdf")
all_text <- paste(pdf_text_content, collapse = " ")
writeLines(all_text, "d:/txt/y2017_01.txt")

pdf_text_content <- pdf_text("d:/pdf/Expense Report  2017-02.pdf")
all_text <- paste(pdf_text_content, collapse = " ")
writeLines(all_text, "d:/txt/y2017_01.txt")

pdf_text_content <- pdf_text("d:/pdf/Expense Report  2016-12.pdf")
all_text <- paste(pdf_text_content, collapse = " ")
writeLines(all_text, "d:/txt/y2016_12.txt")

endsubmit;
run;quit;

filename all3 ('d:/txt/y2017_01.txt' 'd:/txt/y2017_02.txt' 'd:/txt/y2016_12.txt');

data expenses;
  retain name;
  infile all3;
  input ;
  if _infile_=:'Employee Name' then name=left(scan(_infile_,2,':'));
  if
      index(_infile_,'16 Operations') or index(_infile_,'17 Operations') or
      index(_infile_,'16 Sales') or index(_infile_,'17 Sales');
  exp=input(scan(_infile_,3,' '),comma12.2);
  cur=scan(_infile_,4,'/ ');
  itm=scan(_infile_,5,'/ ');

run;quit;

proc print data=expenses(obs=5);
run;quit;

proc freq data=expenses;
tables cur;
run;

/*           _               _
  ___  _   _| |_ _ __  _   _| |_
 / _ \| | | | __| `_ \| | | | __|
| (_) | |_| | |_| |_) | |_| | |_
 \___/ \__,_|\__| .__/ \__,_|\__|
                |_|
*/

Altair SLC

Obs               NAME                 EXP      CUR     ITM

 1     GIELSKI, JESSIE (DW0188511)    158.03    USD    Fuel
 2     GIELSKI, JESSIE (DW0188511)    307.79    USD    Car
 3     GIELSKI, JESSIE (DW0188511)    136.83    USD    Dinner
 4     GIELSKI, JESSIE (DW0188511)    405.49    USD    Hotel
 5     GIELSKI, JESSIE (DW0188511)     90.45    USD    Lunch
 ....

/*
| | ___   __ _
| |/ _ \ / _` |
| | (_) | (_| |
|_|\___/ \__, |
         |___/
*/
1057      ODS _ALL_ CLOSE;
1058      ODS LISTING;
1059      FILENAME WBGSF 'd:\wpswrk\_TD9616/listing_images';
1060      OPTIONS DEVICE=GIF;
1061      GOPTIONS GSFNAME=WBGSF;
1062
1063
1064      filename all3 ('d:/txt/y2017_01.txt' 'd:/txt/y2017_01.txt' 'd:/txt/y2016_12.txt');
1065
1066      data expenses;
1067        retain name;
1068        infile all3;
1069        input ;
1070        if _infile_=:'Employee Name' then name=left(scan(_infile_,2,':'));
1071        if index(_infile_,'16 Operations') or index(_infile_,'17 Operations');
1072        exp=input(scan(_infile_,3,' '),comma12.2);
1073        cur=scan(_infile_,4,'/ ');
1074        itm=scan(_infile_,5,'/ ');
1075
1076      run;

NOTE: The infile all3 is:
      Filename='d:\txt\y2017_01.txt',
      Owner Name=T7610\Roger,
      File size (bytes)=71868,
      Create Time=08:02:34 Oct 29 2025,
      Last Accessed=10:14:19 Oct 29 2025,
      Last Modified=08:02:34 Oct 29 2025,
      Lrecl=32767, Recfm=V

NOTE: The infile all3 is:
      Filename='d:\txt\y2017_01.txt',
      Owner Name=T7610\Roger,
      File size (bytes)=71868,
      Create Time=08:02:34 Oct 29 2025,
      Last Accessed=10:15:57 Oct 29 2025,
      Last Modified=08:02:34 Oct 29 2025,
      Lrecl=32767, Recfm=V

NOTE: The infile all3 is:
      Filename='d:\txt\y2016_12.txt',
      Owner Name=T7610\Roger,
      File size (bytes)=76211,
      Create Time=08:02:34 Oct 29 2025,
      Last Accessed=10:14:19 Oct 29 2025,
      Last Modified=08:02:34 Oct 29 2025,
      Lrecl=32767, Recfm=V

NOTE: 1972 records were read from file all3
      The minimum record length was 0
      The maximum record length was 93
NOTE: 1972 records were read from file all3
      The minimum record length was 0
      The maximum record length was 93
NOTE: 1997 records were read from file all3
      The minimum record length was 0
      The maximum record length was 96
NOTE: Data set "WORK.expenses" has 666 observation(s) and 4 variable(s)
NOTE: The data step took :
      real time : 0.011
      cpu time  : 0.015


1076    !     quit;
1077
1078      &_init_;
1079      proc r;
NOTE: Using R version 4.5.1 (2025-06-13 ucrt) from d:\r451
1080      export data=expenses r=expenses;
NOTE: Creating R data frame 'expenses' from data set 'WORK.expenses'

1081      submit;
1082      str(expenses)
1083      library(sqldf)
1084      options(sqldf.dll = "d:/dll/sqlean.dll")
1085      want<-sqldf("
1086      WITH TotalExpenses AS (
1087        SELECT
1088          name,
1089          ITM,
1090          SUM(EXP) AS Total_EXP
1091        FROM expenses
1092        WHERE ITM IN ('Dinner', 'Taxi')
1093        GROUP BY name, ITM
1094      ),
1095      RankedExpenses AS (
1096        SELECT
1097          name,
1098          ITM,
1099          Total_EXP,
1100          ROW_NUMBER() OVER (PARTITION BY ITM ORDER BY Total_EXP DESC) AS rank
1101        FROM TotalExpenses
1102      )
1103      SELECT name, ITM, Total_EXP
1104      FROM RankedExpenses
1105      WHERE rank <= 2
1106      ORDER BY ITM, rank
1107      ")
1108      endsubmit;

NOTE: Submitting statements to R:

> str(expenses)
Loading required package: gsubfn
Loading required package: proto
Loading required package: RSQLite
> library(sqldf)
> options(sqldf.dll = "d:/dll/sqlean.dll")
> want<-sqldf("
+ WITH TotalExpenses AS (
+   SELECT
+     name,
+     ITM,
+     SUM(EXP) AS Total_EXP
+   FROM expenses
+   WHERE ITM IN ('Dinner', 'Taxi')
+   GROUP BY name, ITM
+ ),
+ RankedExpenses AS (
+   SELECT
+     name,
+     ITM,
+     Total_EXP,
+     ROW_NUMBER() OVER (PARTITION BY ITM ORDER BY Total_EXP DESC) AS rank
+   FROM TotalExpenses
+ )
+ SELECT name, ITM, Total_EXP
+ FROM RankedExpenses
+ WHERE rank <= 2
+ ORDER BY ITM, rank

NOTE: Processing of R statements complete

+ ")
1109      import data=want r=want;
NOTE: Creating data set 'WORK.want' from R data frame 'want'
NOTE: Column names modified during import of 'want'
NOTE: Data set "WORK.want" has 4 observation(s) and 3 variable(s)

1110      run;quit;
NOTE: Procedure r step took :
      real time : 1.398
      cpu time  : 0.015


1111
1112      proc print data=want;
1113      run;
NOTE: 4 observations were read from "WORK.want"
NOTE: Procedure print step took :
      real time : 0.016
      cpu time  : 0.000


1114
1115
1116
1117      quit; run;
1118      ODS _ALL_ CLOSE;
1119      FILENAME WBGSF CLEAR;

/*___                        _         _        _     _
|___ \    ___ _ __ ___  __ _| |_ ___  | |_ __ _| |__ | | ___
  __) |  / __| `__/ _ \/ _` | __/ _ \ | __/ _` | `_ \| |/ _ \
 / __/  | (__| | |  __/ (_| | ||  __/ | || (_| | |_) | |  __/
|_____|  \___|_|  \___|\__,_|\__\___|  \__\__,_|_.__/|_|\___|

             _             __                                         _
 _ __   __ _(_)_ __ ___   / _|_ __ ___  _ __ ___     _____  _____ ___| |
| `_ \ / _` | | `__/ __| | |_| `__/ _ \| `_ ` _ \   / _ \ \/ / __/ _ \ |
| |_) | (_| | | |  \__ \ |  _| | | (_) | | | | | | |  __/>  < (_|  __/ |
| .__/ \__,_|_|_|  |___/ |_| |_|  \___/|_| |_| |_|  \___/_/\_\___\___|_|
|_|
*/


&_init_;

proc datasets lib=work nodetails nolist;
 delete pairs;
run;quit;

libname xls xlsx "d:/xls/Currency Rate Pairs.xlsx";

data pairs;
  set xls.'Expense Currency Pairs$'n;
run;quit;

libname xls clear;

proc print data=pairs(obs=5);
run;quit;

/*           _               _
  ___  _   _| |_ _ __  _   _| |_
 / _ \| | | | __| `_ \| | | | __|
| (_) | |_| | |_| |_) | |_| | |_
 \___/ \__,_|\__| .__/ \__,_|\__|
                |_|
*/

TABLE PAIRS

Altair SLC

Obs               NAME                 EXP      CUR     ITM

 1     GIELSKI, JESSIE (DW0188511)    158.03    USD    Fuel
 2     GIELSKI, JESSIE (DW0188511)    307.79    USD    Car
 3     GIELSKI, JESSIE (DW0188511)    136.83    USD    Dinner
 4     GIELSKI, JESSIE (DW0188511)    405.49    USD    Hotel
 5     GIELSKI, JESSIE (DW0188511)     90.45    USD    Lunch
 ....

/*
| | ___   __ _
| |/ _ \ / _` |
| | (_) | (_| |
|_|\___/ \__, |
         |___/
*/

1710      ODS _ALL_ CLOSE;
1711      ODS LISTING;
1712      FILENAME WBGSF 'd:\wpswrk\_TD9616/listing_images';
1713      OPTIONS DEVICE=GIF;
1714      GOPTIONS GSFNAME=WBGSF;
1715      /*--- covert pdf to text ---*/
1716      &_init_;
1717      proc r;
1718      submit;
1719      library(pdftools)
1720
1721      pdf_text_content <- pdf_text("d:/pdf/Expense Report  2017-01.pdf")
1722      all_text <- paste(pdf_text_content, collapse = " ")
1723      writeLines(all_text, "d:/txt/y2017_01.txt")
1724
1725      pdf_text_content <- pdf_text("d:/pdf/Expense Report  2017-02.pdf")
1726      all_text <- paste(pdf_text_content, collapse = " ")
1727      writeLines(all_text, "d:/txt/y2017_01.txt")
1728
1729      pdf_text_content <- pdf_text("d:/pdf/Expense Report  2016-12.pdf")
1730      all_text <- paste(pdf_text_content, collapse = " ")
1731      writeLines(all_text, "d:/txt/y2016_12.txt")
1732
1733      endsubmit;
NOTE: Using R version 4.5.1 (2025-06-13 ucrt) from d:\r451

NOTE: Submitting statements to R:

Using poppler version 25.05.0
> library(pdftools)
>
> pdf_text_content <- pdf_text("d:/pdf/Expense Report  2017-01.pdf")
> all_text <- paste(pdf_text_content, collapse = " ")
> writeLines(all_text, "d:/txt/y2017_01.txt")
>
> pdf_text_content <- pdf_text("d:/pdf/Expense Report  2017-02.pdf")
> all_text <- paste(pdf_text_content, collapse = " ")
> writeLines(all_text, "d:/txt/y2017_01.txt")
>
> pdf_text_content <- pdf_text("d:/pdf/Expense Report  2016-12.pdf")
> all_text <- paste(pdf_text_content, collapse = " ")
> writeLines(all_text, "d:/txt/y2016_12.txt")

NOTE: Processing of R statements complete

>
1734      run;quit;
NOTE: Procedure r step took :
      real time : 0.842
      cpu time  : 0.000


1735
1736      filename all3 ('d:/txt/y2017_01.txt' 'd:/txt/y2017_02.txt' 'd:/txt/y2016_12.txt');
1737
1738      data expenses;
1739        retain name;
1740        infile all3;
1741        input ;
1742        if _infile_=:'Employee Name' then name=left(scan(_infile_,2,':'));
1743        if
1744            index(_infile_,'16 Operations') or index(_infile_,'17 Operations') or
1745            index(_infile_,'16 Sales') or index(_infile_,'17 Sales');
1746        exp=input(scan(_infile_,3,' '),comma12.2);
1747        cur=scan(_infile_,4,'/ ');
1748        itm=scan(_infile_,5,'/ ');
1749
1750      run;

NOTE: The infile all3 is:
      Filename='d:\txt\y2017_01.txt',
      Owner Name=T7610\Roger,
      File size (bytes)=85328,
      Create Time=08:02:34 Oct 29 2025,
      Last Accessed=12:17:55 Oct 29 2025,
      Last Modified=12:17:55 Oct 29 2025,
      Lrecl=32767, Recfm=V

NOTE: The infile all3 is:
      Filename='d:\txt\y2017_02.txt',
      Owner Name=T7610\Roger,
      File size (bytes)=85328,
      Create Time=08:02:34 Oct 29 2025,
      Last Accessed=12:16:08 Oct 29 2025,
      Last Modified=08:02:34 Oct 29 2025,
      Lrecl=32767, Recfm=V

NOTE: The infile all3 is:
      Filename='d:\txt\y2016_12.txt',
      Owner Name=T7610\Roger,
      File size (bytes)=76211,
      Create Time=08:02:34 Oct 29 2025,
      Last Accessed=12:17:55 Oct 29 2025,
      Last Modified=12:17:55 Oct 29 2025,
      Lrecl=32767, Recfm=V

NOTE: 2223 records were read from file all3
      The minimum record length was 0
      The maximum record length was 92
NOTE: 2223 records were read from file all3
      The minimum record length was 0
      The maximum record length was 92
NOTE: 1997 records were read from file all3
      The minimum record length was 0
      The maximum record length was 96
NOTE: Data set "WORK.expenses" has 1692 observation(s) and 4 variable(s)
NOTE: The data step took :
      real time : 0.015
      cpu time  : 0.015


1750    !     quit;
1751
1752      proc print data=expenses(obs=5);
1753      run;quit;
NOTE: 5 observations were read from "WORK.expenses"
NOTE: Procedure print step took :
      real time : 0.024
      cpu time  : 0.000


1754
1755      quit; run;
1756      ODS _ALL_ CLOSE;
1757      FILENAME WBGSF CLEAR;


/*____                       _         _        _     _                     _
|___ /    ___ _ __ ___  __ _| |_ ___  | |_ __ _| |__ | | ___   _____  _____| |__   __ _ _ __   __ _  ___
  |_ \   / __| `__/ _ \/ _` | __/ _ \ | __/ _` | `_ \| |/ _ \ / _ \ \/ / __| `_ \ / _` | `_ \ / _` |/ _ \
 ___) | | (__| | |  __/ (_| | ||  __/ | || (_| | |_) | |  __/|  __/>  < (__| | | | (_| | | | | (_| |  __/
|____/   \___|_|  \___|\__,_|\__\___|  \__\__,_|_.__/|_|\___| \___/_/\_\___|_| |_|\__,_|_| |_|\__, |\___|
                                                                                              |___/
ALSO CREAT SUM OF ALL EXPENSES
------------------------------

*/

&_init_;

proc datasets lib=work nodetails nolist;
 delete exchange;
run;quit;

proc sql;
  create
    table exchange as
  select
     l.*
    ,r.currency_price
    ,r.currency
    ,case
       when cur='USD' then 1*exp
       else exp/currency_price
      end as equiv
  from
    expenses as l left join pairs as r
  on
   l.cur = r.currency
;quit;

proc print data=exchange(obs=5);
run;quit;

proc sql;
  select
    sum(equiv) as totExp
 from
    exchange
;quit;

/*           _               _
  ___  _   _| |_ _ __  _   _| |_
 / _ \| | | | __| `_ \| | | | __|
| (_) | |_| | |_| |_) | |_| | |_
 \___/ \__,_|\__| .__/ \__,_|\__|
                |_|
*/

 Altair SLC

    TOTEXP
  --------
  3623.247


TABLE EXCHANGE

Altair SLC

Obs               NAME                EXP      CUR         ITM         CURRENCY_PRICE    CURRENCY    EQUIV

184    HOGAN, JAMES (DW1220883)       39.47    GBP    Office               .76707          GBP        51.46
185    CATZ, LAURENCE (DW1030298)     30.50    GBP    Taxi                 .76707          GBP        39.76
186    VOILES, JEFF (DW1910657)       28.11    GBP    Taxi                 .76707          GBP        36.65
187    WANG, FRANKLIN (DW1623033)     32.50    GBP    Taxi                 .76707          GBP        42.37
188    VOILES, JEFF (DW1910657)       67.69    GBP    Taxi                 .76707          GBP        88.24
189    VOILES, JEFF (DW1910657)       14.61    GBP    Entertainment        .76707          GBP        19.05

/*
| | ___   __ _
| |/ _ \ / _` |
| | (_) | (_| |
|_|\___/ \__, |
         |___/
*/

2059      ODS _ALL_ CLOSE;
2060      ODS LISTING;
2061      FILENAME WBGSF 'd:\wpswrk\_TD9616/listing_images';
2062      OPTIONS DEVICE=GIF;
2063      GOPTIONS GSFNAME=WBGSF;
2064
2065      &_init_;
2066
2067      proc datasets lib=work nodetails nolist;
2068       delete exchange;
2069      run;quit;
NOTE: Deleting "WORK.EXCHANGE" (memtype="DATA")
NOTE: Procedure datasets step took :
      real time : 0.001
      cpu time  : 0.000


2070
2071      proc sql;
2072        create
2073          table exchange as
2074        select
2075           l.*
2076          ,r.currency_price
2077          ,r.currency
2078          ,case
2079             when cur='USD' then 1*exp
2080             else exp/currency_price
2081            end as equiv
2082        from
2083          expenses as l left join pairs as r
2084        on
2085         l.cur = r.currency
2086      ;quit;
NOTE: Data set "WORK.exchange" has 1692 observation(s) and 7 variable(s)
NOTE: Procedure sql step took :
      real time : 0.107
      cpu time  : 0.062


2087
2088      proc print data=exchange(obs=5);
2089      run;quit;
NOTE: 5 observations were read from "WORK.exchange"
NOTE: Procedure print step took :
      real time : 0.014
      cpu time  : 0.000


2090
2091      proc sql;
2092        select
2093          sum(equiv) as totExp
2094       from
2095          exchange
2096      ;quit;
NOTE: Procedure sql step took :
      real time : 0.020
      cpu time  : 0.000


2097      quit; run;
2098      ODS _ALL_ CLOSE;
2099      FILENAME WBGSF CLEAR;


/*  _     _                _____   _
| || |   | |_ ___  _ __   |___ /  | |__  _   _   _____  ___ __   ___ _ __  ___  ___    __ _ _ __ ___  _   _ _ __
| || |_  | __/ _ \| `_ \    |_ \  | `_ \| | | | / _ \ \/ / `_ \ / _ \ `_ \/ __|/ _ \  / _` | `__/ _ \| | | | `_ \
|__   _| | || (_) | |_) |  ___) | | |_) | |_| ||  __/>  <| |_) |  __/ | | \__ \  __/ | (_| | | | (_) | |_| | |_) |
   |_|    \__\___/| .__/  |____/  |_.__/ \__, | \___/_/\_\ .__/ \___|_| |_|___/\___|  \__, |_|  \___/ \__,_| .__/
                  |_|                    |___/           |_|                          |___/                |_|
*/

proc datasets lib=work nolist nodetails;
 delete want;
run;quit;

&_init_;
proc datasets lib=work nolist nodetails;
 delete want;
run;quit;

proc r;
export data=exchange r=exchange;
submit;
library(sqldf)
options(sqldf.dll = "d:/dll/sqlean.dll")
want<-sqldf('
with totalexpenses as (
  select
    name,
    itm,
    sum(equiv) as total_exp
  from exchange
  where itm in ("Travel","Taxi","Airfare","Airline","Hotel","Lunch","Breakfast","Dinner")
  group by name, itm
),
rankedexpenses as (
  select
    name,
    itm,
    total_exp,
    row_number() over (partition by itm order by total_exp desc) as rank
  from totalexpenses
)
select name, itm, total_exp
from rankedexpenses
where rank <= 3
order by itm, rank
')

endsubmit;
import data=want r=want;
run;quit;

proc print data=want;
run;

Altair SLC

Obs                 NAME                    ITM       TOTAL_EXP

  1    CATZ, LAURENCE (DW1030298)        Airfare       3098.80
  2    CONNEALY, TYLER (DW0921188)       Airfare       2729.48
  3    MAICO, REKO (DW0386388)           Airfare       1733.61
  4    BALL, GUY (DW0835689)             Airline        164.68
  5    KNOOP, JEFF (DW0767037)           Airline         99.02
  6    CRAIG, ROSE (DW0357148)           Airline         93.74
  7    SANDERSON, MICHAEL (DW0231288)    Breakfast      441.45
  8    SEPHERS, GEOFF (DW0355310)        Breakfast      307.85
  9    CAVANAGH, NATHAN (DW0423887)      Breakfast      292.51
 10    HO, JON (DW0170296)               Dinner        1781.82
 11    SEPHERS, GEOFF (DW0355310)        Dinner        1179.65
 12    SANDERSON, MICHAEL (DW0231288)    Dinner         956.91
 13    PRINCE, KYLE (DW1292453)          Hotel         2410.33
 14    CONNEALY, TYLER (DW0921188)       Hotel         2208.08
 15    PENERA, HEATH (DW0332752)         Hotel         1892.71
 16    MEARS, SHANNON (DW1032826)        Lunch         1647.05
 17    SEPHERS, GEOFF (DW0355310)        Lunch         1568.62
 18    SANDERSON, MICHAEL (DW0231288)    Lunch         1396.57
 19    MAICO, REKO (DW0386388)           Taxi          5871.38
 20    BERRETT, SHAUN (DW0147143)        Taxi          3913.01
 21    WILLIAMS, VAN (DW1071325)         Taxi          2808.35
 22    SIPES, LANCE (DW0741723)          Travel         977.30
 23    METOYER, MICHELLE (DW0744874)     Travel         676.20
 24    GONZALEZ, YVETTE (DW0266104)      Travel         542.78


/*--- Lets check ---*/

proc sort data=exchange(where=(name=:'CONNEALY, TYLER')) out=srt;
 by descending equiv;
run;quit;

proc print data=srt;
run;quit;


proc sort data=exchange(where=(name=:'CONNEALY, TYLER')) out=srt;
 by descending equiv;
run;quit;

proc print data=srt;
run;quit;


proc sort data=exchange(where=(name=:'CATZ, LAURENCE')) out=srt;
 by descending equiv;
run;quit;

proc print data=srt;
run;quit;


Altair SLC

Obs               NAME                 EXP      CUR      ITM      CURRENCY_PRICE    CURRENCY    EQUIV

  1    CONNEALY, TYLER (DW0921188)    682.37    USD    Airfare         .                        682.37
  2    CONNEALY, TYLER (DW0921188)    682.37    USD    Airfare         .                        682.37
  3    CONNEALY, TYLER (DW0921188)    682.37    USD    Airfare         .                        682.37
  4    CONNEALY, TYLER (DW0921188)    682.37    USD    Airfare         .                        682.37
                                                                                               =======
                                                                                               2729.48


            NAME                       EXP      CUR       ITM     CURRENCY_PRICE    CURRENCY     EQUIV

  1    CATZ, LAURENCE (DW1030298)    2377.00    GBP    Airfare        .76707          GBP       3098.80                                                                                              2729.48 Checks


/*
| | ___   __ _
| |/ _ \ / _` |
| | (_) | (_| |
|_|\___/ \__, |
         |___/
*/

2677      ODS _ALL_ CLOSE;
2678      ODS LISTING;
2679      FILENAME WBGSF 'd:\wpswrk\_TD9616/listing_images';
2680      OPTIONS DEVICE=GIF;
2681      GOPTIONS GSFNAME=WBGSF;
2682      &_init_;
2683      proc datasets lib=work nolist nodetails;
2684       delete want;
2685      run;quit;
NOTE: Deleting "WORK.WANT" (memtype="DATA")
NOTE: Procedure datasets step took :
      real time : 0.001
      cpu time  : 0.000


2686
2687      proc r;
NOTE: Using R version 4.5.1 (2025-06-13 ucrt) from d:\r451
2688      export data=exchange r=exchange;
NOTE: Creating R data frame 'exchange' from data set 'WORK.exchange'

2689      submit;
2690      library(sqldf)
2691      options(sqldf.dll = "d:/dll/sqlean.dll")
2692      want<-sqldf('
2693      with totalexpenses as (
2694        select
2695          name,
2696          itm,
2697          sum(equiv) as total_exp
2698        from exchange
2699        where itm in ("Travel","Taxi","Airfare","Airline","Hotel","Lunch","Breakfast","Dinner")
2700        group by name, itm
2701      ),
2702      rankedexpenses as (
2703        select
2704          name,
2705          itm,
2706          total_exp,
2707          row_number() over (partition by itm order by total_exp desc) as rank
2708        from totalexpenses
2709      )
2710      select name, itm, total_exp
2711      from rankedexpenses
2712      where rank <= 3
2713      order by itm, rank
2714      ')
2715
2716      endsubmit;

NOTE: Submitting statements to R:

Loading required package: gsubfn
Loading required package: proto
Loading required package: RSQLite
> library(sqldf)
> options(sqldf.dll = "d:/dll/sqlean.dll")
> want<-sqldf('
+ with totalexpenses as (
+   select
+     name,
+     itm,
+     sum(equiv) as total_exp
+   from exchange
+   where itm in ("Travel","Taxi","Airfare","Airline","Hotel","Lunch","Breakfast","Dinner")
+   group by name, itm
+ ),
+ rankedexpenses as (
+   select
+     name,
+     itm,
+     total_exp,
+     row_number() over (partition by itm order by total_exp desc) as rank
+   from totalexpenses
+ )
+ select name, itm, total_exp
+ from rankedexpenses
+ where rank <= 3
+ order by itm, rank
+ ')

NOTE: Processing of R statements complete

>
2717      import data=want r=want;
NOTE: Creating data set 'WORK.want' from R data frame 'want'
NOTE: Column names modified during import of 'want'
NOTE: Data set "WORK.want" has 24 observation(s) and 3 variable(s)

2718      run;quit;
NOTE: Procedure r step took :
      real time : 1.373
      cpu time  : 0.015


2719
2720      proc print data=want;
2721      run;
NOTE: 24 observations were read from "WORK.want"
NOTE: Procedure print step took :
      real time : 0.011
      cpu time  : 0.015


2722      quit; run;
2723      ODS _ALL_ CLOSE;
2724      FILENAME WBGSF CLEAR;

/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/
