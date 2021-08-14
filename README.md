# test

I. About:
================
  Time-dependent wave packet method for cold atom-diatom scatterings. The
  total wavefunction is decomposed into long-range, asymptotic and interaction
  regions. The split-operator method is used for propagation. Coordinate
  transform RCB method is used for S-matrix collection.

II. Inputs: (atomic units, unless noted specified)
================
  1) split_inputs
       r1, r2        -- The range of r. (r2 > r1)
       nr, nrasy     -- Number of PODVR grids for inteaction and asymptotic
                        regions for r. (nr > nrasy)
       z1, z2        -- The range of z. (z2 > z1)
       nzlon, nz, nzint -- Number of Sinc-DVR grids for long-range, asymptotic
                           and interaction regions for z. (nzlon > nz > nzint)
       jmax, jmaxasy -- The maximum rotational basis used for interaction and
                        asymptotic regions. (jmax >= jmaxasy)
       kcut          -- Cutoff for Kmax used for large Jtot.
       atom          -- The atoms: A, B and C. For A + BC reaction.
       jtotal        -- Total angular momentum quantum number J.
       v0            -- Reactant initial vibrational quantum number.
       j0            -- Reactant initial rotational quantum number.
       l0            -- Oribital angular momentum quantum number.
                        ( j0 + J >= l0 >= |j0 - J| )
       esk0          -- Central collison energy of wave packet (eV).
       delta         -- Width of wave packet.
       z0            -- Central position of wave packet.
       clabs, labs, plabs -- Parameters of damping function for z at long-range
                             region. clabs and plabs are strength of damping
                             function, plabs is starting point.
       czabs, zabs, pzabs -- Parameters of damping function for z at asymptotic
                             region.
       crabs, rabs, prabs -- Parameters of damping function for r at
                             interaction region.
       sotype        -- Type of split-operator method, 22 for seconed order, 46
                        for fourth order, 68 for sixth order.
       dt            -- Time step.
       ttot          -- Total propagation time.
       tprt          -- Print norms of total wavefunction at every tprt times.
       isave         -- If equal to 1, backup TD and TID wavefunctions at every tprt,
                        to restart program.
       iread         -- If equal to 1, read TD and TID wavefunctions and restart
                        program.
       npes          -- Total potential energy surface number.
       ipes          -- Which diabatic PES to construct initial wave packet.

  2) esk_inputs
       neblock       -- Number of collision energy block.
       ecol          -- Collision energy block arrays, ecol(3*neblock) eV. Example:
                          neblock = 2;
                          ecol = 0.1, 0.4, 5,
                                 0.5, 1.0, 10
                        There are 2 energy blocks, 5 energies at [0.1, 0.4] and 10
                        energies at [0.5, 1.0].
       netid         -- If > 0 then collect TID-wavefunction at interaction region
                        at etid(netid) eV.
       etid          -- Collision energy arrays for collecting TID wavefunctions.

  3) reac_inputs
       idreac        -- If equal to 1, calculate flux reaction probabilites.
       rflux         -- Flux postion for r. (r < r2)

  4) rcb_inputs (just set AC+B channel if reactant is homonuclear diatomic molecule)
       idpab         -- If equal to 1, calculate S-matrix for AB + C channel.
       zinfab        -- Smatrix projection plane at product Jacobi coordinate.
       rab1, rab2    -- Range of rab.
       idpac         -- If equal to 1, calculate S-matrix for AC + B channel.
       zinfac        -- Smatrix projection plane at product Jacobi coordinate.
       rac1, rac2    -- Range of rac.

III. Outputs:
================
  1) fort.***
     (a) Split outputs: ( fort.1**, ground is fort.5**, excited is fort.6** )
       101 -- Potential energy curve for cunstructing PODVR.
       102 -- PODVR grids of asymptotic region and corresponding energies.
       103 -- PODVR grids of interaction region and corresponding energies.
       114 -- Initial state wavefunction at r degree.
       118 -- Initial wave packet at z degree.
       119 -- Collision energy coefficient of initial wave packet.
       120 -- Damping function for z at asymptotic region.
       121 -- Damping function for r at interaction region.
       122 -- Damping function for z at long-range region.
       501,601,... -- 2D adiabatic interaction PES as function of r and z.
       502,602,... -- 2D adiabatic asymptotic PES as function of r and z.
       503 or 603 or ... -- 2D adiabatic long-range asymptotic PES.
       510,610,... -- 2D diabatic interaction PES as function of r and z.
       526,626,... -- Norm of time-dependent total wave packet at z degree.
       527,627,... -- Norm of time-dependent total wave packet at r degree.
       528,628,... -- Norm of time-dependent total wave packet for inelastic.

     (b) reac outputs:
       523,623,... -- Total reaction probabilities of AB + C channels at
                      different electronic states.
       524,624,... -- Total reaction probabilities of AC + B channels at
                      different electronic states.
       525,625,... -- Total reaction probabilities at different electronic
                      states.

     (c) RCB outputs: ( AB+C is fort.3**, AC+B is fort.4**)
       301,302,.../401,402,... -- Potential energy curves of AB/AC for PODVR.
       311,312,.../411,412,... -- Integration grids of S-matrix at reactant
                                  Jacobi coordinate. (zbc,rbc,cost,V_ad,V_di)
       321,322,.../421,422,... -- Integration grids of S-matrix at product
                                  Jacobi coordinate. (zinf,rp,costp,Vad,Vdi)
       350/450 -- Norm of product ro-vib WFs at reactant grids. (ip,vp,jp,kp,Evj,an)
       391,392,.../491,492,... -- S-matrix total reaction probabilites of AB/AC
                                  channels at different electronic states.
       (AB/AC)J***L***P*.Smatrix -
               - Parity-adapted state-to-state reaction S-matrix of AB/AC
                 channel. J is total angular momentum, L is oribital angular
                 momentum, P is adiabatic PES.

  2) ***.save: Backups of TD and TID WFs at tt time for restarting program.
    WAVE  -- Total wave functions.
    PHIEr -- Rflux TIWF for flux reaction probabilities.
    PHIEab/PHIEac  -- TIWF for reaction S-matrix of AB/AC channel.

IV. Linking PES:
================
  1) Add PES routine into makefile.

  2) Modify your PES subroutine as:

      subroutine pot0(r1,r2,r3,npes,v)
      real*8,intent(in) :: r1,r2,r3
      integer,intent(in) :: npes
      real*8,intent(out) :: v(npes,npes)

    where r1, r2 and r3 are bond length (bohr) of AB, BC and AC, respectively.
    npes is total number of potential energy surfaces.
    v(npes,npes) returns diabatic potential energy matrix at atomic units.
  3) Need to check if the atom mass is defined correctly in subroutine set_mass
    at sub.f90.

V. Usage: Two commands.
================
  1) make
       Complile the code and one executable file is generated. (scatt.x)

  2) ./scatt.x < inputs
       Running program.

VI. Contact
================
  If you have any feedback, questions, suggestions, do not hestitate to contact: burin201910@gmail.com.

Bayaer Buren
2021.05.29
