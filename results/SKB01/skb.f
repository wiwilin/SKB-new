!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      program skb 
 
!     SIMULATION NUMERIQUE DE LA PROPAGATION DES ONDES SISMO-ACOUSTIQUES
!     EN MILIEU STRATIFIE PLAN
!     Version P-SV-SH du 07/12/99
!     M.D. LGIT Grenoble
!     Langage F77 avec extensions F90 XLF (AIX IBM RS/6000 Version 3.2.4)
 
!     Conventions: 
!     - La source est situee a l'abscisse et a l'ordonnee 0
!     - L'axe des profondeurs est dirige verticalement vers le bas

!=============================================================================

!---- DECLARATIONS GENERALES

!     implicit  none 

!.... Parametres fixes et semi-fixes

      parameter (nkm = 600000, nlm = 200, ntm = 4096, nfm = ntm/2, 
     &           nxm = 1, nym = 1, nzm = 200)

!.... Autres declarations
         character*80 fileou1,fileou2,fileou3,fileou4

      integer   i, ic, idirect, iflag, ifr, ifsurf, inum, iprofil, 
     &          isrc, it, ix, ix0, iy, iy0, iz, iz0, izz, j, k, k1, 
     &          la, lc, lc0dn, lc0up, ls0, lstart, nabove, nbelow, 
     &          nf, nk, nl, nl1, nliq, npos2, nrliq, nt, nx, ny, nz, 
     &          rd_only, ru_only, 
     &          lcount(nlm), lr(nzm), mult(nlm), nkconv(nfm)

      real      a00, a0, a1, a2, a3, a10, a20, a30, ak, ak2, aw, 
     &          df, dfk, dk, dt, eps, ex1, ex2, fimag, fpeak, fref, 
     &          freq, omp2, omref, pi, pi2, rw, t0, tl, xl, 
     &          zbl, zs0, zstart, arg(nkm), bj0(nkm), bj1(nkm),
     &          bbj0(nkm,nxm,nym), bbj1(nkm,nxm,nym),
     &          bjinc(nkm,nxm,nym), cphi(nxm,nym), sphi(nxm,nym),
     &          dens(nlm), qp(nlm), qs(nlm), vp(nlm), vs(nlm), 
     &          zb(0:nlm), rpos(nxm,nym), apos(nxm,nym), zpos(nzm),
     &          zpos2(nzm+1), amp(ntm), y(ntm)

      complex   ai, clnf, om2, omega, q0, q1, q2, q3, q4, 
     &          qq0, qq1, qq2,qq3,      
     &          sk1dn, sk1up, sk2dn, sk2up, sk3dn, sk3up, st1, st2, st3,
     &          r0(nzm,2), rc(nzm,2,2), rmud(nzm,2,2), sfield(nzm,2),  
     &          rcsh(nzm), rmudsh(nzm), sfieldsh(nzm), 
     &          alpha(nlm), ap(nlm), as(nlm), beta(nlm), bulk(nlm), 
     &          emu(nlm), w2(nlm), wa2(nlm), wb2(nlm), 
     &          wza(nlm), wzb(nlm), src(nfm), 
     &          rdi(nlm,2,2), rui(nlm,2,2), tdi(nlm,2,2), tui(nlm,2,2),
     &          rdish(nlm), ruish(nlm), tdish(nlm), tuish(nlm),                 
     &          up(nxm,nym,nzm), ur(nxm,nym,nzm), 
     &          ut(nxm,nym,nzm), uz(nxm,nym,nzm), 
     &          ufp(nfm,nxm,nym,nzm), ufx(nfm,nxm,nym,nzm), 
     &          ufy(nfm,nxm,nym,nzm), ufz(nfm,nxm,nym,nzm), x(nfm), 
     &          rdgen(nzm+1,2,2), rugen(nzm+1,2,2), rdsl(2,2),
     &          tdgen(nzm+1,2,2), tugen(nzm+1,2,2), rufs(2,2),
     &          rdgensh(nzm+1), rugensh(nzm+1), rdslsh,
     &          tdgensh(nzm+1), tugensh(nzm+1), rufssh,
     &          rdrl(nzm,2,2), rdrs(nzm,2,2), rurs(nzm,2,2), 
     &          rufr(nzm,2,2), tdrs(nzm,2,2), turs(nzm,2,2), 
     &          rdrlsh(nzm), rdrssh(nzm), rurssh(nzm), 
     &          rufrsh(nzm), tdrssh(nzm), turssh(nzm), 
     &          reflec(nzm,2,2), reflecsh(nzm),
     &          dw0(nzm), dw1(nzm), dw2(nzm), dw3(nzm)
      complex   ux(nxm,nym,nzm), uy(nxm,nym,nzm)

      common    /mainpar/ idirect, ifsurf, ix0, iy0, iz0, lc0dn, lc0up,
     &                    nabove, nbelow, nk, nl, nliq, nrliq, nt,
     &                    nx, ny, nz, eps, fimag, fpeak, fref, t0, tl, 
     &                    xl, zbl, zs0, ls0, isrc

      data pi, ai /3.1415926535898, (0.,1.)/

 60   format (1x, "Freq =", i4, " /", i4, 5x, "Nombre d'ondes =", i5,
     &        " /", i5)
 61   format (7x, i4, 27x, i5)
 62   format (1x, 14i5)

!---- LECTURE ET VERIFICATION DES PARAMETRES ET DU MODELE
      fileou1='pression'
      fileou2='deplacx'
      fileou3='deplacy'
      fileou4='deplacz'
      open(15,file=fileou1,form='formatted')
      open(16,file=fileou2,form='formatted')
      open(17,file=fileou3,form='formatted')
      open(18,file=fileou4,form='formatted')

      call mread (nkm, nlm, ntm, nxm, nym, nzm, zb, vp, vs, 
     &            dens, qp, qs, mult, lcount, lr, 
     &            rpos, apos, zpos)       
      
!---- INITIALISATIONS

      nl1 = nl - 1
      write (*,*)
      write (*,*) "*** SKB: Sismogrammes synthetiques en milieu stratifi
     &e elastique"
      write (*,*) "         Version P-SV-SH du 07/12/99 -  M.D. LGIT Gre
     &noble"

!.... Parametres temporels et frequentiels

      pi2   = pi + pi
      dt    = tl/real(nt)     ! Echantillonnage en temps
      nf    = nt/2            ! Nombre de frequences
      df    = 1./tl           ! Echantillonnage en frequence
      aw    = -pi2*fimag      ! Partie imaginaire de la frequence
      omp2  = (pi2*fpeak)**2
      omref = pi2*fref

!.... Parametres relatifs a l'integration en nombres d'ondes discrets

      dk  = pi2/xl            ! Echantillonnage en nombre d'onde
      iflag = 0               ! Idem

!.... Fonction d'amplification temporelle (compensation de la partie 
!     imaginaire de la frequence en fin de calcul)

      ex1 = exp(-aw*dt)
      ex2 = 1.
      do it = 1, nt
         amp(it) = ex2
         ex2 = ex2*ex1
      end do

!.... Calcul des fonctions de Bessel (dependance radiale)

      do ix = 1, nx
         do iy = 1, ny
            do k = 1, nk
!              ak = real(k-1)*dk
               ak = (real(k-1) + 0.218)*dk
               arg(k) = ak*rpos(ix,iy)
            end do      
            call sbj0j19 (arg, nk, bj0, bj1)
            do k = 1, nk
               bbj0(k,ix,iy) = bj0(k)
               bbj1(k,ix,iy) = bj1(k)
               if (arg(k) /= 0.) then
                  bjinc(k,ix,iy) = bj1(k)/arg(k)
               else
                  bjinc(k,ix,iy) = 0.5
               end if
            end do
         end do
      end do

!.... Calcul des termes exponentiels (dependance angulaire)

      do ix = 1, nx
         do iy = 1, ny
            cphi(ix,iy) = cos(apos(ix,iy))
            sphi(ix,iy) = sin(apos(ix,iy))
            if (abs(cphi(ix,iy)) == 1.) sphi(ix,iy) = 0.
            if (abs(sphi(ix,iy)) == 1.) cphi(ix,iy) = 0.
         end do
      end do

!.... Attenuation anelastique

      do la = 1, nl
         ap(la) = 1. + ai/(2.*qp(la))
         if (vs(la) /= 0.) then
            as(la) = 1. + ai/(2.*qs(la))
         end if
      end do

!==================== BOUCLE PRINCIPALE SUR LES FREQUENCES ===================
                                                                             !
      do ifr = 1, nf 

!.... Initialisations

      freq  = real(ifr-1)*df
      rw    = pi2*freq
      omega = cmplx(rw,aw)                            ! Pulsation complexe
      om2   = omega**2
      src(ifr) = om2*cexp(-om2/omp2)                  ! Ondelette de Ricker
      src(ifr) = src(ifr)*cexp(ai*omega*t0)           ! Decalage de t0 sec.

      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx
               up(ix,iy,iz) = 0.
               ur(ix,iy,iz) = 0.
               ut(ix,iy,iz) = 0.
               uz(ix,iy,iz) = 0.
            end do
         end do
      end do

!.... Vitesses complexes (attenuation et dispersion)

!      clnf = -clog(omega/omref)/pi
      clnf = -clog(omega/omega)/pi
      do la = 1, nl
         alpha(la) = vp(la)*(ap(la) + clnf/qp(la))
         wa2(la)   = (omega/alpha(la))**2
         if (vs(la) /= 0.) then
            beta(la) = vs(la)*(as(la) + clnf/qs(la))
            wb2(la)  = (omega/beta(la))**2
            emu(la)  = dens(la)*beta(la)**2
         else
            beta(la) = 0.
            wb2(la)  = 0.
            emu(la)  = 0.
         end if
         bulk(la) = dens(la)*(alpha(la)**2 - (4./3.)*beta(la)**2)        
      end do

!------------------- BOUCLE SUR LES NOMBRES D'ONDES -------------------
                                                                      !
      do k = 1, nk

!.... Nombres d'ondes horizontaux et verticaux 

!        ak  = real(k-1)*dk
         ak  = (real(k-1) + 0.218)*dk
         ak2 = ak*ak
         do la = 1, nl
            wza(la) = csqrt(ak2 - wa2(la))
            if (real(wza(la)) < 0.) wza(la) = -wza(la)
            if (vs(la) /= 0.) then
               w2(la)  = wb2(la) - 2.*ak2
               wzb(la) = csqrt(ak2 - wb2(la))
               if (real(wzb(la)) < 0.) wzb(la) = -wzb(la)
            else
               w2(la)  = 0.
               wzb(la) = 0.
            end if
         end do

!.... Calcul des coefficients de reflexion et de transmission 
!     a toutes les interfaces (l-1/l)

         call rtcoef (nl, nlm, nliq, dens, ifsurf, ak, om2, wa2, wb2, 
     &                wza, wzb, w2, emu, beta, rdi, rui, tdi, tui,
     &                rdish, ruish, tdish, tuish)

!---- CALCUL DES COEFFICIENTS DE REFLEXION-TRANSMISSION GENERALISES

         if (idirect < 2) then     ! Non necessaire pour l'onde directe

!.... Propagation de bas en haut: calcul de RdSL, et de RdRL pour tous
!     les recepteurs situes en-dessous de la source

            lstart = nl1
            zstart = zb(nl1)
            lcount(ls0) = lc0dn + 1
            npos2  = nbelow + 1
            zpos2(1) = zs0
            do iz=nabove+1,nz
               izz = iz - nabove + 1
               zpos2(izz) = zpos(iz)
            end do
            rd_only = 1                ! NB: Seul rdgen est calcule
            call prup (lstart, zstart, npos2, nzm, zpos2, lcount,  
     &                 nl, nlm, zb, vs, wza, wzb, rdi, rui, tdi, tui,
     &                 rdish, ruish, tdish, tuish,
     &                 rd_only, rdgen, tugen, rdgensh, tugensh)    
         
            do i=1,2
               do j=1,2
                  rdsl(i,j) = rdgen(1,i,j)
               end do
            end do
            rdslsh = rdgensh(1)
            do iz=nabove+1,nz
               izz = iz - nabove + 1
               do i=1,2
                  do j=1,2
                     rdrl(iz,i,j) = rdgen(izz,i,j)
                  end do
               end do
               rdrlsh(iz) = rdgensh(izz)
            end do

!.... Propagation de haut en bas: calcul de RuFS, et de RuFR pour tous
!     les recepteurs situes au-dessus de la source

            lstart = 1
            zstart = 0.
            lcount(ls0) = lc0up + 1
            npos2  = nabove + 1
            do iz=1,nabove
               zpos2(iz) = zpos(iz)
            end do
            zpos2(npos2) = zs0
            ru_only = 1                ! NB: Seul rugen est calcule
            call prdn (lstart, zstart, npos2, nzm, zpos2, lcount, 
     &                 nl, nlm, zb, vs, wza, wzb, rdi, rui, tdi, tui,
     &                 rdish, ruish, tdish, tuish,
     &                 ru_only, rugen, tdgen, rugensh, tdgensh)    

            do i=1,2
               do j=1,2
                  rufs(i,j) = rugen(npos2,i,j)
               end do      
            end do
            rufssh = rugensh(npos2)
            do iz=1,nabove
               do i=1,2
                  do j=1,2
                     rufr(iz,i,j) = rugen(iz,i,j)
                  end do      
               end do
               rufrsh(iz) = rugensh(iz)
            end do

         end if

!.... Propagation de bas en haut: calcul de RdRS et TuRS pour tous 
!     les recepteurs situes au-dessus de la source

         if (nabove > 0) then
            lstart = ls0
            zstart = zs0
            lcount(ls0) = lc0up
            rd_only = 0
            call prup (lstart, zstart, nabove, nzm, zpos, lcount, 
     &                 nl, nlm, zb, vs, wza, wzb, 
     &                 rdi, rui, tdi, tui, rdish, ruish, tdish, tuish,
     &                 rd_only, rdgen, tugen, rdgensh, tugensh) 
         
            do iz=1,nabove
               do i=1,2
                  do j=1,2
                     rdrs(iz,i,j) = rdgen(iz,i,j)
                     turs(iz,i,j) = tugen(iz,i,j)
                  end do
               end do
               rdrssh(iz) = rdgensh(iz)
               turssh(iz) = tugensh(iz)
            end do
         end if

!.... Propagation de haut en bas: calcul de RuRS et TdRS pour tous
!     les recepteurs situes en-dessous de la source

         if (nbelow > 0) then
            lstart = ls0
            zstart = zs0
            lcount(ls0) = lc0dn
            do iz=nabove+1,nz
               izz = iz - nabove
               zpos2(izz) = zpos(iz)
            end do
            ru_only = 0
            call prdn (lstart, zstart, nbelow, nzm, zpos2, lcount, 
     &                 nl, nlm, zb, vs, wza, wzb, 
     &                 rdi, rui, tdi, tui, rdish, ruish, tdish, tuish,
     &                 ru_only, rugen, tdgen, rugensh, tdgensh)    

            do iz=nabove+1,nz
               izz = iz - nabove 
               do i=1,2
                 do j=1,2
                     rurs(iz,i,j) = rugen(izz,i,j)
                     tdrs(iz,i,j) = tdgen(izz,i,j)
                  end do      
               end do
               rurssh(iz) = rugensh(izz)
               tdrssh(iz) = tdgensh(izz)
            end do
         end if

         if (idirect < 2) then     ! Non necessaire pour l'onde directe

!.... Assemblage des matrices pour les recepteurs situes au-dessus de la source

            if (nabove > 0) call asmup (nabove, nz, nzm, nrliq, 
     &                      rdrs, rufr, turs, rdsl, rufs, 
     &                      rdrssh, rufrsh, turssh, rdslsh, rufssh,
     &                      reflec, reflecsh)

!.... Assemblage des matrices pour les recepteurs situes en-dessous de la source

            if (nbelow > 0) call asmdn (nabove, nz, nzm, nrliq, 
     &                      rurs, rdrl, tdrs, rufs, rdsl, 
     &                      rurssh, rdrlsh, tdrssh, rufssh, rdslsh, 
     &                      reflec, reflecsh)              

         end if

!---- CONTRIBUTION DE LA SOURCE 

!.... Diagramme de radiation

         dfk = ak*df*dk
         q0   = bulk(ls0)/(2.*dens(ls0)*alpha(ls0)**2)
         q1   = 1./(2.*dens(ls0)*om2)
         if (wzb(ls0) /= 0.) then
            q2 = q1*ak/wzb(ls0)
            q3 = q1*ak/wza(ls0)
            q4 = 1./(2.*dens(ls0)*beta(ls0)*wzb(ls0))
         else
            q2 = 0.
            q3 = 0.
            q4 = 0.
         end if
        
         if (isrc == 0) then              ! Source explosive
            sk1up =  dfk*q0/wza(ls0)    
            sk2up =  0.                 
            sk3up =  0.                 
            sk1dn = -sk1up
            sk2dn =  0.  
            sk3dn =  0.  
         else if (isrc == 1) then         ! Force horizontale suivant x
            sk1up =  dfk*q3/2.             
            sk2up =  dfk*q1/2.          
            sk3up = -dfk*q4*ai/2.
            sk1dn = -sk1up
            sk2dn =  sk2up
            sk3dn = -sk3up
         else if (isrc == 2) then         ! Force horizontale suivant y
            sk1up = -dfk*ai*q3/2.    
            sk2up = -dfk*ai*q1/2. 
            sk3up = -dfk*q4/2.
            sk1dn = -sk1up
            sk2dn =  sk2up
            sk3dn = -sk3up
         else if (isrc == 3) then         ! Force verticale suivant z
            sk1up =  dfk*q1           
            sk2up =  dfk*q2          
            sk3up =  0.                 
            sk1dn =  sk1up
            sk2dn = -sk2up
            sk3dn =  0.  
         end if

         if (idirect < 2) then     ! Non necessaire pour l'onde directe

!.... Contribution pour les recepteurs situes au-dessus de la source:
!     RdSL*SIGd - SIGu

            if (nabove > 0) then
               st1 = rdsl(1,1)*sk1dn + rdsl(1,2)*sk2dn - sk1up
               st2 = rdsl(2,1)*sk1dn + rdsl(2,2)*sk2dn - sk2up
               st3 = rdslsh*sk3dn - sk3up

!     Potentiels des ondes P, SV et SH au niveau des recepteurs
               do iz=1,nabove
                  sfield(iz,1) = reflec(iz,1,1)*st1 + reflec(iz,1,2)*st2
                  sfield(iz,2) = reflec(iz,2,1)*st1 + reflec(iz,2,2)*st2        
                  sfieldsh(iz) = reflecsh(iz)*st3
               end do
            end if

!.... Contribution pour les recepteurs situes en-dessous de la source:
!     SIGd - RuFS*SIGu

            if (nbelow > 0) then
               st1 = sk1dn - rufs(1,1)*sk1up - rufs(1,2)*sk2up 
               st2 = sk2dn - rufs(2,1)*sk1up - rufs(2,2)*sk2up 
               st3 = sk3dn - rufssh*sk3up

!     Potentiels des ondes P, SV et SH au niveau des recepteurs
               do iz=nabove+1,nz
                  sfield(iz,1) = reflec(iz,1,1)*st1 + reflec(iz,1,2)*st2
                  sfield(iz,2) = reflec(iz,2,1)*st1 + reflec(iz,2,2)*st2
                  sfieldsh(iz) = reflecsh(iz)*st3
               end do
            end if

         end if

!---- CALCUL DES CHAMPS D'ONDES AUX RECEPTEURS

!.... Conversion des potentiels en deplacement dans la couche contenant
!     le recepteur iz: matrices Mu et Md

         do iz=1,nz
            lc = lr(iz)
            rmud(iz,1,1) = wza(lc)
            rmud(iz,1,2) = ak
            rmud(iz,2,1) = ak
            rmud(iz,2,2) = wzb(lc)
            if (vs(lc) /= 0.) then
               rmudsh(iz)   = 1./beta(lc)
            else
               rmudsh(iz)   = 0.
            end if
         end do

!.... Conversion des potentiels en deplacement pour les recepteurs
!     situes au-dessus de la source: matrice Mu + Md*RuFR

         if (nabove > 0) then
            do iz=1,nabove
               lc = lr(iz)
               r0(iz,1)   = bulk(lc)*wa2(lc)*(1. + rufr(iz,1,1))
               r0(iz,2)   = bulk(lc)*wa2(lc)*rufr(iz,1,2)
               rc(iz,1,1) = -rmud(iz,1,1) + rmud(iz,1,1)*rufr(iz,1,1) 
     &                                    + rmud(iz,1,2)*rufr(iz,2,1)
               rc(iz,1,2) =  rmud(iz,1,2) + rmud(iz,1,1)*rufr(iz,1,2) 
     &                                    + rmud(iz,1,2)*rufr(iz,2,2)
               rc(iz,2,1) =  rmud(iz,2,1) + rmud(iz,2,1)*rufr(iz,1,1) 
     &                                    + rmud(iz,2,2)*rufr(iz,2,1)
               rc(iz,2,2) = -rmud(iz,2,2) + rmud(iz,2,1)*rufr(iz,1,2)         
     &                                    + rmud(iz,2,2)*rufr(iz,2,2)
               rcsh(iz)   =  rmudsh(iz)   + rmudsh(iz)*rufrsh(iz)
            end do
         end if
            
!.... Conversion des potentiels en deplacement pour les recepteurs
!     situes en-dessous de la source: matrice Md + Mu*RdRL

         if (nbelow > 0) then
            do iz=nabove+1,nz
               lc = lr(iz)
               r0(iz,1)   = bulk(lc)*wa2(lc)*(1. + rdrl(iz,1,1))
               r0(iz,2)   = bulk(lc)*wa2(lc)*rdrl(iz,1,2)
               rc(iz,1,1) =  rmud(iz,1,1) - rmud(iz,1,1)*rdrl(iz,1,1) 
     &                                    + rmud(iz,1,2)*rdrl(iz,2,1)
               rc(iz,1,2) =  rmud(iz,1,2) - rmud(iz,1,1)*rdrl(iz,1,2) 
     &                                    + rmud(iz,1,2)*rdrl(iz,2,2)
               rc(iz,2,1) =  rmud(iz,2,1) + rmud(iz,2,1)*rdrl(iz,1,1) 
     &                                    - rmud(iz,2,2)*rdrl(iz,2,1)
               rc(iz,2,2) =  rmud(iz,2,2) + rmud(iz,2,1)*rdrl(iz,1,2) 
     &                                    - rmud(iz,2,2)*rdrl(iz,2,2)
               rcsh(iz)   =  rmudsh(iz)   + rmudsh(iz)*rdrlsh(iz)
            end do
         end if

         if (idirect == 0 .or. idirect == 2) then

!.... Contribution de l'onde directe pour les recepteurs situes 
!     au-dessus de la source: Mu*TuRS*(-SIGu)

            if (nabove > 0) then
               do iz=1,nabove
                  lc = lr(iz)
                  q1 = -turs(iz,1,1)*sk1up - turs(iz,1,2)*sk2up        
                  q2 = -turs(iz,2,1)*sk1up - turs(iz,2,2)*sk2up
                  q3 = -turssh(iz)*sk3up
                  dw0(iz) = q1*bulk(lc)*wa2(lc)
                  dw1(iz) = -rmud(iz,1,1)*q1 + rmud(iz,1,2)*q2
                  dw2(iz) =  rmud(iz,2,1)*q1 - rmud(iz,2,2)*q2       
                  dw3(iz) =  rmudsh(iz)*q3
               end do
            end if

!.... Contribution de l'onde directe pour les recepteurs situes 
!     en-dessous de la source: Md*TdRS*(SIGd)

            if (nbelow > 0) then
               do iz=nabove+1,nz
                  lc = lr(iz)
                  q1 = tdrs(iz,1,1)*sk1dn + tdrs(iz,1,2)*sk2dn        
                  q2 = tdrs(iz,2,1)*sk1dn + tdrs(iz,2,2)*sk2dn 
                  q3 = tdrssh(iz)*sk3dn
                  dw0(iz) = q1*bulk(lc)*wa2(lc)
                  dw1(iz) = rmud(iz,1,1)*q1 + rmud(iz,1,2)*q2
                  dw2(iz) = rmud(iz,2,1)*q1 + rmud(iz,2,2)*q2       
                  dw3(iz) = rmudsh(iz)*q3
               end do
            end if

         end if

         do iz=1,nz            
            do iy = 1, ny
            do ix = 1, nx

!.... Pression

               if (idirect == 0) then          ! Suppression de l'onde directe
                  q0 = r0(iz,1)*sfield(iz,1) + r0(iz,2)*sfield(iz,2) 
     &                 - dw0(iz)
               else if (idirect == 1) then     ! Champ "complet"
                  q0 = r0(iz,1)*sfield(iz,1) + r0(iz,2)*sfield(iz,2)    
               else if (idirect == 2) then     ! Onde directe uniquement
                  q0 = dw0(iz)
               end if        

               if (isrc == 0 .or. isrc == 3) then  ! Ordre azimutal 0
                  qq0 = q0*bbj0(k,ix,iy)      
               else if (isrc == 1) then            ! Ordre azimutal +/- 1 Fx
                  qq0 = 2.*q0*bbj1(k,ix,iy)*cphi(ix,iy)
               else if (isrc == 2) then            ! Ordre azimutal +/- 1 Fy
                  qq0 = 2.*ai*q0*bbj1(k,ix,iy)*sphi(ix,iy)
               end if
               up(ix,iy,iz) = up(ix,iy,iz) + qq0

!.... Deplacement vertical               

               if (idirect == 0) then          ! Suppression de l'onde directe
                  q1 = rc(iz,1,1)*sfield(iz,1) + 
     &                 rc(iz,1,2)*sfield(iz,2) - dw1(iz)
               else if (idirect == 1) then     ! Champ "complet"
                  q1 = rc(iz,1,1)*sfield(iz,1) + 
     &                 rc(iz,1,2)*sfield(iz,2)        
               else if (idirect == 2) then     ! Onde directe uniquement
                  q1 = dw1(iz)
               end if        

               if (isrc == 0 .or. isrc == 3) then  ! Ordre azimutal 0
                  qq1 = q1*bbj0(k,ix,iy)      
               else if (isrc == 1) then            ! Ordre azimutal +/- 1 Fx
                  qq1 = 2.*q1*bbj1(k,ix,iy)*cphi(ix,iy)
               else if (isrc == 2) then            ! Ordre azimutal +/- 1 Fy
                  qq1 = 2.*ai*q1*bbj1(k,ix,iy)*sphi(ix,iy)
               end if
               uz(ix,iy,iz) = uz(ix,iy,iz) + qq1

!.... Deplacements horizontaux radial et tangentiel 

            if (iz > nrliq) then

               if (idirect == 0) then          ! Suppression de l'onde directe
                  q2 = rc(iz,2,1)*sfield(iz,1) + 
     &                 rc(iz,2,2)*sfield(iz,2) - dw2(iz)
                  q3 = rcsh(iz)*sfieldsh(iz)   - dw3(iz)
               else if (idirect == 1) then     ! Champ "complet"
                  q2 = rc(iz,2,1)*sfield(iz,1) + 
     &                 rc(iz,2,2)*sfield(iz,2) 
                  q3 = rcsh(iz)*sfieldsh(iz)
               else if (idirect == 2) then     ! Onde directe uniquement 
                  q2 = dw2(iz) 
                  q3 = dw3(iz) 
               end if       

               if (isrc == 0 .or. isrc == 3) then  ! Ordre azimutal 0
                  qq2 = -q2*bbj1(k,ix,iy)
                  qq3 =  q3*bbj1(k,ix,iy)
               else if (isrc == 1) then            ! Ordre azimutal +/- 1 Fx
                  qq2 = 2.*(q2*bbj0(k,ix,iy) + 
     &                  (ai*q3 - q2)*bjinc(k,ix,iy))*cphi(ix,iy)
                  qq3 = 2.*ai*(-q3*bbj0(k,ix,iy) +
     &                  (q3 + ai*q2)*bjinc(k,ix,iy))*sphi(ix,iy)
               else if (isrc == 2) then            ! Ordre azimutal +/- 1 Fy
                  qq2 = 2.*ai*(q2*bbj0(k,ix,iy) + 
     &                  (ai*q3 - q2)*bjinc(k,ix,iy))*sphi(ix,iy)
                  qq3 = 2.*(-q3*bbj0(k,ix,iy) +
     &                  (q3 + ai*q2)*bjinc(k,ix,iy))*cphi(ix,iy)
               end if
               ur(ix,iy,iz) = ur(ix,iy,iz) + qq2
               ut(ix,iy,iz) = ut(ix,iy,iz) + qq3
 
            else
               ur(ix,iy,iz) = 0.
               ut(ix,iy,iz) = 0.
            end if

!.... Deplacements horizontaux en coordonnees cartesiennes

            ux(ix,iy,iz) = ur(ix,iy,iz)*cphi(ix,iy) - 
     &                     ut(ix,iy,iz)*sphi(ix,iy)       
            uy(ix,iy,iz) = ur(ix,iy,iz)*sphi(ix,iy) + 
     &                     ut(ix,iy,iz)*cphi(ix,iy)       

!.... Test de convergence de la serie en nombres d'ondes 
!     au recepteur (ix0,iy0,iz0)

               if (ix == ix0 .and. iy == iy0 .and. iz == iz0) then
                  a0  = abs(real(qq0)) + abs(aimag(qq0))
                  a1  = abs(real(qq1)) + abs(aimag(qq1))
                  a2  = abs(real(qq2)) + abs(aimag(qq2))
                  a3  = abs(real(qq3)) + abs(aimag(qq3))
                  a00 = (abs( real(up(ix,iy,iz))) + 
     &                   abs(aimag(up(ix,iy,iz))))*eps 
                  a10 = (abs( real(uz(ix,iy,iz))) + 
     &                   abs(aimag(uz(ix,iy,iz))))*eps        
                  a20 = (abs( real(ux(ix,iy,iz))) +
     &                   abs(aimag(ux(ix,iy,iz))))*eps        
                  a30 = (abs( real(uy(ix,iy,iz))) + 
     &                   abs(aimag(uy(ix,iy,iz))))*eps
               end if
            end do
            end do
         end do
         if (k > 1 .and. a1 <= a10 .and. a2 <= a20 .and. a3 <= a30) exit       

      end do                            
                                                                      !
!-------------- FIN DE LA BOUCLE SUR LES NOMBRES D'ONDES --------------

!.... Controle de la convergence

      k1 = k - 1
      nkconv(ifr) = k1
      if (mod(ifr,10) == 1) then
         write (*,*)
         write (*,60) ifr, nf, k1, nk
      else
         write (*,61) ifr, k1
      end if

      if (k1 >= nk) then
         if (iflag == 0 .and. nz > 1) then
            write (*,*)     
            write (*,*) "SKB *** SERIE k TRONQUEE - Essai de correction        
     &***"         
            write (*,*) "Avant correction - ix0, iy0, iz0 : ", 
     &         ix0, iy0, iz0    
            if (iz0 > 1) then
               iz0 = iz0 - 1
            else
               iz0 = iz0 + 1
            end if 
            write (*,*) "Apres correction - ix0, iy0, iz0 : ", 
     &         ix0, iy0, iz0 
            write (*,*)
         else
            stop "SKB *** PAS DE CONVERGENCE ***" 
         end if     
         iflag = 1
      end if

!.... Reponse du milieu stratifie a la frequence Omega

      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx
               ufp(ifr,ix,iy,iz) = up(ix,iy,iz)
               ufx(ifr,ix,iy,iz) = ux(ix,iy,iz)
               ufy(ifr,ix,iy,iz) = uy(ix,iy,iz)
               ufz(ifr,ix,iy,iz) = uz(ix,iy,iz)
            end do
         end do
      end do

      end do                             
                                                                             !
!============== FIN DE LA BOUCLE PRINCIPALE SUR LES FREQUENCES ===============

      write (3,*)
      write (3,*) " Resume - Integration en k:"
      write (3,*)
      write (3,62) (nkconv(ifr), ifr = 1, nf)
      write (3,*)
      write (3,*) " ufp(1,ix0,iy0,iz0) :", ufp(1,ix0,iy0,iz0)
      write (3,*) " ufx(1,ix0,iy0,iz0) :", ufx(1,ix0,iy0,iz0)
      write (3,*) " ufy(1,ix0,iy0,iz0) :", ufy(1,ix0,iy0,iz0)
      write (3,*) " ufz(1,ix0,iy0,iz0) :", ufz(1,ix0,iy0,iz0)
      write (3,*)

!.... Transformee de Fourier Inverse

      do ic = 0, 3
         inum = 0
         do ix = 1, nx
            do iy = 1, ny
               do iz = 1, nz
                  do ifr = 1, nf
                     if (ic == 0) x(ifr) = ufp(ifr,ix,iy,iz)*src(ifr)       
                     if (ic == 1) x(ifr) = ufx(ifr,ix,iy,iz)*src(ifr)
                     if (ic == 2) x(ifr) = ufy(ifr,ix,iy,iz)*src(ifr)
                     if (ic == 3) x(ifr) = ufz(ifr,ix,iy,iz)*src(ifr)
                  end do

                  call crfft9 (x, nt, y, 1)

!.... Compensation de la partie imaginaire de la frequence
                  do it = 1, nt
                     y(it) = y(it)*amp(it)*sqrt(dt*tl*nt)
                     write(15+ic,*)it*dt,y(it)
                  end do

!.... Sortie des resultats (fichiers binaires en access direct)
 
                  inum = inum + 1
!                 write (10+ic, rec=inum) (y(it), it = 1, nt)
                  write (10+ic) (y(it), it = 1, nt)
               end do
            end do
         end do
      end do

      write (3,*) " *** FIN NORMALE *** "
      write (3,*)

      close (1)
      close (2)
      close (3)
      close (10)
      close (11)
      close (12)
      close (13)

      end program skb
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine mread (nkm, nlm, ntm, nxm, nym, nzm, zb, vp, vs, 
     &                  dens, qp, qs, mult, lcount, lr, 
     &                  rpos, apos, zpos)                            
!     LECTURE ET VERIFICATION DES PARAMETRES EN ENTREE DU PROGRAMME SKB
!     Langage F77 avec extensions F90 XLF (AIX IBM RS/6000 Version 3.2.4)

!=============================================================================

      implicit  none
      integer   i, idirect, ifsurf, ix, iy, iz, ix0, iy0, iz0, isrc, la,       
     &          lc0dn, lc0up, ls0, nabove, nbelow, nkm, nk, nl, nlm, 
     &          nrliq, nliq, nt, ntm, nx, ny, nz, nxm, nym, nzm, 
     &          lcount(nlm), lr(nzm), mult(nlm)
      real      dm, dmax, dx, dy, dz, eps, fimag, fnyq, fpeak, fref, pi, 
     &          rmax, t0, tl, xl, zbl, zbnl1, zrc, zs0, vmin, vmax, 
     &          xx, yy, xm, ym, xmin, xmax, ymin, ymax, zmin, zmax, 
     &          zb(0:nlm), vp(nlm), vs(nlm), dens(nlm), qp(nlm),qs(nlm),      
     &          rpos(nxm,nym), apos(nxm,nym), zpos(nzm)      
      character*20 str1, str2, str3, str4
      character*40 datfilp, datfilx, datfily, datfilz, logfile, modfile

      common   /mainpar/ idirect, ifsurf, ix0, iy0, iz0, lc0dn, lc0up,
     &                   nabove, nbelow, nk, nl, nliq, nrliq, nt,
     &                   nx, ny, nz, eps, fimag, fpeak, fref, t0, tl, 
     &                   xl, zbl, zs0, ls0, isrc
      data pi  /3.1415926535898/
      namelist /userpar/ nl, nt, tl, t0, fpeak, zs0, nx, xmin, xmax, 
     &                   ny, ymin, ymax, nz, zmin, zmax 
      namelist /morepar/ isrc, idirect, ifsurf, modfile, logfile, 
     &                   datfilp, datfilx, datfily, datfilz, 
     &                   nk, ix0, iy0, iz0, xl, eps, fimag, fref       

!---- LECTURE DES DONNEES

!.... Parametres generaux du calcul
      
      open (1, file = "skb.nml", form = "formatted")
      read (1, nml  = userpar)
      read (1, nml  = morepar)
      open (2, file = modfile, form = "formatted")
      open (3, file = logfile, form = "formatted")
      open (10,file = datfilp, form = "unformatted") 
!    &      access = "direct", recl = 4*nt)
      open (11,file = datfilx, form = "unformatted")
!    &      access = "direct", recl = 4*nt)
      open (12,file = datfily, form = "unformatted") 
!    &      access = "direct", recl = 4*nt)
      open (13,file = datfilz, form = "unformatted")
!     &      access = "direct", recl = 4*nt)
  
!.... Caracteristiques du modele stratifie

      read (2, *) (zb(la), vp(la), vs(la), dens(la), qp(la), qs(la),        
     &             mult(la), la = 1, nl)           
      zb(0)  = 0.
      zb(nl) = 9999.

!---- VERIFICATION DES TABLEAUX        
      if (nk > nkm .or. nk > nkm .or. nl > nlm .or. nt > ntm .or.
     &   nx > nxm .or. ny > nym .or. nz > nzm) then
         stop  "MREAD *** TABLEAUX MAL DIMENSIONNES ***"
      end if

!---- VERIFICATION DES DONNEES ET INITIALISATIONS                       

!.... Positions des capteurs 

      if (nx == 0 .or. ny == 0 .or. nz == 0) then
         stop  "MREAD *** Parametre nx, ny ou nz incorrect ***"
      end if

      dx = 0.
      dy = 0.
      dz = 0.
      if (nx > 1) dx = (xmax - xmin)/real(nx-1)
      if (ny > 1) dy = (ymax - ymin)/real(ny-1)
      if (nz > 1) dz = (zmax - zmin)/real(nz-1)
      
      xx = xmin
      yy = ymin
      do ix=1,nx
         xx = xmin + real(ix-1)*dx
         do iy=1,ny
            yy = ymin + real(iy-1)*dy
            rpos(ix,iy) = sqrt(xx*xx + yy*yy)
            apos(ix,iy) = atan2(yy,xx)
         end do
      end do

      do iz=1,nz
         zpos(iz) = zmin + real(iz-1)*dz
      end do

!.... Frequence dominante trop elevee ?

      fnyq = real(nt/2)/tl
      if (fpeak == 0.) fpeak = fnyq/3.
      if (fpeak > fnyq/2.) then
         write (*,*) 'fpeak, fnyquist :', fpeak, fnyq
         stop  "MREAD *** FREQUENCE fpeak TROP ELEVEE ***"
      end if

!.... Partie imaginaire de la frequence

      if (fimag == 0.) fimag = 1./tl

!.... Frequence de reference pour la dispersion

      if (fref == 0.) fref = 1.

!.... Parametre eps (test de convergence)

      if (eps == 0.) eps = 1.e-06

!.... Parametre idirect (onde directe)

      if (idirect == 0) then
         str1 = "Non prise en compte"
      else if (idirect == 1) then
         str1 = "Prise en compte"
      else if (idirect == 2) then
         str1 = "Uniquement"
      else
         stop "MREAD *** PARAMETRE idirect MAL DEFINI ***"
      end if

!.... Parametre ifsurf (reflection a la surface libre) 

      if (ifsurf == 0) then
         str2 = "Non prise en compte"
      else if (ifsurf == 1) then 
         str2 = "Prise en compte"
      else
         stop "MREAD *** PARAMETRE ifsurf MAL DEFINI ***"
      end if

!.... Parametre isrc (type de source)

      if (isrc == 0) then
         str3 = "Explosion"
      else if (isrc == 1) then
         str3 = "Force suivant x"
      else if (isrc == 2) then
         str3 = "Force suivant y"
      else if (isrc == 3) then
         str3 = "Force suivant z"
      else
         stop  "MREAD *** PARAMETRE isrc MAL DEFINI ***"
      end if

!.... Parametres ix0, iy0 et iz0 (test de convergence)
!     On prend le recepteur "du milieu" si ix0 = iy0 = iz0 = 0

      if (ix0 > nx .or. iy0 > ny .or. iz0 > nz) then
         stop  "MREAD *** PARAMETRE ix0 ou ix0 ou iz0 MAL DEFINI ***"
      end if
      if (ix0 == 0) ix0 = nx/2 + 1
      if (iy0 == 0) iy0 = ny/2 + 1
      if (iz0 == 0) iz0 = nz/2 + 1

      zbnl1 = zb(nl-1)
      do iz = 1, nz    

!.... Source et recepteur au meme niveau (pas de convergence) ?

         if (zpos(iz) == zs0) then
            write (*,*) 'iz, zpos(iz), zs0 :', iz, zpos(iz), zs0
            stop  "MREAD *** SOURCE ET CAPTEUR AU MEME NIVEAU ***"
         end if

!.... Recepteur dans le demi espace ?

         if (zpos(iz) > zbnl1) then
            write (*,*) 'iz, zpos(iz), zb(nl-1) :', iz, zpos(iz), 
     &                   zbnl1       
            stop  "MREAD *** CAPTEUR(S) DANS LE DEMI-ESPACE ***" 
         end if
     
       end do

!.... Source dans le demi-espace ?

      if (zs0 > zbnl1) then
         write (*,*) 'zs0, zb(nl-1) :', zs0, zbnl1  
         stop  "MREAD *** SOURCE DANS LE DEMI-ESPACE ***"
      end if

      if (vs(1) == 0.) then
         nliq = 1
         zbl  = zb(1)
         str4 = "Oui"
      else
         nliq = 0
         zbl  = 0.
         str4 = "Non"
      end if
      ls0   = 1
      do la = 1, nl

!.... Vitesses incorrectes 

         if (vs(la) >= vp(la) .or. vp(la) <= 0.) then
            write (*,*) 'la, vp(la), vs(la) :', la, vp(la), vs(la)
            stop  "MREAD *** VITESSES INCORRECTES ***"
         end if

!.... Pas d'attenuation si Qp ou Qs est nul

         if (qp(la) == 0.) qp(la) = 100000.
         if (qs(la) == 0.) qs(la) = 100000.

!.... Couches d'epaisseur negative ?

         if (la < nl) then
            if (zb(la) <= zb(la-1)) then
               write (*,*) 'la, zb(la) :', la, zb(la)
               stop  "MREAD *** COUCHE D'EPAISSEUR NEGATIVE ***"
            end if
         end if

!.... Source au niveau d'une interface ?

         if (la < nl .and. zs0 == zb(la)) then
            write (*,*) 'la, zb(la), zs0 :', la, zb(la), zs0
            stop  "MREAD *** SOURCE A UNE INTERFACE ***"
         end if
         if (la < nl .and. zs0 > zb(la-1) .and. zs0 < zb(la)) 
     &       ls0 = la  
       
         lcount(la) = 0
      end do

!.... Recepteur au niveau d'une interface ?
!     Numero de la couche contenant le recepteur (IR), 
!     nombre de recepteurs par couches

      nabove = 0
      nrliq  = 0
      do iz = 1, nz
         lr(iz) = 0
         zrc = zpos(iz)
         if (zrc <= zs0) nabove = nabove + 1
         if (nliq == 1 .and. zrc <= zbl) nrliq  = nrliq  + 1
         do la = 1, nl-1
            if (zrc == zb(la)) then
               write (*,*) 'la, zb(la), iz, zpos(iz) :', 
     &                      la, zb(la), iz, zpos(iz)           
               stop   "MREAD *** CAPTEUR A UNE INTERFACE ***"
            end if
            if (zrc > zb(la-1) .and. zrc < zb(la)) then
               lr(iz) = la 
               lcount(la) = lcount(la) + 1
            end if
         end do
      end do

      nbelow = nz - nabove
      lc0up = 0
      lc0dn = 0
      do iz = 1, nz
         la = lr(iz)
         if (la == ls0 .and. zpos(iz) < zs0) lc0up = lc0up + 1
         if (la == ls0 .and. zpos(iz) > zs0) lc0dn = lc0dn + 1
      end do

!.... Vitesses minimum et maximum

      vmin = 20000.
      vmax = 0.
      do la = 1, nl
         if (vp(la) < vmin)                    vmin = vp(la)
         if (vs(la) < vmin .and. vs(la) /= 0.) vmin = vs(la)  
         if (vp(la) > vmax)                    vmax = vp(la)
      end do

!.... Valeurs calculees de xl et de nk (si non specifiees)

      xm = amax1(abs(xmin),abs(xmax))
      ym = amax1(abs(ymin),abs(ymax))
      rmax = sqrt(xm*xm + ym*ym)
      if (xl == 0.) xl = 1.2*(vmax*(tl + t0) + rmax) 
      if (nk == 0)  nk = int(1.2*fnyq*xl/vmin)

!---- IMPRESSION DES DONNEES 

      write (3,*)
      write (3,*)
      write (3,*) "     ================================================
     &=============="
      write (3,*) "    |      PROPAGATION D'ONDES P-SV-SH EN MILIEU STRA
     &TIFIE PLAN    |"
      write (3,*) "    |                 Algorithme KB (Kennett-Bouchon)           
     &              |"
      write (3,*) "     ================================================
     &=============="         
      write (3,*)
      write (3,*)
      write (3,*) " --- PARAMETRES DU CALCUL --- "
      write (3,*)
      write (3,'("      Nombre total de couches dans le modele    : ",
     &                  i5)') nl
      write (3,*) "     Couche d'eau en surface                   : ",
     &                  str4
      if (nliq == 1) then
      write (3,'("      Epaisseur de la couche d''eau         [km] : ",
     &                  f10.4)') zbl
      end if
      write (3,*) 
      write (3,'("      Duree des sismogrammes                [s] : ",
     &                  f10.4)') tl
      write (3,'("      Decalage en temps                     [s] : ",
     &                  f10.4)') t0
      write (3,'("      Nombre d''echantillons en temps            : ",               
     &                  i5)') nt
      write (3,'("      Frequence de Nyquist du signal       [Hz] : ",
     &                  f10.4)') fnyq          
      write (3,'("      Frequence dominante  du signal       [Hz] : ",
     &                  f10.4)') fpeak          
      write (3,'("      Frequence de ref. pour la dispersion [Hz] : ",
     &                  f10.3)') fref
      write (3,*)                                          
      write (3,'("      Vitesse minimum dans le milieu      [m/s] : ",
     &                  f10.3)') vmin*1000.
      write (3,'("      Vitesse maximum dans le milieu      [m/s] : ",
     &                  f10.3)') vmax*1000.
      write (3,'("      Longueur d''onde minimum (@ fpeak)     [m] : ",
     &                  f10.3)') vmin*1000./fpeak
      write (3,'("      Longueur d''onde maximum (@ fpeak)     [m] : ",
     &                  f10.3)') vmax*1000./fpeak
      write (3,*)
      write (3,'("      Nombre total de positions capteurs        : ",
     &                  i5)') nx*ny*nz
      write (3,'("      Nombre de positions x                     : ",
     &                  i5)') nx
      write (3,'("      Nombre de positions y                     : ",
     &                  i5)') ny
      write (3,'("      Nombre de positions z                     : ",
     &                  i5)') nz
      write (3,'("      Nombre de cotes au-dessus de la source    : ",
     &                  i5)') nabove
      write (3,'("      Nombre de cotes en-dessous de la source   : ",
     &                  i5)') nbelow
      if (nliq == 1) then
      write (3,'("      Nombre de capteurs dans la couche d''eau   : ",
     &                  i5)') nrliq*nx*ny
      end if
      write (3,*)
      write (3,'("      Abscisse minimum   des capteurs      [km] : ",
     &                  f10.4)') xmin         
      write (3,'("      Ordonnee minimum   des capteurs      [km] : ",
     &                  f10.4)') ymin                 
      write (3,'("      Profondeur minimum des capteurs      [km] : ",
     &                  f10.4)') zmin       
      write (3,'("      Couche No                                 : ",
     &                  i5)') lr(1)
      write (3,*)
      write (3,'("      Abscisse maximum   des capteurs      [km] : ",
     &                  f10.4)') xmax         
      write (3,'("      Ordonnee maximum   des capteurs      [km] : ",
     &                  f10.4)') ymax            
      write (3,'("      Profondeur maximum des capteurs      [km] : ",
     &                  f10.4)') zmax       
      write (3,'("      Couche No                                 : ",
     &                  i5)') lr(nz)
      write (3,*)
      write (3,'("      Intervalle en x des capteurs          [m] : ",
     &                  f10.4)') dx*1000.     
      write (3,'("      Intervalle en y des capteurs          [m] : ",
     &                  f10.4)') dy*1000.        
      write (3,'("      Intervalle en z des capteurs          [m] : ",
     &                  f10.4)') dz*1000.   
      write (3,*)      
      write (3,'("      Profondeur de la source              [km] : ",
     &                  f10.4)') zs0
      write (3,'("      Couche No                                 : ",
     &                  i5)') ls0
      write (3,*)
      write (3,*) "     Source ponctuelle                         : ",
     &                  str3 
      write (3,*) "     Onde directe                              : ",
     &                  str1
      write (3,*) "     Reflexion a la surface libre              : ",
     &                  str2 
      write (3,*)                             
      write (3,'("      Distance de discretisation           [km] : ",
     &                  f10.4)') xl
      write (3,'("      Nombre maximum de nombres d''onde          : ",       
     &                  i5)') nkm 
      write (3,'("      Convergence au recepteur (ix0,iy0,iz0)    : ",
     &                  i5,",",i5,",",i5)') ix0, iy0, iz0               
      write (3,*)                             
      write (3,*) " --- MODELE STRATIFIE --- "
      write (3,*)
      write (3,*) "    Zbase        Vp        Vs   Densite      Qp 
     &  Qs   Mult "
      write (3,*) "     (km)    (km/s)    (km/s)   (g/cm3)"
      write (3,*) "_____________________________________________________
     &_____________"
      write (3,*)
      write (3, '(4f10.4, 2f10.1, i5)') (zb(la), vp(la), vs(la), 
     &      dens(la), qp(la), qs(la), mult(la), la = 1, nl)        
      write (3,*)

      end subroutine mread
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine rtcoef (nl, nlm, nliq, dens, ifsurf, ak, om2, 
     &                   wa2, wb2, wza, wzb, w2, emu, beta, 
     &                   rdi, rui, tdi, tui,
     &                   rdish, ruish, tdish, tuish)              

!     CALCUL DES COEFFICIENTS DE REFLEXION ET DE TRANSMISSION A TOUTES
!     LES INTERFACES D'UN MILIEU STRATIFIE PLAN
!     Langage F77 avec extensions F90 XLF (AIX IBM RS/6000 Version 3.2.4)

!     Convention pour les matrices de coefficients de reflexion-transmission: 
!        (1,1): P inc. -> P diffr.       (1,2): S inc. -> P diffr.
!        (2,1): P inc. -> S diffr.       (2,2): S inc. -> S diffr.

!=============================================================================

      implicit  none
      integer   i, ifs, ifsurf, ibot, itop, j, la, nl, nliq, nlm
      real      ak, ak2, dens(nlm)
      complex   c1, c2, c3, c4, ce1, ce2, cf1, cf2, cg, ch, delta,
     &          ee1, ee2, ef1, ef2, eg, eh, om2, q1, q2, q3, ray, 
     &          sh1, sh2, sto,
     &          beta(nlm), emu(nlm), w2(nlm), wa2(nlm), wb2(nlm),
     &          wza(nlm), wzb(nlm), 
     &          rdi(nlm,2,2), rui(nlm,2,2), tdi(nlm,2,2), tui(nlm,2,2),
     &          rdish(nlm), ruish(nlm), tdish(nlm), tuish(nlm)        

!---- INITIALISATIONS

      ak2 = ak*ak
      do la = 1, nl
         do i = 1, 2
            do j = 1, 2
               rdi(la,i,j) = 0.
               rui(la,i,j) = 0.
               tdi(la,i,j) = 0.
               tui(la,i,j) = 0.
            end do
         end do
         rdish(la) = 0.
         ruish(la) = 0.
         tdish(la) = 0.
         tuish(la) = 0.
      end do

!---- COEFFICIENTS DE REFLEXION (1) A LA SURFACE LIBRE

      ifs = 1
      if (ifsurf == 1) then                    ! Reflexions prises en compte
         if (nliq == 0) then                   ! Couche elastique en surface
            q1  = -4.*ak2*wza(ifs)*wzb(ifs)
            q2  = w2(ifs)**2
            ray = q1 + q2
            q3  = -4.*ak*w2(ifs)/ray
            rui(ifs,1,1) = (q1 - q2)/ray
            rui(ifs,1,2) = q3*wzb(ifs)
            rui(ifs,2,1) = q3*wza(ifs)
            rui(ifs,2,2) = rui(ifs,1,1)
            ruish(ifs)   = 1.
         else if (nliq == 1) then              ! Couche d'eau en surface
            rui(ifs,1,1) = -1.
         end if
      end if

!---- COEFFICIENTS DE REFLEXION ET DE TRANSMISSION (L) A L'INTERFACE (L-1/L)

      do la = 2, nl
         itop = la - 1
         ibot = la

!.... Interface Fluide/Solide

         if (la == 2 .and. nliq == 1) then
            c1  = emu(ibot)**2
            c2  = 4.*ak2*wza(ibot)*wzb(ibot)
            c3  = w2(ibot)**2
            c4  = c1*wza(itop)
            ce1 = om2*dens(itop)
            ce2 = om2*dens(ibot)
            ee1 = ce1*ce2*wza(ibot)
            ef1 = c4*(c3 - c2)
            ef2 = c4*(c3 + c2)
            eg  = 2.*emu(ibot)*ce1*wza(itop)
            eh  = 4.*c1*ak
            delta = ee1 + ef1

            rdi(la,1,1) = (-ee1 + ef1)/delta
            rui(la,1,1) = ( ee1 - ef2)/delta
            rui(la,1,2) = -(eh*w2(ibot)*wza(itop)*wzb(ibot))/delta
            rui(la,2,1) = rui(la,1,2)*wza(ibot)/wzb(ibot)
            rui(la,2,2) = (-ee1 - ef2)/delta
            ruish(la)   = 1.

            tdi(la,1,1) = eg*w2(ibot)/delta
            tdi(la,2,1) = 2.*eg*ak*wza(ibot)/delta
            tui(la,1,1) = tdi(la,1,1)*dens(ibot)*wza(ibot)
     &                              /(dens(itop)*wza(itop))
            tui(la,1,2) = tdi(la,2,1)*dens(ibot)*wzb(ibot)
     &                              /(dens(itop)*wza(itop))

            cycle
         end if

!.... Interface Solide/Solide

         c1  = wza(itop)*wzb(itop)
         c2  = wza(ibot)*wzb(ibot)
         c3  = wza(itop)*wzb(ibot)
         c4  = wza(ibot)*wzb(itop)
         ce1 = emu(itop)*wb2(itop)
         ce2 = emu(ibot)*wb2(ibot)
         cf1 = emu(ibot)* w2(ibot)
         cf2 = emu(itop)* w2(itop)
         cg  = cf1 - cf2
         cf1 = cf1 + 2.*ak2*emu(itop)
         cf2 = cf2 + 2.*ak2*emu(ibot)
         ch  = emu(ibot) - emu(itop)
         ee1 = ce1*ce2*(c3 + c4)
         ee2 = ce1*ce2*(c3 - c4)
         ef1 = cf1*cf1*c1
         ef2 = cf2*cf2*c2
         eg  = ak2*cg*cg
         eh  = 4.*ak2*c1*c2*ch*ch

         sto = ee1 + ef1 + ef2 - eg - eh
         q1  =       ef1 - ef2 + eg - eh
         q2  =     - ef1 + ef2 + eg - eh
         q3  = 2.*ce1/sto

         sh1 = emu(itop)*wzb(itop)
         sh2 = emu(ibot)*wzb(ibot) 
            
         rdi(la,1,1) = (q1 + ee2)/sto
         rdi(la,1,2) = 2.*ak*wzb(itop)*(cf1*cg - 2.*c2*cf2*ch)/sto        
         rdi(la,2,1) = rdi(la,1,2)*wza(itop)/wzb(itop)
         rdi(la,2,2) = (q1 - ee2)/sto
         rdish(la)   = (sh1 - sh2)/(sh1 + sh2)
   
         if (la == nl) exit

         rui(la,1,1) = (q2 - ee2)/sto
         rui(la,1,2) = 2.*ak*wzb(ibot)*(cf2*cg - 2.*c1*cf1*ch)/sto
         rui(la,2,1) = rui(la,1,2)*wza(ibot)/wzb(ibot)
         rui(la,2,2) = (q2 + ee2)/sto
         ruish(la)   = -rdish(la)

         tdi(la,1,1) = q3*wza(itop)*(cf1*wzb(itop) + cf2*wzb(ibot))
         tdi(la,1,2) = q3*wzb(itop)*ak*(cg + 2.*c3*ch)
         tdi(la,2,1) = q3*wza(itop)*ak*(cg + 2.*c4*ch)
         tdi(la,2,2) = q3*wzb(itop)*(cf1*wza(itop) + cf2*wza(ibot))
         tdish(la)   = 2.*wzb(itop)*dens(itop)*beta(itop)*beta(ibot)        
     &                 /(sh1 + sh2)

         tui(la,1,1) = tdi(la,1,1)*dens(ibot)*wza(ibot)
     &                           /(dens(itop)*wza(itop))
         tui(la,1,2) = tdi(la,2,1)*dens(ibot)*wzb(ibot)
     &                           /(dens(itop)*wza(itop))
         tui(la,2,1) = tdi(la,1,2)*dens(ibot)*wza(ibot)
     &                           /(dens(itop)*wzb(itop))
         tui(la,2,2) = tdi(la,2,2)*dens(ibot)*wzb(ibot)
     &                           /(dens(itop)*wzb(itop))
         tuish(la)   =   tdish(la)*dens(ibot)*wzb(ibot)
     &                           /(dens(itop)*wzb(itop))

      end do 

      end subroutine rtcoef  
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine prup (lstart, zstart, nz, nzm, zz, lcount, 
     &                 nl, nlm, zb, vs, wza, wzb, 
     &                 rdi, rui, tdi, tui, rdish, ruish, tdish, tuish,       
     &                 rd_only, rdgen, tugen, rdgensh, tugensh)              

!     PROPAGATION VERS LE HAUT D'UNE ONDE PLANE DANS UN MILIEU
!     STRATIFIE PLAN
!     Langage F77 avec extensions F90 XLF (AIX IBM RS/6000 Version 3.2.4)

!     Convention pour les matrices de coefficients de reflexion-transmission: 
!        (1,1): P inc. -> P diffr.       (1,2): S inc. -> P diffr.
!        (2,1): P inc. -> S diffr.       (2,2): S inc. -> S diffr.

!     RdI, RuI, TdI, TuI : coefficients de reflexion-transmission 
!                          a l'interface (L-1/L)
!     RdG, TuG           : coefficients de reflexion-transmission generalises
!                          apres passage d'une interface
!     Rd , Tu            : coefficients de reflexion-transmission generalises
!                          apres dephasage

!=============================================================================

      implicit none
      integer  i, ip, iz, j, l1, la, lc, lstart, nl, nlm, nz, nzm, 
     &         rd_only, lcount(nlm)
      real     dzp, zref, zstart,
     &         vs(nlm), zb(0:nlm), zz(nzm+1) 
      complex  det, pha, phb, q1, q2, 
     &         ash, bsh, rdsh, tush, tugsh, rdgsh,
     &         a(2,2), b(2,2), rd(2,2), tu(2,2), rdg(2,2), tug(2,2), 
     &         rdi(nlm,2,2), rui(nlm,2,2), tdi(nlm,2,2), tui(nlm,2,2),
     &         rdish(nlm), ruish(nlm), tdish(nlm), tuish(nlm),         
     &         wza(nlm), wzb(nlm), 
     &         rdgen(nzm+1,2,2), tugen(nzm+1,2,2), 
     &         rdgensh(nzm+1), tugensh(nzm+1)                     

!---- INITIALISATIONS

      do i = 1, 2
         do j = 1, 2
            rdg(i,j) = 0.
            tug(i,j) = 0.
         end do
      end do
      tug(1,1) = 1.
      tug(2,2) = 1.
      rdgsh = 0.
      tugsh = 1.      

      l1  = lstart + 1
      if (zstart == zb(lstart)) then
         rdg(1,1) = rdi(l1,1,1)
         rdg(1,2) = rdi(l1,1,2)
         rdg(2,1) = rdi(l1,2,1)
         rdg(2,2) = rdi(l1,2,2)
         rdgsh    = rdish(l1)
      end if

!==== APPLICATION DES RELATIONS ITERATIVES
!     1) Dephasage des coefficients jusqu'au sommet de la couche (L)
!     2) Passage de l'interface (L-1/L)

      iz = nz
      do la = lstart,1,-1              ! Propagation du bas vers le haut
         if (la == lstart) then
            zref = zstart
         else
            zref = zb(la)
         end if

!---- 1A) Dephasages pour les recepteurs se trouvant dans la couche courante
             
         lc = lcount(la)
         if (lc /= 0) then
            do ip = 1, lc              
               dzp = zref - zz(iz)
               pha = -wza(la)*dzp
               q1  = cexp(pha)
               if (vs(la) /= 0.) then
                  phb = -wzb(la)*dzp
                  q2  = cexp(phb)
               end if
               
!.... Rd = E*RdG*E
               rdgen(iz,1,1) = rdg(1,1)*q1*q1
               rdgen(iz,1,2) = rdg(1,2)*q1*q2
               rdgen(iz,2,1) = rdg(2,1)*q1*q2
               rdgen(iz,2,2) = rdg(2,2)*q2*q2
               rdgensh(iz)   = rdgsh*q2*q2

!.... Tu = E*TuG
               if (rd_only /= 1) then               
                  tugen(iz,1,1) = tug(1,1)*q1
                  tugen(iz,1,2) = tug(1,2)*q1
                  tugen(iz,2,1) = tug(2,1)*q2
                  tugen(iz,2,2) = tug(2,2)*q2
                  tugensh(iz)   = tugsh*q2
               end if

               iz = iz - 1
               if (iz < 1) return
            end do
         end if

!---- 1B) Dephasage dans la couche courante 

         dzp = zref - zb(la-1)
         pha = -wza(la)*dzp
         q1  = cexp(pha)
         if (vs(la) /= 0.) then
            phb = -wzb(la)*dzp
            q2  = cexp(phb)
         end if
               
!.... Rd = E*RdG*E
         rd(1,1) = rdg(1,1)*q1*q1
         rd(1,2) = rdg(1,2)*q1*q2
         rd(2,1) = rdg(2,1)*q1*q2
         rd(2,2) = rdg(2,2)*q2*q2
         rdsh   =  rdgsh*q2*q2

!.... Tu = E*TuG
         if (rd_only /= 1) then       
            tu(1,1) = tug(1,1)*q1
            tu(1,2) = tug(1,2)*q1
            tu(2,1) = tug(2,1)*q2
            tu(2,2) = tug(2,2)*q2
            tush    = tugsh*q2
         end if

!---- 2A) Construction de la matrice Rd 

!.... RuI*Rd                                                      
         a(1,1) = rui(la,1,1)*rd(1,1) + rui(la,1,2)*rd(2,1)     
         a(1,2) = rui(la,1,1)*rd(1,2) + rui(la,1,2)*rd(2,2)
         a(2,1) = rui(la,2,1)*rd(1,1) + rui(la,2,2)*rd(2,1)
         a(2,2) = rui(la,2,1)*rd(1,2) + rui(la,2,2)*rd(2,2)
         ash    = ruish(la)*rdsh 

!.... Inv{I - RuI*Rd}
         q1  = 1. - a(1,1)                                
         q2  = 1. - a(2,2)
         det = q1*q2 - a(1,2)*a(2,1)
         a(1,1) = q2/det
         a(1,2) = a(1,2)/det
         a(2,1) = a(2,1)/det
         a(2,2) = q1/det
         ash    = 1./(1. - ash)

!.... Inv{I - RuI*Rd}*TdI
         b(1,1) = a(1,1)*tdi(la,1,1) + a(1,2)*tdi(la,2,1)
         b(1,2) = a(1,1)*tdi(la,1,2) + a(1,2)*tdi(la,2,2)
         b(2,1) = a(2,1)*tdi(la,1,1) + a(2,2)*tdi(la,2,1)
         b(2,2) = a(2,1)*tdi(la,1,2) + a(2,2)*tdi(la,2,2)
         bsh    = ash*tdish(la)

!.... Rd*Inv{I - RuI*Rd}*TdI
         a(1,1) = rd(1,1)*b(1,1)  + rd(1,2)*b(2,1)
         a(1,2) = rd(1,1)*b(1,2)  + rd(1,2)*b(2,2)
         a(2,1) = rd(2,1)*b(1,1)  + rd(2,2)*b(2,1)
         a(2,2) = rd(2,1)*b(1,2)  + rd(2,2)*b(2,2) 
         ash    = rdsh*bsh

!.... RdG = RdI + TuI*Rd*Inv{I - RuI*Rd}*TdI
         rdg(1,1) = rdi(la,1,1) + 
     &              tui(la,1,1)*a(1,1) + tui(la,1,2)*a(2,1)
         rdg(1,2) = rdi(la,1,2) + 
     &              tui(la,1,1)*a(1,2) + tui(la,1,2)*a(2,2)
         rdg(2,1) = rdi(la,2,1) + 
     &              tui(la,2,1)*a(1,1) + tui(la,2,2)*a(2,1)
         rdg(2,2) = rdi(la,2,2) + 
     &              tui(la,2,1)*a(1,2) + tui(la,2,2)*a(2,2)
         rdgsh    = rdish(la) + ash

         if (rd_only == 1) cycle

!---- 2B) Construction de la matrice Tu

!.... Rd*RuI 
         a(1,1) = rd(1,1)*rui(la,1,1) + rd(1,2)*rui(la,2,1)     
         a(1,2) = rd(1,1)*rui(la,1,2) + rd(1,2)*rui(la,2,2)
         a(2,1) = rd(2,1)*rui(la,1,1) + rd(2,2)*rui(la,2,1)
         a(2,2) = rd(2,1)*rui(la,1,2) + rd(2,2)*rui(la,2,2)
         ash    = rdsh*ruish(la)

!.... Inv{I - Rd*RuI}
         q1  = 1. - a(1,1)                                
         q2  = 1. - a(2,2)
         det = q1*q2 - a(1,2)*a(2,1)
         a(1,1) = q2/det
         a(1,2) = a(1,2)/det
         a(2,1) = a(2,1)/det
         a(2,2) = q1/det
         ash    = 1./(1. - ash)

!.... Inv{I - Rd*RuI}*Tu
         b(1,1) = a(1,1)*tu(1,1)  + a(1,2)*tu(2,1)
         b(1,2) = a(1,1)*tu(1,2)  + a(1,2)*tu(2,2)
         b(2,1) = a(2,1)*tu(1,1)  + a(2,2)*tu(2,1)
         b(2,2) = a(2,1)*tu(1,2)  + a(2,2)*tu(2,2)
         bsh    = ash*tush

!.... TuG = TuI*Inv{I - Rd*RuI}*Tu
         tug(1,1) = tui(la,1,1)*b(1,1) + tui(la,1,2)*b(2,1)
         tug(1,2) = tui(la,1,1)*b(1,2) + tui(la,1,2)*b(2,2)
         tug(2,1) = tui(la,2,1)*b(1,1) + tui(la,2,2)*b(2,1)
         tug(2,2) = tui(la,2,1)*b(1,2) + tui(la,2,2)*b(2,2)
         tugsh    = tuish(la)*bsh

      end do

      end subroutine prup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine prdn (lstart, zstart, nz, nzm, zz, lcount, 
     &                 nl, nlm, zb, vs, wza, wzb, 
     &                 rdi, rui, tdi, tui, rdish, ruish, tdish, tuish,
     &                 ru_only, rugen, tdgen, rugensh, tdgensh)

!     PROPAGATION VERS LE BAS D'UNE ONDE PLANE DANS UN MILIEU
!     STRATIFIE PLAN
!     Langage F77 avec extensions F90 XLF (AIX IBM RS/6000 Version 3.2.4)

!     Convention pour les matrices de coefficients de reflexion-transmission: 
!        (1,1): P inc. -> P diffr.       (1,2): S inc. -> P diffr.
!        (2,1): P inc. -> S diffr.       (2,2): S inc. -> S diffr.

!     RdI, RuI, TdI, TuI : coefficients de reflexion-transmission 
!                          a l'interface (L-1/L)
!     RuG, TdG           : coefficients de reflexion-transmission generalises
!                          apres passage d'une interface
!     Ru , Td            : coefficients de reflexion-transmission generalises
!                          apres dephasage

!=============================================================================

      implicit none
      integer  i, ip, iz, j, l1, la, lc, lstart, nl, nlm, nz, nzm,        
     &         ru_only, lcount(nlm)
      real     dzp, zref, zstart,
     &         vs(nlm), zb(0:nlm), zz(nzm+1) 
      complex  det, pha, phb, q1, q2, 
     &         ash, bsh, rush, tdsh, rugsh, tdgsh,
     &         a(2,2), b(2,2), ru(2,2), td(2,2), rug(2,2), tdg(2,2),
     &         rdi(nlm,2,2), rui(nlm,2,2), tdi(nlm,2,2), tui(nlm,2,2),
     &         rdish(nlm), ruish(nlm), tdish(nlm), tuish(nlm),        
     &         wza(nlm), wzb(nlm),
     &         rugen(nzm+1,2,2), tdgen(nzm+1,2,2),
     &         rugensh(nzm+1), tdgensh(nzm+1)                

!---- INITIALISATIONS

      do i = 1, 2
         do j = 1, 2
            rug(i,j) = 0.
            tdg(i,j) = 0.
         end do
      end do
      tdg(1,1) = 1.
      tdg(2,2) = 1.
      rugsh = 0.
      tdgsh = 1.

      if (zstart == zb(lstart-1)) then
         rug(1,1) = rui(lstart,1,1)
         rug(1,2) = rui(lstart,1,2)
         rug(2,1) = rui(lstart,2,1)
         rug(2,2) = rui(lstart,2,2)
         rugsh    = ruish(lstart)
      end if

!==== APPLICATION DES RELATIONS ITERATIVES
!     1) Dephasage des coefficients jusqu'a la base de la couche (L)
!     2) Passage de l'interface (L/L+1)

      iz = 1
      do la = lstart, nl-1              ! Propagation du haut vers le bas
         if (la == lstart) then
            zref = zstart
         else
            zref = zb(la-1)
         end if

!---- 1A) Dephasages pour les recepteurs se trouvant dans la couche courante
             
         lc = lcount(la)
         if (lc /= 0) then
            do ip = 1, lc              
               dzp = zz(iz) - zref
               pha = -wza(la)*dzp
               q1  = cexp(pha)
               if (vs(la) /= 0.) then
                  phb = -wzb(la)*dzp
                  q2  = cexp(phb)
               end if
               
!.... Ru = E*RuG*E
               rugen(iz,1,1) = rug(1,1)*q1*q1
               rugen(iz,1,2) = rug(1,2)*q1*q2
               rugen(iz,2,1) = rug(2,1)*q1*q2
               rugen(iz,2,2) = rug(2,2)*q2*q2
               rugensh(iz)   = rugsh*q2*q2

!.... Td = E*TdG
               if (ru_only /= 1) then       
                  tdgen(iz,1,1) = tdg(1,1)*q1
                  tdgen(iz,1,2) = tdg(1,2)*q1
                  tdgen(iz,2,1) = tdg(2,1)*q2
                  tdgen(iz,2,2) = tdg(2,2)*q2
                  tdgensh(iz)   = tdgsh*q2
!!!               if (vs(la) == 0.) then
!!!                  tdgen(iz,2,2) = 0.
!!!                  tdgensh(iz)   = 0.
!!!               end if
               end if

               iz = iz + 1
               if (iz > nz) exit
            end do
         end if

!---- 1B) Dephasage dans la couche courante 

         dzp = zb(la) - zref
         pha = -wza(la)*dzp
         q1  = cexp(pha)
         if (vs(la) /= 0.) then
            phb = -wzb(la)*dzp
            q2  = cexp(phb)
         end if
               
!.... Ru = E*RuG*E
         ru(1,1) = rug(1,1)*q1*q1
         ru(1,2) = rug(1,2)*q1*q2
         ru(2,1) = rug(2,1)*q1*q2
         ru(2,2) = rug(2,2)*q2*q2
         rush    = rugsh*q2*q2

!.... Td = E*TdG 
         if (ru_only /= 1) then      
            td(1,1) = tdg(1,1)*q1
            td(1,2) = tdg(1,2)*q1
            td(2,1) = tdg(2,1)*q2
            td(2,2) = tdg(2,2)*q2
            tdsh    = tdgsh*q2
!!!         if (vs(la) == 0.) then
!!!            td(2,2) = 0.
!!!            tdsh    = 0.
!!!         end if
         end if

         l1 = la + 1

!---- 2A) Construction de la matrice Ru

!.... RdI*Ru 
         a(1,1) = rdi(l1,1,1)*ru(1,1) + rdi(l1,1,2)*ru(2,1)     
         a(1,2) = rdi(l1,1,1)*ru(1,2) + rdi(l1,1,2)*ru(2,2)
         a(2,1) = rdi(l1,2,1)*ru(1,1) + rdi(l1,2,2)*ru(2,1)
         a(2,2) = rdi(l1,2,1)*ru(1,2) + rdi(l1,2,2)*ru(2,2)
         ash    = rdish(l1)*rush

!.... Inv{I - RdI*Ru}
         q1  = 1. - a(1,1)                                
         q2  = 1. - a(2,2)
         det = q1*q2 - a(1,2)*a(2,1)
         a(1,1) = q2/det
         a(1,2) = a(1,2)/det
         a(2,1) = a(2,1)/det
         a(2,2) = q1/det
         ash    = 1./(1. - ash)

!.... Inv{I - RdI*Ru}*TuI
         b(1,1) = a(1,1)*tui(l1,1,1) + a(1,2)*tui(l1,2,1)
         b(1,2) = a(1,1)*tui(l1,1,2) + a(1,2)*tui(l1,2,2)
         b(2,1) = a(2,1)*tui(l1,1,1) + a(2,2)*tui(l1,2,1)
         b(2,2) = a(2,1)*tui(l1,1,2) + a(2,2)*tui(l1,2,2)
         bsh    = ash*tuish(l1)

!.... Ru*Inv{I - RdI*Ru}*TuI
         a(1,1) = ru(1,1)*b(1,1) + ru(1,2)*b(2,1)
         a(1,2) = ru(1,1)*b(1,2) + ru(1,2)*b(2,2)
         a(2,1) = ru(2,1)*b(1,1) + ru(2,2)*b(2,1)
         a(2,2) = ru(2,1)*b(1,2) + ru(2,2)*b(2,2)
         ash    = rush*bsh 

!.... RuG = RuI + TdI*Ru*Inv{I - RdI*Ru}*TuI
         rug(1,1) = rui(l1,1,1) + 
     &              tdi(l1,1,1)*a(1,1) + tdi(l1,1,2)*a(2,1)       
         rug(1,2) = rui(l1,1,2) + 
     &              tdi(l1,1,1)*a(1,2) + tdi(l1,1,2)*a(2,2)
         rug(2,1) = rui(l1,2,1) + 
     &              tdi(l1,2,1)*a(1,1) + tdi(l1,2,2)*a(2,1)
         rug(2,2) = rui(l1,2,2) + 
     &              tdi(l1,2,1)*a(1,2) + tdi(l1,2,2)*a(2,2)
         rugsh    = ruish(l1) + tdish(l1)*ash
      
         if (ru_only == 1) cycle

!---- 2B) Construction de la matrice Td

!.... Ru*RdI                                                      
         a(1,1) = ru(1,1)*rdi(l1,1,1) + ru(1,2)*rdi(l1,2,1)     
         a(1,2) = ru(1,1)*rdi(l1,1,2) + ru(1,2)*rdi(l1,2,2)
         a(2,1) = ru(2,1)*rdi(l1,1,1) + ru(2,2)*rdi(l1,2,1)
         a(2,2) = ru(2,1)*rdi(l1,1,2) + ru(2,2)*rdi(l1,2,2)
         ash    = rush*rdish(l1)

!.... Inv{I - Ru*RdI}
         q1  = 1. - a(1,1)                                
         q2  = 1. - a(2,2)
         det = q1*q2 - a(1,2)*a(2,1)
         a(1,1) = q2/det
         a(1,2) = a(1,2)/det
         a(2,1) = a(2,1)/det
         a(2,2) = q1/det
         ash    = 1./(1. - ash)

!.... Inv{I - Ru*RdI}*Td
         b(1,1) = a(1,1)*td(1,1) + a(1,2)*td(2,1)
         b(1,2) = a(1,1)*td(1,2) + a(1,2)*td(2,2)
         b(2,1) = a(2,1)*td(1,1) + a(2,2)*td(2,1)
         b(2,2) = a(2,1)*td(1,2) + a(2,2)*td(2,2)
         bsh    = ash*tdsh

!.... TdG = TdI*Inv{I - Ru*RdI}*Td
         tdg(1,1) = tdi(l1,1,1)*b(1,1) + tdi(l1,1,2)*b(2,1)
         tdg(1,2) = tdi(l1,1,1)*b(1,2) + tdi(l1,1,2)*b(2,2)
         tdg(2,1) = tdi(l1,2,1)*b(1,1) + tdi(l1,2,2)*b(2,1)
         tdg(2,2) = tdi(l1,2,1)*b(1,2) + tdi(l1,2,2)*b(2,2)
         tdgsh    = tdish(l1)*bsh

      end do

      end subroutine prdn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine asmup (nabove, nz, nzm, nrliq, 
     &                  rdrs, rufr, turs, rdsl, rufs, 
     &                  rdrssh, rufrsh, turssh, rdslsh, rufssh, 
     &                  reflec, reflecsh)              

!     CALCUL DE LA REFLECTIVITE GLOBALE D'UN MILIEU STRATIFIE PLAN
!     (Assemblage des matrices de reflexion-transmission generalisees)
!     Langage F77 avec extensions F90 XLF (AIX IBM RS/6000 Version 3.2.4)

!     Convention pour les matrices de coefficients de reflexion-transmission: 
!        (1,1): P inc. -> P diffr.       (1,2): S inc. -> P diffr.
!        (2,1): P inc. -> S diffr.       (2,2): S inc. -> S diffr.

!=============================================================================

      implicit none
      integer  i, ifs, iz, j, nabove, nrliq, nz, nzm
      complex  det, q1, q2, ash, bsh, rufssh, rdslsh,
     &         a(2,2), b(2,2), rufs(2,2), rdsl(2,2),
     &         rdrs(nzm,2,2), rufr(nzm,2,2), turs(nzm,2,2),
     &         rdrssh(nzm), rufrsh(nzm), turssh(nzm), 
     &         reflec(nzm,2,2), reflecsh(nzm)
 
      do iz = 1, nabove

      if (iz > nrliq) then

!.... RdRS*RuFR      
         a(1,1) = rdrs(iz,1,1)*rufr(iz,1,1) + rdrs(iz,1,2)*rufr(iz,2,1)        
         a(1,2) = rdrs(iz,1,1)*rufr(iz,1,2) + rdrs(iz,1,2)*rufr(iz,2,2)
         a(2,1) = rdrs(iz,2,1)*rufr(iz,1,1) + rdrs(iz,2,2)*rufr(iz,2,1)
         a(2,2) = rdrs(iz,2,1)*rufr(iz,1,2) + rdrs(iz,2,2)*rufr(iz,2,2)
         ash    = rdrssh(iz)*rufrsh(iz)

!.... Inv[I - RdRS*RuFR]
         q1  = 1. - a(1,1)
         q2  = 1. - a(2,2)
         det = q1*q2 - a(1,2)*a(2,1)
         a(1,1) = q2/det
         a(1,2) = a(1,2)/det
         a(2,1) = a(2,1)/det
         a(2,2) = q1/det
         ash    = 1./(1. - ash)

!.... Inv[I - RdRS*Rt]*TuRS
         b(1,1) = a(1,1)*turs(iz,1,1) + a(1,2)*turs(iz,2,1)
         b(1,2) = a(1,1)*turs(iz,1,2) + a(1,2)*turs(iz,2,2)
         b(2,1) = a(2,1)*turs(iz,1,1) + a(2,2)*turs(iz,2,1)
         b(2,2) = a(2,1)*turs(iz,1,2) + a(2,2)*turs(iz,2,2)
         bsh    = ash*turssh(iz)  

!.... RdSL*RuFS
         a(1,1) = rdsl(1,1)*rufs(1,1) + rdsl(1,2)*rufs(2,1)
         a(1,2) = rdsl(1,1)*rufs(1,2) + rdsl(1,2)*rufs(2,2)
         a(2,1) = rdsl(2,1)*rufs(1,1) + rdsl(2,2)*rufs(2,1)
         a(2,2) = rdsl(2,1)*rufs(1,2) + rdsl(2,2)*rufs(2,2)
         ash    = rdslsh*rufssh

!.... Inv[I - RdSL*RuFS]
         q1  = 1. - a(1,1)
         q2  = 1. - a(2,2)
         det = q1*q2 - a(1,2)*a(2,1)
         a(1,1) = q2/det
         a(1,2) = a(1,2)/det
         a(2,1) = a(2,1)/det
         a(2,2) = q1/det
         ash    = 1./(1. - ash)

!.... Reflec = Inv[I - RdRS*RuFR]*TuRS*Inv[I - RdSL*RuFS]
         reflec(iz,1,1) = b(1,1)*a(1,1) + b(1,2)*a(2,1)
         reflec(iz,1,2) = b(1,1)*a(1,2) + b(1,2)*a(2,2)
         reflec(iz,2,1) = b(2,1)*a(1,1) + b(2,2)*a(2,1)
         reflec(iz,2,2) = b(2,1)*a(1,2) + b(2,2)*a(2,2)
         reflecsh(iz)   = bsh*ash

      else
         q1 = 1./(1. - rdrs(iz,1,1)*rufr(iz,1,1))
         q2 = 1./(1. - rdsl(1,1)*rufs(1,1))
         reflec(iz,1,1) = turs(iz,1,1)*q1*q2
         reflec(iz,1,2) = 0.
         reflec(iz,2,1) = 0.
         reflec(iz,2,2) = 0.
         reflecsh(iz)   = 0.
      end if

      end do

      end subroutine asmup
   
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine asmdn (nabove, nz, nzm, nrliq, 
     &                  rurs, rdrl, tdrs, rufs, rdsl, 
     &                  rurssh, rdrlsh, tdrssh, rufssh, rdslsh, 
     &                  reflec, reflecsh)              

!     CALCUL DE LA REFLECTIVITE GLOBALE D'UN MILIEU STRATIFIE PLAN
!     (Assemblage des matrices de reflexion-transmission generalisees)
!     Langage F77 avec extensions F90 XLF (AIX IBM RS/6000 Version 3.2.4)

!     Convention pour les matrices de coefficients de reflexion-transmission: 
!        (1,1): P inc. -> P diffr.       (1,2): S inc. -> P diffr.
!        (2,1): P inc. -> S diffr.       (2,2): S inc. -> S diffr.

!=============================================================================

      implicit none
      integer  i, ifs, iz, j, nabove, nrliq, nz, nzm
      complex  det, q1, q2, ash, bsh, rufssh, rdslsh,
     &         a(2,2), b(2,2), rufs(2,2), rdsl(2,2),
     &         rurs(nzm,2,2), rdrl(nzm,2,2), tdrs(nzm,2,2), 
     &         rurssh(nzm), rdrlsh(nzm), tdrssh(nzm), 
     &         reflec(nzm,2,2), reflecsh(nzm)

      do iz = nabove+1, nz

      if (iz > nrliq) then

!.... RuRS*RdRL      
         a(1,1) = rurs(iz,1,1)*rdrl(iz,1,1) + rurs(iz,1,2)*rdrl(iz,2,1)        
         a(1,2) = rurs(iz,1,1)*rdrl(iz,1,2) + rurs(iz,1,2)*rdrl(iz,2,2)
         a(2,1) = rurs(iz,2,1)*rdrl(iz,1,1) + rurs(iz,2,2)*rdrl(iz,2,1)
         a(2,2) = rurs(iz,2,1)*rdrl(iz,1,2) + rurs(iz,2,2)*rdrl(iz,2,2)
         ash    = rurssh(iz)*rdrlsh(iz)

!.... Inv[I - RuRS*RdRL]
         q1  = 1. - a(1,1)
         q2  = 1. - a(2,2)
         det = q1*q2 - a(1,2)*a(2,1)
         a(1,1) = q2/det
         a(1,2) = a(1,2)/det
         a(2,1) = a(2,1)/det
         a(2,2) = q1/det
         ash    = 1./(1. - ash)

!.... Inv[I - RuRS*RdRL]*TdRS
         b(1,1) = a(1,1)*tdrs(iz,1,1) + a(1,2)*tdrs(iz,2,1)
         b(1,2) = a(1,1)*tdrs(iz,1,2) + a(1,2)*tdrs(iz,2,2)
         b(2,1) = a(2,1)*tdrs(iz,1,1) + a(2,2)*tdrs(iz,2,1)
         b(2,2) = a(2,1)*tdrs(iz,1,2) + a(2,2)*tdrs(iz,2,2)
         bsh    = ash*tdrssh(iz)

!.... RuFS*RdSL
         a(1,1) = rufs(1,1)*rdsl(1,1) + rufs(1,2)*rdsl(2,1)
         a(1,2) = rufs(1,1)*rdsl(1,2) + rufs(1,2)*rdsl(2,2)
         a(2,1) = rufs(2,1)*rdsl(1,1) + rufs(2,2)*rdsl(2,1)
         a(2,2) = rufs(2,1)*rdsl(1,2) + rufs(2,2)*rdsl(2,2)
         ash    = rufssh*rdslsh

!.... Inv[I - RuFS*RdSL]
         q1  = 1. - a(1,1)
         q2  = 1. - a(2,2)
         det = q1*q2 - a(1,2)*a(2,1)
         a(1,1) = q2/det
         a(1,2) = a(1,2)/det
         a(2,1) = a(2,1)/det
         a(2,2) = q1/det
         ash    = 1./(1. - ash)

!.... Reflec = Inv[I - RuRS*RdRL]*TdRS*Inv[I - RuFS*RdSL]
         reflec(iz,1,1) = b(1,1)*a(1,1) + b(1,2)*a(2,1)
         reflec(iz,1,2) = b(1,1)*a(1,2) + b(1,2)*a(2,2)
         reflec(iz,2,1) = b(2,1)*a(1,1) + b(2,2)*a(2,1)
         reflec(iz,2,2) = b(2,1)*a(1,2) + b(2,2)*a(2,2)
         reflecsh(iz)   = bsh*ash

      else
         q1 = 1./(1. - rurs(iz,1,1)*rdrl(iz,1,1))
         q2 = 1./(1. - rufs(1,1)*rdsl(1,1))
         reflec(iz,1,1) = tdrs(iz,1,1)*q1*q2
         reflec(iz,1,2) = 0.
         reflec(iz,2,1) = 0.
         reflec(iz,2,2) = 0.
         reflecsh(iz)   = 0.
      end if

      end do

      end subroutine asmdn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine sbj0j19 (arg, n, bj0, bj1)

! CALCULE LA VALEUR DES FONCTIONS DE BESSEL J0 ET J1 POUR UN ENSEMBLE DE VALEURS
! REELLES arg(n) EN UTILISANT LES APPROXIMATIONS POLYNOMIALES D'ABRAMOVITZ
! ET STEGUN (1964, pp. 369-370).

! Entree:        arg(i), i = 1, n : n valeurs reelles rangees par
!                                   ordre croissant

! Sortie:        bj0(i), i = 1, n : Valeurs J0 et J1 correspondant aux 
!                bj1(i), i = 1, n   arguments arg(i), i = 1, n

! Utilisation:   program main
!                integer n
!                real    arg(n), bj0(n), bj1(n)
!                call sbj0j19 (arg, n, bj0, bj1)

! Langage F77 avec extensions F90 xlf IBM RS/6000 (AIX Version 3.2.3)
! M. Dietrich LGIT 04/93

!==============================================================================

      implicit none
      integer  i, iac, iac1, n
      real     a1, a2, f0, f1, t, z0, z1
      real     arg(n), bj0(n), bj1(n)

! Test arg > 3 ?

      iac = n
      do i = 1, n
         a1 = arg(i)
         if (a1 <= 3.) cycle
         iac = i - 1
         exit
      end do

! Petits arguments 

      do i = 1, iac
         a1 = arg(i)
         t  = (a1/3.)**2
         bj0(i) = 1. - t*(2.2499997 - t*(1.2656208 - t*(.3163866
     &               - t*(.0444479 - t*(.0039444 - t*.00021)))))
         bj1(i) = a1*(.5 - t*(.56249985 - t*(.21093573 - t*(.03954289
     &               - t*(.00443319 - t*(.00031761 - t*.00001109))))))        
      end do

! Grands arguments

      iac1 = iac + 1
      do i = iac1, n
         a1 = arg(i)
         a2 = abs(a1)
         t  = 3./a2
         f0 = .79788456 + t*(-.00000077 + t*(-.0055274 + t*(-.00009512
     &        + t*(.00137237 + t*(-.00072805 + t*.00014476)))))
         z0 = .78539816 - t*(-.04166397 + t*(-.00003954 + t*(.00262573  
     &        + t*(-.00054125 + t*(-.00029333 + t*.00013558)))))
         bj0(i) = f0*cos(a2 - z0)/sqrt(a2)

         f1 = .79788456 + t*(.00000156 + t*(.01659667 + t*(.00017105
     &        + t*(-.00249511 + t*(.00113653 - t*.00020033)))))
         z1 = -2.35619449 + t*(.12499612 + t*(.0000565 + t*(-.00637879
     &        + t*(.00074348 + t*(.00079824 - t*.00029166)))))
         bj1(i) = sign(f1,a1)*cos(a2 + z1)/sqrt(a2)
      end do

      end subroutine sbj0j19

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine crfft9 (sigf, nt, sigt, isign) 

! FFT COMPLEXE A REEL - UTILISE LA SUBROUTINE FOUR1 DES NUMERICAL RECIPES

! Entree:        sigf(i), i = 1, nt/2 : signal complexe de nt/2 frequences
!                                       positives

! Sortie:        sigt(i), i = 1, nt   : signal reel de nt echantillons 

! Utilisation:   program main
!                integer nt
!                real    sigt(nt)
!                complex sigf(nt/2)
!                call crfft9 (sigf, nt, sigt, isign)

! Langage F77 avec extensions F90 xlf IBM RS/6000 (AIX Version 3.2.3)
! M. Dietrich LGIT 04/93  

!==============================================================================

!     implicit  none
      parameter (ipm = 13, ntm = 2**ipm)        ! ici nt = 8192 pts max
      integer   i, i1, i2, ifr, isign, ip, ip2,
     &          it, itest, n2, n2m, n3, nf, nt 
      real      sigt(nt), work(2*ntm)
      complex   sigf(nt/2)

      nf  = nt/2            ! Nombre de frequences positives
      n2  = 2*nt            ! Dimension utile du tableau (reel) WORK
      n2m = 2*ntm           ! Dimension maximale du tableau WORK 
      n3  = n2 + 3

! Teste si nt n'excede pas la valeur maximale autorisee 
!    dans cette mise en oeuvre. La longueur maximum du signal
!    reel en sortie est arbitrairement fixee a nt = ntm = 2**ipm

      if (nt > ntm) then
         print *, 'nt, ntm :', nt, ntm
         stop     "CRFFT9 *** TABLEAU work SOUS DIMENSIONNE ***"
      end if

! Teste si nt est une puissance de 2

      itest = 0
      do ip = 1, ipm
         ip2 = 2**ip
         if (mod(nt,ip2) == 0) itest = 1
      end do
      if (itest == 0) then
         print *, 'nt =', nt
         stop     "CRFFT9 *** PARAMETRE nt MAL DEFINI (/= 2**p) ***"        
      end if

! Initialisation et chargement du tableau WORK 

      do i = 1, n2m
         work(i) = 0.
      end do

      do ifr = 1, nf
         i2 = 2*ifr    
         i1 = i2 - 1
         work(i1) =  real(sigf(ifr))
         work(i2) = aimag(sigf(ifr))
         if (ifr > 1) then
            work(n3 - i2) =  work(i1)
            work(n3 - i1) = -work(i2)
         end if
      end do

! Transformee de Fourier utilisant la subroutine FOUR1 des Numerical Recipes
!    (Algorithme original de Cooley-Tukey)

      call nrfour1 (work, nt, isign)
 
! Recuperation du signal reel SIGT

      do it = 1, nt  
         sigt(it) = work(2*it - 1)/real(nt)
      end do

      end subroutine crfft9

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      SUBROUTINE NRFOUR1(DATA,NN,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

