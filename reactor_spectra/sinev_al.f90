! Program by Valery Sinev.
! Compute reactor spectra for average composition of isotopes during reactor cycle.
! Isotope composition corresponds to Rovno reactor. 
! V. I. Kopeikin, L. A. Mikaelyan, and V. V.  Sinev 
! Yad. Fiz. 60, 230 (1997) [Phys. At.  Nucl. 60, 172 (1997)].
    subroutine sourc0
       implicit real*8 (a-h,o-z)
!    Calculation of antineutrino spectrum
    common/neutrino/ en(1200),yn(1200),ep(1200),yp(1200)
    common/step/ step,nspe,nspec
       a1=5.09
       a2=0.648
       a3=0.0273
       a4=1.411
!
     step=0.01
     nspe=1200
!
    ymax=9.25
!
    do i=1,nspe
     en(i)=1.8+step*(i-1)
     ec=en(i)
    if(ec.gt.ymax) then
       nspe=i
       goto 5
    endif
!
       yn(i)=a1*exp(-a2*ec-a3*ec**2-a4*(0.125*ec)**10)
       print '(f5.2,a,e12.5)', ec, ' ', yn(i)
       ! print '(e7.2,a,e12.5)', ec, ',', yn(i)
       ! How to print left justified numbers in Fortran:
       ! https://jblevins.org/log/leftjust
       ! print '(a7,f7.2,a2,e12.5)', 'E_nu=', ec, ' ', yn(i)
!
    enddo
 5    continue
!   print '(a7,f7.2,a2,i5)',' E_max=',ymax,'  ',nspe
    return
    end 

PROGRAM MAIN
print "(a)", "# E_nu vs number of antineutrino per MeV per fission"
CALL sourc0
END

