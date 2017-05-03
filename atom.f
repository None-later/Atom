C   April 26 2017 Wed. Can calculate n-Gaussian single-electron atom
C   April 26 2017 Wed. Can calculate single-Gaussian single-electron atom
C   April 23 2017 Sunday adding Coulomb (e-nucl) and Kinetic Enegry integrals
C   Now we read the name of a basis file from the command argumnet

      module commondat
         implicit none
         real, parameter :: Pi = 3.1415927
         integer el, cntr
         integer n(7), Z(7)
         character(2) symbol, orb
         character(20) string
C    *********************************************************************************
C    alf(8,7,40): maximum 8 elements; 7 types of orbitals: s,p,d,f,g,h,i; 60 gaussians
C    n(o): curent contraction position for orbital o; contr(el,o,i) = contraction
C    *********************************************************************************
         real(8) alf(8,7,60), coef(8,7,60), coef1(8,7,60), KE1, PE1
         integer contr(8,7,62)
         interface
      end interface
      end module commondat
 
      program atom
      use commondat
      implicit none

      interface
         subroutine get(o, p)
            integer o, i, shift
            integer, optional :: p
         end subroutine get
         real(8) function norm(e,o,i)
            integer e, o, i, k, l
            real(8) add
         end function norm
         real(8) function KE(e,o,i)
            integer e, o, i, k, l
            real(8) add
         end function KE
         real(8) function PE(e,o,i)
            integer e, o, i, k
         end function PE
      end interface

      integer o, i, e, k
      real(8) s
      el = 0
      call getarg(1, string)
      open(unit=7,status="old",file=string)
      do 
         READ(7,'(A2)',end=100) symbol
         if (symbol .eq. '**') exit
      enddo

 20   READ(7,'(A2)') symbol
         el = el + 1
         select case (symbol)
         case ('H');  Z(el) = 1
         case ('He'); Z(el) = 2
         case ('Li'); Z(el) = 3
         case ('Be'); Z(el) = 4
         case ('B');  Z(el) = 5
         case ('C');  Z(el) = 6
         case ('N');  Z(el) = 7
         case ('O');  Z(el) = 8
         case default; goto 100
         end select
      do o = 1, 7
         n(o) = 1
         contr(el,o,1) = 1
      enddo

 10   READ(7,'(A6)') string
      if (string(1:1) .ne. "*" ) then 
         READ(string,'(A2,I5)') orb, cntr
         select case (orb)
         case ('S');  call get(1) 
         case ('SP'); call get(1,2)
         case ('P');  call get(2) 
         case ('D');  call get(3)
         case ('F');  call get(4)
         case ('G');  call get(5) 
         case ('H');  call get(6) 
         case ('I');  call get(7) 
         end select
         goto 10
      else
         goto 20
      endif
 
100   close(7)

      do e = 1, el-1
         print '(A8,I2)', 'Element', e 
         do o = 1, 7
            if (contr(e,o,2) .ne. 0) then 
               print '(I2)', o
               do i = 1, 60
                  do k = contr(e,o,i), contr(e,o,i+1)-1
                     print '(2F20.9)', alf(e,o,k), coef(e,o,k)
                     coef(e,o,k)=coef(e,o,k)*((2*alf(e,o,k)/Pi)**0.75) 
                  enddo
                  KE1 = KE(e,o,i)/norm(e,o,i) 
                  PE1 = PE(e,o,i)/norm(e,o,i) 
                  print 300, 'NORM:', norm(e,o,i), 'KE:', KE1, 
     +'PE:', PE1, 'SCF:', KE1-PE1
                  if (contr(e,o,i+2) .eq. 0) goto 200
               enddo
            endif
200      enddo
      enddo

300   format (4(A8,F11.7))
      
      end program atom


      subroutine get(o, p)
      use commondat
      implicit none  
      integer o, i, shift
      integer, optional :: p

      contr(el,o,n(o)+1) = contr(el,o,n(o)) + cntr
      contr(el,o,n(o)+2) = 0
      if (present(p)) then 
         contr(el,p,n(p)+1) = contr(el,p,n(p)) + cntr
         contr(el,p,n(p)+2) = 0
         shift = contr(el,p,n(p)) - contr(el,o,n(o))
         do  i = contr(el,o,n(o)), contr(el,o,n(o)+1)-1
            read(7,'(3F23.9)'),  alf(el,o,i), coef(el,o,i), 
     +                           coef(el,p,i+shift)
            alf(el,p,i+shift) = alf(el,o,i)
         enddo
         n(p) = n(p) + 1
      else 
         do  i = contr(el,o,n(o)), contr(el,o,n(o)+1)-1
            read(7,'(2F23.9)'),  alf(el,o,i), coef(el,o,i)
         enddo
      endif
      n(o) = n(o) + 1
      end subroutine get
       
      real(8) function norm(e,o,i)
         use commondat
         implicit none
         integer e, o, i, k, l
         real(8) :: add = 0.0
         norm = 0.0
         do  k = contr(e,o,i), contr(e,o,i+1)-1
            do  l = contr(e,o,i), contr(e,o,i+1)-1
               add = (Pi/(alf(e,o,k)+alf(e,o,l)))**1.5
               add = add*coef(e,o,k)*coef(e,o,l)
               norm = norm + add
            enddo
         enddo
         return
      end function norm

      real(8) function PE(e,o,i)
         use commondat
         implicit none
         integer e, o, i, k, l
         real(8) :: add = 0.0
         PE = 0.0
         do  k = contr(e,o,i), contr(e,o,i+1)-1
            do  l = contr(e,o,i), contr(e,o,i+1)-1
               add = 2*Pi/(alf(e,o,k)+alf(e,o,l))
               add = add*coef(e,o,k)*coef(e,o,l)
               PE = PE + add
            enddo
         enddo
         PE = PE*Z(e)
         return
      end function PE

      real(8) function KE(e,o,i)
         use commondat
         implicit none
         integer e, o, i, k, l
         real(8) :: add = 0.0
         KE = 0.0
         do  k = contr(e,o,i), contr(e,o,i+1)-1
            do  l = contr(e,o,i), contr(e,o,i+1)-1
               add=3*alf(e,o,k)*alf(e,o,l)/(alf(e,o,k)+alf(e,o,l))**2.5
               add = add*coef(e,o,k)*coef(e,o,l)*Pi**1.5
               KE = KE + add
            enddo
         enddo
         return
      end function KE

