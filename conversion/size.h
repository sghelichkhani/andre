!...  Number of grid intervals along icosahedral diamond edge:
      integer, parameter :: mt=512
 
!...  Number of grid intervals along edge of local subdomain:
      integer, parameter :: nt=32
 
!...  Number of diamonds mapped to a local process (5 or 10):
      integer, parameter :: nd=5

!...  Note: Number of processes is given by (mt/nt)**2*10/nd.
 
!...  Number of radial layers:
      integer, parameter :: nr= mt/2
 
!...  Number of diamonds dimension for opr array (1, 5 or 10):
      integer, parameter :: ndo= nd
!...  Note:  If there is no lateral viscosity variation, ndo may
!...         be set to 1; otherwise, it must be set equal to nd.
 
!...  Number of grid points (vertices) in local subdomain:
      integer, parameter :: nv= (nt+1)**2*nd*(nr+1)

 
!...  Scaling of plate velocities according to RMS surface velocities
!...  obtained from a free slip calculation (2.8 for viscosity increase
!...  of 32, 3.7 for 100 and 4.3 for 200 with visc==1e21 and 3500K t_icb)
!...  with t_icb==4200 (35% Q_CMB): 3.7*1.215~=4.5 for 100 and visc=1e21
!...  with t_icb==2900 ( 5% Q_CMB): 3.7*0.785~=2.9 for 100 and visc=1e21
!...  asthenosphere down to 600km: 5.27 for 32, 4.37 for 100 and 6.51 for 200 with visc==5e20 and 3500K t_icb
!...  asthenosphere down to 400km: 6.43 for 32, 5.97 for 100 and 4.39 for 200 with visc==5e20 and 3500K t_icb)
!      real, parameter :: plate_scale=3.7
