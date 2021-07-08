cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            THIS FILE CONTAINS IN FORTRAN FORM ALL THE PARAMETERS
c                 THAT ARE GOING TO BE USED FOR THE CODE TO RUN
c
c            EXCEPT FROM PROBLEM PARAMETERS IT ALSO HAS TIME AND OTHER STUFF
c


      DOUBLE PRECISION CLup,kon,koff,kdeg
      DOUBLE PRECISION Vp,Ve,DoseA,C0
      DOUBLE PRECISION tend
C==============================================================C
C================      Kinetic Parameters     =================C
C==============================================================C
      parameter(CLup=0.167d0)               ! L/h   
      parameter(kon=0.559d0)                ! nM/h
      parameter(koff=23.9d0)          ! 1/h
      parameter(kdeg=25.0d0)                ! 1/h
C==============================================================C
C================     Compartment volume      =================C
C==============================================================C
      parameter(Vp=3.1d0)                   ! L
      parameter(Ve=0.34d0)                  ! L
C==============================================================C
C================   Baseline concentrations   =================C
C==============================================================C
      parameter(DoseA=24.0d+4)         ! nmol
      parameter(C0=4.98d+4)            ! nM   
C==============================================================C
C================    Simulation paremeters    =================C
C==============================================================C
      parameter(tend=5.1d+3)        ! day

