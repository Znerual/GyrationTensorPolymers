    !  GyrationTensorPolymers.f90
    !
    !  FUNCTIONS:
    !  GyrationTensorPolymers - Entry point of console application.
    ! Authors: Marina Bishara and Laurenz Ruzicka
    !
    !****************************************************************************
    !
    !  PROGRAM: GyrationTensorPolymers
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !  COMPILE1: gfortran -c LinearSystems.f90
    !  COMPILE: gfortran LinearSystems.o GyrationTensorPolymers.f90 -o gyration -L//lib -llapack -lblas -cpp
    !  RUN: ./gyration
    !  PLOT: gnuplot plot.gp
    !****************************************************************************

    program GyrationTensorPolymers
    use LinearSystems
#if defined(__INTEL_COMPILER__)
    use LAPACK95
#endif
    implicit none


    INTEGER :: i,j, k
    INTEGER, PARAMETER :: m=100, n=100 ! 2501
    
    DOUBLE PRECISION :: rg_av, rg_pos_av, asphericity_av, prolateness_av
    DOUBLE PRECISION :: lambda1,lambda2,lambda3
    CHARACTER(2) :: ncstring
    !Zum automatischen Fileinlesen
    character(20) :: filename
    character(2) :: stiffness
    character(5) :: form
    integer :: stat, length, num, conf=100
    !berechnen einer konfiguration
    double precision, dimension(:,:), allocatable :: position
    double precision, dimension(3) :: schwerpunkt = 0d0
    double precision, dimension(3, 3) :: gyrationTensor, gyrationTensorLAPACK
    double precision :: sxx = 0d0, syy= 0d0, sxy= 0d0, sxz= 0d0, syz= 0d0, szz= 0d0
    !berechnen der housholdertransformation
    double precision, dimension(3) :: u1
    double precision, dimension(2) :: u2
    
    double precision, dimension(3, 3) :: H1, H2, Q, R
    double precision, dimension(9) :: valI = [1d0,0d0, 0d0,0d0, 1d0, 0d0, 0d0, 0d0, 1d0]
    double precision, dimension(3,3) :: einheit
    
    double precision, dimension(3) :: eigenwert
    !LAPACK Eigenwerte
    double precision, dimension(12) :: lwork
    double precision, dimension(3) :: wr, wi
    double precision, dimension(3,3) :: vl, vr
    integer :: info
    
    !Sortierte EG und berechnung der physikalischen Größen
    double precision, dimension(3) :: eigenwert_sort
    DOUBLE PRECISION :: rg,rg_pos,asphericity,prolateness
    
    
    einheit =  reshape(valI, [3,3])
    ! Speichern der Namen der .dat Dateien in FILES.txt
#if defined(__GFORTRAN__)
    call system('ls -l | grep *.dat > FILES.txt') !Linux
#else
    call system('DIR *.DAT /B > FILES.txt') !Windows
#endif

    open(unit= 11,file="FILES.txt",action="read")
    open(unit =17, file="Average-chain.txt", action="write")
    open(unit =18, file="Average-ring.txt", action="write")
    open(unit =19, file="Average-star.txt", action="write")
    open(unit=15, file="Eigenvalues.txt", action="write")
    do
        read(11,*, iostat=stat) filename
        filename = trim(filename)
        length = LEN(trim(filename))
        if (stat /= 0) exit

        ! Form des Polymers
        form = filename(6:length-7)
        ! Anzahl an Monomeren festlegen
        if (form == "chain" .or. form == "ring") then
            num = 100
        else if (form == "star") then
            num = 2501
        endif

        !bestimmen der größen
        allocate(position(num, 3))

        ! Stiffness des Polymers
        stiffness = filename(length-5:length-4)

        ! Test
        print*, filename, LEN(trim(filename)), stiffness, ' ' , form, ' ', num
        
        ! Aufrufen einer der .dat Dateien
        open(unit=12, file=filename, action="read")
        
        open(unit=13, file=trim("data" // filename(1:length-3) // "txt"),action="write")
        ! Loop über die verschiedenen Konfigurationen, die in einer Datei gespeichert sind
        do i = 1,conf
            schwerpunkt = 0d0
            do j = 1, num
                read(12, *) position(j,1), position(j,2), position(j,3)
                schwerpunkt = schwerpunkt + position(j,1:3)
            end do
            
            
            schwerpunkt = schwerpunkt / dble(num)
            
            rg_pos = 0d0
            sxx = 0d0
            syy = 0d0
            szz = 0d0
            sxy = 0d0
            sxz = 0d0
            syz = 0d0
            do j = 1, num
                sxx = sxx + (position(j,1) - schwerpunkt(1))**2
                syy = syy + (position(j,2) - schwerpunkt(2))**2
                szz = szz + (position(j,3) - schwerpunkt(3))**2
                sxy = sxy + (position(j,1) - schwerpunkt(1))*(position(j,2) - schwerpunkt(2))
                sxz = sxz + (position(j,1) - schwerpunkt(1))*(position(j,3) - schwerpunkt(3))
                syz = syz + (position(j,2) - schwerpunkt(2))*(position(j,3) - schwerpunkt(3))
                rg_pos = rg_pos +  sum((position(j, 1:3) - schwerpunkt(1:3))**2)

            end do
            
            
            
            rg_pos = sqrt(rg_pos / dble(num))
            
            gyrationTensor(1,1:3) = [sxx, sxy, sxz]
            gyrationTensor(2,1:3) = [sxy, syy, syz]
            gyrationTensor(3,1:3) = [sxz, syz, szz]
            gyrationTensor = gyrationTensor / dble(num)
            
            
            
            gyrationTensorLAPACK = gyrationTensor
            !call showMatrix(gyrationTensor, "Gyration Tensor")
            call dgeev('N','N',3, gyrationTensorLAPACK, 3,wr,wi,vl,3,vr,3, lwork,12,info)
            
            
            k = 0
            do while(.true.)
                if (k >= 100) exit
                k = k + 1
                call QRFaktorisierung3d(gyrationTensor, Q,R)
                gyrationTensor = matmul(R,Q)
                if (abs(gyrationTensor(2,1)) < 1d-10 .and. &
                    abs(gyrationTensor(3,1)) < 1d-10 .and. abs(gyrationTensor(3,2)) < 1d-10) then
                    exit
                end if
            end do
            do k= 1, 3
                eigenwert(k) = gyrationTensor(k,k)
            end do
           write(15, '(A,1x,A,1x,I3, 1x,F16.6, 1x,F16.6, 1x,F16.6, 1x,F16.6, 1x,F16.6, 1x,F16.6, 1x)') &
               form, stiffness, i, eigenwert(1), eigenwert(2), eigenwert(3), wr(1), wr(2), wr(3)
            !Sortieren der Eigenwerte
            eigenwert_sort(1) = MINVAL(wr)
            wr(MINLOC(wr)) = MAXVAL(wr) + 1d0
            eigenwert_sort(2) = MINVAL(wr)
            wr(MINLOC(wr)) = MAXVAL(wr) + 1d0
            eigenwert_sort(3) = MINVAL(wr)
            
            
            wr = 0d0
            rg =0d0
            asphericity = 0d0
            prolateness = 0d0
            rg = sum(eigenwert_sort)
            
            asphericity = eigenwert_sort(3) - 1d0/2d0 * (eigenwert_sort(1) + eigenwert_sort(2))
            
            prolateness = ((3 * eigenwert_sort(1) - rg)*(3*eigenwert_sort(2) - rg)*&
                ( 3*eigenwert_sort(3) - rg))/ rg**3
            
            rg = sqrt(rg)
            
            write(13,'(I3,1x, F16.6, 1x,F16.6, 1x,F16.6,1x, F16.6, 1x,F16.6,1x, F16.6,1x, F16.6)') i, rg,rg_pos, &
                asphericity, prolateness, eigenwert_sort(1), eigenwert_sort(2),eigenwert_sort(3) 
        end do
        close(12)
        close(13)
        deallocate(position)
        
        OPEN(33,file=trim("data" // filename(1:length-3) // "txt"),status='old', action="read")
        rg_av=0.d0
        rg_pos_av=0.d0
        asphericity_av=0.d0
        prolateness_av=0.d0
        DO j=1,conf
            READ(33,*)i,rg,rg_pos,asphericity,prolateness,lambda1,lambda2,lambda3
            rg_av=rg_av+rg
            rg_pos_av=rg_pos_av+rg_pos
            asphericity_av=asphericity_av+asphericity
            prolateness_av=prolateness_av+prolateness
        END DO
        CLOSE(33)
        rg_av=rg_av/dble(conf)
        rg_pos_av=rg_pos_av/dble(conf)
        asphericity_av=asphericity_av/dble(conf)
        prolateness_av=prolateness_av/dble(conf)
        if (form == "chain") then
            write(17, '(A, 1x, F16.6, 1x, F16.6, 1x, F16.6, 1x, F16.6)') stiffness, rg_av, rg_pos_av,&
                asphericity_av, prolateness_av    
        else if(form == "ring") then
            write(18, '(A, 1x, F16.6, 1x, F16.6, 1x, F16.6, 1x, F16.6)') stiffness, rg_av, rg_pos_av, &
                asphericity_av, prolateness_av    
        else if(form=="star") then
            write(19, '(A, 1x, F16.6, 1x, F16.6, 1x, F16.6, 1x, F16.6)') stiffness, rg_av, rg_pos_av, &
                asphericity_av, prolateness_av     
        else
            print*, "Unexpected Form!"
            exit
        end if
        
        
    end do    
    close(11)
    close(14)
    close(15)
    close(17)
    close(18)
    close(19)

  
    call system('gnuplot plot.gp') 
    

    end program GyrationTensorPolymers

