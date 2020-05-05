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
    !****************************************************************************

    program GyrationTensorPolymers
    use LinearSystems
    use LAPACK95
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
    call system('DIR *.DAT /B > FILES.txt') !Windows
    !call system('ls -l *.dat./inFiles > FILES.txt') !Linux

    open(unit= 11,file="FILES.txt",action="read")
    open(unit =14, file="Average.txt", action="write")
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
            
            do j = 1, num
                read(12, *) position(j,1), position(j,2), position(j,3)
                schwerpunkt = schwerpunkt + position(j,1:3)
            end do
            
            
            schwerpunkt = schwerpunkt / dble(num)
            rg_pos = 0d0
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
                if (abs(gyrationTensor(2,1)) < 1d-10 .and. abs(gyrationTensor(3,1)) < 1d-10 .and. abs(gyrationTensor(3,2)) < 1d-10) then
                    exit
                end if
            end do
            do k= 1, 3
                eigenwert(k) = gyrationTensor(k,k)
            end do
           write(15, '(A,1x,A,1x,I3, 1x,F16.6, 1x,F16.6, 1x,F16.6, 1x,F16.6, 1x,F16.6, 1x,F16.6, 1x)') form, stiffness, i, eigenwert(1), eigenwert(2), eigenwert(3), wr(1), wr(2), wr(3)
            !Sortieren der Eigenwerte
            eigenwert_sort(1) = MINVAL(wr)
            wr(MINLOC(wr)) = MAXVAL(wr) + 1d0
            eigenwert_sort(2) = MINVAL(wr)
            wr(MINLOC(wr)) = MAXVAL(wr) + 1d0
            eigenwert_sort(3) = MINVAL(wr)
            
            
            wr = 0d0
            
           
            rg = sum(eigenwert_sort)
            
            asphericity = eigenwert_sort(3) - 1d0/2d0 * (eigenwert_sort(1) + eigenwert_sort(2))
            
            prolateness = ((3 * eigenwert_sort(1) - rg)*(3*eigenwert_sort(2) - rg)*( 3*eigenwert_sort(3) - rg))/ rg**3
            
            rg = sqrt(rg)
            
            write(13,'(I3,1x, F16.6, 1x,F16.6, 1x,F16.6,1x, F16.6, 1x,F16.6,1x, F16.6,1x, F16.6)') i, rg,rg_pos, asphericity, prolateness, eigenwert_sort(1), eigenwert_sort(2),eigenwert_sort(3) 
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
        rg_av=rg_av/dble(num)
        rg_pos_av=rg_pos_av/dble(num)
        asphericity_av=asphericity_av/dble(num)
        prolateness_av=prolateness_av/dble(num)
        write(14, '(A, 1x, A, 1x, F16.6, 1x, F16.6, 1x, F16.6, 1x, F16.6)') form, stiffness, rg_av, rg_pos_av, asphericity_av, prolateness_av
        
        
    end do    
    close(11)
    close(14)
    close(15)


    !! Store the data for each configuration j
    !!OPEN(33,file='data.dat',status='unknown')
    !!
    !!!DO j=1,m
    !!
    !!    ! Compute the center of mass of the polymer
    !!
    !!    ! Compute the gyration tensor A
    !!
    !!    ! Computer the eigenvalues of the gyration tensor
    !!
    !!    ! ... with home-made QR-algorithm
    !!    ! Compute  u1
    !!    ! Compute  H1
    !!    ! Compute the product H1A
    !!
    !!    ! Compute  u2
    !!    ! Compute  H2
    !!    ! Compute the product H2A
    !!
    !!    ! Store the three eigenvalues
    !!
    !!    ! ... with LAPACK
    !!
    !!    ! Store the three eigenvalues
    !!
    !!    ! Compute the gyration radius, the asphericity and the prolateness
    !!    ! from the invariants of the gyration tensor
    !!
    !!    ! Compute the gyration radius from the monomers positions directly
    !!
    !!!END DO
    !!
    !!CLOSE(30)
    !!CLOSE(33)
    !!
    !!! Average over the m independent configurations
    !!OPEN(33,file='data.dat',status='old')
    !!rg_av=0.d0
    !!rg_pos_av=0.d0
    !!asphericity_av=0.d0
    !!prolateness_av=0.d0
    !!DO j=1,m
    !!    READ(33,*)i,rg,rg_pos,asphericity,prolateness,lambda1,lambda2,lambda3
    !!    rg_av=rg_av+rg
    !!    rg_pos_av=rg_pos_av+rg_pos
    !!    asphericity_av=asphericity_av+asphericity
    !!    prolateness_av=prolateness_av+prolateness
    !!END DO
    !!CLOSE(33)
    !!rg_av=rg_av/dble(m)
    !!rg_pos_av=rg_pos_av/dble(m)
    !!asphericity_av=asphericity_av/dble(m)
    !!prolateness_av=prolateness_av/dble(m)
    !!
    !!! Save the averaged data
    !!OPEN(30,file='data-averaged.dat',status='unknown', access='append')
    !!WRITE(30,*)ncstring,sqrt(rg_av),sqrt(rg_pos_av),asphericity_av/rg_av,prolateness_av
    !!CLOSE(30)
    !!print *, 'Hello World'

    end program GyrationTensorPolymers

