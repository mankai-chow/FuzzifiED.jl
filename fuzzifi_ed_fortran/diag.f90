module diag 

contains 

subroutine vec_prod(dim_d, dim_f, sym_q, nel, colptr, rowid, elval, st_d, st_f)
    use omp_lib 
    implicit none 
    integer(8), intent(in) :: dim_d, dim_f, sym_q, nel
    integer(8), intent(in) :: colptr(dim_d + 1), rowid(nel)
    complex(8), intent(in) :: elval(nel), st_d(dim_d) 
    complex(8) :: st_f(dim_f)
    complex(8) :: val
    complex(8), allocatable :: st_f1(:)
    integer(8) :: i, j, i1

    st_f = 0
    !$omp parallel shared(dim_d, dim_f, sym_q, nel, colptr, rowid, elval, st_d, st_f) private(st_f1, i, j, i1, val)
    allocate(st_f1(dim_d))
    st_f1 = 0 
    !$omp do 
    do i = 1, dim_d
        do j = colptr(i), colptr(i + 1) - 1
            i1 = rowid(j)
            val = elval(j)
            st_f1(i1) = st_f1(i1) + val * st_d(i)
            if (sym_q == 0 .or. i == i1) cycle 
            if (sym_q == 1) val = conjg(val)
            st_f1(i) = st_f1(i) + val * st_d(i1)
        end do
    end do
    !$omp end do
    !$omp critical 
    st_f = st_f + st_f1 
    !$omp end critical
    deallocate(st_f1)
    !$omp end parallel
end subroutine

subroutine scal_prod(dim_d, dim_f, sym_q, nel, colptr, rowid, elval, st_d, st_f, prod)
    use omp_lib 
    implicit none 
    integer(8), intent(in) :: dim_d, dim_f, sym_q, nel
    integer(8), intent(in) :: colptr(dim_d + 1), rowid(nel)
    complex(8), intent(in) :: elval(nel), st_d(dim_d), st_f(dim_f)
    complex(8), intent(out) :: prod
    complex(8) :: prod1, val
    integer(8) :: i, j, i1

    prod = 0
    !$omp parallel shared(dim_d, dim_f, sym_q, nel, colptr, rowid, elval, st_d, st_f, prod) private(prod1, i, j, i1, val)
    prod1 = 0 
    !$omp do 
    do i = 1, dim_d
        do j = colptr(i), colptr(i + 1) - 1
            i1 = rowid(j)
            val = elval(j)
            prod1 = prod1 + conjg(st_f(i1)) * val * st_d(i)
            if (sym_q == 0 .or. i == i1) cycle 
            if (sym_q == 1) val = conjg(val)
            prod1 = prod1 + conjg(st_f(i)) * val * st_d(i1)
        end do
    end do
    !$omp end do
    !$omp critical 
    prod = prod + prod1
    !$omp end critical
    !$omp end parallel
end subroutine

subroutine diagonalisation(dim, sym_q, nel, colptr, rowid, elval, nst, tol, eigval, eigvec)
    implicit none
    integer(8), intent(in) :: dim, sym_q, nel, nst
    integer(8), intent(in) :: colptr(dim + 1), rowid(nel)
    complex(8), intent(in) :: elval(nel)
    real(8), intent(in) :: tol
    complex(8), intent(out) :: eigval(nst + 1), eigvec(dim, nst + 1)

    integer(4) :: nit

    ! znaupd type_basis
    character :: bmat * 1, which * 2
    integer(4) :: ido, n, ncv, nev, ldv, iparam(11), ipntr(14), lworkl, info
    complex(8), allocatable :: resid(:), v(:,:), workd(:), workl(:)
    real(8), allocatable :: rwork(:)

    ! zneupd type_basis
    logical :: rvec
    integer(4) :: ldz
    character :: howmny * 1
    complex(8) :: sigma
    logical, allocatable :: select(:)
    complex(8), allocatable :: workev(:)

    ido = 0
    bmat = 'I'
    n = dim
    which = 'SR'
    allocate(resid(n))
    nev = min(nst, n - 1)
    ncv = min(max(2 * nev, nev + 10), n)
    allocate(v(n, ncv))
    ldv = n
    iparam(1) = 1
    iparam(3) = 300
    iparam(7) = 1
    allocate(workd(3 * n))
    lworkl = 3 * ncv ** 2 + 5 * ncv
    allocate(workl(lworkl))
    allocate(rwork(ncv))
    info = 0

    print *, 'Diagonisation begins.'
    call znaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, &
        ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)
    nit = 0
    do while (ido == 1 .or. ido == -1)
        call vec_prod(dim, dim, sym_q, nel, colptr, rowid, elval, workd(ipntr(1)),  workd(ipntr(2)))
        call znaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, &
            ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)
        nit = nit + 1
        if (mod(nit, 100) == 0) print *, 'Diagonisation, iteration : ', nit
    end do
    if (info < 0 .or. ido /= 99) print *, 'Errors in znaupd, info =', info

    rvec = .true.
    howmny = 'A'
    allocate(select(ncv))
    ldz = n
    allocate(workev(2 * ncv))

    call zneupd(rvec, howmny, select, eigval, eigvec, ldz, sigma, workev, bmat, n, which, nev, &
        tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)

    if (info == 0) print *, 'Diagonisation successful, total iteration : ', nit

    deallocate(resid)
    deallocate(v)
    deallocate(workd)
    deallocate(workl)
    deallocate(rwork)
    deallocate(select)
    deallocate(workev)
end subroutine

end module 