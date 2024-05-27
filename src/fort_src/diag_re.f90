module diag_re

contains 

subroutine vec_prod_re(dim_d, dim_f, sym_q, nel, colptr, rowid, elval, st_d, st_f)
    use omp_lib 
    implicit none 
    integer(8), intent(in) :: dim_d, dim_f, sym_q, nel
    integer(8), intent(in) :: colptr(dim_d + 1), rowid(nel)
    real(8), intent(in) :: elval(nel), st_d(dim_d) 
    real(8) :: st_f(dim_f)
    real(8) :: val
    real(8), allocatable :: st_f1(:)
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

subroutine scal_prod_re(dim_d, dim_f, sym_q, nel, colptr, rowid, elval, st_d, st_f, prod)
    use omp_lib 
    implicit none 
    integer(8), intent(in) :: dim_d, dim_f, sym_q, nel
    integer(8), intent(in) :: colptr(dim_d + 1), rowid(nel)
    real(8), intent(in) :: elval(nel), st_d(dim_d), st_f(dim_f)
    real(8), intent(out) :: prod
    real(8) :: prod1, val
    integer(8) :: i, j, i1

    prod = 0
    !$omp parallel shared(dim_d, dim_f, sym_q, nel, colptr, rowid, elval, st_d, st_f, prod) private(prod1, i, j, i1, val)
    prod1 = 0 
    !$omp do 
    do i = 1, dim_d
        do j = colptr(i), colptr(i + 1) - 1
            i1 = rowid(j)
            val = elval(j)
            prod1 = prod1 + st_f(i1) * val * st_d(i)
            if (sym_q == 0 .or. i == i1) cycle 
            prod1 = prod1 + st_f(i) * val * st_d(i1)
        end do
    end do
    !$omp end do
    !$omp critical 
    prod = prod + prod1
    !$omp end critical
    !$omp end parallel
end subroutine

subroutine diagonalisation_re(dim, sym_q, nel, colptr, rowid, elval, nst, tol, ncv_in, eigval, eigvec)
    implicit none
    integer(8), intent(in) :: dim, sym_q, nel, nst, ncv_in
    integer(8), intent(in) :: colptr(dim + 1), rowid(nel)
    real(8), intent(in) :: elval(nel)
    real(8), intent(in) :: tol
    real(8), intent(out) :: eigval(nst + 1), eigvec(dim, nst + 1)

    integer(4) :: nit

    ! dnaupd type_basis
    character :: bmat * 1, which * 2
    integer(4) :: ido, n, ncv, nev, ldv, iparam(11), ipntr(14), lworkl, info
    real(8) :: eigval_im(nst + 1)
    real(8), allocatable :: resid(:), v(:,:), workd(:), workl(:)
    real(8), allocatable :: rwork(:)

    ! dneupd type_basis
    logical :: rvec
    integer(4) :: ldz
    character :: howmny * 1
    real(8) :: sigmar, sigmai
    logical, allocatable :: select(:)
    real(8), allocatable :: workev(:)

    ido = 0
    bmat = 'I'
    n = dim
    which = 'SR'
    allocate(resid(n))
    nev = min(nst, n - 1)
    ncv = ncv_in
    allocate(v(n, ncv))
    ldv = n
    iparam(1) = 1
    iparam(3) = 300
    iparam(7) = 1
    allocate(workd(3 * n))
    lworkl = 3 * ncv ** 2 + 6 * ncv
    allocate(workl(lworkl))
    allocate(rwork(ncv))
    info = 0

    print *, 'Diagonisation begins.'
    call dnaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, &
        ldv, iparam, ipntr, workd, workl, lworkl, info)
    nit = 0
    do while (ido == 1 .or. ido == -1)
        call vec_prod_re(dim, dim, sym_q, nel, colptr, rowid, elval, workd(ipntr(1)),  workd(ipntr(2)))
        call dnaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, &
            ldv, iparam, ipntr, workd, workl, lworkl, info)
        nit = nit + 1
        if (mod(nit, 100) == 0) print *, 'Diagonisation, iteration : ', nit
    end do
    if (info < 0 .or. ido /= 99) print *, 'Errors in dnaupd, info =', info

    rvec = .true.
    howmny = 'A'
    allocate(select(ncv))
    ldz = n
    allocate(workev(3 * ncv))

    call dneupd(rvec, howmny, select, eigval, eigval_im, eigvec, ldz, sigmar, sigmai, workev, bmat, n, which, nev, &
        tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)

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