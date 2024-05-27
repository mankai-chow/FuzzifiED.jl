! ifort -shared -fPIC -larpack -qopenmp -O3 -o lib_fuzzifi_ed.so ./fort_src/*.f90

module cfs
    
contains

subroutine count_cfs(no, nor, nqnu, qnu_s, qnu_o, modul, ncf, lid)
    implicit none
    integer(8), intent(in) :: no, nor, nqnu
    integer(8), intent(in) :: qnu_s(nqnu), qnu_o(no, nqnu), modul(nqnu)
    integer(8), intent(out) :: lid(kibset(0_8, no - nor) + 1)
    integer(8), intent(out) :: ncf
    integer(8) :: qnu_1(nqnu), i, j

    !$omp parallel shared(no, nor, nqnu, qnu_s, qnu_o, modul, ncf, lid) private(qnu_1, i, j)
    !$omp do
    do i = 0, kibset(0_8, no - nor) - 1
        qnu_1 = qnu_s
        do j = nor + 1, no
            if (kibits(i, j - nor - 1, 1_8) == 0) cycle
            qnu_1 = qnu_1 - qnu_o(j, :)
        end do
        lid(i + 2) = 0
        call count_cfs_rec(no, nor, nqnu, qnu_1, qnu_o, modul, lid(i + 2))
    end do
    !$omp end do
    !$omp end parallel
    lid(1) = 0
    do i = 1, kibset(0_8, no - nor)
        lid(i + 1) = lid(i + 1) + lid(i)
    end do
    ncf = lid(kibset(0_8, no - nor) + 1)
end subroutine

recursive subroutine count_cfs_rec(no, nom, nqnu, qnu_1, qnu_o, modul, ct)
    implicit none
    integer(8), intent(in) :: no, nom, nqnu
    integer(8), intent(in) :: qnu_1(nqnu), qnu_o(no, nqnu), modul(nqnu)
    integer(8), intent(inout) :: ct
    logical :: flag
    integer(8) :: i, nom1, cyc1

    flag = .true. 
    do i = 1, nqnu
        cyc1 = modul(i)
        if (cyc1 == 1 .and. qnu_1(i) == 0) cycle
        if (cyc1 > 1 .and. mod(qnu_1(i), cyc1) == 0) cycle
        flag = .false.
        if (cyc1 == 1 .and. qnu_1(i) < 0) return
    end do
    if (flag) then 
        ct = ct + 1 
        return
    end if
    do nom1 = 1, nom
        call count_cfs_rec(no, nom1 - 1, nqnu, qnu_1 - qnu_o(nom1, :), qnu_o, modul, ct)
    end do
end subroutine

subroutine generate_cfs(no, nor, nqnu, qnu_s, qnu_o, modul, ncf, lid, rid, conf)
    implicit none
    integer(8), intent(in) :: no, nor, nqnu, ncf
    integer(8), intent(in) :: qnu_s(nqnu), qnu_o(no, nqnu), modul(nqnu)
    integer(8), intent(in) :: lid(kibset(0_8, no - nor) + 1)
    integer(8), intent(out) :: conf(ncf), rid(kibset(0_8, nor))
    integer(8) :: qnu_1(nqnu), i, j, ct

    rid = 0
    !$omp parallel shared(no, nor, nqnu, qnu_s, qnu_o, modul, ncf, lid, rid, conf) private(qnu_1, i, j, ct)
    !$omp do
    do i = 0, kibset(0_8, no - nor) - 1
        if (lid(i + 1) == lid(i + 2)) cycle 
        qnu_1 = qnu_s
        do j = nor + 1, no
            if (kibits(i, j - nor - 1, 1_8) == 0) cycle
            qnu_1 = qnu_1 - qnu_o(j, :)
        end do
        ct = lid(i + 1)
        call generate_cfs_rec(no, nor, nor, nqnu, qnu_1, qnu_o, modul, ncf, ct, kishft(i, nor), lid, rid, conf)
    end do
    !$omp end do
    !$omp end parallel
end subroutine

recursive subroutine generate_cfs_rec(no, nor, nom, nqnu, qnu_1, qnu_o, modul, ncf, ct, tmp, lid, rid, conf)
    implicit none
    integer(8), intent(in) :: no, nor, nom, nqnu, tmp, ncf
    integer(8), intent(in) :: qnu_1(nqnu), qnu_o(no, nqnu), modul(nqnu)
    integer(8), intent(in) :: lid(kibset(0_8, no - nor) + 1)
    integer(8), intent(inout) :: ct, conf(ncf), rid(kibset(0_8, nor))
    logical :: flag
    integer(8) :: i, nom1, li, ri, cyc1

    flag = .true. 
    do i = 1, nqnu
        cyc1 = modul(i)
        if (cyc1 == 1 .and. qnu_1(i) == 0) cycle
        if (cyc1 > 1 .and. mod(qnu_1(i), cyc1) == 0) cycle
        flag = .false.
        if (cyc1 == 1 .and. qnu_1(i) < 0) return
    end do
    if (flag) then 
        ct = ct + 1 
        conf(ct) = tmp 
        ri = kibits(tmp, 0_8, nor)
        if (rid(ri + 1) > 0) return
        li = kibits(tmp, nor, no - nor)
        rid(ri + 1) = ct - lid(li + 1)
        return
    end if
    do nom1 = 1, nom
        call generate_cfs_rec(no, nor, nom1 - 1, nqnu, qnu_1 - qnu_o(nom1, :), qnu_o, modul, ncf, ct, kibset(tmp, nom1 - 1), lid, rid, conf)
    end do
end subroutine

function search_conf(no, nor, lid, rid, cf) result(i)
    implicit none
    integer(8), intent(in) :: no, nor
    integer(8), intent(in) :: lid(kibset(0_8, no - nor) + 1), rid(kibset(0_8, nor))
    integer(8), intent(in) :: cf
    integer(8):: i
    i = lid(kibits(cf, nor, no - nor) + 1) + rid(kibits(cf, 0_8, nor) + 1)
end function

end module
