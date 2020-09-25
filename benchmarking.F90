program benchmarking
    !
    !  Description:
    !    Running all benchmarks to test code performance, and reproducing all the test
    !
    use benchmarks
    implicit none
    
    call afterglow_syn_lc
    write(*, *) '=======  FINISHED  ======='
end program benchmarking