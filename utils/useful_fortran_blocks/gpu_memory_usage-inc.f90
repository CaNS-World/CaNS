#if defined(_OPENACC)
      block
        use openacc
        real :: mem_use,mem_tot,mem_free
        integer :: mydev
        integer(acc_device_kind) :: dev_type
        dev_type = acc_get_device_type()
        mydev    = acc_get_device_num(dev_type)
        mem_tot  = acc_get_property(mydev,dev_type,acc_property_memory     )/(1.*1024**3)
        mem_free = acc_get_property(mydev,dev_type,acc_property_free_memory)/(1.*1024**3)
        mem_use  = mem_tot-mem_free
        if(myid == 0) print*,'Device memory usage (Gb): ', mem_use, ' / ', mem_tot
#if defined(__NVCOMPILER) && 0
        block
          use accel_lib
          !
          ! alternatively to what is done above
          !
          !mem_tot = acc_get_memory()/(1.*1024**3)
          !mem_use = (acc_get_memory()-acc_get_free_memory())/(1.*1024**3)
          if(myid == 0) print*,'Device memory allocated via OpenACC data or compute regions (Gb): ', acc_bytesalloc()/(1.*1024**3)
          if(myid == 0) call acc_present_dump ! print all data present on the device corresponding to myid == 0
        end block
#endif
      end block
#endif
