#!/usr/bin/env python
def test_ldc():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    data_one_step,xp,yp,zp,xu,yv,zw = read_single_field_binary("data-one-step/","vey_fld_0001500.bin",np.array([1,1,1]))
    data_two_step,xp,yp,zp,xu,yv,zw = read_single_field_binary("data-two-step/","vey_fld_0001500.bin",np.array([1,1,1]))
    np.testing.assert_allclose(data_one_step[:,:,:], data_two_step[:,:,:], rtol=1e-10, atol=0)
if __name__ == "__main__":
    test_ldc()
    print("Passed!")
