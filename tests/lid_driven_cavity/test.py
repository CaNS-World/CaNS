#!/usr/bin/env python
def test_ldc():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("vey_fld_0001500.bin",np.array([1,1,1]))
    data_ref = np.loadtxt("data_ldc_re1000.txt")
    islice = int(np.size(data[0,0,:])/2)
    np.testing.assert_allclose(data[0,islice,:], data_ref[:,1], rtol=1e-7, atol=0)
if __name__ == "__main__":
    test_ldc()
    print("Passed!")
