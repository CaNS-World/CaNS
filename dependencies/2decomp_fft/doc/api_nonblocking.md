## Non-blocking API for Overlap of Communication and Computation

Transpose-based parallelisation is inherently communication intensive. For large-scale applications, it is not unusual that communication accounts for more than half of the total cost. Application performance may be significantly improved if algorithms can be redesigned to allow overlap of communication and computation. From version 1.4, 2DECOMP&FFT provides a low-level communication API to facilitate such effort.

The API is based on ideas of non-blocking MPI collectives (such as MPI_IALLTOALL and MPI_IALLTOALLV) introduced in MPI version 3.

[Old users of 2DECOMP&FFT may recall the use of third-party library libNBC, which implemented the non-blocking MPI collectives using existing MPI 1 functions, to support such features. Using third party libraries is no longer necessary.]

### The API

Each of the four transposition routines in the base [decomposition library](api_decomposition.md) contains three key elements: algorithm to pack the MPI send buffers, MPI_ALLTOALL(V) communication, and algorithms to unpack the MPI receive buffers. When the non-blocking version of the MPI_ALLTOALL(V) is used, these routines are broken into smaller routines. For example, when transposing from X pencils to Y pencils, the blocking version of the communication routine is:
```
	call transpose_x_to_y(in, out, decomp)
```
The corresponding non-blocking routines are:
```
	call transpose_x_to_y_start(handle, in, out, sbuf, rbuf, decomp)
	call transpose_x_to_y_wait(handle, in, out, sbuf, rbuf, decomp)
```
The *start* routine packs the MPI send buffer, starts the non-blocking MPI_ALLTOALL(V) communication, and returns immediately. Later, a call to the corresponding *wait* routine ensures the communication is completed and then unpacks the MPI receive buffer. The first parameter `handle` is used to uniquely identify each communication session. Because several non-blocking communications may be ongoing at the same time, each has to define its own send buffer *sbuf* and receive buffer *rbuf*<a href="#note1" id="note1ref"><sup>1</sup></a>. It is up to the applications to supply (and if possible, reuse) these buffers, the size and shape of which should match the corresponding input array in and output array out. Between a *start* call and the corresponding *wait* call, the content of *sbuf* should not be modified and the content of *out* should not be referenced, to avoid unpredictable results. Other unrelated computations may be carried out while the communication is ongoing.

There are similar *start/wait* routines defined to all other transposition routines.

These routines are useful on systems with dedicated networking hardware to process the communication stack. On systems without such hardware, one has to call `MPI_TEST` explicitly from the user thread to progress the non-blocking communication. A utility routine is provided for this purpose:
```
	call transpose_test(handle)
```
This needs to be called from time to time from the computational part of application, in order to progress the communication identified by `handle`. Of course, the practical difficulty is where and how frequently this should be called, a matter that is entirely application dependent. 

Currently, the author is not aware of any stable and high-quality software implementation that progresses all-to-all type of communication asynchronously<a href="#note2" id="note2ref"><sup>2</sup></a>.

#### A Sample Application

To demonstrate the use of this API, a sample application (non_blocking) is provided to compute multiple independent FFTs, using both the blocking and non-blocking versions of the communication library. The idea of overlapping the communication of one 3D FFT and the computation of another, as described by Kandalla et al.[1], is implemented. The algorithm's pseudo-code looks like:
```
      1D FFT in X for V_1
      call transpose_x_to_y for V_1 (blocking)
      1D FFT in Y for V_1
      call transpose_y_z_start for V_1
      do k=2,N
        1D FFT in X for V_k
        call transpose_x_to_y for V_k (blocking)
        1D FFT in Y for V_k
        call transpose_y_to_z_start for V_k
        call transpose_y_to_z_wait for V_(k-1)
        1D FFT in Z for V_(k-1)
      end do
      call transpose_y_to_z_wait for V_N to complete
      1D FFT in Z for V_N
```

This algorithm compute multiple independent 3D FFTs on dataset *V<sub>k</sub> (k=1,N)*. As can be seen, the Y=>Z transpose for dataset *k* and the computation of 1D FFT in Z for dataset *k-1* are overlapped. Note that in the sample application the computations are done using loops of 1D FFTs, rather than with FFTW's advanced interface that allows multiple 1D FFTs to be done in one go. This design is to allow `MPI_TEST` calls to be inserted to progress the communication.

It is up to the application developers to identify opportunities in their algorithms that may benefit from this non-blocking API.

#### References

[1] K. Kandalla, H. Subramoni, K. Tomko, D. Pekurovsky, S. Sur and D.K. Panda, "High-performance and scalable non-blocking all-to-all with collective offload on InfiniBand clusters: a study with parallel 3D FFT", *Computer Science - Research and Development*, vol. 26(3-4):237-246, 2011.


---

<a id="note1" href="#note1ref"><sup>1</sup></a>The blocking version also needs to define send/recv buffers. But because there is only one communication at any time, the buffers are temporarily allocated as required by the library, or for performance reason defined globally and shared by multiple communication calls.

<a id="note2" href="#note2ref"><sup>2</sup></a>There are *asynchronous progress control* in Intel MPI library. However, the only supported non-blocking collective calls are *Ibcast*, *Ireduce* and *Iallreduce*.