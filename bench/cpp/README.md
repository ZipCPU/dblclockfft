Here are the bench tests for the pipelined FFT.  In general, there's a
`*_tb.cpp` file corresponding to every unit within the FFT.  Feel free to
try them.

Be aware, however, the [fft_tb](fft_tb.cpp) doesn't truly
check for success--I just haven't gotten to the point of verifying that
the FFT result is *close enough* to the right answer in spite of actually
calculating the right answer.  Instead, it creates a data file that can be
read in Octave via [fft_tb.m](fft_tb.m).  That will show the first test output.
The second and subsequent outputs can be read via `k=k+1;` followed by calling
[plottst](plottst.m).

As another note (before I clean things up more), you'll need the `*.hex` files
in the same directory as the one you call [fft_tb](fft_tb.cpp) or
[fftstage_tb](fftstage_tb.cpp) from.

I expect the IFFT will work: it's just an FFT with conjugate twiddle factors,
although I haven't fully tested it yet.
