## VolumeServer: How it works

This document provides a high level overview of how the DensityServer works.

## Overview

- Data is stored in using block layout to reduce the number of disk seeks/reads each query requires.
- Data is downsampled by ``1/2``, ``1/4``, ``1/8``, ... depending on the size of the input.
- To keep the server response time/size small, each query is satisfied using the appropriate downsampling level.
- The server response is encoded using the [BinaryCIF](https://github.com/dsehnal/BinaryCIF) format.
- The contour level is preserved using relative instead of absolute values.

## Data Layout

To enable efficient access to the 3D data, the density values are stored in a "block level" format. 
This means that the data is split into ``NxNxN`` blocks (by default ``N=96``, which corresponds to ``96^3 * 4 bytes = 3.375MB`` disk read 
per block access and provides good size/performance ratio).  This data layout 
enables to access the data from a hard drive using a bounded number of disk seeks/reads which
greatly reduces the server latency.

## Downsampling 

- The input is density data with ``[H,K,L]`` number of samples along each axis (i.e. the ``extent`` field in the CCP4 header).
- To downsample, use the kernel ``C = [1,4,6,4,1]`` (customizable on the source code level) along each axis, because it is "separable":

    ```
    downsampled[i] = C[0] * source[2 * i - 2] + ... + C[4] * source[2 * i + 2]
    ```

    The downsampling step is applied in 3 steps:

    ```
    [H,K,L] => [H/2, K, L] => [H/2, K/2, L] => [H/2, K/2, L/2]
    ```

    (if the dimension is odd, the value ``(D+1)/2`` is used instead).

- Apply the downsampling step iteratively until the number of samples along the largest dimension is smaller than "block size" (or the smallest dimension has >2 samples).

## Satisfying the query

When the server receives a query for a 3D region, it chooses the the appropriate downsampling level based on the required details so that 
the number of voxels in the response is small enough. This enables sub-second response time even for the largest of entries.

### Encoding the response

The [BinaryCIF](https://github.com/dsehnal/BinaryCIF) format is used to encode the response. Floating point data are quantized into 1 byte values (256 levels) before being
sent back to the client. This quantization is performed by splitting the data interval into uniform pieces.

## Preserving the contour level

Downsampling the data results in changing of absolute contour levels. To mitigate this effect, relative values are always used when displaying the data.

- Imagine the input data points are ``A = [-0.3, 2, 0.1, 6, 3, -0.4]``: 
- Downsampling using every other value results in ``B = [-0.3, 0.1, 3]``.
- The "range" of the data went from (-0.4, 6) to (-0.3,3).
- Attempting to use the same absolute contour level on both "data sets" will likely yield very different results.
- The effect is similar if instead of skipping values they are averaged (or weighted averaged in the case of the ``[1 4 6 4 1]`` kernel) only not as severe.
- As a result, the "absolute range" of the data changes, some outlier values are lost, but the mean and relative proportions (i.e. deviation ``X`` from mean in ``Y = mean + sigma * X``) are preserved. 

----------------------

## Compression Analysis

- Downsampling: ``i-th`` level (starting from zero) reduces the size by approximate factor ``1/[(2^i)^3]`` (i.e. "cubic" of the frequency).
- BinaryCIF: CCP4 mode 2 (32 bit floats) is reduced by factor of 4, CCP4 mode 1 (16bit integers) by factor of 2, CCP4 mode 0 (just bytes) is not reduced. This is done by single byte quantization, but smarter than CCP4 mode 0
- Gzip, from observation:
  - Gzipping BinaryCIF reduces the size by factor ~2 - ~7 (2 for "dense" data such as x-ray density, 7 for sparse data such such an envelope of a virus)
  - Gzipping CCP4 reduces the size by 10-25% (be it mode 2 or 0)
- Applying the downsampling kernel helps with the compression ratios because it smooths out the values.

### Toy example:

```
Start with 3.5GB compressed density data in the CCP4 mode 2 format (32-bit float for each value)
    => ~4GB uncompressed CCP4
    => Downsample by 1/4 => 4GB * (1/4)^3 = 62MB
    => Convert to BinaryCIF => 62MB / 4 = ~16MB
    => Gzip: 2 - 8 MB depending on the "density" of the data 
        (e.g. a viral shell data will be smaller because it is "empty" inside)
```