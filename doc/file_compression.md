# File compression

netCDF compression has been optimised and made available in `nco` tools

## Command example
```bash
ncks -7 -L 5 --baa=4 --ppc default=3 in.nc out.nc
```
where

- `-7` 
- `-L 5` indicates the _deflate level_ (_Lempel-Ziv_ deflation/compression), which should be an integer from 0 (no compression) to 9 (high compression). Compression takes longer when using a larger level.
- `--baa=4` sets the _bit-adjustment-algorithm_. See [here](https://nco.sourceforge.net/nco.html#Quantization-Algorithms-1) for details.
- `--ppc default=N` sets the default level of the quantization algorithm.

## Results

| Deflation level | File size |  
| :-----------: | ----------- |
| 0  | 1.1G |
| 1  | 147M |
| 2  | 146M |
| 3  | 143M |
| 4  | 139M |
| 5  | 137M |
| 6  | 135M |
| 7  | 134M |
| 8  | 131M |
| 9  | 131M |

Whatever the deflation level, the file size is approximately divided by 10.      
Following the recommendations, we will continue the test with `L=5`, as it offers a good compromise.

We then use a `bit-adjustment-algorithm`:

| Deflation level | File size |  
| :-----------: | ----------- |
Granular BitRound (NSD) | 64M |
BitGroom (NSD) | 82M |
BitShave (NSD) | 76M |
BitSet (NSD) | 76M |
DigitRound (NSD) | 71M | 
Granular BitRound (NSD) | 64M |  
BitGroomRound (NSD) | 76M | 
HalfShave (NSB) | 65M | 
BruteForce (NSD) | 70M | 
BitRound (NSB) | 65M |

The choice of the algorithm can further reduce the file size. The best results are obtained with the _Granular BitRound_ algorithm.

## References

- https://gmd.copernicus.org/articles/9/3199/2016/
- https://www.nature.com/articles/s43588-021-00156-2.pdf
- https://github.com/observingClouds/xbitinfo
- http://climate-cms.wikis.unsw.edu.au/NetCDF_Compression_Tools

```bash
for ii in $(seq 1 9); do 
    ncks -L ${ii} test01.nc test01_deflation-${ii}.nc 
done
```
    

Granular BitRound (NSD)
BitGroom (NSD)
BitShave (NSD)
BitSet (NSD)
DigitRound (NSD)
Granular BitRound (NSD)
BitGroomRound (NSD)
HalfShave (NSB)
BruteForce (NSD)
BitRound (NSB)