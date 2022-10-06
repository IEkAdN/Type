# Type
## Usage
In order to run Type, execute the following command:  

```
Type [R1.sam] [R2.sam] [reference.fasta]
```

## Inputs
As shown in the usage section, `Type` requires following files:
* Two SAM files containing mapping results of Illumina paired-end reads
  * We verified `Type` using SAM files of BWA-mem (version 0.7.15)
* A FASTA file containing reference sequences against which the reads were mapped

## Outputs
Homology and coverage are output to the standard output as TSV. An example of the output is as follows:  

|  Reference IDs  |  Homology  |  Coverage  |  Coverage*  |
| ---- | ---- | ---- | ---- |
|  Stx2c_AJ605767.1  |  0.999539  |  1  |  1  |
|  Stx2c_AY633461.1  |  0.999539  |  1  |  1  |
|  Stx2c_AF461167.1  |  0.999461  |  0.995165  |  0.995165  |
|  Stx2d_AF500189.1  |  0.9996  |  0.840451  |  0.871072  |

*`Type` calculates homology and coverage (the 3rd column) considering only top-hit reads. However, as supplementary information, it also calculates coverage (the 4th columns) considering all hit reads.  
`Type` calculates this homology by following 2 steps:  
1. Calculate homology at each reference position by averaging homology between the reference sequence and top-hit reads mapped to the position
1. Average the homology (calculated in the 1st step) across the reference sequence covered by top-hit reads (the averaged values are output)

## Installation
If you cannot run the Linux binary `Type` in this repository, you can compile `Type` as follows:

```
git clone https://github.com/IEkAdN/Type.git
cd Type/
make
```
