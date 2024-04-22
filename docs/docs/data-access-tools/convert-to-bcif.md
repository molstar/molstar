# Convert CIF to BinaryCIF
BinaryCIF is an efficient, binary flavor of the CIF format. See [specification](https://github.com/molstar/BinaryCIF) and [publication](https://doi.org/10.1371/journal.pcbi.1008247) for further details.

This script reads data in CIF format and converts it lossless to a BinaryCIF file that can be read by Mol* or other 
applications.

## Example
```sh
node lib/commonjs/cli/cif2bcif/index.js file.cif file.bcif
```

## Usage
| Argument | Description |
| --- | --- |
| `src` | Source CIF to convert (can be gzipped) |
| `out` | Generated BinaryCIF output path |
| `-c` | Path to optional config file |
| `-f` | Path to optional filter file |

```sh
index.js [-h] [-c CONFIG] [-f FILTER] src out
```

### Config file
Controls how certain columns will be encoded. This is a JSON array of instructions:
```ts
interface EncodingStrategyHint {
    categoryName: string,
    columnName: string,
    encoding: 'pack' | 'rle' | 'delta' | 'delta-rle',
    precision?: number
}
```
Identify a particular CIF columns by its name and override the encoding by Integer Packing, Run-Length Encoding, Delta 
Encoding, or Delta & Run-Length Encoding. You can optionally control the precision if dealing with float values.

### Filter file
Specifies which categories and columns will be written. This is a plain text file, each line represents one entry. 
You can specify explicitly which categories or columns to include by adding `category_name` or 
`category_name.field_name`. You can also choose to ignore some categories or columns by adding `!category_name` or 
`!category_name.field_name`.