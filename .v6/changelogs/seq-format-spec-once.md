## Core: format specs parsed once

- `formatValue` parses each format spec string once, which cuts the allocations of `formatGenomic` (a genomic ruler's labels while panning).
