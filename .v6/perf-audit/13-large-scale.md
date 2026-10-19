# Slice 13: large-scale structures

The pipeline of [00-profile.md](00-profile.md) at 1–5M atoms and 1.7k–11k units, with and without symmetry sharing, plus
the interactive paths that scale with units. Measured 2026-10-19 on master (`a0652c96a`, the audit fixes #1986–#1992
merged), Node 25.6, Apple M4 (10 cores, 16 GB), headless, no GPU. Read-only on `packages/`.

**Machine load.** Other sessions ran benches at the same time: load average 2–14, swap 10–15 GB used. Absolute times are
noisy: the same 186k-atom stage differed up to 5× between runs (`flat-32` gaussian overpaint 330 vs 1,782 ms in the
first pass). The tables below are a second pass at load 2–8; allocated bytes, retained bytes and render-buffer sizes are
reproducible. 1M-vs-5M ratios compare separate processes, so they are indicative; super-linear verdicts below rest on
the isolated reruns (`-iso`: one process, only the stage and its prerequisites) and on the profiles.

## Fixtures

All from `1tup.bcif` (5,828 atoms, 9 units incl. water; C 1 2 1, 4 operators per cell with the centering), expanded by
its crystal operators over a block of cells (`buildSymmetryRange`), then written two ways by
`tests/src/pipeline/fixtures.ts`:

| Id         | Cells | Copies |     Atoms |  Units | Symmetry groups | File (local, `tests/build/pipeline/fixtures/`)                             |
| ---------- | ----- | -----: | --------: | -----: | --------------: | -------------------------------------------------------------------------- |
| `asm-32`   | 2×2×2 |     32 |   186,496 |    288 |               9 | `1tup-asm-32.cif` (607 KB): asymmetric unit + assembly '1' of 32 operators |
| `asm-192`  | 4×4×3 |    192 | 1,118,976 |  1,728 |               9 | `1tup-asm-192.cif` (621 KB)                                                |
| `asm-864`  | 6×6×6 |    864 | 5,035,392 |  7,776 |               9 | `1tup-asm-864.cif` (679 KB)                                                |
| `flat-32`  | 2×2×2 |     32 |   186,496 |    416 |             416 | `1tup-flat-32.bcif` (3.4 MB): one model; = 00-profile's `big-bcif`         |
| `flat-192` | 4×4×3 |    192 | 1,118,976 |  2,496 |           2,496 | `1tup-flat-192.bcif` (20 MB); `.cif` text 124 MB                           |
| `flat-864` | 6×6×6 |    864 | 5,035,392 | 11,232 |          11,232 | `1tup-flat-864.bcif` (90 MB)                                               |

`asm-*` is the capsid / large-assembly case: the file holds one copy and `pdbx_struct_oper_list` the operators, the
structure comes from `buildAssembly`, and all copies share 9 symmetry groups (instanced geometry). `flat-*` is a large
deposited entry or a cellPACK-style model: every unit is its own group, nothing is shared. Water splits into more units
when written as one model, hence 11,232 vs 7,776. `label_asym_id` A–M repeats in every copy of a `flat-*` file
(non-standard; Mol\* still splits units per copy).

## How to rerun

The bench is in the repo (`tests/src/pipeline/`, run through `scripts/pipeline-bench.mjs`). The script bundles the TS
sources with esbuild under the `molstar-src` condition, so no package build is needed. **Node 25+** is required for the
Alloc column (`v8.getHeapStatistics().total_allocated_bytes`; older Node prints a warning and NaN).

```sh
node scripts/pipeline-bench.mjs fixtures                                                 # once, ~10 s, all three sizes
node scripts/pipeline-bench.mjs run --fixture asm-192 --runs 2 --out out/asm-192.json
node scripts/pipeline-bench.mjs run --fixture flat-864 --skip 'repr:molecular-surface*,repr:interactions*'
node scripts/pipeline-bench.mjs run --fixture asm-864 --quality low --only 'repr:*,overpaint*,transparency*,release:*'
node scripts/pipeline-bench.mjs run --fixture flat-192 --profile prof/flat-192           # CPU, one .cpuprofile per stage
node scripts/pipeline-bench.mjs run --fixture flat-32 --profile prof/a --alloc --no-warmup # sampled JS-heap allocations
node scripts/pipeline-bench.mjs analyze prof/flat-192 10                                 # self / total per stage, TS lines
```

Fixtures: also `7qpd` and `1tup-2x2x2` (00-profile), or a path (`--fixture x.bcif [--assembly 1]`). `--heap <MB>` sets
V8's old space (default 14000). Stages print as they finish: ms, MB allocated, MB retained, peak RSS, GC pause ms.

Method changes against 00-profile:

- Each representation runs its theme change and interactive stages right after it is built, then `release:<repr>`
  destroys it. The release stage shows the memory given back as negative retained (−10.5 GB for ball-and-stick on
  `flat-864`), so peak memory holds one representation, not all.
- GC: full `last-resort` collections with a short wait before and after each stage. Freed array buffers are swept
  concurrently, and without the wait their release was credited to the next stage.
- `--only`/`--skip` keep the prerequisites (parse, model, structure, lookup, bonds, secondary structure, component
  producers) running unmeasured. Valence model, interactions and ASA are not prerequisites: a representation that needs
  them computes them inside its own stage, so the interactions representation rows below were run with the
  `interactions` stage on.
- Profiles carry their own source map (`run --profile` copies it), and `analyze` names GC `(garbage collector)`.
- One run per size, two for the 32-copy sizes: the 5M stages take minutes.

**Visual quality.** Representations default to `quality: 'auto'`, which picks the quality from the root structure's
element count (`DefaultQualityThresholds`, [repr/util.ts:91](../../packages/graphics/src/repr/util.ts#L91)):

- 186k atoms → `low`.
- ≥1M atoms → `lowest`: 8 Å surfaces, 1 linear and 2 radial segments, no color smoothing (resolution ≥ 3).

The default run therefore measures cheaper geometry at 1M and 5M than at 186k. A second table pins `--quality low` (what
186k gets) to show how the same geometry scales.

## Stage tables

ms / MB allocated (JS heap) / MB retained (heap + array buffers); negative retained = memory given back. `–`: not run.

### Load and properties

| Stage                          |         asm-32 |           asm-192 |            asm-864 |         flat-32 |          flat-192 |             flat-864 |
| ------------------------------ | -------------: | ----------------: | -----------------: | --------------: | ----------------: | -------------------: |
| `parse`                        |      9 / 2 / 2 |         7 / 2 / 2 |         10 / 3 / 2 |       7 / 3 / 3 |       20 / 3 / 19 |          89 / 5 / 86 |
| `model`                        |    14 / 15 / 1 |       15 / 24 / 1 |        16 / 28 / 1 |  139 / 134 / 58 |   567 / 616 / 343 |   2933 / 2696 / 1535 |
| `structure`                    |    6 / 5 / 0.5 |        8 / 23 / 2 |        16 / 43 / 8 |      13 / 9 / 2 |       36 / 47 / 9 |       261 / 236 / 40 |
| `lookup+boundary`              |  27 / 26 / 0.3 |     40 / 83 / 0.4 |    126 / 245 / 0.5 |    101 / 58 / 8 |    348 / 233 / 45 |     1409 / 888 / 202 |
| `bonds-intra`                  |  18 / 11 / 0.5 |     14 / 19 / 0.6 |      17 / 20 / 0.6 |    230 / 41 / 7 |   1200 / 256 / 39 |    7305 / 1152 / 175 |
| `bonds-inter`                  |   73 / 103 / 2 |     388 / 746 / 9 |   2132 / 3648 / 38 |   55 / 63 / 0.1 |   318 / 412 / 0.3 |    2959 / 1990 / 0.7 |
| `secondary-structure`          |  0.8 / 0.5 / 0 |     0.9 / 0.6 / 0 |        1 / 0.4 / 0 |   3 / 0.8 / 0.3 |        9 / 10 / 2 |          51 / 17 / 7 |
| `valence-model`                |  66 / 94 / -12 |     228 / 491 / 6 |   1198 / 2132 / 28 |   69 / 98 / -30 |     264 / 496 / 7 |     1186 / 2156 / 31 |
| `interactions`                 | 551 / 491 / 20 | 2229 / 1242 / 130 | 12998 / 5160 / 579 | 1005 / 710 / 11 | 4837 / 3926 / 255 | 30793 / 17654 / 1143 |
| `asa`                          |  994 / 112 / 2 |   5633 / 626 / 14 |  31379 / 2818 / 61 |   1166 / 16 / 2 |    5611 / 74 / 14 |     28350 / 314 / 61 |
| `theme-eval:element-symbol`    |      9 / 2 / 0 |      26 / 1 / 0.1 |      121 / 3 / 0.1 |      10 / 2 / 0 |      31 / 2 / 0.1 |        274 / 3 / 0.1 |
| `theme-eval:residue-name`      |     10 / 5 / 0 |       44 / 27 / 0 |      237 / 117 / 0 |     8 / 0.4 / 0 |      30 / 0.8 / 0 |          154 / 1 / 0 |
| `component:chain`              |  0.6 / 0.4 / 0 |     0.5 / 0.3 / 0 |        1 / 0.3 / 0 |   0.7 / 0.3 / 0 |       1 / 0.2 / 0 |          4 / 0.8 / 0 |
| `component:chain-derived`      |  0.5 / 0.4 / 0 |       1 / 0.7 / 0 |        1 / 0.6 / 0 |   0.9 / 0.6 / 0 |     0.9 / 0.7 / 0 |          2 / 0.8 / 0 |
| `component:every-copy`         |  0.5 / 0.3 / 0 |     0.7 / 0.1 / 0 |      2 / 0.8 / 0.1 |     1 / 0.4 / 0 |       4 / 0.6 / 0 |         23 / 2 / 0.1 |
| `component:every-copy-derived` |      6 / 6 / 0 |       28 / 43 / 0 |      131 / 186 / 0 |       8 / 5 / 0 |     46 / 34 / 0.2 |      204 / 119 / 0.6 |
| `component:polymer`            |      5 / 7 / 0 |     14 / 32 / 0.1 |     78 / 140 / 0.4 |       8 / 6 / 0 |    30 / 33 / -0.1 |     457 / 145 / -0.2 |
| `component:polymer-derived`    |   71 / 101 / 1 |     414 / 749 / 6 |   2033 / 3756 / 30 |   73 / 55 / 0.2 |   426 / 398 / 0.7 |      1954 / 1976 / 3 |
| peak RSS (MB)                  |           1047 |              1202 |               2531 |            1185 |              2766 |                 4117 |

Components: `chain` = a bundle of one chain under one operator (plugin `StructureComponent` 'bundle'), `every-copy` =
that chain in every copy, `polymer` = a static polymer query; `-derived` = what a representation on the component
computes first (boundary, lookup, inter-unit bonds, symmetry groups).

Text mmCIF (`flat-192.cif`, 124 MB, first pass): parse 1.0–1.2 s with 320 MB retained token buffers, model 0.7–0.8 s.

### Representations, default quality (`auto`: `low` at 186k, `lowest` at 1M and 5M)

| Stage                                         |          asm-32 |         asm-192 |            asm-864 |           flat-32 |           flat-192 |               flat-864 |
| --------------------------------------------- | --------------: | --------------: | -----------------: | ----------------: | -----------------: | ---------------------: |
| `repr:cartoon`                                |     23 / 37 / 4 |     25 / 37 / 3 |        69 / 43 / 6 |   376 / 374 / 126 |  1285 / 1377 / 379 |    11382 / 5988 / 1697 |
| `theme-apply:cartoon:secondary-structure`     |     2 / 0.4 / 0 |   2 / 0.2 / 0.1 |      4 / 0.5 / 0.1 |        27 / 7 / 0 |    112 / 38 / -0.5 |         951 / 132 / -1 |
| `repr:ball-and-stick`                         |    42 / 33 / 10 |    39 / 34 / 11 |       83 / 31 / 24 |   811 / 414 / 390 | 6250 / 2272 / 2334 | 191182 / 10128 / 10497 |
| `theme-apply:ball-and-stick:chain-id`         |       6 / 2 / 0 |     7 / 0.3 / 0 |        8 / 0.3 / 0 |    84 / 11 / -0.1 |    584 / 53 / -0.6 |        3780 / 190 / -1 |
| `repr:spacefill`                              |       9 / 6 / 1 |      13 / 9 / 2 |        20 / 10 / 7 |    194 / 126 / 42 |    818 / 691 / 253 |       2636 / 565 / 197 |
| `theme-apply:spacefill:chain-id`              |     3 / 0.9 / 0 |     2 / 0.3 / 0 |        4 / 0.6 / 0 |        29 / 6 / 0 |    158 / 31 / -0.2 |          443 / 0.8 / 0 |
| `repr:gaussian-surface`                       |     41 / 25 / 4 |     18 / 16 / 3 |        52 / 39 / 8 |   840 / 158 / 195 |    930 / 544 / 324 |     6668 / 2384 / 1459 |
| `theme-apply:gaussian-surface:element-symbol` |     32 / 42 / 4 |     3 / 0.8 / 0 |        3 / 0.7 / 0 |   392 / 151 / -72 |    153 / 31 / -0.1 |       672 / 117 / -0.5 |
| `repr:molecular-surface`                      |   311 / 144 / 3 |    171 / 96 / 3 |      228 / 138 / 8 | 6856 / 1529 / 171 | 21306 / 8413 / 323 |                      – |
| `theme-apply:molecular-surface:residue-name`  |   16 / 10 / 0.1 |     4 / 0.9 / 0 |        3 / 0.9 / 0 |   267 / 121 / -64 |    147 / 30 / -0.1 |                      – |
| `repr:interactions`                           | 205 / 101 / 166 | 807 / 423 / 442 | 5292 / 1761 / 1988 |   456 / 178 / 473 | 3311 / 1046 / 1348 |                      – |
| peak RSS (MB)                                 |            1047 |            1202 |               2531 |              1185 |               2766 |                   4117 |

### Representations at pinned `low` quality

| Stage                                         |     flat-32-low |      asm-192-low |        asm-864-low |        flat-192-low |         flat-864-low |
| --------------------------------------------- | --------------: | ---------------: | -----------------: | ------------------: | -------------------: |
| `repr:cartoon`                                | 225 / 418 / 127 |      30 / 50 / 5 |        51 / 53 / 8 |   1378 / 2080 / 783 |  13248 / 9024 / 3512 |
| `theme-apply:cartoon:secondary-structure`     |    15 / 8 / 0.2 |    2 / 0.8 / 0.1 |      3 / 0.2 / 0.1 |      88 / 40 / -0.3 |     612 / 135 / -0.8 |
| `repr:ball-and-stick`                         | 539 / 454 / 395 |     29 / 30 / 11 |       57 / 30 / 24 |  3883 / 2457 / 2364 |                    – |
| `theme-apply:ball-and-stick:chain-id`         |   53 / 9 / -0.1 |      7 / 0.5 / 0 |        4 / 0.3 / 0 |     335 / 48 / -0.6 |                    – |
| `repr:spacefill`                              |  129 / 127 / 43 |       10 / 8 / 2 |        31 / 10 / 7 |     826 / 692 / 253 |                    – |
| `theme-apply:spacefill:chain-id`              |   18 / 7 / -0.1 |      2 / 0.3 / 0 |        4 / 0.5 / 0 |     143 / 30 / -0.2 |                    – |
| `repr:gaussian-surface`                       | 600 / 160 / 195 |      33 / 37 / 5 |       67 / 46 / 10 |   5741 / 894 / 1170 |  42432 / 3632 / 5265 |
| `theme-apply:gaussian-surface:element-symbol` | 326 / 161 / -71 |      61 / 99 / 8 |        58 / 64 / 8 |   2116 / 727 / -440 | 10385 / 2998 / -1989 |
| `repr:molecular-surface`                      |               – |    275 / 168 / 5 |     341 / 214 / 10 | 45534 / 8703 / 1024 |                    – |
| `theme-apply:molecular-surface:residue-name`  |               – |    22 / 19 / 0.1 |        51 / 73 / 0 |   1854 / 626 / -383 |                    – |
| `repr:interactions`                           |               – | 1309 / 416 / 980 | 6775 / 1767 / 4416 |                   – |                    – |
| peak RSS (MB)                                 |             823 |             3179 |               3499 |                3413 |                 5122 |

`flat-864-low` ran cartoon and gaussian surface only.

### Isolated reruns (flat; one process each, model, intra bonds and three representations)

| Stage                        | flat-192           | flat-864                      | ×(5M / 1M) |  GC pause 1M / 5M |
| ---------------------------- | ------------------ | ----------------------------- | ---------: | ----------------: |
| `model`                      | 697 / 615 / 342    | 3688 / 2674 / 1535            |        5.3 |       89 / 745 ms |
| `bonds-intra`                | 1416 / 254 / 39    | 6412 / 1150 / 175             |        4.5 |        12 / 44 ms |
| `repr:cartoon` (`lowest`)    | 1100 / 1394 / 380  | 5925 / 6022 / 1698            |        5.4 |      114 / 629 ms |
| `theme-apply:cartoon`        | 124 / 40 / 0       | 449 / 137 / −1                |        3.6 |        11 / 49 ms |
| `repr:ball-and-stick`        | 4581 / 2454 / 2364 | **41780** / 10887 / **10633** |    **9.1** | 378 / **5396** ms |
| `theme-apply:ball-and-stick` | 398 / 46 / −1      | 2020 / 189 / −1               |        5.1 |        12 / 66 ms |
| `release:ball-and-stick`     | 4 / 0 / −2326      | 26 / 0 / −10463               |            |                   |
| `repr:spacefill`             | 767 / 692 / 253    | 2637 / 552 / 236              |        3.4 |       67 / 103 ms |

Peak RSS stayed at 4.8 GB while ball-and-stick retained 10.6 GB, so most of it was not resident: paged out or compressed
by macOS. In the full runs (props and components held, machine swapping) the same 5M ball-and-stick took 191 s (16.9 s
of GC pauses).

### Render buffers (sum of typed arrays in the render objects' values)

| Repr                                | asm-864              | flat-192    | flat-864          | Largest                                        |
| ----------------------------------- | -------------------- | ----------- | ----------------- | ---------------------------------------------- |
| cartoon (`lowest` / `low`)          | 2.2 / 3.5 MB         | 60 / 211 MB | 271 / **948 MB**  | `elements`, `aPosition`, `aNormal`             |
| ball-and-stick (impostor)           | 21.6 MB              | 991 MB      | **4,459 MB**      | `aMapping`, `aStart`, `aEnd` 850 MB each at 5M |
| spacefill (impostor)                | 6.5 MB               | 27 MB       | 197 MB            | `tPositionGroup`                               |
| gaussian surface (`lowest` / `low`) | 6.7 / 9.2 MB         | 12 / 437 MB | 52 / **1,964 MB** | `elements`, `aPosition`, `aNormal`             |
| interactions (`auto` / `low`)       | 1,986 / **4,414 MB** | 832 MB      | –                 | dashed mesh cylinders, `aPosition`/`aNormal`   |

Symmetric assemblies upload one copy of the geometry plus per-instance data. At 864 copies `tMarker` (group × instance)
is the largest buffer: 14 MB for ball-and-stick. Unshared structures upload everything. Retained vs buffers:
ball-and-stick 2.4× (the merged copy, 07/P7); cartoon 6.3× (`lowest`: small buffers, the retained also holds per-visual
geometry and trace data).

### Interactive paths, one chain (ms)

`hover` = highlight + remove highlight; `select` one chain; `select-every-copy` the same chain in all copies (192 or 864
units); `overpaint` / `transparency` one chain (a bundle layer, as the plugin's overpaint transform).

| Stage                                | asm-32 | asm-864 | flat-32 | flat-192 | flat-864 | flat-192-iso | flat-864-iso |
| ------------------------------------ | -----: | ------: | ------: | -------: | -------: | -----------: | -----------: |
| `hover:cartoon`                      |    0.6 |       1 |       3 |       13 |      191 |           18 |           65 |
| `select:cartoon`                     |    0.4 |     0.7 |       1 |        7 |       36 |            8 |           34 |
| `select-every-copy:cartoon`          |    0.5 |       2 |       3 |       39 |  **685** |           24 |      **200** |
| `overpaint:cartoon`                  |    0.5 |       2 |       4 |       44 |      152 |           24 |          106 |
| `transparency:cartoon`               |    0.8 |      29 |       5 |       39 |      184 |            – |            – |
| `hover:ball-and-stick`               |      1 |       2 |       2 |       15 |      108 |           15 |           68 |
| `select:ball-and-stick`              |    0.6 |     0.8 |       1 |        7 |       54 |            6 |           38 |
| `select-every-copy:ball-and-stick`   |      3 |      56 |       7 |       55 |  **506** |           40 |      **360** |
| `overpaint:ball-and-stick`           |      1 |       7 |       5 |       32 |      255 |           34 |          665 |
| `transparency:ball-and-stick`        |      4 |     102 |       9 |       44 |      317 |            – |            – |
| `hover:gaussian-surface`             |    0.4 |     0.5 |       1 |        8 |       35 |            – |            – |
| `select-every-copy:gaussian-surface` |    0.2 |     0.8 |       1 |       18 |      196 |            – |            – |
| `overpaint:gaussian-surface`         |    401 |       2 |     301 |       18 |       86 |            – |            – |
| `transparency:gaussian-surface`      |    265 |      19 |     208 |       19 |      197 |            – |            – |

At pinned `low` quality (color smoothing on); ms / MB allocated / peak RSS after (MB):

| Stage                            | flat-32         | asm-192               | asm-864                 | flat-192          | flat-864                |
| -------------------------------- | --------------- | --------------------- | ----------------------- | ----------------- | ----------------------- |
| `overpaint:gaussian-surface`     | 238 / 200 / 823 | **2421** / 914 / 2610 | **12564** / 4006 / 3311 | 1690 / 896 / 3413 | **11673** / 3907 / 5122 |
| `transparency:gaussian-surface`  | 169 / 157 / 823 | 1551 / 865 / 2610     | 6725 / 3704 / 3499      | 998 / 887 / 3413  | 6072 / 3901 / 5122      |
| `overpaint:molecular-surface`    | –               | 1841 / 83 / 3179      | **11697** / 3279 / 3499 | 1444 / 727 / 3413 | –                       |
| `transparency:molecular-surface` | –               | 1043 / 3 / 3179       | 6382 / 477 / 3499       | 1090 / 734 / 3413 | –                       |

On asm-192 peak RSS went from 426 MB before the surface overpaint to 2,610 MB after it. The clears (`overpaint-clear`,
`transparency-clear`) cost ≤ 0.35 s everywhere. Spacefill and molecular-surface marks are ≤ 20 ms at every size.

## Scaling and super-linear stages

Ratio of the 5M to the 1M time (×4.5 atoms and units; linear ≈ 4.5).

| Stage                                          |                 asm (full run) |      flat (full run) | flat (isolated) | Verdict                                                               |
| ---------------------------------------------- | -----------------------------: | -------------------: | --------------: | --------------------------------------------------------------------- |
| `model`                                        |                              – |                  5.2 |             5.3 | ~linear; GC 520–745 ms at 5M                                          |
| `structure`, `lookup+boundary`                 |                       2.0, 3.2 |             7.3, 4.1 |               – | linear (structure is small)                                           |
| `bonds-intra`                                  |                       (shared) |                  6.1 |             4.5 | linear isolated                                                       |
| `bonds-inter`                                  |                            5.5 | 9.3 (first pass 6.2) |               – | inconclusive on flat; 3.6 GB transient on asm (05/P4)                 |
| `valence-model`, `asa`                         |                       5.3, 5.6 |             4.5, 5.1 |               – | linear, but 1.2–31 s                                                  |
| `interactions`                                 |                            5.8 |                  6.4 |               – | mildly; 13–31 s                                                       |
| `component:polymer-derived`                    | 4.9 (component-only rerun 3.9) |                  4.6 |               – | linear (the first pass's 8.1 was noise)                               |
| `repr:cartoon`                                 |                    (instanced) |                  8.9 |             5.4 | linear isolated; 8.9 with everything held                             |
| `repr:ball-and-stick`                          |                    (instanced) |                   31 |         **9.1** | **super-linear**: 10.6 GB retained, 5.4 s GC, most of it not resident |
| `repr:gaussian-surface`                        |                              – |      7.2 (`low` 7.4) |               – | not isolated; 1.5–5.3 GB retained                                     |
| `repr:interactions`                            |                6.6 (`low` 5.2) |                    – |               – | ~linear; 2–4.4 GB of buffers                                          |
| `theme-apply` (unshared)                       |                              – |              2.8–8.5 |         3.6–5.1 | linear isolated                                                       |
| `hover`, `select` (unshared)                   |                              – |               4.4–15 |         3.6–6.3 | ~linear per visual: every visual is visited (07/P3)                   |
| `select-every-copy` (unshared)                 |                              – |               9.2–18 |     **8.3–9.0** | **super-linear**: visuals × loci elements, both ∝ copies (07/P3)      |
| `overpaint:ball-and-stick` (unshared)          |                              – |                  8.0 |              20 | runs while ball-and-stick's 10 GB is live                             |
| one-chain overpaint, instanced surface (`low`) |                            5.2 |                    – |               – | grid and splat over every instance (13/P1)                            |

The kernels stay linear. Two stages are super-linear by construction: ball-and-stick on unshared structures (its working
set passes physical memory) and marking a selection that spans many units (07/P3). In the full runs, where a session
holds the structure, its properties and components, the other unshared stages reach 6–9×. GC pauses explain only part of
that: excluding them, cartoon is still 8.0× and intra bonds 6.1×. The remainder coincides with the machine swapping
(14.8 of 16 GB). The cause is not established, but it means a 5M unshared session on a 16 GB machine is memory-bound
before it is CPU-bound.

## What breaks first

In order of the size at which each one stops a workflow (16 GB machine, `auto` quality unless noted):

1. **Overpaint or transparency on a symmetric surface**, from 186k atoms at `low` quality and finer:
   - one chain colored on a 32-copy assembly: 0.4 s;
   - 192 copies at `low`: 2.4 s, peak RSS 0.4 → 2.6 GB;
   - 864 copies at `low`: 11.7–12.6 s, transient allocation 3.3–4 GB (13/P1).

   At ≥1M atoms the `auto` → `lowest` switch hides it by turning color smoothing off.

2. **Ball-and-stick on unshared structures**:
   - 1M atoms: 3.9–6.3 s, 1 GB render buffers, 2.3 GB retained;
   - 5M atoms: 42 s isolated, 191 s in a full session; 4.5 GB render buffers, 10.6 GB retained. That is beyond what a
     browser tab or a GPU upload can hold (08/P7, 07/P7).
3. **Interactive latency at ~2.5k+ unshared units**:
   - hover 13–18 ms at 2.5k units, 65–191 ms at 11k;
   - select every copy 0.2–0.7 s;
   - one-chain overpaint 0.1–0.7 s on cartoon and ball-and-stick;
   - recolor 0.45–3.8 s (07/P3, 07/P1, 13/P4, 13/P5).
4. **Unshared load at 5M**: model 2.9–3.7 s, intra bonds 6.4–7.3 s, cartoon 5.9–11.4 s; peak RSS ≤ 3.6 GB until
   ball-and-stick (13/P3).
5. **Text mmCIF past ~4.8M atoms** fails outright (13/B1). Text costs ~111 bytes per atom, and V8's maximum string
   length is 2²⁹ − 24 characters (512 MB).
6. **On-demand properties** run on the main thread (06/P3, 09/P1, 13/P6):
   - interactions: 13 s on asm-864, 31 s on flat-864;
   - ASA: 28–31 s;
   - the interactions representation on assemblies: another 5.3–6.8 s and 2–4.4 GB of buffers.
7. **Molecular surface on unshared**: 21 s at 1M (`lowest`), 46 s at `low`; not run at 5M (minutes). On assemblies it
   stays at 0.2–0.35 s for any copy count.

Not broken: everything on symmetric assemblies up to 5M atoms except the above.

- Structure and assembly 16 ms; lookup 0.13 s; intra bonds 17 ms (shared).
- Cartoon, ball-and-stick, spacefill and surfaces 20–341 ms.
- Marks ≤ 2 ms per chain; peak RSS 2.5 GB including props.
- Theme evaluation 25–55 ns per element (121–274 ms at 5M).

## Hotspots at scale

From the first-pass profiles (load 6–14). Profiled times are shorter than run times, since the profile run did not hold
the properties.

| Stage [fixture]                                   | Function (file:line)                                                                                                                             |         Share | Note                                                                                                                                                                                                                                         |
| ------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ | ------------: | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| repr:ball-and-stick [flat-864, 110 s profiled]    | `ChunkedArray.allocateNext` [chunked-array.ts:38](../../packages/core/src/data/util/chunked-array.ts#L38)                                        |     22 % self | (garbage collector) 24 %, `copyVertexBlock` `gl/renderable/merged.ts:521` 8 %; `createLinkCylinderImpostors` 42 % total, merge 13 % total. `allocateNext` creates fresh typed-array chunks; that it is slow from page faults is a hypothesis |
| repr:ball-and-stick [flat-192, 7 s]               | `allocateNext` 11 %, `copyVertexBlock` 7 %                                                                                                       |               | `createIntraUnitBondCylinderImpostors` 48 % total, `mergeRenderObjectsByType` 16 % total                                                                                                                                                     |
| repr:cartoon [flat-864, 11 s]                     | (garbage collector) 20 %, `Vec3.dot`/`sub` 11 %                                                                                                  |               | `createPolymerTraceMesh` 29 % total, `createMergedValues` (`rebuild` `merged.ts:635`) 19 % total                                                                                                                                             |
| model [flat-864, 5 s]                             | (garbage collector) 18 %, `arrayView` [column.ts:485](../../packages/core/src/data/db/column.ts#L485) 16 %                                       |               | `addMissingResiduesFromEntityPolySeq` 13 %, BCIF `stringArray`/`runLength` 13 %, `sortAtomSite` 22 % total (04/P5, 04/P6)                                                                                                                    |
| bonds-intra [flat-192, 1.7 s]                     | `GridLookup3D.query` [grid.ts:315](../../packages/core/src/math/geometry/lookup3d/grid.ts#L315)                                                  |     49 % self | `findBonds` 18 %; 57 % in `find` total (02/P3, 05/P7)                                                                                                                                                                                        |
| bonds-inter [asm-864, 3 s, 3.6 GB alloc]          | `GridLookup3D.find` `grid.ts:45` 23 %, `query` `grid.ts:315` 22 %                                                                                |               | `transformMat4` 12 %, `findPairBonds` 11 % (05/P4)                                                                                                                                                                                           |
| component:polymer-derived [asm-864, 2.8 s]        | same as bonds-inter                                                                                                                              |               | `findBonds2` 96 % total: the child's inter-unit bonds from scratch (05/P5)                                                                                                                                                                   |
| interactions [asm-864, 15 s]                      | `query` `grid.ts:315` 31 %, `addStructureContacts` [contacts.ts:197](../../packages/model/src/props/computed/interactions/contacts.ts#L197) 22 % |               | `findInterUnitContacts` 93 % total via `eachUnitPair` (06/P3, 06/N1)                                                                                                                                                                         |
| asa [asm-864, 42 s]                               | `computeRange` [area.ts:74](../../packages/model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts#L74) 42 %                      |               | `query` 37 %, structure lookup `add` 4 % (06/P7)                                                                                                                                                                                             |
| valence-model [asm-864, 2.2 s]                    | `getEdgeIndices` [inter-unit-graph.ts:57](../../packages/core/src/math/graph/inter-unit-graph.ts#L57)                                            |          34 % | 02/P9                                                                                                                                                                                                                                        |
| overpaint:gaussian-surface [asm-192 `low`, 3.6 s] | `calcMeshColorSmoothing` [color-smoothing.ts:151](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L151)                         | **84 %** self | 95 % total, all 9 visuals, every instance (13/P1)                                                                                                                                                                                            |
| theme-apply:cartoon [flat-864, 1.4 s]             | `forEachCell` → `add` [render-object.ts:191](../../packages/graphics/src/gl/render-object.ts#L191) 18 %                                          |               | `hasForeignResources` 35 % total, `prepareUpdate` 44 % total (13/P4)                                                                                                                                                                         |
| theme-apply:ball-and-stick [flat-864, 2.9 s]      | `createGroupColor` [color-data.ts:268](../../packages/graphics/src/geo/geometry/color-data.ts#L268) 53 % total                                   |               | `LocationIterator.move` 24 %, `bond.ts:174` `getLocation` 9 %                                                                                                                                                                                |
| hover:cartoon [flat-864, 48 ms]                   | `mark` [units-visual.ts:525](../../packages/graphics/src/repr/structure/units-visual.ts#L525)                                                    |     39 % self | `visuals.forEach` 88 % total: per-visual `invariantId` scan (07/P3)                                                                                                                                                                          |
| select-every-copy:cartoon [flat-864, 186 ms]      | `mark` units-visual.ts:525 56 %                                                                                                                  |               | `lociApply` 37 %: 11k visuals × 864 loci elements (07/P3, 07/P5)                                                                                                                                                                             |
| overpaint:cartoon [flat-864, 190 ms]              | `ValueCell.set` 21 %, `setOverpaint` [visual.ts:260](../../packages/graphics/src/repr/visual.ts#L260) 14 %                                       |               | 11k visuals each allocate and clear an overpaint texture (13/P5)                                                                                                                                                                             |

Sampled JS-heap allocation (`--alloc`, flat-32):

- **Cartoon:** 373 MB for 186k atoms (`addTube` 11 %, `getBoundary` 6 %, `interpolateNormals`, `slerp`, `normalize`,
  `orthogonalize`).
- **Ball-and-stick:** 537 MB (`createLinkCylinderImpostors` 17 % self, the per-bond `position` callback 7 %,
  `calculateShiftDir` 6 %, `Sphere3D` 5 %).
- **Model:** 99 MB (`arrayView` 33 %, `stringArray` 22 %, `addMissingResiduesFromEntityPolySeq` 21 %).

The geometry loops allocate nothing visible per element. Plausibly most of it is boxed doubles returned from callbacks
V8 does not inline (`radius`, `position`, `size`), plus temporaries in the trace iterator. Not verified further.

## Findings

New ids are `13/…`; existing ids are re-measured at scale. Estimates are from the measurements above unless marked.

### 13/P1 Instanced surfaces: overpaint and transparency smoothing span every instance (high)

- **Mechanism.** On instanced geometry the overpaint / transparency type is `groupInstance`, so `calcMeshColorSmoothing`
  ([color-smoothing.ts:151](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L151)) works across every
  instance:
  - it builds one grid over the bounding box of _all_ instances (`getVertexBox` transforms the vertex box by every
    instance, :129);
  - it splats `vertexCount / stride` samples per instance (:193);
  - it does so for every visual, touched or not (07/P1).

  The grid grows with the assembly's volume and is capped only by `maxTextureSize`.

- **Evidence.** One chain overpainted on the gaussian surface:
  - asm-32: 0.4 s;
  - asm-192 at `low`: 2.4 s, 914 MB allocated, peak RSS 0.4 → 2.6 GB;
  - asm-864 at `low`: 12.6 s, 4 GB allocated; the molecular surface 11.7 s.

  Profile: `calcMeshColorSmoothing` 84 % self, 95 % total.

- **Fix.** Combine with 07/P1: skip visuals the layers do not touch. Within a touched visual, splat only instances that
  have a non-zero layer value (one chain → one instance), and bound the grid to their box. Untouched instances must then
  read zero outside the grid: check the shader's clamp, or keep the full-box grid but splat only touched instances (time
  win only). Smoothing in the invariant frame once per touched instance would need a shader layout change: a design
  question for the WebGPU owners.
- **Estimate.** One chain on asm-864 goes from ~12 s to milliseconds (one instance of one visual). Effort S (skip
  untouched instances) or M (bounded grid). Not frozen (`geo/`), but the layout reaches the shaders. Risk: the
  merged-object define rule of 07/P1.

### 13/P2 `quality: 'auto'` is set by total atoms, not by geometry per group (question)

- `getStructureQuality` scores `structure.root.elementCount`
  ([repr/util.ts:148](../../packages/graphics/src/repr/util.ts#L148)). An 864-copy assembly of 5.8k atoms gets `lowest`
  (8 Å surfaces, 1 linear segment), although each instanced visual holds 5.8k atoms; at `low` the same assembly builds
  cartoon in 51 ms and surfaces in 67 ms.
- Conversely, `lowest` does not bound the unshared cases that hurt: impostors ignore it, and cartoon at `lowest` still
  holds 1.7 GB at 5M.
- **Question.** Score instanced representations by the largest symmetry group's element count (or atoms ÷ instances)? It
  changes the look of every large assembly, and it re-exposes 13/P1, so fix 13/P1 first.

### 13/P3 Unshared 5M: allocation and live heap (medium, broad)

- **Observation.**
  - Allocation per stage is linear: model 2.7 GB, cartoon 6 GB, ball-and-stick 11 GB of JS heap at 5M.
  - Isolated, model and cartoon scale 5.3–5.4× for 4.5× the atoms, with 0.6–0.75 s of GC pauses.
  - With the session's properties held they reach 5.2–8.9×, and GC pauses rise (cartoon 2.3 s).
  - GC explains only part of the excess. The rest coincides with swapping and is not attributed.
- **Fix direction.** Reduce the per-element allocators that `--alloc` shows:
  - `arrayView` copies of all columns (04/P5) and template-string keys (04/P6);
  - polymer-trace temporaries, and possibly boxed doubles from per-bond callbacks (07/P9, the 02/P4 pattern);
  - builder over-allocation (08/P3, 08/P8).

  Each is S; together they target a 1–2 s GC share at 5M.

- **Measure.** `--alloc` per stage and the GC column on an idle machine; the aim is less allocation, not a faster
  kernel.

### 13/P4 Recolor pays a per-visual resource check (medium)

- `UnitsVisual.prepareUpdate` calls `hasForeignResources(renderObject, graphics)` on every update
  ([units-visual.ts:216](../../packages/graphics/src/repr/structure/units-visual.ts#L216)).
- That walks every value cell twice (`forEachCell`,
  [render-object.ts:190](../../packages/graphics/src/gl/render-object.ts#L190)).
- It is 35 % of a cartoon recolor at 11k visuals (474 of 1,375 ms profiled), although the context changes only on a
  backend switch.
- **Fix.** Check only when `graphics` differs from the context the visual last checked against (keep it in the visual's
  closure). The call site is in `repr/` (not frozen); `render-object.ts` stays as is. The same pattern is in
  `complex-visual.ts:175`, `volume/visual.ts:198` and `particles/visual.ts:197`. Effort S.

### 13/P5 Overpaint and transparency touch every visual (medium; extends 07/P1)

- Without smoothing, `setOverpaint` still allocates or clears the overpaint texture of every visual and runs `lociApply`
  on each (`visual.ts:260`). One chain costs 106–152 ms on cartoon and 255–665 ms on ball-and-stick at 11k visuals.
- **Fix.** Apply the 07/P1 rule (skip visuals the layers do not touch) before `createOverpaint`. Use 07/P3's per-loci
  unit index to find the touched visuals without scanning. Keep defines uniform within merged groups (07/P1 risk).
  Effort S on top of 07/P1/P3.

### 13/P6 ASA ignores symmetry (low priority, on demand)

- `asa` costs about the same on asm and flat (31.4 vs 28.4 s at 5M). Shrake-Rupley runs over every atom of the structure
  through the structure lookup (06/P7), so the 864 copies of an assembly pay 864×.
- **Question.** Symmetry-equivalent units differ only at their interfaces. Computing per group with the neighbours of
  the group's first unit is exact for unit interiors but not at the copies' edges, so this is a spec decision (as for
  06/P3).
- Until then ASA is computed on demand only, which is already the case: the plugin behavior has `autoAttach: false`, and
  the ASA color theme, its tooltip and MolScript queries attach it.

### 13/B1 Text mmCIF stops at ~4.8M atoms (limit)

- `CIF.parse` takes a string. Text costs ~111 bytes per atom (`1tup-flat-192.cif`: 124 MB for 1.12M atoms), so a 5M-atom
  file is ~560 MB. That is past V8's maximum string length (536,870,888 characters,
  `buffer.constants.MAX_STRING_LENGTH`), and reading it as text fails before parsing. Chromium has the same V8 limit;
  Firefox and Safari have their own, not checked.
- **Fix.** A tokenizer over `Uint8Array` (ASCII) for text CIF, effort L; or document BCIF for >4M atoms.

### Existing ids at scale

| Id                     | At scale                                                                                                                                  | Re-ranked                        |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| 08/P7 (frozen) + 07/P7 | ball-and-stick buffers 0.99 GB at 1M, 4.46 GB at 5M (`aMapping`/`aStart`/`aEnd` 850 MB each); retained 2.3 / 10.6 GB with the merged copy | **top memory item for unshared** |
| 07/P3, 07/P5           | select every copy 0.2–0.7 s at 11k visuals, super-linear (8.3–9×); hover 65–191 ms                                                        | medium → high for unshared       |
| 07/P1                  | one-chain overpaint on unshared `low` surfaces 1.7 s (1M) and 11.7 s (5M)                                                                 | high                             |
| 06/P3, 06/N1           | interactions 13 s on asm-864: every one of the symmetry-equivalent unit pairs                                                             | high for assemblies              |
| 09/P1                  | interactions representation 2–4.4 GB of mesh buffers, 5.3–6.8 s on asm-864                                                                | high for assemblies              |
| 05/P5                  | polymer component on asm-864 1.6–2 s, 3.8 GB allocated (linear)                                                                           | medium (components are routine)  |
| 05/P4                  | inter-unit bonds 2.1 s and 3.6 GB transient on asm-864                                                                                    | medium → high                    |
| 05/P3                  | every-copy bundle component derived data 131 ms on asm-864                                                                                | medium                           |
| 02/P3, 05/P7           | intra bonds 6.4–7.3 s at 5M unshared, `query` ~50 %                                                                                       | medium                           |
| 02/P9                  | valence model 1.2 s on asm-864, `getEdgeIndices` 34 %                                                                                     | low                              |
| 05/P2                  | lookup+boundary 1.4 s at 5M unshared                                                                                                      | low                              |
| 04/P5, 04/P6           | model 2.9–3.7 s at 5M unshared (13/P3)                                                                                                    | medium                           |
| 08/P3, 08/P8           | builder over-allocation feeds 13/P3                                                                                                       | medium                           |

## Top wins at large scale, ranked

Rank = impact at 1–10M atoms × how common the workflow is ÷ effort. "asm" = symmetric assemblies, "flat" = unshared.
Estimates marked ~ are extrapolated, not measured.

| #   | Ids                                | What                                                                                | Case                    | Estimate                                                                                                                        | Effort         | Frozen            |
| --- | ---------------------------------- | ----------------------------------------------------------------------------------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------------------- | -------------- | ----------------- |
| 1   | 13/P1 + 07/P1                      | Smooth only touched visuals and instances for overpaint / transparency              | asm (and flat surfaces) | one chain on 864 copies: 12 s and 4 GB allocated → ~ms; flat 1M 1.7 s → ~visual-local                                           | S–M            | no (shader check) |
| 2   | 08/P7 + 07/P7                      | Cylinder impostor data once per cylinder; release or skip the merged CPU copy       | flat                    | ball-and-stick 5M: buffers 4.46 → ~1 GB, retained 10.6 → ~2–5 GB, back under physical memory (42–191 s → ~20 s, linear from 1M) | L / S (policy) | **yes** (gl) / no |
| 3   | 07/P3 + 07/P5 + 13/P5              | Index loci by unit once per mark; single-pass `Visual.mark`; skip untouched visuals | flat                    | select every copy 0.2–0.7 s → ~tens of ms; hover 65–191 → ~few ms; one-chain overpaint 0.1–0.7 s → ~ms at 11k units             | M              | no                |
| 4   | 06/P3 + 06/N1                      | Interactions per symmetry-equivalent unit pair class                                | asm                     | 13 s on asm-864 → ~10× less (06/P3: 441 classes for 4,282 pairs at 32 copies); spec question                                    | M              | no                |
| 5   | 09/P1                              | Interactions as cylinder impostors                                                  | asm                     | buffers ~5× smaller (09/P1: 422 → ~80 MB at 186k), so 2–4.4 GB → ~0.4–0.9 GB; the build ~7× faster                              | M              | no                |
| 6   | 05/P5 + 05/P4                      | Child inter-unit bonds from the parent; `findPairBonds` via A's grid                | asm (components)        | polymer component 1.6–2 s → < 0.1 s; root inter bonds 2.1 s / 3.6 GB → ~1 s / < 0.1 GB                                          | M              | no                |
| 7   | 13/P4                              | Skip the foreign-resource walk unless the context changed                           | flat                    | recolor −35 % at 11k visuals (cartoon 0.45–0.95 s)                                                                              | S              | no                |
| 8   | 13/P3 (04/P5, 04/P6, 08/P3, 08/P8) | Cut per-element JS allocation in model and geometry builders                        | flat                    | GC pauses −0.5–2 s per stage at 5M; less memory pressure                                                                        | S each         | no                |
| 9   | 02/P3 + 05/P7                      | Finer grid cells and 3.1 Å query for intra bonds                                    | flat                    | query −15–40 % (02/P3) and bonds −15–25 % (05/P7): intra bonds 6.4–7.3 s → ~5 s at 5M                                           | S–M            | no                |
| 10  | 13/B1                              | Byte-level text CIF tokenizer                                                       | flat (text)             | > 4.8M atoms loadable as text                                                                                                   | L              | no                |
| 11  | 13/P2                              | Quality by atoms per group for instanced representations                            | asm                     | look, not time; do after #1                                                                                                     | S              | no                |
| 12  | 13/P6, 06/P7                       | ASA per group (spec) and world grid                                                 | both                    | 28–31 s → ~1/copies on asm                                                                                                      | M–L            | no                |

Frozen items go to the WebGPU owners with these numbers: the #2 impostor layout and member release in
`gl/renderable/merged.ts`. The non-frozen half of #2 is the merge policy (07/P7): do not merge large members, or merge
without keeping member arrays. The bench measures it as retained ≈ buffers instead of ≈ 2.4× buffers.

## Caveats

- One run per 1M/5M stage on a loaded machine: compare the isolated ratios and the profiles, not single numbers. Rerun
  the rows a fix targets with `--only` on an idle machine before and after.
- Headless CPU paths only: no GPU upload, frame time or GPU memory. GPU surfaces (`tryUseGpu`) are not measured; the CPU
  gaussian surface is what headless and fallback paths pay.
- 1tup copies are crystal packing: dense, so inter-unit work (bonds, contacts) is near its worst case per unit pair. A
  capsid shell has fewer neighbours per unit; a cellPACK model has more distinct entities (more groups even when
  instanced).
- `flat-*` written by `to_mmCIF` keeps one entity set and repeats chain ids per copy; real large entries also carry
  thousands of `struct_conn` and `chem_comp` rows, not measured here.
- The source fixture `1tup.bcif` is read from Prism's test fixtures by path (`tests/src/pipeline/common.ts`); moving it
  to `data/` would decouple the bench.
