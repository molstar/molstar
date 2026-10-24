# Profile v3: where pipeline time and memory go now

Third run of the pipeline profile ([v2](../perf-audit-2/00-profile.md), [v1](../perf-audit/00-profile.md)), after the
PRs #2024–#2048 (`git log 0944823d4..16277d822`). Read-only: nothing under `packages/` was changed. Raw JSON, logs,
profiles and scripts are in the session scratchpad (`pa3/00/`), not in the repo.

## Machine and conditions

| Item | Value |
| --- | --- |
| CPU / RAM | Apple M4, 10 cores, 16 GB |
| Node | **v26.12.1** (v2: v25.6.1; `total_allocated_bytes` available) |
| HEAD | `16277d822` (= master, 2026-10-25) |
| Base | v2 HEAD `0944823d4` (`git archive` of `packages/`, same `bench.ts`) |
| Load | 00:28 start 7.0; **31–36 from 00:31 to 00:45** (other auditors); 4–12 from 00:50 on; 5–9 for the reruns, profiles, `many`, playback |
| Swap | 16.0–17.2 GB of 16.4–18.4 GB in use during the first session |

**Read the ms columns as ranges, not values.**

- The first session ran 7qpd, 3k7f, 8zzh and 1tup-2x2x2 at load 31–36. Its times are 3–5× v2's (7qpd gaussian surface
  2.0 s, v2 0.38 s), and runs of one stage swing up to 6×. The tables below use a **rerun of those four fixtures at
  load 5–12** (`*-r2`), in which `now` ran **before** `base`. The first session ran base first, so the two orders
  cancel: both show the same direction on every bold row.
- The large fixtures (flat-32 … flat-192) ran base first at falling load (base at 11 → 4, now at 6 → 4). That order
  favours `now` by maybe 10–30 %. Trust their alloc and retained columns.
- **Allocation is not load-proof either (L1).** A JIT-sensitive kernel allocates boxed doubles until it is optimized,
  and under load optimization comes late. 7qpd `repr:molecular-surface` allocated 115–478 MB per run at load 31, and
  146–181 MB at load 6. Medians of 5–6 runs at load ≤ 12 reproduce within ±15 %. Single high-load runs do not.

## Method

- **Bench**: `tests/src/pipeline/{bench,many,playback}.ts` of HEAD, bundled as `scripts/pipeline-bench.mjs` does
  (esbuild, `molstar-src` condition, `NODE_ENV=production`). The trees are `pa3/00/base` (v2 HEAD) and `pa3/00/now`
  (HEAD); `now` is a `git archive` too, so both are bundled the same way.
- **Base check**: the bundle's metafile resolves 572 of 572 Mol\* source inputs from the archive and 0 from the repo
  (`build.mjs` prints it per bundle). `nodePaths` holds third-party packages only. `@molstar/*` resolves through the
  archive's own `node_modules` links.
- **API adaptation**: none needed. `Structure.transform` exists at the base. The `transform:*` stages (added by #2025)
  run on both trees.
- **8zzh (12/B1)**: still crashes the repo bench (`bench.ts:324`, `first` undefined without an atomic polymer unit).
  Both trees ran a scratch copy that falls back to `s.units[0]`.
- **Runs**: median of 5 (7qpd, 3k7f, 8zzh, 1tup-2x2x2), 3 (flat-32), 2 (asm-32, asm-192) and 1 (flat-192). flat-192 skipped
  `repr:ball-and-stick`, `repr:molecular-surface` and `repr:interactions` (swap was full).
- **Columns**: **ms** = median wall time. **alloc** = V8 `total_allocated_bytes` delta (JS heap; typed-array stores not
  counted). **retained** = heap + array buffers after GC, minus before. **buf** = typed-array MB of the render objects.
  Bold: beyond ±20 % (floor 3 ms / 3 MB). ÷n / ×n for factors of 2 and more.

## Fixtures

As in v2: `7qpd` (11,521 atoms), `3k7f` (13,904, qFit alt locs, path), `8zzh` (3,528 coarse, patched bench),
`1tup-2x2x2` (186,496, 288 units in 9 groups), `flat-32` (186,496, 416 groups), `asm-32` / `asm-192` (186k / 1.12M, 9
groups), `flat-192` (1.12M, 2,496 groups), `many` (500 entries, 1.2k–11.5k atoms each).

## Headlines (base `0944823d4` → HEAD; ms at load ≤ 12)

| Stage | Fixture | ms | MB alloc | Why |
| --- | --- | ---: | ---: | --- |
| `interactions` | flat-32 / asm-192 / flat-192 | 906 → 364 / 2489 → 1065 / 5473 → 2091 | 579 → 194 / 1478 → 448 / 3023 → 893 | #2025 (06/P4 + P10), #2047 (06/P12) |
| `repr:molecular-surface` | flat-32 | 8.7 → 6.2 s | **1509 → 204** | #2025 (02/P10, the Vec3 contagion) |
| `repr:molecular-surface` | 7qpd / 8zzh | 1747 → 1337 / 1829 → 632 | 181 → 173 / 199 → 211 | #2033 (02/P11) |
| `repr:molecular-surface` | 3k7f | 2326 → 2089 | **114 → 223 (×2)** | **R4** |
| `transform:bonds` (rigid move) | asm-192 / flat-192 | 295 → 0.7 / 388 → 2.3 | 157 → 0.3 | #2025 (05/P14) |
| `transform:polymer-bonds` | asm-192 | 234 → 12 | 124 → 8 | #2025 (05/P14) |
| `model` | flat-192 / flat-32 | 563 → 384 / 135 → 76 | 520 → 152 / 94 → 32; retained 343 → 133 | #2036 (lazy BCIF) |
| `parse` (text CIF) | 7qpd / 8zzh | 23 → 9.8 / 9.1 → 3.8 | ≈ | #2036 (byte tokenizer) |
| `lookup+boundary` | asm-192 / flat-192 | 49 → 15 / 376 → 216 | 81 → 21 / 240 → 73 | #2025 (02/P10) |
| `bonds-inter` | asm-192 / flat-192 | 288 → 168 / 272 → 144 | 107 → 44 / 74 → 39 | #2045 (05/P16, P17) |
| `asa` | asm-192 | 6.0 → 6.1 s | 613 → 121 | **still the largest stage**: all 192 copies computed (13/P6) |
| `repr:interactions` | flat-32 / asm-192 | 259 → 96 / 386 → 243 | 145 → 89 / 274 → 130 | #2025 (07/P11) |
| `repr:ball-and-stick` | flat-32 | 577 → 398 | 399 → 190 | #2033 (08/P16), #2025 |
| `repr:spacefill` | flat-192 | 1098 → 718 | 700 → 363 | #2034 (08/P11) |
| `repr:cartoon` | 7qpd | 79 → 48 | **79 → 79** | R1 still open (alloc ×2 of v1) |
| `repr:gaussian-surface` | 7qpd / 1tup-2x2x2 | 658 → 517 / 58 → 60 | 69 → 71 / **29 → 43** | R2 still open; **R5** |
| `theme-apply:gaussian-surface` | flat-32 | 365 → 283 | **123 → 11** | #2025 |
| `represet` (many, N = 500) | – | 4147 → 235 | 3095 → 173 | #2025 (03/P1) |
| `repr:auto` (many, N = 500) | – | 12.8 → 9.3 s | 4478 → 2917 | per structure 25.7 → 18.5 ms |
| startup (import PluginContext + DefaultPluginSpec) | – | 290–370 → 128–176 | heap 35.4 → 13.5 MB | #2026 (12/P1) |

## Stage tables (base → now)

ms:

| Stage | 7qpd | 3k7f | 8zzh | 1tup-2x2x2 |
| --- | ---: | ---: | ---: | ---: |
| `parse` | **23 → 9.8 ÷2** | **7.7 → 4.4 −43%** | **9.1 → 3.8 ÷2** | **13 → 8.0 −39%** |
| `model` | **15 → 7.8 −48%** | 15 → 16 | **16 → 9.9 −39%** | **10 → 6.5 −35%** |
| `lookup+boundary` | 4.5 → 3.4 | 4.6 → 4.9 | 3.7 → 2.4 | **15 → 5.2 ÷3** |
| `bonds-intra` | **28 → 18 −35%** | 34 → 34 | 0.0 → 0.0 | **17 → 11 −34%** |
| `bonds-inter` | 5.5 → 3.8 | 2.6 → 3.2 | 0.7 → 0.6 | **64 → 39 −39%** |
| `valence-model` | **11 → 7.8 −32%** | 18 → 16 | 0.3 → 0.2 | **103 → 82 −21%** |
| `interactions` | **63 → 30 ÷2** | **188 → 87 ÷2** | 1.9 → 1.4 | **724 → 243 ÷3** |
| `asa` | **89 → 59 −34%** | 153 → 127 | **throws → 0.8** | **1900 → 1345 −29%** |
| `component:polymer-derived` | 0.7 → 0.6 | 0.7 → 0.6 | 0.9 → 0.6 | **8.2 → 3.2 ÷3** |
| `transform:bonds` | **6.2 → 0.2 ÷31** | **3.4 → 0.2 ÷17** | 0.2 → 0.1 | **66 → 0.1 ÷658** |
| `transform:polymer-bonds` | **5.2 → 0.3 ÷17** | 2.8 → 0.3 | 0.1 → 0.1 | **54 → 1.7 ÷32** |
| `repr:cartoon` | **79 → 48 −39%** | **93 → 42 ÷2** | **60 → 39 −35%** | 23 → 20 |
| `repr:ball-and-stick` | **58 → 39 −32%** | **81 → 55 −33%** | 11 → 9.4 | **57 → 33 −42%** |
| `repr:spacefill` | **15 → 11 −27%** | **22 → 16 −27%** | **11 → 6.5 −38%** | 12 → 8.9 |
| `repr:gaussian-surface` | **658 → 517 −22%** | 616 → 566 | **1024 → 649 −37%** | 58 → 60 |
| `theme-apply:gaussian-surface:element-symbol` | **160 → 123 −23%** | 134 → 127 | **93 → 73 −22%** | 24 → 20 |
| `repr:molecular-surface` | **1747 → 1337 −23%** | 2326 → 2089 | **1829 → 632 ÷3** | 429 → 378 |
| `theme-apply:molecular-surface:residue-name` | **111 → 88 −20%** | 95 → 110 | **149 → 68 ÷2** | **22 → 16 −26%** |
| `repr:interactions` | **17 → 11 −36%** | **39 → 28 −27%** | **12 → 6.2 −49%** | **140 → 49 ÷3** |

| Stage | flat-32 | asm-32 | asm-192 | flat-192 |
| --- | ---: | ---: | ---: | ---: |
| `parse` | 8.7 → 6.0 | 7.5 → 5.0 | **13 → 8.1 −36%** | 23 → 26 |
| `model` | **135 → 76 −44%** | 11 → 10 | **15 → 12 −21%** | **563 → 384 −32%** |
| `lookup+boundary` | **83 → 36 ÷2** | **11 → 6.4 −42%** | **49 → 15 ÷3** | **376 → 216 −43%** |
| `bonds-intra` | 176 → 144 | 12 → 11 | **19 → 30 +59%** ¹ | 1022 → 823 |
| `bonds-inter` | **45 → 30 −34%** | **38 → 26 −31%** | **288 → 168 −42%** | **272 → 144 −47%** |
| `valence-model` | **78 → 56 −29%** | 63 → 55 | 332 → 312 | **301 → 222 −26%** |
| `interactions` | **906 → 364 ÷2** | **448 → 178 ÷3** | **2489 → 1065 ÷2** | **5473 → 2091 ÷3** |
| `asa` | 932 → 906 | 996 → 954 | 5967 → 6122 | **7412 → 5750 −22%** |
| `component:polymer` | **5.7 → 0.8 ÷7** | **4.4 → 0.7 ÷6** | **20 → 4.0 ÷5** | **28 → 6.1 ÷5** |
| `component:polymer-derived` | 1.5 → 1.1 | **5.8 → 1.9 ÷3** | **35 → 15 ÷2** | 30 → 25 |
| `transform:bonds` | **55 → 0.1 ÷545** | **47 → 0.1 ÷471** | **295 → 0.7 ÷421** | **388 → 2.3 ÷169** |
| `transform:polymer-bonds` | **44 → 0.2 ÷219** | **38 → 1.2 ÷32** | **234 → 12 ÷20** | **275 → 1.2 ÷229** |
| `repr:cartoon` | **197 → 130 −34%** | **21 → 16 −26%** | 24 → 24 | 1026 → 850 |
| `theme-apply:cartoon:secondary-structure` | **17 → 11 −34%** | 2.6 → 1.4 | **7.1 → 2.6 ÷3** | **120 → 83 −31%** |
| `repr:ball-and-stick` | **577 → 398 −31%** | **41 → 23 −42%** | 37 → 38 | – |
| `repr:spacefill` | **169 → 113 −33%** | 8.4 → 6.7 | **15 → 12 −22%** | **1098 → 718 −35%** |
| `repr:gaussian-surface` | 858 → 700 | 39 → 37 | 18 → 17 | 892 → 940 |
| `theme-apply:gaussian-surface:element-symbol` | **365 → 283 −23%** | 19 → 17 | 3.6 → 3.4 | 177 → 164 |
| `repr:molecular-surface` | **8716 → 6224 −29%** | 302 → 317 | 173 → 169 | – |
| `theme-apply:molecular-surface:residue-name` | **313 → 219 −30%** | 16 → 17 | 2.7 → 2.8 | – |
| `repr:interactions` | **259 → 96 ÷3** | **88 → 58 −34%** | **386 → 243 −37%** | – |

¹ +11 ms on 2 runs, allocation equal (14 → 14 MB): read as noise.

MB allocated:

| Stage | 7qpd | 3k7f | 8zzh | 1tup-2x2x2 |
| --- | ---: | ---: | ---: | ---: |
| `parse` | 1.5 → 2.1 | 1.9 → 0.9 | 1.0 → 1.1 | **6.7 → 3.7 −45%** |
| `lookup+boundary` | 1.0 → 1.8 | 1.5 → 1.0 | 1.8 → 1.9 | **14 → 5.4 ÷3** |
| `bonds-inter` | 2.0 → 2.0 | 0.9 → 1.1 | 0.5 → 0.5 | **16 → 11 −32%** |
| `interactions` | **36 → 19 −46%** | **98 → 39 ÷2** | 0.7 → 0.7 | **297 → 108 ÷3** |
| `asa` | 6.8 → 6.1 | 7.1 → 6.9 | throws → 0.6 | **101 → 41 ÷2** |
| `transform:bonds` | 1.8 → 0.5 | 1.2 → 0.5 | 0.5 → 0.5 | **26 → 0.4 ÷65** |
| `repr:cartoon` | 79 → 79 | 66 → 66 | **54 → 67 +24%** | 23 → 22 |
| `repr:ball-and-stick` | 18 → 17 | 19 → 20 | 2.4 → 2.8 | **24 → 17 −31%** |
| `repr:gaussian-surface` | 69 → 71 | 65 → 67 | 81 → 79 | **29 → 43 +47%** |
| `repr:molecular-surface` | 181 → 173 | **114 → 223 +95%** | 199 → 211 | 133 → 140 |
| `repr:interactions` | 5.1 → 4.1 | 9.0 → 8.8 | 1.6 → 1.7 | **52 → 27 −49%** |

| Stage | flat-32 | asm-32 | asm-192 | flat-192 |
| --- | ---: | ---: | ---: | ---: |
| `model` | **94 → 32 ÷3** | 17 → 15 | 16 → 14 | **520 → 152 ÷3** |
| `lookup+boundary` | **55 → 11 ÷5** | **25 → 13 −49%** | **81 → 21 ÷4** | **240 → 73 ÷3** |
| `bonds-intra` | 41 → 37 | **15 → 11 −30%** | 14 → 14 | 199 → 192 |
| `bonds-inter` | 9.4 → 9.6 | 23 → 19 | **107 → 44 ÷2** | **74 → 39 −47%** |
| `valence-model` | 92 → 87 | 92 → 91 | 488 → 436 | **498 → 346 −31%** |
| `interactions` | **579 → 194 ÷3** | **333 → 121 ÷3** | **1478 → 448 ÷3** | **3023 → 893 ÷3** |
| `asa` | **16 → 12 −26%** | **128 → 45 ÷3** | **613 → 121 ÷5** | 70 → 74 |
| `component:polymer` | **5.8 → 0.7 ÷8** | **6.3 → 1.2 ÷5** | **32 → 1.9 ÷17** | **34 → 2.6 ÷13** |
| `transform:bonds` | **24 → 0.4 ÷60** | **28 → 0.4 ÷69** | **157 → 0.3 ÷524** | **152 → 0.7 ÷217** |
| `transform:polymer-bonds` | **15 → 0.4 ÷38** | **21 → 1.7 ÷13** | **124 → 8.1 ÷15** | **89 → 0.7 ÷128** |
| `repr:cartoon` | **220 → 140 −36%** | 38 → 35 | 27 → 28 | **1250 → 738 −41%** |
| `theme-apply:cartoon:secondary-structure` | 5.8 → 6.8 | 0.8 → 0.6 | 0.6 → 0.4 | **33 → 46 +38%** |
| `repr:ball-and-stick` | **399 → 190 ÷2** | **30 → 24 −22%** | 25 → 20 | – |
| `repr:spacefill` | **129 → 76 −41%** | 5.8 → 5.1 | 5.6 → 4.0 | **700 → 363 −48%** |
| `repr:gaussian-surface` | 153 → 138 | 28 → 32 | 14 → 12 | **539 → 384 −29%** |
| `theme-apply:gaussian-surface:element-symbol` | **123 → 11 ÷11** | **12 → 8.0 −31%** | 0.8 → 0.9 | 29 → 33 |
| `repr:molecular-surface` | **1509 → 204 ÷7** | 146 → 159 | **99 → 61 −39%** | – |
| `theme-apply:molecular-surface:residue-name` | **101 → 11 ÷9** | 16 → 13 | 0.9 → 0.6 | – |
| `repr:interactions` | **145 → 89 −39%** | **53 → 34 −36%** | **274 → 130 ÷2** | – |

MB retained (rows that differ; every other load, property and representation row is within ±1 MB):

| Stage | flat-32 | flat-192 |
| --- | ---: | ---: |
| `model` | **57 → 22 ÷3** | **343 → 133 ÷3** |
| `lookup+boundary` | 6.7 → 6.0 | **45 → 33 −26%** |
| `bonds-inter` | **0 → 3.0** | **0.3 → 17** (R6) |
| `interactions` | **−21 → +11** (R6) | 256 → 257 |
| `repr:spacefill` | 43 → 37 | 254 → 220 |
| `repr:gaussian-surface` | 176 → 170 | 203 → 169 |

Render buffers (MB): **identical on every fixture and representation** (cartoon 14 / 12 / 9.6 / 1.2 / 35 / 0.7 / 60;
ball-and-stick 165 on flat-32; interactions 2.4 / 9.8 / 45 / 83 / 265 on 7qpd / 3k7f / 1tup-2x2x2 / flat-32 / asm-192;
gaussian surface 30 / 73 / 12 on 7qpd / flat-32 / flat-192). No PR since v2 changed what reaches the GPU.

Interactive stages (flat-32 / asm-192 / flat-192): `hover:cartoon` 0.5 / 1.3 / 5.1 ms, `select-every-copy:cartoon` 0.4 / 1.7
/ 5.7 ms, `overpaint:cartoon` 3.4 / 0.8 / 34 ms, `overpaint:gaussian-surface` 2.3 / 0.6 / 22 ms. All equal to base
within noise; flat-192 is one run.

## Many structures and playback

`many --transaction stage`, one run each (ms / MB alloc; per structure in brackets):

| Stage | N = 100 base → now | N = 500 base → now |
| --- | ---: | ---: |
| `load` | 623 / 527 → 493 / 373 | 3422 / 2356 → 2421 / **1468** (4.7 → 2.9 MB each) |
| `bonds` | 1072 / 313 → 910 / 233 | 4582 / 1284 → 3425 / 859 |
| `repr:auto` | 1858 / 1044 → 1471 / 735 (18.6 → 14.7) | 12832 / 4478 → 9265 / **2917** (25.7 → 18.5) |
| `theme:one` | 12 → 17 | 96 → 80 |
| `hover` (per event, 21 events) | 1.2 → 1.2 | 8.7 → 8.5 |
| `select` (per event, 10 events) | 1.1 → 1.3 | 8.8 → 11.8 |
| `represet` | 219 / 183 → 72 / 45 | **4147 / 3095 → 235 / 173** |
| `remove:all` | 13 → 17 | 237 → 73 |

- `repr:auto` is still super-linear: 14.7 ms per structure at N = 100, 18.5 at N = 500. Per-event select and hover grow
  ~8× for 5× N (v2's P7; the Canvas3D marks are frozen).
- `load` allocates 38 % less (#2035, #2036, #2039: downloads, text CIF read as bytes).

Playback (defaults): every counter is equal base → now. That covers ticks, applied frames, cells updated per tick
(4.39–6.83), marks (108, 468), `setSnapshot` (31 for 3 applied frames) and the changed arrays per applied frame (14.2 MB
for the synth preset, 31 MB for the gaussian surface). The ms per applied frame are equal within noise (gaussian
surface 527 → 416 ms). **v2's P9 is fixed** (#2034): mutative's "return value does not contain any draft" warning was printed 6×
on base and 0× now.

## CPU and allocation hot spots now (HEAD; 100 µs sampling after a warm-up run; % of the stage)

| Stage [fixture, wall, alloc] | Top self (file:line) | Note |
| --- | --- | --- |
| repr:molecular-surface [flat-32, 6.2 s, 204 MB] | `projectTorus` `core/src/math/geometry/molecular-surface.ts:269` 44 %, `obscured` :132 24 %, `projectPointsRange` :181 14 %, `query` `lookup3d/grid.ts:315` 8 % | alloc: MC `interpolate` `graphics/src/geo/util/marching-cubes/algorithm.ts:189` 57 MB, `projectPointsRange` 32 MB, `projectTorus` 24 MB |
| repr:molecular-surface [3k7f, 2.1 s, 223 MB] | – | alloc: `projectPointsRange` **120 MB (base 52)**, `projectTorus` 39 (base 46), `interpolate` 36 (R4) |
| asa [asm-192, 6.1 s, 121 MB] | `computeRange` `model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts:74` 42 %, `query` grid.ts:315 41 %, `add` `lookup3d/common.ts:18` 6 % | all 192 copies computed (13/P6); alloc: `computeRange` 66 MB, `SymmetryOperator` x/y/z `core/src/math/geometry/symmetry-operator.ts:366–376` 19 MB |
| asa [flat-32, 0.9 s] | `computeRange` 47 %, `query` 39 % | |
| interactions [asm-192, 1.1 s, 448 MB] | `query` grid.ts:315 41 %, `addStructureContacts` `model/src/props/computed/interactions/contacts.ts:330` 4 % (69 % total via `eachUnitPair` `structure.ts:1333`) | alloc: `(set)` 61 MB, `refineInteractions` `refine.ts:30` 50 MB, `InteractionsInterContacts` `common.ts:102` 35 MB |
| interactions [flat-32, 0.36 s, 194 MB] | `query` 31 % | inter-unit 44 % total |
| valence-model [asm-192, 0.31 s, **436–530 MB**] | `getEdgeIndices` `core/src/math/graph/inter-unit-graph.ts:57` 16 %, `advance` `unit/bonds.ts:230` 13 %, `get bonds` `unit.ts:459` 12 % | alloc: `calculateHydrogensCharge` `valence-model.ts:104` 302 MB, `get interUnitBonds` `structure.ts:242` 210 MB self (inlined callees) |
| repr:gaussian-surface [flat-32, 0.70 s, 138 MB] | `accumulateRange` `core/src/math/geometry/gaussian-density/cpu.ts:70` 31 %, `stridedRow` algorithm.ts:342 17 %, `interpolate` :189 13 % | alloc: `interpolate` 45 MB (R2) |
| theme-apply on surfaces [flat-32, 0.28 s] | `splatVertex` `graphics/src/geo/geometry/mesh/color-smoothing.ts:332` 69 %, `calcMeshColorSmoothing` :414 20 % | 08/P1 (partly fixed in v2) |
| repr:ball-and-stick [flat-32, 0.40 s, 190 MB] | GC 8 %, `CylindersBuilder.add` `cylinders-builder.ts:108` 5 %, `LocationIterator.move` 5 %, `size4` `theme/size/physical.ts:47` 4 %, `copyVertexBlock` `gl/renderable/merged.ts:521` 4 % (frozen) | alloc: `createLinkCylinderImpostors` `repr/structure/visual/util/link.ts:361` 41 MB, `size4` 13 MB, `ValueCell.create` 12 MB |
| repr:interactions [asm-192, 0.24 s, 130 MB] | `CylindersBuilder.add` 19 %, `fillMappingAndIndices2` :270 9 %, inter-unit `position` `interactions-inter-unit-cylinder.ts:72` 9 % | alloc: `addFixedCountDashes` `cylinders-builder.ts:184` 69 MB, inter-unit `position` 42 MB |
| bonds-intra [flat-32, 0.14 s] | `query` grid.ts:315 45 %, `findBonds` `intra-compute.ts:164` 20 % | alloc: `findBonds` 28 MB |
| bonds-inter [asm-192, 0.17 s] | `query` 54 %, `findPairBonds` `inter-compute.ts:91` 13 % | |
| model [flat-32, 76 ms, 32 MB] | `fusedNumbers` `io/src/common/binary-cif/decoder.ts:368` 10 %, `addMissingResiduesFromEntityPolySeq` `basic/properties.ts:64` 6 %, `getAtomicIndex` 5 % | alloc: `(set)` 6 MB, `areBackboneConnected` 4 MB |
| repr:cartoon [3k7f, 42 ms, 66 MB] | – | alloc: `addTube` `mesh/builder/tube.ts:52` 36 MB, `interpolateNormals` `polymer/curve-segment.ts:164` 14 MB, the same as base (R1) |

### `--trace-deopt` on 7qpd (cartoon, gaussian, molecular surface; 5 runs each, base → now)

| Event | base | now | Reading |
| --- | ---: | ---: | --- |
| `accumulateRange` "exit from OSR'd inner loop" | 2,713 | 2,980 | 02/P13 (v2's P4) **still open** |
| `projectPointsRange` OSR exits | 668 | 685 | still open |
| `projectToriiRange` OSR exits | 1,627 | **0** | gone (#2033, 02/P11 restructured it) |
| "embedded weak objects cleared": MC `interpolate` / `processCell` / `stridedRow` / `slice` / `get3dOffsetFromEdgeInfo` | 10 each | 10–11 each | **R2 / 08/P15 still open**: optimized MC code dropped at every major GC, ~2 per run |
| same, `interpolateNormals` / `addTube` / `addSheet` | 4 / 4 / 2 | 4 / 3 / 2 | **R1 / 08/P12 still open** |
| all bailouts / weak-object markings | 5,302 / 185 | 3,942 / 178 | |

R1 and R2 are not fixed: no PR since v2 touched `mesh-builder.ts` or the MC state map (`algorithm.ts:117`). 7qpd cartoon
allocates 79 MB on both trees (v1 base: 38 MB).

## Startup and bundle

| Measure | base | now | Note |
| --- | ---: | ---: | --- |
| import `PluginContext` + `DefaultPluginSpec` (Node, bundled, load 6, 7 runs) | 290–370 ms | 128–176 ms | #2026 (12/P1: lazy ElementSymbol table) |
| JS heap after import (after GC) | 35.4 MB | 13.5 MB | −21.9 MB, as #2026 measured |
| `new PluginContext` / `init()` | 24–31 / 16–30 ms | 11–24 / 15–24 ms | |
| Prism start closure (`index.js` + static chunks, production) | 6.31 MB raw / 1.96 MB gzip (v2) | **5.29 MB / 1.59 MB** (90 of 136 files) | v2 after #2026: 5.16 / 1.55; **+0.13 MB raw since**; budget 5.5 MB (96 % used) |
| Prism single file (`prism-viewer.js`) | – | 7.24 MB raw / 2.22 MB gzip | budget 8 MB |

- **What is left of startup** (CPU profile, now): 101 ms of `compileSourceTextModule` (the bundle itself), 53 ms of
  module `__init` bodies, GC 21 ms. Shader source normalization regexes run at import (~13 ms;
  `graphics/src/gl/webgl/programs.ts:80`, frozen). `zip/util.ts:70` (CRC table) takes 1.5 ms.
- **Start closure by package (kB)**: graphics 1,713, prism/viewer 731, model 559, plugin/core 459, core 455,
  prism/ui 442, io 179, **query-language 131** (Jmol `properties.ts` alone 23 kB), mvs runtime + builder 188.
  `sources/scene.ts` is the largest Prism input (51 kB). Sequences: 0 (P9 holds).

## v2 → now (this slice's items)

| Id | v2 status | Now | PR | Note |
| --- | --- | --- | --- | --- |
| R1 (= 08/P12) cartoon small entries | regressed | **still open** | – | 7qpd 79 → 79 MB; weak-object deopts unchanged |
| R2 (= 08/P15) surface alloc small entries | regressed | **still open** | – | 7qpd gaussian 69 → 71 MB; MC weak deopts unchanged |
| R3 3k7f `lookup+boundary` ×2 | (not reproduced warm in v2) | **dropped** | – | 4.6 → 4.9 ms, alloc 1.5 → 1.0 MB |
| P3 `projectTorus` 7.4 GB on flat-192 | open | **fixed** | #2025 (02/P10), #2033 (02/P11) | flat-32 molsurf alloc 1509 → 204 MB |
| P4 (= 02/P13) OSR exits | open | **partly**: torii gone, density and points not | #2033 | 2,980 / 685 exits per 5 runs |
| P5 nucleotide ring mesh 14 % of flat-192 cartoon | open | not re-profiled (flat-192 cartoon 1.03 → 0.85 s, alloc −41 %) | | |
| P6 `lookup+boundary` EPOS | open | **partly** | #2025 | ÷2–4 alloc, −43 % to ÷3 time; per-unit boundary still the cost |
| P7 many: super-linear `repr:auto`, select | open | **still open** | | 14.7 → 18.5 ms per structure; frozen marks |
| P8 interactions buffers 1.2 GB asm-864 | open | **still open** (frozen schema) | | asm-192 265 MB, identical |
| P9 focus `produce` warning | open | **fixed** | #2034 (as 12/P7) | 6 → 0 warnings in playback |
| 06/B3 ASA on coarse | open | **fixed** | #2032 | 8zzh `asa` throws → 0.8 ms |
| 12/B1 bench crash on coarse | new | **still open** | | `bench.ts:324` |
| 12/B4 errors read as fast stages | new | **still open** | | `stage.ts:99`, `bench.ts:209` |
| 01/B1 = 04/B1 `memoize1` pins `BasicData` | open | **still open** (by code) | | `properties.ts:189`, `memoize.ts:30` |

Counts: fixed 4 (P3, P9, 06/B3; R3 dropped), partly 2 (P4, P6), still open 7 (R1, R2, P7, P8, 12/B1, 12/B4, 01/B1), not
rechecked 1 (P5).

## Regressions and new costs (A/B)

- **R4 (new) `repr:molecular-surface` allocation on 3k7f ×2.**
  - Isolated runs, 2 × 6 per tree: base 113–222 MB (median ~140), now 230–354 MB (median ~320).
  - Time is not worse: 2.2–3.8 → 2.1–2.9 s.
  - Allocation profile: `projectPointsRange` (`molecular-surface.ts:181`) self 52 → 120 MB; `projectTorus` 46 → 39.
  - Not seen on 7qpd (181 → 173 MB) and not on flat-32 (÷7 better).
  - 3k7f is the qFit entry: many alt locs, so dense neighbourhoods.
  - Candidates in the slice-02/07 diff: #2034 `updateChunk` (02/B2) changes how many atoms each OSR'd call covers;
    #2034 07/P6 compact per-unit `Float64Array` coordinates (`repr/structure/visual/util/molecular-surface.ts:47`);
    #2033's scalar-local `projectTorus` (02/P11).
  - Slice 02 should bisect.
- **R5 (new) `repr:gaussian-surface` on 1tup-2x2x2: +50 % allocation, +15 % time.**
  - 28–33 → 43–51 MB, and 60–70 → 68–88 ms (isolated, 2 × 6 runs).
  - Allocation profile: MC `interpolate` 13 → 22 MB, `accumulateRange` 2 → 8 MB.
  - Only the small-unit assembly is hit (7qpd, flat-32 and flat-192 are equal or better).
  - Candidate: #2034's `updateChunk` in `gaussian-density/cpu.ts:68` (02/B2). It changes the calls per unit of the OSR'd
    `accumulateRange` (P4). Slices 02 and 08.
- **R6 (new, small) retained memory after bond and interaction stages.**
  - `bonds-inter` retains 0 → 3 MB on flat-32 and 0.3 → 17 MB on flat-192.
  - `interactions` on flat-32 retains +11 MB where base gave 21 MB back.
  - Candidates: #2045's struct_conn partner cache (`inter-compute.ts:457`, a WeakMap keyed by `elements`) or the bond
    providers (#2016). Slices 05 and 06.
- **Weak**: 8zzh `repr:cartoon` alloc 54 → 67 MB (+24 %; the first session gave 60 → 65); flat-192
  `theme-apply:cartoon:secondary-structure` 33 → 46 MB (one run).

## Hot spots by slice (what each slice should look at)

| Slice | Lead | Evidence |
| --- | --- | --- |
| 01 core data/util | `ChunkedArray.create`/`_compact`, `ValueCell.create` and builtin `Set` in builders | `create2` 17 ms of flat-32 `repr:interactions`; `ValueCell.create` 12 MB in flat-32 ball-and-stick; `(set)` 61 MB of asm-192 `interactions` |
| 02 core math | **R4, R5** (bisect #2033/#2034); 02/P13 OSR exits still 2,980 / 685; `GridLookup3D.query` is 31–54 % of `interactions`, `asa`, `bonds-*`; molsurf `projectTorus` 44 % + `obscured` 24 % of the largest flat-32 stage (6.2 s); `SymmetryOperator` x/y/z box 19 MB in asm-192 ASA | tables above |
| 03 core state/task | snapshots playback: 31 `setSnapshot` for 3 applied frames (28 no-op restores); `many` select per event 1.3 → 11.8 ms (N 100 → 500) | playback, many |
| 04 model formats | `model` on flat-32: `fusedNumbers` 10 %, `addMissingResiduesFromEntityPolySeq` 6 %, `(set)` 6 MB; 01/B1 = 04/B1 still pins `BasicData` | profile |
| 05 model structure | **R6** (`bonds-inter` retained 17 MB at flat-192); `get interUnitBonds` self-attributed 210 MB (valence) + 29 MB (interactions) on asm-192: find which inlined accessor allocates; `bonds-intra` `query` 45 % | alloc profile |
| 06 model props/script | **`asa` on asm-192 6.1 s** (largest stage, all copies, 13/P6); `valence-model` 436–530 MB on asm-192 (`calculateHydrogensCharge` 302 MB); `interactions` `refineInteractions` 50 MB, `InteractionsInterContacts` 35 MB; R6's `interactions` retained +32 MB on flat-32 | profiles |
| 07 graphics repr | R1 (`interpolateNormals`, `addTube` dropped at GC); `repr:interactions` `addFixedCountDashes` 69 MB + inter-unit `position` 42 MB on asm-192; ball-and-stick `createLinkCylinderImpostors` 41 MB on flat-32 | alloc profiles |
| 08 graphics geo | R2 (MC state map in optimized code, 10–11 drops per 5 runs); R5; `splatVertex` 69 % of surface recolor (283 ms flat-32); `CylindersBuilder.add` 19 % of asm-192 interactions repr | deopt trace, profiles |
| 09 graphics theme/props | `size4` `physical.ts:47` 13 MB per flat-32 ball-and-stick build; `theme-apply:cartoon:secondary-structure` alloc +38 % on flat-192 (one run); `LocationIterator.move` 5 % of ball-and-stick | profile |
| 10 prism viewer | start closure 5.29 MB, +0.13 MB since #2026 and 96 % of its budget; query-language 131 kB and MVS 188 kB load at start | metafile |
| 11 prism ui/assistant | `prism/ui` 442 kB static (generated lucide icons 20 kB, `canvas.module.css.ts` 19 kB); `ai/src` 27 kB static | metafile |
| 12 app level | 12/B1, 12/B4 still open; L1 (alloc counters vs load) for H10's guard; startup 128–176 ms is now mostly bundle compile; `many` `repr:auto` per structure still super-linear | above |
| 13 io | text CIF parse ÷2 (#2036); BCIF `parse` on flat-192 23 → 26 ms (one run); `many` `load` alloc −38 % | tables |

## Leads no slice owns

- **L1 Allocation counters are not load-proof.**
  - JIT-sensitive kernels allocate boxed doubles until they are optimized, and under load concurrent compilation is
    late.
  - Seen here: 7qpd molsurf allocated 115–478 MB per run at load 31, and 146–181 MB at load 6. 3k7f cartoon 40–73 MB.
  - A clock-free regression guard (v2 H10) on `total_allocated_bytes` needs a median of ≥ 5 runs, a pinned
    `--single-threaded`-style flag, or both. Untested.
  - The v2 claim "allocated MB reproduce" holds only at low load.
- **L2 Node moved from 25.6 to 26.12 between v2 and v3.** Absolute numbers are not comparable with v2's report; only
  in-session A/B is.
- **L3 Profiling cost on tight loops.** `--profile` (100 µs sampling) inflates flat-32 molsurf from 6.2 to 11 s and
  asm-192 ASA from 6 to 11 s. The percentages hold; the walls in profile runs do not.
- **L4 The bench has no option for symmetry ranges on a path fixture** (v2 noted it for 1tup-cell). It also has no
  `--seed` or order control for A/B. A built-in `--base <dir>` that interleaves two bundles per run would remove the
  order bias seen here.

## Not run / not measured

- flat-192: `repr:ball-and-stick`, `repr:molecular-surface` and `repr:interactions` were skipped (swap 16–17 GB in use);
  the flat-192 rows are one run.
- asm-864, flat-864, flat-192-cif: not run (brief: too large for the shared machine).
- Pinned `--quality low` on asm-192 (v2's 13/P1 table): not rerun.
- The Classic Viewer bundle, its WebGL variant and browser timing: not built (only Prism's closure).
- P5 (nucleotide ring mesh on flat-192): not re-profiled.
- Base-side startup was measured, but base-side Prism closure was not (base needs the pre-#2027 `extensions/` layout).
  v2's 6.31 / 1.96 MB stands in.

## How to rerun

Scripts in the scratchpad `pa3/00/bench/`:

- `build.mjs <root> <bench|many|playback|startup> <out.mjs>` bundles against a tree. For a base tree it links
  `node_modules/@molstar/*` to the archive and prints how many Mol\* inputs came from where.
- `queue.sh`, `queue2.sh`: the A/B runs, profiles, deopt traces, `many`, playback.
- `compare.mjs <outdir> ms|alloc|ret|buf a b …`: the tables.
- `deopt.mjs`: the bailout counts.
- `ms3k7f.sh`, `alloc3k7f.sh`: R4 and R5.
- `startup.sh`, `startup/app.mjs` + `closure.mjs`: startup and the Prism build. This is a copy of
  `scripts/esbuild/app.mjs` with its output in the scratchpad.

Base tree:

```sh
git archive -o base.tar 0944823d4 packages scripts/workspace/inventory.json tests/tsconfig.json tsconfig.json tsconfig.base.json
tar -x -C base -f base.tar
cp -R tests/src/pipeline base/tests/src/   # HEAD's bench; 8zzh needs the `?? s.units[0]` fallback at bench.ts:324
ln -s <repo>/data base/data
ln -s <repo>/tests/build/pipeline/fixtures base/tests/build/pipeline/fixtures
node build.mjs base bench base/tests/build/pipeline/bench.mjs
```
