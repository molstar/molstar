# Profile v2: where pipeline time and memory go now

Second run of [v1's 00-profile](../perf-audit/00-profile.md), after the ~30 perf and bug-fix PRs from `60fe4c2f1` to
`0944823d4` (2026-10-22). Read-only: nothing under `packages/` was changed. Raw JSON, logs and profiles are in the
session scratchpad (`perf-audit-2/profile/`), not in the repo.

## Machine

| Item | Value |
| --- | --- |
| CPU / cores / RAM | Apple M4, 10 cores, 16 GB (`sysctl machdep.cpu.brand_string hw.ncpu hw.memsize`) |
| Node | v25.6.1 (`total_allocated_bytes` available: Alloc column valid) |
| HEAD | `0944823d4` (branch `v6-perf-audit-v2` = master) |
| Load | 13:54 before: 2.06 / 2.10 / 2.05; during: 2.7–4.7; 14:29 after: 4.44 / 3.94 / 3.58 |
| Swap | 9.2 GB used at the start, 11.1 GB from 14:22 on (other sessions; the machine is shared) |

**Load is moderate, not idle.** Absolute times carry ±10–20 % noise; one row moved +30 % and was noise on a rerun (¹).
Allocated and retained MB reproduce. Every v1 → now comparison below is therefore an **A/B in this session**, not a
comparison with v1's published numbers (next section).

## Method: what is compared with what

- **Bench**: the repo's `scripts/pipeline-bench.mjs` (bundles `tests/src/pipeline/*.ts` with esbuild, `molstar-src`
  condition, `NODE_ENV=production`).
- **v1 → now = A/B on the same machine, same session, same `bench.ts`.**
  - *base*: the v1 audit base `60fe4c2f1`. Its `packages/` were extracted with `git archive` into the scratchpad, and
    the current `bench.ts` was bundled against them.
  - *now*: HEAD.
  - The base bundle was checked to resolve only the archived sources: `compactTrimmed` (#2007) is absent, and ASA
    allocates 7–36 GB there as in v1.
- **Why not v1's published numbers.**
  - 00-profile ran an earlier scratch `bench.ts`. Stage definitions differ:
    - no `component:*`, `hover`/`select`/`overpaint`/`transparency` or `release:*` stages;
    - `theme-eval:accessible-surface-area` and `theme-apply:cartoon:accessible-surface-area` existed, now gone;
    - `theme-apply:ball-and-stick:chain-id` is new;
    - every representation stayed alive to the end of the run (now each is destroyed in `release:<repr>`);
    - plain `gc()` between stages (now two `last-resort` GCs and a 20 ms wait, so freed array buffers are not credited
      to the next stage);
    - no `NODE_ENV=production` define.
  - So 00's numbers are off by up to 2× on small stages: 7qpd `repr:spacefill` reads 4 ms in 00 against 7.5 ms on the
    base code with today's bench.
  - 13-large-scale used the repo bench (the same stages, minus the theme-eval rows `chain-id` and
    `secondary-structure`), but at load 2–14 with 10–15 GB swap.
  - v1-report → now cross-checks are given only where no A/B was run (asm-864, flat-864).
- **Not reproducible.**
  - `1tup-cell`: the bench has no symmetry-range option for a path fixture, so it was not run.
  - `8zzh`: the repo bench crashes on it (B1). It ran with a scratch copy of `bench.ts` that falls back to
    `s.units[0]` in `selections()`.
  - 3k7f ran by path.
  - flat-864 ran without `repr:ball-and-stick`, `repr:molecular-surface` and `repr:interactions`. 13 measured 10.6 GB
    retained for ball-and-stick, and swap on this shared machine was already 11 GB. Its v1 baseline is 13's
    full-session run.
- Runs: median of 5 (7qpd, 3k7f, 8zzh, 1tup-2x2x2), 3 (flat-32), 2 (asm-32, asm-192), 1 (flat-192, 864s). Columns as
  in v1: **ms** = median wall; **alloc** = V8 `total_allocated_bytes` delta (JS heap; typed-array stores not counted);
  **retained** = heap + array buffers after GC minus before. Bold: change beyond ±20 % (floor 3 ms / 3 MB).

## Fixtures

| Id | File | Atoms | Units (groups) |
| --- | --- | ---: | ---: |
| `7qpd` | `data/examples/7qpd.fw2.cif` (1.8 MB text, no `chem_comp`) | 11,521 | 7 (7) |
| `3k7f` | `packages/prism/viewer/src/_stories/dynamic-pdb/fixtures/3k7f-qfit.bcif` (path) | 13,904 | 8 (8) |
| `8zzh` | `data/examples/ihm-coarse-8zzh.cif` (coarse only; patched bench, B1) | 3,528 | 4 (4) |
| `1tup-2x2x2` | 1tup.bcif, `buildSymmetryRange` 000..111 | 186,496 | 288 (9) |
| `flat-32` (= v1 `big-bcif`) | `tests/build/pipeline/fixtures/1tup-flat-32.bcif` (3.4 MB) | 186,496 | 416 (416) |
| `asm-32` / `asm-192` / `asm-864` | `1tup-asm-{32,192,864}.cif` (assembly '1') | 186k / 1.12M / 5.04M | 288 / 1,728 / 7,776 (9) |
| `flat-192` / `flat-864` | `1tup-flat-{192,864}.bcif` (20 / 90 MB) | 1.12M / 5.04M | 2,496 / 11,232 (unshared) |
| `flat-192-cif` | `1tup-flat-192.cif` (124 MB text) | 1.12M | 2,496 |
| many | `tests/build/pipeline/fixtures/many/` (500 files, `many-fixtures --count 500`) | 1.2k–11.5k each | – |

## How to rerun

```sh
node scripts/pipeline-bench.mjs fixtures --sizes 4x4x3 --text      # asm/flat-192 + flat-192.cif (3 s)
node scripts/pipeline-bench.mjs fixtures --sizes 6x6x6             # asm/flat-864 (5 s); 2x2x2 likewise
node scripts/pipeline-bench.mjs run --fixture 7qpd --runs 5 --out out/7qpd.json          # 1tup-2x2x2 likewise
node scripts/pipeline-bench.mjs run --fixture packages/prism/viewer/src/_stories/dynamic-pdb/fixtures/3k7f-qfit.bcif --runs 5 --out out/3k7f.json
node scripts/pipeline-bench.mjs run --fixture flat-32 --runs 3 --out out/flat-32.json
node scripts/pipeline-bench.mjs run --fixture asm-192 --runs 2 --out out/asm-192.json    # asm-32 likewise
node scripts/pipeline-bench.mjs run --fixture flat-192 --out out/flat-192.json
node scripts/pipeline-bench.mjs run --fixture flat-192-cif --runs 2 --only parse,model
node scripts/pipeline-bench.mjs run --fixture asm-192 --quality low --only 'repr:*,theme-apply*,overpaint*,transparency*,release:*'
node scripts/pipeline-bench.mjs run --fixture asm-864 --out out/asm-864.json
node scripts/pipeline-bench.mjs run --fixture asm-864 --quality low --only 'repr:*,overpaint*,transparency*,release:*'
node scripts/pipeline-bench.mjs run --fixture flat-864 --skip 'repr:molecular-surface*,repr:interactions*,repr:ball-and-stick*'
node scripts/pipeline-bench.mjs run --fixture 7qpd --profile prof/7qpd && node scripts/pipeline-bench.mjs analyze prof/7qpd 8
node scripts/pipeline-bench.mjs run --fixture 7qpd --profile prof/a7qpd --alloc --only 'repr:*'
node scripts/pipeline-bench.mjs many-fixtures --count 500
node scripts/pipeline-bench.mjs many --count 300 --transaction stage --runs 2      # also --count 100, 500; --count 200 --transaction none
node scripts/pipeline-bench.mjs playback --out out/playback.json
```

The A/B base takes four steps:

1. Extract the old sources:
   `git archive 60fe4c2f1 packages scripts/workspace/inventory.json tests/tsconfig.json tsconfig.json tsconfig.base.json | tar -x -C <dir>`.
2. Link the packages: `<dir>/node_modules/@molstar/<name>` → `<dir>/<path>` for every `packages/` entry of
   `inventory.json`.
3. Copy `tests/src/pipeline/*.ts` into `<dir>/tests/src/pipeline`, and import `AssetManager` as `DefaultAssetManager`
   (the rename of #2019).
4. Bundle it as `scripts/pipeline-bench.mjs` does. Set `absWorkingDir: <dir>`, and set `nodePaths` to the repo's
   `packages/*/node_modules`, `packages/*/*/node_modules`, `tests/node_modules` and `node_modules/.pnpm/node_modules`
   (third-party only). Write the bundle to `<dir>/tests/build/pipeline/bench.mjs`, and symlink `<dir>/data` and the
   fixtures folder. Run it with `node --expose-gc --max-old-space-size=14000 --enable-source-maps`.

The bench never runs two at a time: a sequential queue script ran everything (13:55–14:29).

## Stage tables: 00-profile's fixtures (base `60fe4c2f1` → now, same bench)

ms:

| Stage | 7qpd | 3k7f | 8zzh | 1tup-2x2x2 | flat-32 | Likely PR |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `parse` | 12 → 12 | 3.9 → 4.0 | 4.2 → 4.2 | 6.5 → 6.5 | 6.3 → 6.4 |  |
| `model` | 7.1 → 7.0 | 8.0 → 8.1 | 7.0 → 6.6 | 4.9 → 4.9 | 103 → 95 | #2007 (04/P6), #2017 |
| `structure` | 1.0 → 1.1 | 1.2 → 1.3 | 0.6 → 0.6 | 2.0 → 2.1 | 9.8 → 10 |  |
| `lookup+boundary` | 4.2 → 4.3 | 2.4 → 4.8 **×2.0** | 1.6 → 1.5 | 12 → 7.2 **−38%** | 69 → 66 | #1994 |
| `bonds-intra` | 22 → 13 **−39%** | 30 → 16 **−46%** | 0 → 0 | 11 → 8.7 **−24%** | 248 → 130 **−47%** | #2007 (05/P7), #1992 |
| `bonds-inter` | 3.2 → 2.9 | 1.2 → 1.2 | 0.4 → 0.4 | 60 → 32 **−46%** | 54 → 32 **−40%** | #1994 (05/P4), #1992 |
| `secondary-structure` | 0.5 → 0.5 | 3.6 → 3.5 | 0.2 → 0.2 | 0.5 → 0.5 | 2.0 → 1.9 |  |
| `valence-model` | 6.1 → 6.0 | 8.1 → 8.0 | 0.1 → 0.1 | 51 → 55 | 52 → 52 |  |
| `interactions` | 38 → 35 | 98 → 88 | 0.9 → 1.1 | 508 → 390 **−23%** | 944 → 804 | #1996 (06/N1), #1992 |
| `asa` | 114 → 49 **−57%** | 185 → 75 **−59%** | throws → throws | 1873 → 937 **−50%** | 1813 → 863 **−52%** | #1986 (06/P1) |
| `theme-eval:element-symbol` | 1.9 → 1.9 | 2.1 → 2.0 | 0.5 → 0.5 | 8.6 → 9.5 | 8.8 → 9.7 |  |
| `theme-eval:residue-name` | 1.1 → 1.0 | 0.9 → 0.9 | 0.5 → 0.5 | 4.5 → 4.6 | 6.0 → 5.5 |  |
| `component:every-copy-derived` | 0.4 → 0.3 | 0.5 → 0.3 | 0.4 → 0.3 | 4.4 → 1.5 **−66%** | 3.0 → 1.0 **−67%** | #1994 (05/P5) |
| `component:polymer-derived` | 3.2 → 0.5 **÷6** | 1.1 → 0.4 | 0.5 → 0.4 | 61 → 5.4 **÷11** | 51 → 1.3 **÷39** | #1994 (05/P5) |
| `repr:cartoon` | 32 → 39 **+24%** | 26 → 30 | 22 → 27 **+27%** | 13 → 14 | 208 → 168 | #1993, #2007 (08/P3) — see R1 |
| `theme-apply:cartoon:secondary-structure` | 1.5 → 1.6 | 1.4 → 1.5 | 2.1 → 2.1 | 1.3 → 1.4 | 13 → 14 |  |
| `repr:ball-and-stick` | 29 → 28 | 30 → 31 | 5.6 → 5.7 | 29 → 32 | 485 → 515 |  |
| `theme-apply:ball-and-stick:chain-id` | 4.6 → 4.4 | 4.8 → 4.9 | 2.4 → 2.5 | 4.1 → 4.0 | 57 → 53 |  |
| `repr:spacefill` | 7.5 → 7.7 | 8.4 → 8.5 | 4.6 → 5.3 | 6.6 → 6.7 | 123 → 126 |  |
| `repr:gaussian-surface` | 626 → 375 **−40%** | 562 → 377 **−33%** | 588 → 330 **−44%** | 44 → 33 **−24%** | 769 → 610 **−21%** | #1987, #1990 |
| `theme-apply:gaussian-surface:element-symbol` | 261 → 82 **−68%** | 264 → 80 **−70%** | 180 → 45 **−75%** | 43 → 13 **−69%** | 684 → 319 **−53%** | #1991 (08/P1) |
| `repr:molecular-surface` | 1732 → 929 **−46%** | 2204 → 1266 **−43%** | 574 → 364 **−37%** | 429 → 243 **−43%** | 12437 → 6512 **−48%** | #1988 (02/P1), #1987 |
| `theme-apply:molecular-surface:residue-name` | 217 → 66 **−69%** | 221 → 64 **−71%** | 191 → 50 **−74%** | 44 → 12 **−74%** | 670 → 273 **−59%** | #1991 (08/P1) |
| `repr:interactions` | 30 → 7.9 **−73%** | 106 → 21 **−80%** | 3.7 → 4.3 | 451 → 74 **÷6** | 996 → 184 **÷5** | #2004 (09/P1) |

MB allocated:

| Stage | 7qpd | 3k7f | 8zzh | 1tup-2x2x2 | flat-32 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `parse` | 1.1 → 1.6 | 1.4 → 1.6 | 0.5 → 0.7 | 6.8 → 6.2 | 1.6 → 1.9 |
| `model` | 4.7 → 6.1 **+30%** | 6.5 → 6.6 | 11 → 11 | 3.3 → 4.5 **+36%** | 112 → 91 |
| `structure` | 0.6 → 0.6 | 0.2 → 0.8 | 0.4 → 0.4 | 2.2 → 2.5 | 6.2 → 5.3 |
| `lookup+boundary` | 0.8 → 1.8 | 1.0 → 1.1 | 1.8 → 1.5 | 12 → 13 | 54 → 53 |
| `bonds-intra` | 6.4 → 7.0 | 7.2 → 6.0 | 0.5 → 0.4 | 3.5 → 6.1 **+74%** | 32 → 41 **+27%** |
| `bonds-inter` | 2.4 → 1.4 | 1.0 → 0.7 | 0.5 → 0.5 | 94 → 16 **÷6** | 53 → 9.7 **÷5** |
| `secondary-structure` | 0.3 → 0.3 | 1.8 → 2.0 | 0.5 → 0.5 | 0.2 → 0.2 | 0.7 → 0.7 |
| `valence-model` | 9.3 → 8.6 | 10 → 10 | 0.5 → 0.5 | 80 → 80 | 93 → 91 |
| `interactions` | 31 → 29 | 90 → 100 | 0.5 → 0.6 | 283 → 266 | 684 → 580 |
| `asa` | 773 → 6.8 **÷114** | 851 → 6.1 **÷140** | throws → throws | 7300 → 83 **÷88** | 7447 → 16 **÷471** |
| `theme-eval:element-symbol` | 1.1 → 1.1 | 0.9 → 1.0 | 1.0 → 0.8 | 1.3 → 1.5 | 1.5 → 1.5 |
| `theme-eval:residue-name` | 0.7 → 0.6 | 0.4 → 0.3 | 0.5 → 0.4 | 0.6 → 0.6 | 0.5 → 0.4 |
| `component:every-copy-derived` | 0.5 → 0.3 | 0.5 → 0.3 | 0.4 → 0.4 | 5.5 → 1.4 **−75%** | 2.6 → 0.3 |
| `component:polymer-derived` | 2.2 → 0.3 | 1.0 → 0.3 | 0.4 → 0.3 | 99 → 7.7 **÷13** | 54 → 0.5 **÷107** |
| `repr:cartoon` | 38 → 78 **×2.1** | 31 → 67 **×2.2** | 26 → 56 **×2.1** | 18 → 27 **+52%** | 354 → 220 **−38%** |
| `theme-apply:cartoon:secondary-structure` | 0.4 → 0.4 | 0.5 → 0.4 | 0.7 → 0.5 | 0.2 → 0.3 | 5.4 → 5.7 |
| `repr:ball-and-stick` | 20 → 18 | 18 → 18 | 2.6 → 2.7 | 25 → 23 | 407 → 394 |
| `theme-apply:ball-and-stick:chain-id` | 1.6 → 1.7 | 0.9 → 1.3 | 1.4 → 1.2 | 1.6 → 2.2 | 10 → 11 |
| `repr:spacefill` | 2.9 → 3.1 | 6.8 → 7.1 | 2.1 → 2.1 | 5.4 → 5.9 | 124 → 124 |
| `repr:gaussian-surface` | 52 → 77 **+49%** | 50 → 77 **+53%** | 72 → 81 | 38 → 29 **−24%** | 217 → 156 **−28%** |
| `theme-apply:gaussian-surface:element-symbol` | 16 → 7.8 **−50%** | 34 → 12 **−63%** | 14 → 3.8 **−73%** | 31 → 10 **−67%** | 146 → 116 **−21%** |
| `repr:molecular-surface` | 110 → 154 **+40%** | 92 → 122 **+33%** | 150 → 189 **+25%** | 116 → 139 | 1570 → 1510 |
| `theme-apply:molecular-surface:residue-name` | 14 → 5.6 **−59%** | 27 → 6.9 **−75%** | 13 → 8.4 **−37%** | 31 → 9.9 **−68%** | 148 → 101 **−32%** |
| `repr:interactions` | 6.3 → 4.5 **−29%** | 18 → 9.8 **−45%** | 1.0 → 1.0 | 973 → 53 **÷18** | 1851 → 146 **÷13** |

## Stage tables: 13-large-scale's fixtures (base → now, same bench)

ms (`auto` quality: `low` at 186k, `lowest` at 1.12M):

| Stage | asm-32 | asm-192 | flat-192 | Likely PR |
| --- | ---: | ---: | ---: | --- |
| `parse` | 6.3 → 6.4 | 8.9 → 8.4 | 17 → 17 |  |
| `model` | 9.1 → 8.9 | 9.4 → 8.9 | 531 → 483 | #2007 (04/P6), #2017 |
| `structure` | 3.5 → 3.5 | 5.5 → 5.8 | 34 → 34 |  |
| `lookup+boundary` | 17 → 9.9 **−43%** | 32 → 32 | 338 → 304 | #1994 |
| `bonds-intra` | 12 → 10 | 14 → 12 | 1490 → 821 **−45%** | #2007 (05/P7), #1992 |
| `bonds-inter` | 59 → 32 **−46%** | 420 → 182 **−57%** | 359 → 204 **−43%** | #1994 (05/P4), #1992 |
| `secondary-structure` | 0.6 → 0.7 | 0.6 → 0.7 | 7.2 → 7.2 |  |
| `valence-model` | 49 → 51 | 275 → 358 ¹ | 210 → 217 |  |
| `interactions` | 519 → 379 **−27%** | 2836 → 2220 **−22%** | 5692 → 4627 | #1996 (06/N1), #1992 |
| `asa` | 1856 → 901 **−51%** | 10894 → 5648 **−48%** | 11126 → 5500 **−51%** | #1986 (06/P1) |
| `theme-eval:element-symbol` | 7.7 → 7.7 | 28 → 29 | 32 → 32 |  |
| `theme-eval:residue-name` | 8.0 → 8.0 | 38 → 37 | 24 → 23 |  |
| `component:every-copy-derived` | 4.7 → 2.8 **−40%** | 25 → 8.3 **−67%** | 29 → 6.6 **−77%** | #1994 (05/P5) |
| `component:polymer-derived` | 61 → 5.2 **÷12** | 444 → 26 **÷17** | 363 → 20 **÷18** | #1994 (05/P5) |
| `repr:cartoon` | 18 → 17 | 17 → 19 | 897 → 865 | #1993, #2007 (08/P3) — see R1 |
| `theme-apply:cartoon:secondary-structure` | 1.7 → 1.9 | 1.7 → 1.8 | 86 → 81 |  |
| `repr:ball-and-stick` | 33 → 31 | 28 → 27 | 4326 → 4658 |  |
| `theme-apply:ball-and-stick:chain-id` | 4.2 → 4.3 | 5.3 → 5.8 | 348 → 355 |  |
| `repr:spacefill` | 6.8 → 7.1 | 9.4 → 9.2 | 693 → 696 |  |
| `repr:gaussian-surface` | 41 → 29 **−28%** | 33 → 13 **−61%** | 793 → 678 | #1987, #1990 |
| `theme-apply:gaussian-surface:element-symbol` | 42 → 16 **−61%** | 2.2 → 2.2 | 137 → 124 | #1991 (08/P1) |
| `repr:molecular-surface` | 441 → 246 **−44%** | 204 → 131 **−36%** | 30478 → 19908 **−35%** | #1988 (02/P1), #1987 |
| `theme-apply:molecular-surface:residue-name` | 42 → 14 **−67%** | 1.9 → 1.9 | 107 → 109 | #1991 (08/P1) |
| `repr:interactions` | 454 → 70 **÷6** | 1354 → 327 **−76%** | 3440 → 1423 **−59%** | #2004 (09/P1) |

¹ noise: an isolated rerun (`--only valence-model --runs 5`) gives 178 → 184 ms, allocation 470 → 471 MB.

MB allocated:

| Stage | asm-32 | asm-192 | flat-192 |
| --- | ---: | ---: | ---: |
| `parse` | 2.0 → 1.5 | 1.2 → 1.8 | 3.4 → 2.8 |
| `model` | 13 → 14 | 17 → 14 | 641 → 500 **−22%** |
| `structure` | 5.3 → 7.2 **+36%** | 13 → 14 | 46 → 45 |
| `lookup+boundary` | 25 → 25 | 66 → 66 | 236 → 252 |
| `bonds-intra` | 11 → 15 **+33%** | 11 → 15 **+30%** | 255 → 199 **−22%** |
| `bonds-inter` | 104 → 23 **−78%** | 735 → 104 **÷7** | 412 → 80 **÷5** |
| `secondary-structure` | 0.3 → 0.3 | 0.4 → 0.4 | 5.2 → 7.6 **+46%** |
| `valence-model` | 95 → 96 | 489 → 488 | 497 → 499 |
| `interactions` | 496 → 336 **−32%** | 1226 → 1351 | 3915 → 3069 **−22%** |
| `asa` | 8581 → 130 **÷66** | 35913 → 615 **÷58** | 34587 → 89 **÷388** |
| `theme-eval:element-symbol` | 1.7 → 2.2 | 1.4 → 1.0 | 1.5 → 1.6 |
| `theme-eval:residue-name` | 4.9 → 4.8 | 27 → 27 | 1.1 → 1.1 |
| `component:every-copy-derived` | 6.2 → 2.9 **−53%** | 39 → 7.6 **÷5** | 34 → 7.6 **−77%** |
| `component:polymer-derived` | 100 → 8.9 **÷11** | 746 → 42 **÷18** | 397 → 20 **÷20** |
| `repr:cartoon` | 35 → 37 | 28 → 27 | 1421 → 1267 |
| `theme-apply:cartoon:secondary-structure` | 0.4 → 0.3 | 0.6 → 0.2 | 35 → 39 |
| `repr:ball-and-stick` | 31 → 32 | 27 → 27 | 2271 → 2264 |
| `theme-apply:ball-and-stick:chain-id` | 1.9 → 1.8 | 1.2 → 0.6 | 49 → 54 |
| `repr:spacefill` | 4.7 → 5.9 **+26%** | 7.9 → 7.7 | 689 → 692 |
| `repr:gaussian-surface` | 37 → 26 **−31%** | 14 → 21 **+49%** | 554 → 549 |
| `theme-apply:gaussian-surface:element-symbol` | 28 → 11 **−60%** | 0.7 → 0.7 | 30 → 29 |
| `repr:molecular-surface` | 124 → 145 | 120 → 105 | 8392 → 8422 |
| `theme-apply:molecular-surface:residue-name` | 33 → 17 **−50%** | 0.8 → 0.7 | 34 → 33 |
| `repr:interactions` | 978 → 55 **÷18** | 2649 → 277 **÷10** | 5294 → 850 **÷6** |

MB retained (only rows that are not identical; every load and property row retains the same MB on both):

| Stage | 7qpd | 3k7f | 1tup-2x2x2 | flat-32 | asm-32 | asm-192 | flat-192 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `valence-model` | -0.8 → -0.8 | -3.1 → -3.1 | -25 → -25 | -61 → -61 | -12 → -12 | -73 → -73 | 7.1 → 7.1 |
| `interactions` | -3.8 → -3.8 | -1.3 → -1.3 | 17 → 17 | -22 → -21 | 20 → 20 | 128 → 128 | 255 → 257 |
| `repr:cartoon` | 24 → 15 **−40%** | 17 → 12 **−27%** | 3.1 → 2.6 | 125 → 105 | 3.6 → 3.0 | 3.0 → 2.8 | 379 → 350 |
| `repr:gaussian-surface` | 34 → 31 | 31 → 30 | 3.9 → 3.2 | 195 → 176 | 4.0 → 3.3 | 2.7 → 2.2 | 324 → 203 **−37%** |
| `repr:molecular-surface` | 28 → 24 | 26 → 26 | 3.4 → 2.9 | 171 → 157 | 3.4 → 2.9 | 2.6 → 2.1 | 323 → 201 **−38%** |
| `repr:interactions` | 30 → 3.1 **÷10** | 121 → 10 **÷12** | 228 → 45 **÷5** | 643 → 144 **−78%** | 228 → 46 **÷5** | 442 → 266 **−40%** | 1348 → 863 **−36%** |
| `repr:ball-and-stick` | 12 → 12 | 12 → 12 | 9.9 → 9.9 | 389 → 390 | 10 → 10 | 11 → 11 | 2334 → 2335 |
| `repr:spacefill` | 1.0 → 1.0 | 1.1 → 1.1 | 1.2 → 1.2 | 42 → 42 | 1.2 → 1.2 | 2.4 → 2.4 | 254 → 254 |

### 864 copies (v1 = 13-large-scale's report → now; not A/B)

13 ran these at load 2–14 with swap; flat-864's v1 row is from the full session that also held ball-and-stick (10.6 GB).

| Stage | asm-864 | flat-864 |
| --- | ---: | ---: |
| `model` | 16 → 11 **−29%** | 2933 → 2221 **−24%** |
| `lookup+boundary` | 126 → 171 **+36%** | 1409 → 1256 |
| `bonds-intra` | 17 → 11 **−34%** | 7305 → 3564 **−51%** |
| `bonds-inter` | 2132 → 784 **−63%** | 2959 → 934 **−68%** |
| `valence-model` | 1198 → 941 **−21%** | 1186 → 922 **−22%** |
| `interactions` | 12998 → 9800 **−25%** | 30793 → 21209 **−31%** |
| `asa` | 31379 → 26778 | 28350 → 26200 |
| `component:polymer-derived` | 2033 → 190 **÷11** | 1954 → 82 **÷24** |
| `repr:cartoon` | 69 → 52 **−24%** | 11382 → 5082 **−55%** |
| `theme-apply:cartoon:secondary-structure` | 4.0 → 2.5 **−38%** | 951 → 400 **−58%** |
| `repr:ball-and-stick` | 83 → 54 **−35%** | 191182 → not run |
| `repr:spacefill` | 20 → 21 | 2636 → 1525 **−42%** |
| `repr:gaussian-surface` | 52 → 38 **−26%** | 6668 → 3376 **−49%** |
| `repr:molecular-surface` | 228 → 166 **−27%** | – |
| `repr:interactions` | 5292 → 1321 **−75%** | – |

MB allocated:

| Stage | asm-864 | flat-864 |
| --- | ---: | ---: |
| `model` | 28 → 28 | 2696 → 2120 **−21%** |
| `lookup+boundary` | 245 → 300 **+22%** | 888 → 913 |
| `bonds-intra` | 20 → 25 **+24%** | 1152 → 838 **−27%** |
| `bonds-inter` | 3648 → 443 **÷8** | 1990 → 289 **÷7** |
| `valence-model` | 2132 → 2131 | 2156 → 2155 |
| `interactions` | 5160 → 6223 **+21%** | 17654 → 13442 **−24%** |
| `asa` | 2818 → 2849 | 314 → 420 **+34%** |
| `component:polymer-derived` | 3756 → 217 **÷17** | 1976 → 85 **÷23** |
| `repr:cartoon` | 43 → 50 | 5988 → 5501 |
| `theme-apply:cartoon:secondary-structure` | 0.5 → 0.9 | 132 → 135 |
| `repr:ball-and-stick` | 31 → 31 | 10128 → not run |
| `repr:spacefill` | 10 → 9.7 | 565 → 724 **+28%** |
| `repr:gaussian-surface` | 39 → 40 | 2384 → 2304 |
| `repr:molecular-surface` | 138 → 125 | – |
| `repr:interactions` | 1761 → 1214 **−31%** | – |

### Text mmCIF (`flat-192-cif`, 124 MB)

| Stage | v1 (13, first pass) | now (2 runs) |
| --- | ---: | ---: |
| `parse` | 1.0–1.2 s, 320 MB retained | 576–698 ms, 319–320 MB retained |
| `model` | 0.7–0.8 s | 389–391 ms / 269–291 MB alloc |

### Interactive paths (ms; base → now)

| Stage | flat-32 | asm-192 | flat-192 | flat-864 (v1 13 → now) |
| --- | ---: | ---: | ---: | ---: |
| `hover:cartoon` | 1.1 → 0.5 | 0.5 → 0.8 | 11 → 3.2 | 191 → 19 |
| `select:cartoon` | 0.8 → 0.1 | 0.3 → 0.4 | 6.5 → 0.6 | 36 → 0.6 |
| `select-every-copy:cartoon` | 1.6 → 0.3 | 0.5 → 1.0 | 30 → 2.6 | 685 → 13 |
| `select-every-copy:ball-and-stick` | 3.9 → 0.9 | 9.8 → 3.7 | 40 → 4.8 | 506 → – |
| `overpaint:cartoon` | 2.3 → 2.3 | 0.5 → 0.4 | 24 → 29 | 152 → 140 |
| `overpaint:ball-and-stick` | 3.2 → 2.9 | 1.2 → 1.5 | 27 → 29 | 255 → – |
| `overpaint:gaussian-surface` | **628 → 2.2** | 0.5 → 0.4 | 12 → 11 | 86 → 78 |
| `transparency:gaussian-surface` | **420 → 2.5** | 2.8 → 2.8 | 12 → 12 | 197 → 66 |
| `overpaint:molecular-surface` | **559 → 2.3** | 0.5 → 0.5 | 11 → 12 | – |

At pinned `low` quality (where 13/P1 bit; base → now, ms / MB alloc / peak RSS MB after):

| Stage | asm-192-low | asm-864-low (v1 13 → now) |
| --- | ---: | ---: |
| `overpaint:gaussian-surface` | **4692 / 896 / 2534 → 20 / 10 / 334** | **12564 / 4006 / 3311 → 31 / 12 / 463** |
| `transparency:gaussian-surface` | **3341 / 856 → 14 / 21** | 6725 / 3704 → 32 / 16 |
| `overpaint:molecular-surface` | **4091 / 732 → 15 / 8** | **11697 / 3279 → 44 / 25** |
| `transparency:molecular-surface` | **2946 / 722 → 7.6 / 5** | 6382 / 477 → 36 / 17 |
| `theme-apply:gaussian-surface:element-symbol` | 42 → 19 | 58 → – |
| `repr:molecular-surface` | 436 → 248 | 341 → 299 |

`repr:interactions` in `--only 'repr:*'` runs computes the interactions property inside its own stage (2.7 s
asm-192-low, 12.8 s asm-864-low). 13's rows had the `interactions` stage on, so those cells are not comparable.

## Render buffers (MB of typed arrays in render-object values; base → now)

| Repr | 7qpd | 3k7f | 8zzh | 1tup-2x2x2 | flat-32 | asm-192 | flat-192 | asm-864 / flat-864 (now) | Largest now |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| cartoon | **23.4 → 13.8** | **16.1 → 11.6** | **13.2 → 9.6** | 1.9 → 1.2 | 36.2 → 35.1 | 0.8 → 0.7 | 60.3 → 60.3 | 2.0 / 271 | `aPosition`, `aNormal` (flat-192 21 MB each) |
| ball-and-stick | 10.9 → 10.9 | 11.3 → 11.3 | 0.0 | 8.3 → 8.3 | 165 → 165 | 8.9 → 8.9 | 991 → 991 | 21.6 / – | `aMapping`/`aStart`/`aEnd` 189 MB each at 1M (frozen) |
| spacefill | 0.5 | 0.5 | 0.1 | 0.5 | 4.5 | 1.6 | 27.2 | 6.5 / 197 | `tPositionGroup` |
| gaussian-surface | 33.6 → 30.3 | 30.2 → 30.0 | **19.0 → 10.7** | 3.2 → 2.5 | 72.8 → 72.8 | 1.9 → 1.4 | 11.6 → 11.6 | 6.3 / 52 | `elements`, `aPosition` |
| molecular-surface | 27.0 → 23.9 | 25.5 → 25.4 | **19.4 → 11.6** | 2.7 → 2.2 | 62.8 → 62.8 | 1.9 → 1.4 | 10.6 → 10.6 | 6.3 / – | |
| interactions | **29.1 → 2.4** | **121 → 9.8** | 0 | **227 → 45** | **422 → 83** | **441 → 265** | **832 → 500** | 1,195 / – | cylinder impostors (#2004): `aMapping`/`aStart`/`aEnd` 6× per cylinder (frozen schema) |

The trimmed copies of #2007 (08/P3) reach the cartoon buffers of small entries (−40 %) and the surfaces of 8zzh.
Interactions are now impostors: −80 to −92 % of buffers, but still 500 MB at 1M unshared and 1.2 GB at asm-864.

## Many structures (14) and playback

### `many --transaction stage` (ms / MB alloc; per-structure in brackets)

| Stage | N = 100 | N = 300 (2 runs) | N = 500 | N = 200, `none` |
| --- | ---: | ---: | ---: | ---: |
| `load` | 442 / 534 (4.4) | 1268 / 1388 (4.2) | 2109 / 2324 (4.2) | 918 / 1163 (4.6) |
| `model+structure` | 39 / 94 (0.39) | 151 / 177 (0.50) | 323 / 373 (0.65) | 704 / 1496 (3.5) |
| `bonds` | 669 / 321 (6.7) | 1818 / 789 (6.1) | 2984 / 1230 (6.0) | 1251 / 540 (6.3) |
| `repr:auto` | 1172 / 1040 (11.7) | 4228 / 2639 (14.1) | 9549 / 4486 (**19.1**) | 7349 / **11401** (36.7) |
| `theme:all` | 81 / 77 | 228 / 184 | 431 / 302 | 150 / 138 |
| `theme:one` | 7 | 19 | **39** | 9 |
| `hover` (per event) | 0.77 | 2.1 | **3.5** | 1.5 |
| `select` (per event) | 0.54 | 1.5 | **4.9** | 1.2 |
| `deselect-all` | 3 | 11 | 22 | 8 |
| `overpaint:one` / `overpaint:all` | 6 / 23 | 21 / 76 | 30 / 132 | 10 / 52 |
| `remove:one` / `remove:all` | 6 / 12 | 22 / 54 | 32 / 110 | 11 / 39 |

Against 14's "after" (F1–F8):

| Metric | Before (14) | 14 after | Now |
| --- | ---: | ---: | ---: |
| Load to full scene (load + model + bonds + secondary structure + `repr:auto`), N = 100 / 300 / 500 | – | 2.8 / 8.8 / 16.3 s | 2.3 / 7.5 / 15.0 s |
| `repr:auto` per structure, N = 100 / 300 / 500 | – | 14.7 / 14.9–17.9 / 20.8 ms | 11.7 / 14.1 / 19.1 ms |
| `repr:auto` allocated, N = 300 | – | 2.65 GB | 2.64 GB |
| overpaint of all, N = 100 / 300 / 500 | – | 31 / 92 / 187 ms | 23 / 76 / 132 ms |
| theme / overpaint / removal of one, N = 500 | – | 45 / 46 / 46 ms | 39 / 30 / 32 ms |
| hover / select / deselect all, N = 500 | – | 5.1 / 5.4 / 24 ms | 3.5 / 4.9 / 22 ms |
| `none` mode, N = 200: `repr:auto` | 28.6 s / 51 GB | 8.6 s / 11 GB | 7.3 s / 11.4 GB |
| `none` mode, N = 200: model and structure | 3.4 s | 0.9 s | 0.70 s |

Still super-linear in N:

- `repr:auto` per structure: 11.7 → 19.1 ms from N = 100 to 500.
- `model+structure` per structure: 0.39 → 0.65 ms.
- Per-event select: 0.54 → 4.9 ms (×9 for ×5 N).
- Single-structure actions grow linearly with N: `theme:one` 7 → 39 ms.

These are 14's "not fixed" items 1 and 2 (the canvas marks every representation, and the hierarchy is rebuilt by a
whole-state walk).

### `playback` (defaults; no v1 baseline: this is the baseline)

| Scene | Ticks / applied | ms/tick med / max | ms per applied tick | KB alloc/tick | cells updated/tick | marks | setSnapshot | changed arrays KB/applied |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `snapshots:plain` | 181 / 3 | 0 / 1.7 | 1.22 | 102 | 0.05 | 0 | 31 | 1 |
| `snapshots:effects` | 181 / 3 | 0 / 1 | 0.93 | 77 | 0.05 | 0 | 31 | 0 |
| `traj:nmr:preset` | 20 / 9 | 0.11 / 4.7 | 2.78 | 2881 | 4.5 | 0 | 0 | 236 |
| `traj:nmr:preset+focus+select` | 20 / 9 | 0.03 / 10.4 | 6.58 | 5278 | 6.75 | 108 | 0 | 302 |
| `traj:nmr:preset+focus` | 20 / 9 | 0.03 / 10.4 | 5.07 | 3316 | 6.75 | 8 | 0 | 302 |
| `traj:synth:preset` | 80 / 39 | 0.03 / 28.6 | 9.39 | 5320 | 4.39 | 0 | 0 | 14176 |
| `traj:synth:preset+focus+select` | 80 / 39 | 0.03 / 18.3 | 11.88 | 5864 | 6.83 | 468 | 0 | 14235 |
| `traj:synth:preset+overpaint` | 80 / 39 | 0.04 / 15.2 | 9.13 | 4426 | 5.85 | 0 | 0 | 14177 |
| `traj:synth:cartoon` | 80 / 39 | 0.03 / 13.3 | 8.12 | 3710 | 2.92 | 0 | 0 | 14127 |
| `traj:synth:ball-and-stick` | 80 / 39 | 0.03 / 9.9 | 4.97 | 1489 | 2.92 | 0 | 0 | 4831 |
| `traj:synth:gaussian-surface` | 80 / 39 | 0.04 / 335.7 | 289.93 | 10970 | 2.92 | 0 | 0 | 31152 |

- In `traj:synth:*` (7qpd topology, 20 jittered frames) every applied frame rebuilds geometry. That is 14 MB of
  changed render arrays per frame for the preset or cartoon, 4.8 MB for ball-and-stick, and 31 MB plus 290 ms for the
  CPU gaussian surface (there is no GPU in the stub).
- `+focus+select` adds 12 marks per applied frame and 2.5 ms.
- The two `snapshots:*` scenes call the whole-tree `setSnapshot` 31× for 3 applied frames (28 no-op restores, < 1 ms
  each).
- Focus scenes print mutative's warning "The return value does not contain any draft" (P9).

## Hotspots now (CPU, 100 µs sampling, after a warm-up run; self unless "total"; % of the stage)

| Stage [fixture, wall] | Function (file:line) | Share | Note |
| --- | --- | ---: | --- |
| repr:molecular-surface [flat-192, 20.6 s, **8.4 GB alloc**] | `projectTorus` `core/src/math/geometry/molecular-surface.ts:282` | 50 % self, 63 % total | **7.4 GB of the 8.3 GB sampled allocation is `projectTorus` self** (P3); `query` grid.ts:315 12 %, `obscured` :145 12 %, `projectPointsRange` :194 6 %; mesh 23 % total |
| repr:molecular-surface [7qpd, 0.96 s] | `projectTorus` :282 43 %, `obscured` :145 17 %, `projectPointsRange` :194 15 % | | marching cubes `stridedRow` `graphics/src/geo/util/marching-cubes/algorithm.ts:342` 8 % |
| asa [asm-192 5.9 s; flat-192 5.7 s; asm-864 26.8 s] | `computeRange` `model/src/props/computed/accessible-surface-area/shrake-rupley/area.ts:74` | 44–45 % | `GridLookup3D.query` `core/src/math/geometry/lookup3d/grid.ts:315` 35–42 %, `add` `lookup3d/common.ts:18` 6–8 %; asm allocates 0.6 GB through the structure lookup (`structure/util/lookup3d.ts:24`); same cost on asm as on flat (13/P6) |
| interactions [flat-192 4.9 s / 3.2 GB; asm-192 2.3 s / 1.4 GB] | `query` grid.ts:315 | 32 % | `addStructureContacts` `model/src/props/computed/interactions/contacts.ts:242` 8–16 % self, 50–83 % total via `eachUnitPair` `structure/structure.ts:1252`; intra-unit 33 % on flat; `Vec3.dot` vec3.ts:372 8 % |
| repr:ball-and-stick [flat-192, 4.2 s, 2.3 GB] | GC 10 %, `copyVertexBlock` `graphics/src/gl/renderable/merged.ts:521` 9 % (frozen), `CylindersBuilder.add` `geo/geometry/cylinders/cylinders-builder.ts:107` 8 % | | `createLinkCylinderImpostors` `repr/structure/visual/util/link.ts:361` 39 % total; `createUnitsRenderObject` 22 %, merge 19 % total (07/P7, frozen) |
| repr:interactions [flat-192, 1.5 s, 844 MB] | GC 16 %, `ChunkedArray._compact` `core/src/data/util/chunked-array.ts:122` 7 %, `addFixedCountDashes` cylinders-builder.ts:130 7 % (23 % total) | | `createLinkCylinderImpostors` 54 % total; intra-unit 36 %, inter-unit 24 % |
| repr:cartoon [flat-192, 0.86 s, 1.2 GB] | GC 15 %, `Vec3.sub`/`dot` 13 % | | `createPolymerTraceMesh` `repr/structure/visual/polymer-trace-mesh.ts:62` 31 % total; `createNucleotideRingMesh` `nucleotide-ring-mesh.ts:78` **14 %** total |
| repr:cartoon [7qpd, 43 ms, 74 MB] | `createPolymerTraceMesh` :62 17 %, GC 16 %, `addTube` `geo/geometry/mesh/builder/tube.ts:52` 10 % | | alloc: `addTube` 36 MB (base 20), `interpolateNormals` `polymer/curve-segment.ts:164` 14 MB (base 1) (R1) |
| repr:gaussian-surface [7qpd, 0.40 s] | `accumulateRange` `core/src/math/geometry/gaussian-density/cpu.ts:67` | **36 %** | `stridedRow` algorithm.ts:342 23 %, GC 10 %, `interpolate` algorithm.ts:189 8 % (42 MB alloc, R2); `accumulateRange` OSR-exits 1,850× in 3 runs (P4) |
| repr:gaussian-surface [flat-192, 0.69 s] | `accumulateRange` 17 %, `Vec3.dot` 12 % (boundaries) | | `createMergedValues` `gl/renderable/merged.ts:308` 23 % total (frozen) |
| repr:spacefill [flat-192, 0.68 s, 674 MB] | `Vec3.dot` 13 %, GC 11 %, `computeExtrema` `core/src/math/geometry/boundary-helper.ts:30` 5 % | | `createValues` `geo/geometry/spheres/spheres.ts:301` 36 % total; merge 26 % total (frozen) |
| bonds-intra [flat-192, 0.76 s] | `query` grid.ts:315 | 55 % | `findBonds` `model/src/model/structure/structure/unit/bonds/intra-compute.ts:176` 19 % self, 95 % total |
| lookup+boundary [flat-192, 0.36 s, 303 MB] | `Vec3.dot` vec3.ts:372 37 %, `computeSphereExtrema` boundary-helper.ts:43 18 % | | per-unit `getBoundary` `core/src/math/geometry/boundary.ts:33` 83 % total |
| model [flat-192, 0.58 s, 468 MB] | `arrayView` `core/src/data/db/column.ts:485` 17 %, GC 13 %, BCIF `stringArray` `io/src/common/binary-cif/decoder.ts:219` 9 %, `_makeBuckets` `core/src/data/util/buckets.ts:20` 7 % | | `sortAtomSite` `formats/structure/basic/sort.ts:17` 29 % total |
| valence-model [asm-192 0.32 s, 484 MB] | `getEdgeIndices` `core/src/math/graph/inter-unit-graph.ts:57` 23 % | | `calculateHydrogensCharge` `props/computed/chemistry/valence-model.ts:104` 92 % total, `isConjugated` :41 40 % total; flat: `advance` `unit/bonds.ts:230` 14 % |
| bonds-inter [asm-192, 0.19 s] | `query` grid.ts:315 45 %, `findPairBonds` `unit/bonds/inter-compute.ts:90` 14 % | | `hasCommonStructConnRecord` inter-compute.ts:449 9 % |
| theme-apply on surfaces [7qpd, 66–82 ms] | `splatVertex` `graphics/src/geo/geometry/mesh/color-smoothing.ts:332` | 66 % | `calcMeshColorSmoothing` :414 20 % |
| theme-apply:ball-and-stick [flat-192, 0.32 s] | `LocationIterator.move` `geo/util/location-iterator.ts:85` 15 %, `createGroupColor` `geo/geometry/color-data.ts:268` 11 % (67 % total) | | `chain-id` `color` `theme/color/chain-id.ts:97` 9 %, `bond.ts:174` `getLocation` 6 % |

## Observations

**Improvements** (A/B, same bench):

1. **ASA** −50 to −59 % time everywhere; allocation 7–36 GB → 6–615 MB (#1986).
2. **Molecular surface** −35 to −48 % (flat-32 12.4 → 6.5 s, 7qpd 1.73 → 0.93 s; #1988, #1987).
3. **Surface recolor** −53 to −75 % (#1991).
4. **One-chain overpaint or transparency on surfaces**: from ½–12 s to 2–44 ms on every fixture (#2005, #2012). At
   asm-192-low the peak RSS after the stage drops from 2.5 GB to 0.33 GB.
5. **Interactions representation** ÷4–6 time, ÷6–18 allocation, buffers −40 to −92 % (#2004).
6. **Component derivations** (inter-unit bond subsets) ÷11–39 (#1994).
7. **`bonds-inter`** −40 to −68 %, allocation ÷5–8 (#1994).
8. **`bonds-intra`** −24 to −51 % (#2007, #1992).
9. **Marks on unshared structures**: hover −70 %, select and select-every-copy ÷8–17 on flat-192; flat-864
   select-every-copy 685 → 13 ms (#1995).
10. **Gaussian surface (CPU)** −21 to −61 % (#1987, #1990).
11. **Large unshared loads at 5M**: model −24 %, `bonds-intra` −51 %, `repr:cartoon` −55 %. Part of that is 13's
    swapping baseline.

**Regressions** (A/B):

- **R1 `repr:cartoon` on small entries.**
  - 7qpd: 32 → 39 ms (+24 %), allocation 38 → 78 MB (×2.1). 3k7f and 8zzh are the same, ×2.1–2.2 allocation.
  - Only from the second build on: the first build is 32 ms / 57 MB on both. With cartoon alone, 10 runs: base
    31 ms / 37 MB, now 37 ms / 80 MB.
  - Allocation profile: `addTube` `tube.ts:52` 36 MB (base 20 MB) and `interpolateNormals` `curve-segment.ts:164`
    14 MB (base 1 MB).
  - `--trace-deopt`: `interpolateNormals`' Maglev code is deoptimized after each full GC ("embedded weak objects
    cleared"). The scalar-local rewrite (#1993) evidently runs in a tier that boxes doubles.
  - flat-32 got faster (208 → 168 ms, −38 % allocation), so the cost is per build of a small structure. Interactive
    cartoon rebuilds (playback: 8 ms per frame on 7qpd) are this case.
  - Slice 08 / repr should look at it.
- **R2 Surface allocation on small entries.**
  - 7qpd `repr:gaussian-surface` 52 → 77 MB (+49 %) and `repr:molecular-surface` 110 → 154 MB (+40 %), although
    both got 40–46 % faster.
  - Marching cubes `interpolate` (`algorithm.ts:189`) allocates 42 MB, against 10 MB at `algorithm.ts:163` in base.
    Likely boxed doubles in the strided branch of #1987.
- **R3 `lookup+boundary` on 3k7f** 2.4 → 4.8 ms, steady over runs 3–5. Small in absolute terms. Entry with qFit alt
  locs; check #1992 / #1994.
- **Noise, not regressions.** `valence-model` asm-192 (+30 %, a rerun showed it equal). `repr:ball-and-stick`
  flat-192 4.3 → 4.7 s (+8 %, allocation equal). The +20–34 % allocation of `bonds-intra` on 186k assemblies is from
  #2007's per-residue dictionary walk (time −24 to −47 %).

**Biggest remaining stages**, in seconds:

| Stage | Cost | Notes |
| --- | --- | --- |
| `repr:molecular-surface` unshared | 6.5 s at 186k, 19.9 s at 1M (`lowest`), 8.4 GB allocated | P3 |
| `asa` | 0.9 s at 186k, 5.5–5.6 s at 1M, 26–27 s at 5M | same on assemblies, 13/P6 |
| `interactions` | 0.4–0.8 s at 186k, 2.2–4.6 s at 1M, 9.8 s (asm) and 21 s (flat, 13.4 GB allocated) at 5M | |
| `repr:ball-and-stick` unshared | 4.7 s and 2.3 GB retained at 1M | 13: 42–191 s at 5M; not rerun |
| unshared 5M load | model 2.2 s, `bonds-intra` 3.6 s, `repr:cartoon` 5.1 s | peak RSS 4.4 GB without ball-and-stick |
| CPU gaussian surface | 0.38 s per build on 7qpd | 290 ms per playback frame |

**Perf leads for the slice auditors** (new; not yet verified beyond the measurement):

| Id | Lead | Where |
| --- | --- | --- |
| P1 | = R1 (cartoon, small entries) | `tube.ts:52`, `curve-segment.ts:164` |
| P2 | = R2 (marching-cubes `interpolate` allocation) | `algorithm.ts:189` |
| P3 | `projectTorus` allocates 7.4 GB self on flat-192 (89 % of the stage's sampled bytes; 42 MB on 7qpd, base 20 MB). Look for doubles boxed per probe or voxel (closure-held `lastClip`/`tCount`, the `px` shadowing, Maglev). | `molecular-surface.ts:282` |
| P4 | `accumulateRange`, the top self function of the CPU gaussian surface (36 % on 7qpd), deoptimizes 1,850× in 3 runs with "exit from OSR'd inner loop", as does `projectPointsRange` (475×). Each call starts unoptimized and OSRs: split the per-atom kernel out so that it gets a regular optimized version. | `gaussian-density/cpu.ts:67`, `molecular-surface.ts:194` |
| P5 | `nucleotide-ring-mesh` is 14 % of the cartoon on flat-192 (1tup's DNA) | `nucleotide-ring-mesh.ts:78` |
| P6 | `lookup+boundary`: `Vec3.dot` in per-unit `computeSphereExtrema` is 55 % of 0.36 s at 1M unshared. The EPOS boundary, v1 item, still open. | `boundary-helper.ts:43` |
| P7 | many: `repr:auto` per structure 11.7 → 19.1 ms (N 100 → 500); `select` per event ×9 for ×5 N | 14 "not fixed" 1–2 (Canvas3D marks: frozen) |
| P8 | Interactions buffers still 1.2 GB on asm-864 (inter-unit contacts × 864 copies; the 6× per-cylinder impostor schema is frozen) | `interactions-inter-unit-cylinder.ts:218` |
| P9 | Focus updates return new objects from `produce` recipes (`(old) => ({ ...old, bundle })`). On every focus change mutative walks the returned value recursively (`handleReturnValue`) and prints a `console.warn` (seen in playback). | `plugin/core/src/behavior/dynamic/selection/structure-focus-representation.ts:212,217,271,272` → `core/src/state/state/builder.ts:341,345` |

## Bugs (re-verified or new)

| Id | Status | Where | Failure scenario |
| --- | --- | --- | --- |
| 06/B3 | **still open** | `model/src/props/computed/accessible-surface-area.ts:62` (`isApplicable: () => true`) | `asa` on 8zzh: "Property only available for atomic models." (base and now) |
| 01/B1, 04/B1 | **still open** (by code) | `model/src/formats/structure/basic/properties.ts:189` `getUniqueComponentNames = memoize1(...)`; `core/src/util/memoize.ts:30` `new Map` | every load without `chem_comp` (7qpd) pins its `BasicData` and source text |
| B1 (new) | bench | `tests/src/pipeline/bench.ts:324` `selections()` | `run --fixture data/examples/ihm-coarse-8zzh.cif`: `TypeError: Cannot read properties of undefined (reading 'elements')` after `theme-eval`. `first` is `undefined` without an atomic polymer unit, and the component stages are required, so no `--skip` avoids it. Fall back to `s.units[0]`, or skip the component and interactive stages. |
