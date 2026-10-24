# Verify A: slices 01, 02, 03, 13 (adversarial pass)

HEAD `16277d822`, v2 base `0944823d4`. Load fell to ~4 during this pass, so the re-runs below are less inflated than the
slice reports' (load 25–30). Probes of this pass: `scratchpad/pa3/verA/` (`gzrss.ts`, `mp-ab.ts` + `decode-slice.ts`);
re-runs of `01/sniff4.mjs`, `01/vpal.mjs`, `13/msgpack.mjs`, `13/ascii.mjs`.

| Id | Verdict | Corrected severity / estimate | Reason |
| --- | --- | --- | --- |
| 01/P13 = 12/P11 | CORRECTED (duplicate pair; keep one) | ~0.27 s per `auto` open of flat-192.cif at load 4 (5 × ~55 ms native scans), +0.11 s over the string search; not 0.6 s, not 2.1–2.5 s | mechanism right; 12's timing includes a 124 MB `new Uint8Array(bytes)` copy and round-1 JIT; 01's "one pass" JS table is no faster (277 vs 267 ms in its own probe) |
| 01/P14 | CONFIRMED | ~0.2 s and 91 MB per 1.1M-value palette at load 4 (loop 4–17 ms) | min/max over `toNumber` gives the same [min, max]; reached by Seq* `makePaletteFunction` with every row |
| 01/P15 | CONFIRMED, fix CORRECTED | one extra transient copy per held download on `.molx` save, as stated | "cache the made `File`" would retain the bytes twice; give the serializer `bytes` when held (0 copies vs v2's 1) |
| 01/B17 | CORRECTED | low (not medium) | allocations real (358 MiB declared → `arrayBuffers` +363 MiB) but pages untouched: RSS +25 MiB; the file fails either way; the harm is the empty `TypeError` and a RangeError risk on huge declared sizes |
| 01/B18 | CONFIRMED | low | storage compared before the transform; only caller (`atomic.ts:158`, same schema) gets false negatives only (lost sharing); API also gives a false positive for columns whose transforms differ (not reached) |
| 02/B6 | CONFIRMED | medium (structure-level surface visuals and very large units only) | formula gives 93.5k atoms/chunk at res 8 vs 6.7k before; update-count ratios 4.3×/8.9×/14× at res 1/2/8 match the probe's 442/104, 222/26, 56/4 |
| 02/B7 | CONFIRMED | medium | hand check: d = (0, 0, 4.5) in the L = 10 dodecahedron → (−5, −5, −2.57), \|m\| 7.52; only caller is centre pinning |
| 02/P14 | CONFIRMED | small: −6 % at res 8 only | early exit is safe by the triangle inequality (atom `a` can never satisfy it); no writes are skipped; `lastClip` reset by the next call |
| 03/P18 | PLAUSIBLE, not reachable in one plugin | low | each plugin has one data `State`; its `updateQueue` runs transforms one at a time; the out-of-state parses (`.run()` with no observer) are synchronous; needs ≥ 2 plugins or app code running observable parses in parallel |
| 13/P1 | CONFIRMED, fix CORRECTED | whole `decodeMsgPack` of flat-192.bcif 33 → 3 ms (median of 7, same process, load 4): −30 ms, ~80 % of `parseBinary` (37 ms) | output identical; write `Uint8Array.prototype.slice.call(buf, o, o + n)`: `Buffer#slice` is a view, so a Node `Buffer` input would break `view(data.buffer)` |
| 13/P4 | CONFIRMED | low: 5–7 bytes 49–87 → 28–34 ms per 1M, 48–96 → 24 B transient | retained strings are the same (flat below 13 chars); reached only by unique 5–7-byte `str` values (CSV, names), not the mmCIF build |
| 13/P5 | PLAUSIBLE, fix CORRECTED | ~−10–20 ms (estimated, not measured idle) | ~10 grows, intermediates sum to ~40 % of the final size; "uncapped after 4096 rows" brings back the over-reserve the cap prevents (a 10k-row loop near the top extrapolated to the whole file); use a higher cap (8–16×) |
| 13/B3 | CONFIRMED | low | the cache key is the requested name, the lookup is case-insensitive; reached only by mixed-case requests (MVS/user field names) |
| 13/B4 | REFUTED | — | GROMACS writes and reads `.gro` columns by C byte widths; counting bytes matches the writer, and v2's char count was the one off for multi-byte names |
| 13/B5 | CONFIRMED | very low | `parseFloatBytes` stops at the first bad byte (`abc` → 0, `1,5` → 1, also `nan` → 0); the other text parsers are lenient the same way |
| 13/B1 vs 01/B10 | DUPLICATE, fix CORRECTED | low (keep 01/B10, add 13's bgzip case) | same root cause; 13's fix is wrong: Node's `DecompressionStream('gzip')` also throws "Trailing junk…" on two members (checked); loop over the members with the JS inflater (it knows its end position) |

## Evidence

- **01/P13 / 12/P11.** `asciiIndexOf` anchors: `v`, `r`, `r`, `a`, `a`; a plain mmCIF has none of the needles, so each
  of the 5 needles scans the whole file. Re-run of `sniff4.mjs` (min of 7, load 4): asciiIndexOf ×5 267 ms,
  `string.indexOf` ×5 153 ms, one `Uint8Array.indexOf` 57 ms, one-pass JS table 277 ms, `Buffer.indexOf` 3 ms. The two
  reports agree on the ratio (1.4–2.5×) but not on the absolute time: 12's `sniff.ts` copies the file
  (`new Uint8Array(bytes)`) inside the timed region and reports rounds 1–3 raw; 01 reports a min. A cheaper exact
  fix neither report gives: the needles share prefixes (`_refln.pdbx_` ×2, `atom_site` + `_`/`.` ×2), so 3 scans give
  the same answers (−40 %). 12's "server variants from the header only" and "decide after `ParseCif`" remain the
  bigger wins.
- **01/P14.** `vpal.mjs` re-run: 188–240 ms, 91–95 MB vs loop 3–17 ms. `parseFloat(\`${n}\`)` is `n` for every finite
  number and ±Infinity. The only difference is −0 → "0" → +0, which leaves `(x − min)` unchanged.
- **01/P15.** v2 kept a `File` per download (`base/.../assets.ts:236`); HEAD keeps `bytes`, and `get file()` makes a
  `File` (a copy) on each access. `snapshots.ts:437` is the only reader of `entry.file`.
- **01/B17.** `gzrss.ts`: a 5 MB text gzipped, cut by 6 bytes, declares 358 MiB. `ungzip` throws `TypeError ""`;
  `arrayBuffers` 6 → 369 MiB, RSS +25 MiB. The report's RSS +2.4–6.9 GB is not reproducible: `new Uint8Array(n)`
  pages stay unmapped until written, and inflating writes ≤ 5 MB. The two allocations (native path, then
  `ungzipBlocks`) are real. Callers `data-source.ts:164`, `fetch.ts:240` and `host.ts:504` skip `MaxDeflateRatio`, as
  stated.
- **01/B18.** `areColumnsEqual` → `areArraysEqual(sa.array, sb.array)`: raw storage, while `value` lowercases or
  uppercases. Only `valueType` is compared, not `transform`.
- **02/B6.** `molecular-surface.ts:489`. At maxRadius 1.7 and probe 1.4, the chunk is 2e6 / (6.2·scale + 2)³ atoms
  (3.6k / 15k / 91k at res 1 / 2 / 8) against v2's 1e5 / (1.7⁹·scale) (843 / 1.7k / 6.7k). `shouldUpdate` is checked
  only between chunks (`:252`, `:437`). Unit visuals (~720 atoms) are one chunk either way.
- **02/B7.** `periodic-box.ts:47` rounds each fractional coordinate on its own. The dodecahedron case above was
  computed by hand: f = (−0.318, −0.318, 0.636) → (−0.318, −0.318, −0.364). The doc comment's claim of exactness is
  false. Fix as stated (27 shifts).
- **02/P14.** Exit condition `|c − mid| < r_c − rInt(1+1e-6) − 1e-4` ⇒ every probe point is strictly inside `c`, so
  `obscured` is true for all of them and no voxel is written. For `c = a`: `dmp > rA − rInt` always, so `a` never
  exits. Not re-timed: the auditor's interleaved CPU figures stand.
- **03/P18.** The mechanism holds even with 2 tasks: each yields only after a full interval of its own chunks, so the
  other's resumed first chunk already exceeds 250 ms. A shared ctx (child tasks with `runInContext` of one update
  task) does not collapse: its `lastUpdatedTime` is reset by whichever chunk updates. `Task.run()` without an observer
  is `SyncRuntimeContext` (`task.ts:39`), so the parses in MVS `annotation-prop.ts:560` and Prism `source.ts:56` take
  `f(MAX)`. A single long main-thread stall costs one 4× shrink, recovered in the next interval.
- **13/P1.** `mp-ab.ts`: HEAD decode vs a copy with `Uint8Array.prototype.slice.call`; deep-equal with byteOffset 0
  and an exact-size buffer for both Uint8Array and Buffer inputs. The auditor's "loop 82–152 ms" came from a
  standalone replica, larger than the whole decode (62 ms in their own run). `13/msgpack.mjs` re-run: replica loop
  37–237 ms, slice 1.5–9 ms.
- **13/P4.** `ascii.mjs` (`--expose-gc`): now / switch: 5 B 49 / 28 ms, 6 B 74 / 33, 7 B 87 / 34; 8 B equal; 12–24 B
  `apply(subarray)` worse (203–740 ms, 248–352 B), as the report says.
- **13/P5.** `parser.ts:686`: capacity is in rows (column-major), steps 32 → 52 → 84 → ×4 … → estimate. Not
  re-profiled (the 11 % self time was measured while swapping).
- **13/B3.** `binary/parser.ts:38`: `cache[name]` is keyed by the requested name, after a case-insensitive
  `normalizedMap` lookup.
- **13/B5.** `xyz/parser.ts` → `parseFloatBytes`. v2 used `+fields[k]`.
- **13/B1 / 01/B10.** In Node, `new Blob([a.gz, b.gz]).stream().pipeThrough(new DecompressionStream('gzip'))` gives
  "Trailing junk found after the end of the compressed stream".
