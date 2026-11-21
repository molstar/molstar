# Lane FX5: the author's look and the decided gaps

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-5`. Model:
fast. The fix lane's rules ([fx-fix-lane.md](fx-fix-lane.md)) hold: one commit an item, a quick-lane test failing on
the parent for each fix, every modified file's header updated (`@author David Sehnal <david.sehnal@gmail.com>` is the
author line already there; the year range to 2026), hot files (`package.json`s, the lockfile, `biome.json`, `.seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`) never staged: list what you need in the report.
Fetches are allowed for the items below (the author, 2026-11-21: "download what you need for fixtures etc"): each
fixture recorded once into the repo with its source URL, date and licence in a `README.md` or the data file's header
beside it; tests never reach the network; keep fixtures small (trim to what the story needs).

## Items, in order

1. **The extension's sites show the card but no hover emphasis** (the author's look:
   `Features/Extensibility/Toggled At Run Time`, the `demo` extension registered). Hovering a lollipop (S6, S15, …)
   shows the hover card, but the lollipop is not emphasised as a built-in item is (ux-foundations.md, Emphasis: hover
   looks the same everywhere). Find why: the extension layer type's draw path likely emits no hover/selection emphasis
   (or no hit ids the emphasis pass reads). Fix it in the extension contract so every extension layer type gets
   emphasis the way built-ins do (prefer the generic path: the emphasis pass drawing over an item's extent/shape, not a
   per-extension hook, unless the design says otherwise: read `.seq/design/extension-openings.md` and core's
   `_spec/README.md` on extensions and emphasis). Selection too, if it has the same gap. A test that hovers an
   extension item and finds the emphasis in the draw list, failing on the parent. Look in light and dark.
2. **The PAE in AlphaFold DB's colours, on real data** (the author: "for PAE maps, can we use the same color theme as
   alphafolddb does? and perhaps real data too, so it's easier to verify visually?"). Stories: `Examples/Complex/PAE and
   contacts` (`Pae`, `Over The Complex`), `packages/seq/examples/src/complex/{pae,complex-grid,data}.ts`.
   - **Colours**: AlphaFold DB's PAE plot palette (a green ramp: confident dark green at 0 Å to white/light at 30 Å).
     Take the exact stops from AlphaFold DB's own viewer source or EBI's published PAE viewer (fetch allowed) and name
     the source in a comment. A theme palette entry (e.g. `alphafold-pae`) beside `viridis` if themes have a palette
     registry; otherwise explicit colours in the examples. Keep the contacts readable on it (they are white now: pick
     a colour that reads on dark green and on white, in both themes).
   - **Real data**: AlphaFold DB has CDK2 (P24941) and cyclin A2 (P20248) as monomers with their PAE JSON
     (`AF-<acc>-F1-predicted_aligned_error_v<N>.json`, CC-BY 4.0; the confidence files in `fixtures/complex/` are
     already v6). The complex's cross-chain PAE is not in AlphaFold DB: search ModelArchive (CC-BY) for a CDK2–cyclin A2
     model with PAE, and use it if one exists. Rule: **never mix real and synthetic values in one plot without the
     label saying which blocks are which**. If no real complex PAE exists: make the `Pae` story (or a new one) show a
     real monomer PAE (CDK2's, whole chain, in AlphaFold DB's colours so it can be compared with the AFDB page side by
     side) and keep the complex grid synthetic with its "(synthetic)" label; say in the report which you found. The
     fixtures as compact as the formats allow (the PAE as the `seq` format's grid or a trimmed JSON: 298² values for
     CDK2 is fine; 730² is fine too).
3. **T-DM1's CCD ids** (fixes.md, open): DM1 and the SMCC/MCC linker are `UNL` placeholders in
   `packages/seq/examples/src/antibody/adc.ts`. Look them up in the PDB's Chemical Component Dictionary (fetch
   allowed: rcsb or PDBe pages), use the real ids and names, record the source.
4. **The hover card on the structure stories' viewers and compact legends** (novice F1, F22; decided):
   `.seq/execution/fourth-slice/novice-findings.md` says where a first user looked for names and found none (residues in
   the structure link's stories; a compact legend's swatches). Turn the viewer's hover card on in those stories, and
   give compact legend swatches their names on hover (a `title` or the card, whichever `@molstar/ui`'s legend part
   already supports; no new public API past `packages/seq` or `@molstar/ui`).
5. **If time allows**: the viewer's quick-lane folder (`packages/seq/viewer`, 2.78 s against 2 s): fewer jsdom files
   (tests that need no DOM moved off `// @vitest-environment jsdom`). Measure before and after on a quiet machine
   (`uptime` under 3); stop if no gain.

## Checks before the report

`pnpm test:unit packages/seq/`, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`,
`pnpm -s -C packages/seq/{core,viewer,examples} typecheck:test`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`,
`pnpm format:check`, `node .seq/_checks/{ec-ids,closed-work,file-size}.ts`, `pnpm -s -C packages/seq/core size` (bounds:
230 / 830 / 1330 / 940 kB). The look: the stories you touched in a Storybook in your worktree (port 6009:
`pnpm --filter @molstar/storybook exec storybook dev -p 6009 --ci --no-open`, run in the background and stopped after),
light and dark (`&globals=theme:dark`); describe what you saw, the orchestrator looks again.

**Report**: under 300 words: commits per item; checks and results; tests failing on the parent; public surface added;
the fixtures with their sources and licences; what is left open and why; the `.v6/changelogs/` entry's text.
