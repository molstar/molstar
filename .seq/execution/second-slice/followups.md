# Followups

- (M) `maxBytes` on `file://` reads in Node (`fetchData`'s file path ignores the four options).
- (M) The asset conformance suite sees a manager through `fetch`; a manager that downloads by other means is outside it.
- (L) The real p53 fixtures round-tripped JSON ↔ MessagePack in a seq package test (io cannot import them).
- (L) `repr.codec.msgpack.syntax` cataloged in core's `CodecIssues` beside the JSON and S-expression codes.
- (L) Duplicate map keys in MessagePack keep the last silently (JSON reports them): an option in `common/msgpack` later.
- (L) Payload alignment is by position: pad in `common/msgpack` if zero-copy decoding ever matters more.
- (L) `checkCodec`'s trees in core gain a leading-U+FEFF string, key and category entry; the codec's issue codes
  (`syntax`, `duplicate_key`, `key_type`) move into core's `CodecIssues`; the Codec's documentation says `encode` throws on
  a tree no reader would give back (a lone surrogate).
- (L) Encoder padding for dependable zero-copy; a canonical-category fast path; key-order determinism (with JSON).
- (K) Derive the next revision's ref index at the end of each apply: a structural edit on a 20,000-ref spec costs
  about 10 ms, mostly rebuilding the old root's index (0.7 ms with it cached); the 16 ms edit target.
- (K) `render/export.ts` still requests data for hidden layers.
- (K) `repr/tree/address.ts` (the child index) is reachable through core's `./*`; io imports `repr/tree/index` by
  path, so closing `./repr/tree/*` needs io moved to `./repr/tree` first.
- (R) Core's other `_test` folders have 73 type errors; `tsconfig.test.json` covers `data/immutable` only.
- The Complex dataset (examples.md, the author 2026-11-18): CDK2 with cyclin A, its interface, PAE and ligands with
  AlphaFill's transplants. Its fixtures and views go to the slice that builds pairs across chains, dot plots and ligands
  (the third slice's draft or later); the entries chosen then and checked against their sources.
- (K2 review) Many moves under one parent are quadratic in `applyPatch` (reversing 16k children: 105 ms alone, 229 ms
  in a merge); inserts and diffs are linear.
- (K2 review) `invalid` discards every change of theirs from `partial`, an unrelated one included (by design, N-8);
  lane N's UI should say so.
- (R+K) Law 3 fails with stored patches in a tree: the patch fuzz's derived mode, seed 3, case 572 (`fuzz(600, 3, true)`),
  the inverse refused as `repr.patch.dangling` with `inverseNeeds`. Skipped there with a comment; the orchestrator fixes it.
- (R+K) `RefSite.path` holds child indices, so an insert in a wide parent re-creates its later siblings' sites; nothing
  reads `path`: drop it or make it lazy.
- (O) An item's coordinate system as a prefix of location forms (`(item :id k)` gives elements): a `(system :item k)`
  form recommended. Refs in data nodes' params and in a section's `of`/`follow` are refused for now.
- (O) Declared numberings: the way back prints `ordinal` or `index` only until the language declares others.
- The frame lane's dropped-tick leak is fixed in `@molstar/ui/dom` (`frame.ts`, the scheduler remembered); lane O's
  twins can drop their per-test `DomHost.requestAnimationFrame` swap.
- `packages/sequences`' `stream/tiles/server.test.ts` fails under `--sequence.shuffle` (seed 9), on the parent too.
