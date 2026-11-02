# Bond orders

A [Prism](../../packages/prism/AGENTS.md) viewer with the
[bond-orders extension](../../packages/extensions/bond-orders/README.md): bond-order perception from the 3D coordinates,
registered as a bond provider beside the built-in bonds.

```sh
pnpm dev:examples -- bond-orders   # http://localhost:1338/examples/bond-orders/build/, watching the sources
pnpm -C examples/bond-orders build # esbuild → build/
pnpm -C examples/bond-orders typecheck
```

Built like the other examples (`scripts/esbuild/app.mjs`), as an ES module with lazy chunks, from the workspace sources,
so nothing needs building first.

## The page

[`src/spec.ts`](src/spec.ts) is the default Prism viewer plus the extension's `BondOrders` registry entry (and
ball-and-stick, so the bonds are what you see). [`src/index.ts`](src/index.ts) loads a PDB entry in the PDB format (no
bond orders, CONECT records at most; the mmCIF files state the ligands' orders, which `auto` keeps), 1NCO by default
(the neocarzinostatin chromophore CHR: double, triple and aromatic bonds, a non-aromatic five-membered ring), or 5KRE
(6WG, a covalent inhibitor bonded to SER A 124). The header switches the structure's bonds through the `structure.bonds`
operation, as the Bonds field of the structure's inspector does (Structure panel, pick the structure):

| Button                    | Choice                                          |
| ------------------------- | ----------------------------------------------- |
| Default                   | the built-in bonds                              |
| Perceived (auto)          | `bond-order-perception`, `mode: 'auto'`         |
| Perceived (force compute) | `bond-order-perception`, `mode: 'forceCompute'` |
