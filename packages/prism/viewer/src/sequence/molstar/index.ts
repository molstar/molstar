/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Mol* adapter (`@molstar/prism-viewer/sequence/molstar`; seqspec:libraries#the-mol-adapter, sequence.md §1):
// structures read headless and turned into Seq* complexes (`molstar:structure`), the structure formats and the `pdbe`
// and `afdb` services, the structure providers (pLDDT, PAE, secondary structure, SIFTS coverage), the link between a
// structure and a Seq* set (`StructureLink`, registered as external systems on a `Translator`: core's `SystemRegistry`,
// a host's own with no viewer or a viewer's `ctx.systems`), the viewer-to-3D link (`linkViews`) and SeqQL over Mol* (`molql`). Prism's sequence is made on it.

export * from './chains.js';
export * from './complex.js';
export * from './coverage.js';
export * from './link.js';
export * from './link-views.js';
export * from './pae.js';
export * from './plddt.js';
export * from './query.js';
export * from './registry.js';
export * from './secondary.js';
export * from './seq-units.js';
export * from './sifts.js';
export * from './source.js';
export * from './structure-ref.js';
