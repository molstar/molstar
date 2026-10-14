/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The numeric codes of the resolved tables' coded columns: what a column holds, named once, so the resolver that writes
// them and a viewer that reads them say the same.

export const SectionAxis = { blocks: 0, columns: 1 } as const;
export const Layout = { tracks: 0, text: 1 } as const;
export const Arrangement = { side: 0, stacked: 1 } as const;
export const LigandForm = { position: 0, atoms: 1 } as const;
export const Height = { content: 0, box: 1 } as const;
export const RowSizing = { uniform: 0, natural: 1 } as const;
export const Sizing = { row: 0, free: 1 } as const;
