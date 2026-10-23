/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Glyphs for molecular concepts that no generic icon set covers. Drawn on the same 24x24 grid and in the same
// stroke style as Lucide (2px round strokes, no fill; `Icon` sets the stroke attributes on the root), so they sit
// next to Lucide icons. Keep them simple enough to read at 16-20 px.

export const MolecularIcons = {
  /** Two balls joined by a bond line, bent at a third. */
  'ball-and-stick':
    '<circle cx="5" cy="18" r="2.5"/><circle cx="12" cy="6" r="3"/><circle cx="19" cy="17" r="2.5"/><path d="m6.3 15.8 4.2-7.2"/><path d="m13.6 8.5 4 6.4"/>',
  /** A helix: two phase-shifted strands. */
  cartoon: '<path d="M3 15c3-10 6-10 9 0s6 10 9 0"/><path d="M3 9c3 10 6 10 9 0s6-10 9 0"/>',
  /** A hexagonal ring with a substituent. */
  ligand: '<path d="M10 7l5.2 3v6L10 19l-5.2-3v-6z"/><path d="m15.2 10 3.5-2"/><circle cx="20" cy="7" r="1.5"/>',
  /** Snapshots: a frame holding a molecule (a ring with a substituent), another frame behind it. */
  snapshots:
    '<rect x="7" y="3" width="14" height="14" rx="2"/><path d="M3 7v12a2 2 0 0 0 2 2h12"/><path d="M13 6.4l2.6 1.5v3L13 12.4l-2.6-1.5v-3z"/><path d="m15.6 10.9 1.6 1"/><circle cx="18.1" cy="12.5" r="1.1"/>',
  /** Two distance arrows with end bars. */
  'measure-distance': '<path d="M3 6v12M21 6v12M3 12h18M7 9l-4 3 4 3M17 9l4 3-4 3"/>',
  /** Two rays from a vertex with an arc between them. */
  'measure-angle': '<path d="M4 20h17M4 20 17 6M12 20a8 8 0 0 0-2.54-5.85"/>',
  /** A row of residue blocks with a position ruler above. */
  sequence:
    '<rect x="2" y="9" width="5" height="8" rx="1"/><rect x="9.5" y="9" width="5" height="8" rx="1"/><rect x="17" y="9" width="5" height="8" rx="1"/><path d="M4.5 4v2M12 4v2M19.5 4v2"/>',
  /** A kinemage: two atoms and the contact dots between them, as Probe writes them (the Kinemage tab). */
  kinemage:
    '<circle cx="6.5" cy="15.5" r="4"/><circle cx="17.5" cy="8.5" r="4"/><path d="M10 9h.01M12.5 11.5h.01M15 14h.01M9 12h.01M14.5 17h.01M11.5 15h.01"/>',
  /** An assembly: three subunits, overlapping (the Structure tab). */
  structure: '<circle cx="12" cy="8" r="5"/><circle cx="7.5" cy="15.5" r="5"/><circle cx="16.5" cy="15.5" r="5"/>',
  /** Overlapping, filled-looking spheres. */
  spacefill: '<circle cx="9" cy="10" r="6"/><circle cx="16" cy="15" r="5"/>',
  /** A blob outline with a fold. */
  surface:
    '<path d="M12 3c3 0 4 2 7 3s3 4 2 7-1 5-4 6-5 3-8 2-6-3-6-6 0-5 2-7 4-5 7-5z"/><path d="M9 11.5c1.5-2 4-2 6 0"/>',
  /** Two offset outlines. */
  superpose:
    '<path d="M9.5 3.5l5.63 3.25v6.5L9.5 16.5l-5.63-3.25v-6.5z"/><path d="M14.5 7.5l5.63 3.25v6.5L14.5 20.5l-5.63-3.25v-6.5z"/>',
  /** Stacked frames. */
  trajectory:
    '<rect x="3" y="10" width="12" height="11" rx="2"/><path d="M6 7h10a2 2 0 0 1 2 2v9"/><path d="M9 4h10a2 2 0 0 1 2 2v9"/>',
  /** A box with a molecule crossing its face: a simulation's periodic box. */
  'periodic-box':
    '<path d="M4 7h12v12H4z"/><path d="M4 7l4-4h12v12l-4 4M16 7l4-4"/><circle cx="16" cy="13" r="1.6"/><circle cx="19.5" cy="11" r="1.2"/>',
  /** A wireframe cube around a density contour. */
  volume:
    '<path d="M12 2.5l8.5 4.75v9.5L12 21.5l-8.5-4.75v-9.5z"/><path d="M8.5 12.5c0-2.5 2-4 3.5-4s3.5 1.5 3.5 4-1.5 3-3.5 3-3.5-.5-3.5-3z"/>',
  /**
   * Prism's mark as a delta (Home): a triangle and, inside it, a smaller one in a thicker line, a core (the wordmark's
   * mark turned to point up).
   */
  prism:
    '<path stroke-linejoin="miter" d="M12 3 21.5 19.5h-19z"/><path stroke-width="2.6" stroke-linejoin="miter" d="M12 11.2 14.6 15.7H9.4z"/>',
  /** The Mol* mark: the outline of an eight-armed star (four bars through the centre), sharp corners as in the logo. */
  molstar:
    '<path stroke-linejoin="miter" d="M9.7 1.4h4.6v5.05l3.57-3.57 3.25 3.25-3.57 3.57h5.05v4.6h-5.05l3.57 3.57-3.25 3.25-3.57-3.57v5.05H9.7v-5.05l-3.57 3.57-3.25-3.25 3.57-3.57H1.4V9.7h5.05L2.88 6.13l3.25-3.25 3.57 3.57z"/>',
  // Picking levels (granularity), after stories-untold's pick-level glyphs: one mark per level, simple at 15 px.
  /** One atom. */
  'pick-atom': '<circle cx="12" cy="12" r="4"/>',
  /** A residue: one ring. */
  'pick-residue': '<path d="M12 4.5l6.5 3.75v7.5L12 19.5l-6.5-3.75v-7.5z"/>',
  /** A chain: beads on a string. */
  'pick-chain':
    '<circle cx="5" cy="15" r="2.5"/><circle cx="12" cy="9" r="2.5"/><circle cx="19" cy="15" r="2.5"/><path d="m7 13.3 3-2.6M14 10.7l3 2.6"/>',
  /** An entity: the same chain, many copies. */
  'pick-entity':
    '<circle cx="6" cy="8" r="2"/><circle cx="12" cy="8" r="2"/><circle cx="18" cy="8" r="2"/><circle cx="6" cy="16" r="2"/><circle cx="12" cy="16" r="2"/><circle cx="18" cy="16" r="2"/>',
  /** A whole model: a filled-looking outline. */
  'pick-model': '<rect x="4" y="4" width="16" height="16" rx="4"/><circle cx="12" cy="12" r="3"/>',
} as const;
