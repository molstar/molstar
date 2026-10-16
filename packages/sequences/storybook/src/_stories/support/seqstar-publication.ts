/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text view's publication tracks (`text-view.stories.tsx`): p53 (UniProt P04637, the committed fixture) as text,
// with tracks drawn in full under every line (seqspec:viewer#the-text-view): its Kyte–Doolittle hydrophobicity as a line
// (computed from the sequence), a pLDDT-like confidence as an area in AlphaFold DB's bands (illustrative), its zinc
// sites as glyphs, its regions as labeled ranges, and the DNA-binding domain's secondary structure (illustrative, after
// the core domain's crystal structures), in the `figure` preset (tracks in full, the axes and the labels on). No React
// here (the smoke test imports it).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { DocNode, SeqViewSpecDoc } from '@molstar/seq-core/format';
import p53Json from '../../../../core/_test/fixtures/p53.json?raw';

interface P53 {
  readonly sequence: string;
  /** [type, first, last (1-based, inclusive), description]. */
  readonly features: readonly (readonly [string, number, number, string])[];
}

/** Kyte and Doolittle's scale (1982). */
const KD: Readonly<Record<string, number>> = {
  A: 1.8,
  R: -4.5,
  N: -3.5,
  D: -3.5,
  C: 2.5,
  Q: -3.5,
  E: -3.5,
  G: -0.4,
  H: -3.2,
  I: 4.5,
  L: 3.8,
  K: -3.9,
  M: 1.9,
  F: 2.8,
  P: -1.6,
  S: -0.8,
  T: -0.7,
  W: -0.9,
  Y: -1.3,
  V: 4.2,
};

/** Each position's mean over a centered window of `w` residues (fewer at the ends). */
export function hydrophobicity(seq: string, w = 9): number[] {
  const h = Math.floor(w / 2);
  return Array.from(seq, (_, i) => {
    let sum = 0;
    let n = 0;
    for (let j = Math.max(0, i - h); j <= Math.min(seq.length - 1, i + h); j++) {
      sum += KD[seq[j]] ?? 0;
      n++;
    }
    return Math.round((sum / n) * 100) / 100;
  });
}

/** Illustrative confidence: the structured DNA-binding domain and tetramerization region high, the tails low. */
function confidence(n: number): number[] {
  return Array.from({ length: n }, (_, i) => {
    const p = i + 1;
    const core = p >= 95 && p <= 292 ? 92 : p >= 323 && p <= 356 ? 82 : p <= 60 ? 38 : 55;
    return Math.round(core + 6 * Math.sin(p / 7));
  });
}

/** AlphaFold DB's pLDDT bands. */
const AFDB_BANDS = {
  kind: 'discrete',
  mode: 'absolute',
  colors: [
    ['#FF7D45', 0],
    ['#FFDB13', 50],
    ['#65CBF3', 70],
    ['#0053D6', 90],
  ],
};

/** The DNA-binding domain's main strands and helices (illustrative, 1-based inclusive). */
const SECONDARY: readonly (readonly [number, number, string])[] = [
  [110, 112, 'strand'],
  [124, 127, 'strand'],
  [132, 135, 'strand'],
  [141, 146, 'strand'],
  [156, 164, 'strand'],
  [177, 180, 'helix'],
  [195, 199, 'strand'],
  [204, 207, 'strand'],
  [214, 219, 'strand'],
  [228, 236, 'strand'],
  [251, 258, 'strand'],
  [264, 274, 'strand'],
  [278, 286, 'helix'],
  [326, 333, 'strand'],
  [335, 354, 'helix'],
];

/** The section of the publication document: the text view and the tracks view are two presets of it. */
export const PUBLICATION_SECTION = 'text.main';
/**
 * The annotations the story's panel adds to: two bookmark sets of the section, each declared empty with a ranges track
 * over it, so an annotation is one `set` of its rows (one undo step) and is drawn at once, as an underline lane in the
 * text and as a ranges row in the tracks, each in the row's own `color`. Below is the section's bookmarks set
 * (`<section>.bookmarks`, what `bookmark.add` writes without a `data`), its track after the letters; above is another
 * (`bookmark.add`'s `data`), its track before them: the band order of the section puts it over the letters' line.
 */
export const PUBLICATION_ANNOTATIONS = `${PUBLICATION_SECTION}.bookmarks`;
export const PUBLICATION_ANNOTATIONS_ABOVE = `${PUBLICATION_ANNOTATIONS}.above`;
/**
 * Two more sets: frames (seqspec:viewer#frames, a `frame` track whose span runs from the sequence's track to the
 * secondary structure's, its titles over the letters) and markers above the sequence (a `points` track, a glyph each).
 */
export const PUBLICATION_FRAMES = `${PUBLICATION_SECTION}.frames`;
export const PUBLICATION_MARKERS = `${PUBLICATION_SECTION}.markers`;
/** The Annotations panel's choice (its `places`): the sets by name and kind, the first the panel starts on. */
export const PUBLICATION_PLACES = [
  { data: PUBLICATION_ANNOTATIONS_ABOVE, label: 'Above', kind: 'ranges' },
  { data: PUBLICATION_ANNOTATIONS, label: 'Below', kind: 'ranges' },
  { data: PUBLICATION_FRAMES, label: 'Frame', kind: 'frame' },
  { data: PUBLICATION_MARKERS, label: 'Markers', kind: 'points' },
] as const;

/** The annotations' line style: each row's own `line_style` column (what the panel's picker writes; none is solid). */
const LINE_STYLE = { theme: 'field_line_style', params: { field: 'line_style' } };

/** The section's two presets: `null` is the document's own (`figure`, the text view). */
export type PublicationView = 'text' | 'tracks';
export const PUBLICATION_PRESETS: Readonly<Record<PublicationView, string | null>> = { text: null, tracks: 'tracks' };

/** The rows of an annotations set (default: below) as the shown document holds them. */
export function annotationRows(
  doc: SeqViewSpecDoc | null,
  ref: string = PUBLICATION_ANNOTATIONS,
): readonly { id: string }[] {
  const find = (n: DocNode): DocNode | undefined =>
    n.ref === ref ? n : n.children?.map(find).find((x) => x !== undefined);
  const value = doc ? find(doc.root)?.params?.value : undefined;
  return Array.isArray(value) ? value : [];
}

/**
 * p53 as text with its tracks drawn in full, the axes and the labels on; the same section as tracks is the preset
 * `tracks` (`ctx.ui.setStyle`), and the two annotations tracks (empty until the panel adds to them) are drawn in both.
 */
export function publicationTracksDoc(): SeqViewSpecDoc {
  const p53 = JSON.parse(p53Json) as P53;
  const n = p53.sequence.length;
  const b = createSeqViewSpec({ title: 'p53: publication tracks in the text view' });
  const a = b.data.sequence({ ref: 'p53', label: 'p53 (P04637)', alphabet: 'protein', value: p53.sequence });
  const hyd = b.data.annotation({
    ref: 'hydrophobicity',
    geometry: 'values',
    target: a,
    value: { value: hydrophobicity(p53.sequence) },
  });
  const conf = b.data.annotation({ ref: 'confidence', geometry: 'values', target: a, value: { value: confidence(n) } });
  const zn = p53.features.filter((f) => f[0] === 'Zn');
  const sites = b.data.annotation({
    ref: 'zinc',
    geometry: 'ranges',
    target: a,
    value: { start: zn.map((f) => f[1] - 1), label: zn.map((f) => f[3]) },
  } as never);
  const regions = p53.features.filter((f) => f[0] === 'Region' || f[0] === 'DNA binding' || f[0] === 'Motif');
  const ranges = b.data.annotation({
    ref: 'regions',
    geometry: 'ranges',
    target: a,
    value: {
      start: regions.map((f) => f[1] - 1),
      stop: regions.map((f) => f[2]),
      label: regions.map((f) => f[3]),
      category: regions.map((f) => f[0]),
    },
  } as never);
  const ss = b.data.annotation({
    ref: 'secondary',
    geometry: 'ranges',
    target: a,
    value: {
      start: SECONDARY.map((s) => s[0] - 1),
      stop: SECONDARY.map((s) => s[1]),
      category: SECONDARY.map((s) => s[2]),
    },
  } as never);
  const notes = b.data.annotation({
    ref: PUBLICATION_ANNOTATIONS,
    label: 'Annotations below',
    geometry: 'ranges',
    value: [],
  } as never);
  const above = b.data.annotation({
    ref: PUBLICATION_ANNOTATIONS_ABOVE,
    label: 'Annotations above',
    geometry: 'ranges',
    value: [],
  } as never);
  const frames = b.data.annotation({
    ref: PUBLICATION_FRAMES,
    label: 'Frames',
    geometry: 'ranges',
    value: [],
  } as never);
  const markers = b.data.annotation({
    ref: PUBLICATION_MARKERS,
    label: 'Markers',
    geometry: 'ranges',
    value: [],
  } as never);
  const view = b.view({ ref: 'text', label: 'Text', layout: 'text' });
  const sec = view.section({
    ref: PUBLICATION_SECTION,
    style: 'figure',
  } as never);
  sec.track({ ref: 'frames', label: 'Frames' }).frame({
    data: frames,
    span: { from: 'letters', to: 'ss' },
    color: { theme: 'field_color' },
    line_style: LINE_STYLE,
  } as never);
  sec.track({ ref: 'markers', label: 'Markers', height: 1 }).points({
    data: markers,
    glyph: { theme: 'field_glyph', params: { field: 'glyph' } },
    color: { theme: 'field_color' },
  } as never);
  sec
    .track({ ref: 'annotations-above', label: 'Annotations (above)' })
    .ranges({ data: above, labels: true, color: { theme: 'field_color' }, line_style: LINE_STYLE } as never);
  sec.track({ ref: 'letters', label: 'Sequence' }).letters();
  sec
    .track({ ref: 'annotations-track', label: 'Annotations (below)' })
    .ranges({ data: notes, labels: true, color: { theme: 'field_color' }, line_style: LINE_STYLE } as never);
  sec
    .track({ ref: 'hyd', label: 'Hydrophobicity (Kyte–Doolittle, 9)', height: 2 })
    .values({ data: hyd, encoding: 'line', domain: [-4.5, 4.5] } as never);
  sec.track({ ref: 'conf', label: 'Confidence (pLDDT-like, illustrative)', height: 1.5 }).values({
    data: conf,
    encoding: 'area',
    domain: [0, 100],
    color: { theme: 'field_palette', params: { field: 'value', palette: AFDB_BANDS } },
  } as never);
  sec
    .track({ ref: 'zn', label: 'Zn²⁺ binding', height: 1 })
    .points({ data: sites, color: { theme: 'uniform', params: { color: '#7a5195' } } } as never);
  sec.track({ ref: 'regions-track', label: 'Regions', height: 1 }).ranges({
    data: ranges,
    labels: true,
    color: { theme: 'field_palette', params: { field: 'category', palette: { kind: 'categorical', colors: 'Set2' } } },
  } as never);
  sec.track({ ref: 'ss', label: 'Secondary structure (illustrative)', height: 1 }).secondary({ data: ss });
  return b.getDocument();
}
