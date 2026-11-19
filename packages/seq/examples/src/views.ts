/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import * as complexBlocks from './complex/blocks.js';
import * as navigators from './navigators.js';
import * as complexPae from './complex/pae.js';
import * as trna from './trna/trna.js';
import { generatedRegistry } from './globins/data.js';
import * as globins from './globins/family.js';
import * as plasmidRing from './plasmid/ring.js';
import * as addressing from './p53/addressing.js';
import * as axis from './p53/axis.js';
import * as bookmarks from './p53/bookmarks.js';
import * as geneByRegion from './p53/gene-by-region.js';
import * as layers from './p53/layers.js';
import * as proteinOnGene from './p53/protein-on-gene.js';
import * as query from './p53/query.js';
import * as story from './p53/story.js';
import * as themes from './p53/themes.js';
import * as text from './p53/text.js';
import * as tracks from './p53/tracks.js';
import { AlignmentFileRuntime } from '@molstar/seqx-core/data/alignment-readers';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { type DomainFile, PlasmidGffFile } from './files/data.js';
import * as p53Entry from './files/p53-entry.js';
import * as plasmidGff from './files/plasmid-gff.js';
import * as tables from './files/tables.js';

// The examples' views by dataset: what a dataset's index lists and a test walks.

export interface ExampleView {
  readonly dataset: string;
  /** The module's name under the dataset (`tracks`). */
  readonly id: string;
  readonly title: string;
  /** What the view shows, a sentence. */
  readonly summary: string;
  build(): Node;
  /**
   * What a viewer of the view needs besides the default registry (a runtime of a format its sources name); answers
   * delayed by `delayMs` (the view's own pacing when absent; 0 in tests: no timer).
   */
  registry?(options?: { readonly delayMs?: number }): readonly SeqRegistryEntry[];
}

export const P53Views: readonly ExampleView[] = [
  {
    dataset: 'p53',
    id: 'tracks',
    title: 'Tracks',
    summary: 'The protein’s letters, its domains in lanes and a conservation signal with two residues missing.',
    build: tracks.build,
  },
  {
    dataset: 'p53',
    id: 'text',
    title: 'Text',
    summary:
      'The protein as text: lines of whole groups of ten wrapped to the width, a number a line and a group, its secondary structure and domains as thin bands under the letters, each left out of a line it is empty on.',
    build: text.build,
  },
  {
    dataset: 'p53',
    id: 'text-centered',
    title: 'Text, centered',
    summary:
      'The text figure centered in at most 640 px, its group numbers inline in the gaps, a number at each end of a line.',
    build: text.centered,
  },
  {
    dataset: 'p53',
    id: 'text-under-tracks',
    title: 'Text under its tracks',
    summary:
      'A tracks overview of the domains and the protein as text under it, in one view: the text section’s own layout.',
    build: text.underTracks,
  },
  {
    dataset: 'p53',
    id: 'protein-on-gene',
    title: 'Protein on its gene',
    summary:
      'The gene’s section follows the protein’s through the codon map (minus strand, three exons), reading along it (`direction: follow`) on its coding pieces: each intron a hairline, three nucleotides a residue, both moving the same way.',
    build: proteinOnGene.build,
  },
  {
    dataset: 'p53',
    id: 'one-axis',
    title: 'Protein and gene on one axis',
    summary:
      'One gene section with the protein’s letters a codon wide, its domains and the gene’s letters: it reads along the protein (`{ along: p53 }`) on its coding pieces; a codon split by an intron is drawn in each part; R175H, the gene’s one variant, is marked as a conflict.',
    build: axis.oneAxis,
  },
  {
    dataset: 'p53',
    id: 'genomic-follower',
    title: 'Genomic follower',
    summary:
      'For contrast: the gene follows its protein in genomic order (`direction: genomic`, introns shown), so on the minus strand it moves against the protein.',
    build: axis.genomicFollower,
  },
  {
    dataset: 'p53',
    id: 'reversed-gene',
    title: 'Reversed gene',
    summary:
      'The gene alone, reversed: positions count down left to right, its letters read on the other strand (the complement), the domains through the codon map; Shift+C collapses its introns.',
    build: axis.reversedGene,
  },
  {
    dataset: 'p53',
    id: 'gene-by-region',
    title: 'Gene by region',
    summary:
      'The gene’s region read region by region: transcripts arrive as dynamic nodes, a track each, and TP53’s bring their domains, shown through each transcript’s codon map.',
    build: geneByRegion.build,
    registry: geneByRegion.registry,
  },
  {
    dataset: 'p53',
    id: 'gene-region-ideogram',
    title: 'Gene region on chr17’s ideogram',
    summary:
      'The gene’s region with chromosome 17’s ideogram above it (GRCh38’s bands, the centromere two triangles): the region, a sliver of the chromosome, marked on it, and a zoom strip under it over the region alone, its window dragged and resized there. Hover a band for its name, the mark for the region.',
    build: navigators.geneRegionIdeogram,
    registry: geneByRegion.registry,
  },
  {
    dataset: 'p53',
    id: 'tracks-strip',
    title: 'Tracks with a strip',
    summary:
      'p53’s tracks with a strip of them small above: drag the window to pan, its grips to zoom, click outside it to centre there, drag outside it to sweep a region; Alt+Left goes back.',
    build: navigators.tracksStrip,
  },
  {
    dataset: 'p53',
    id: 'protein-gene-strip',
    title: 'Protein strip with its gene’s window',
    summary:
      'The protein on its gene with a strip of the protein: its window and its gene follower’s, carried along the codon map; a drag on the gene’s window moves the protein, which leads.',
    build: navigators.proteinGeneStrip,
  },
  {
    dataset: 'p53',
    id: 'themes',
    title: 'Themes',
    summary:
      'The protein’s letters under each scheme that fits it; domains by category with their legend; a confidence as pLDDT’s four bands and as a ramp; two signals on one shared scale; variants by consequence, sized by a log-scaled frequency; secondary structure. Its export draws the legends below.',
    build: themes.build,
  },
  {
    dataset: 'p53',
    id: 'story',
    title: 'A story',
    summary:
      'The tracks as a story of four steps: the whole protein, its DNA-binding domain (the conservation taller), the R175H hotspot (marked by a track of its own, the domains recolored) and the tetramerization domain (the conservation gone); each a view and a patch over the spec.',
    build: story.build,
  },
];

/** Views made for a feature's stories (`Features/…`), over the same datasets. */
export const FeatureViews: readonly ExampleView[] = [
  {
    dataset: 'p53',
    id: 'addressing',
    title: 'Addressing',
    summary: 'The tracks view with refs naming its domains and hotspots, and a text whose links act on addresses.',
    build: addressing.build,
  },
  {
    dataset: 'p53',
    id: 'bookmarks',
    title: 'Bookmarks',
    summary: 'Gene by region with bookmarks: refs on a range and on transcripts that arrive with their region.',
    build: bookmarks.build,
    registry: bookmarks.registry,
  },
  {
    dataset: 'p53',
    id: 'query',
    title: 'Queried',
    summary:
      'The protein with its domains and conservation and the gene reading along it, for the query stories: an address typed or written, its result selected or highlighted in the view.',
    build: query.build,
  },
  {
    dataset: 'globins',
    id: 'query-alignment',
    title: 'Queried alignment',
    summary: 'The globins family for the query stories: rows and columns addressed, a result selected in the view.',
    build: () => globins.build(),
  },
  {
    dataset: 'p53',
    id: 'themes-editing',
    title: 'Editing themes',
    summary:
      'The themes view edited: a scheme previewed while hovered and picked in one patch; the shared scale’s domain set for both signals in one step and undone in one.',
    build: themes.build,
  },
];

/** Views made for the small layers' stories (`Features/Layers`), over p53. */
export const LayerViews: readonly ExampleView[] = [
  {
    dataset: 'p53',
    id: 'ruler',
    title: 'Ruler',
    summary:
      'The ruler on the protein and, reversed on its coding pieces, on its gene: the labels are steps of 1, 2 or 5 × 10ⁿ that never overlap, and the gene’s count down.',
    build: layers.ruler,
  },
  {
    dataset: 'p53',
    id: 'secondary',
    title: 'Secondary structure',
    summary:
      'Helices and strands as blocks, a turn thin; each strand ends in an arrowhead on the piece holding its C-terminal end, on the protein and through the codon map on the gene reading along it.',
    build: layers.secondary,
  },
  {
    dataset: 'p53',
    id: 'variants',
    title: 'Variants',
    summary:
      'A stem and a head a site, the heads of one site stacked (R248Q and R248W, R273H and R273C), on the protein and on the gene; zoom out and they become a count strip.',
    build: layers.variants,
  },
  {
    dataset: 'p53',
    id: 'values-line',
    title: 'Values as a line and an area',
    summary:
      'The conservation as bars, a line (held over half cells at its ends and at the two missing residues, broken there) and an area under it.',
    build: layers.valuesLine,
  },
  {
    dataset: 'p53',
    id: 'ranges-lines',
    title: 'Ranges as lines and arrows',
    summary:
      'The domains in each line style (the dash and the wave run from each domain’s start, so they move with it as you pan) and with arrowheads on their ends.',
    build: layers.rangesLines,
  },
];

/** The Plasmid (a stand-in until its dataset is recorded): the ring with its linear detail, the rotated linear view. */
export const PlasmidViews: readonly ExampleView[] = [
  {
    dataset: 'plasmid',
    id: 'ring',
    title: 'Ring with a linear detail',
    summary:
      'The construct as rings (each strand’s features, the sites, a ruler; labels outside, the name in the hole) and a linear detail following the ring’s window; bla and a primer site cross the origin, each one arc. Drag the window, or pan the detail: the other follows.',
    build: () => plasmidRing.build(),
  },
  {
    dataset: 'plasmid',
    id: 'rotated-linear',
    title: 'Rotated linear view',
    summary:
      'The construct in tracks from inside bla over the origin to past lacZα: a `range` past the origin, two pieces meeting at the seam with no break.',
    build: plasmidRing.rotatedLinear,
  },
  {
    dataset: 'plasmid',
    id: 'detail-strip',
    title: 'Linear detail with a strip',
    summary:
      'The ring and its linear detail with a strip of the whole construct over the detail: drag the ring’s window over the origin and the strip’s window wraps in two pieces.',
    build: navigators.plasmidStrip,
  },
];

/** The Globins (a generated family until the Pfam seed and its tree are recorded): rows, trees, folds and scale. */
export const GlobinsViews: readonly ExampleView[] = [
  {
    dataset: 'globins',
    id: 'family',
    title: 'Family, logo and tree',
    summary:
      'The family’s rows by its tree beside them, its conservation logo above: hover a cell and its column lights in every row; click a tree node to fold its clade (session state); Cmd/Ctrl+click a node selects its rows; the wheel scrolls the rows, `]` and `[` zoom them.',
    build: () => globins.build(),
  },
  {
    dataset: 'globins',
    id: 'folded',
    title: 'Folded clades',
    summary:
      'Two clades folded when the view opens (the layer’s `collapsed`): each one summary row, its consensus opaque by its share, labeled “first +n”, a wedge in the tree. Click a node to unfold it; `clades.keep` writes the folds shown back.',
    build: globins.folded,
  },
  {
    dataset: 'globins',
    id: 'figure',
    title: 'Alignment figure',
    summary: 'The family as an ESPript-style text figure: groups of ten, every row on every line, a ruler above.',
    build: globins.figure,
  },
  {
    dataset: 'globins',
    id: 'scale',
    title: '10,000 rows',
    summary:
      'A generated variant of 10,000 rows × 300 columns, held, by a balanced tree: zoomed out the rows draw from the shared reduction as image tiles; `]` and `[` zoom the rows, the wheel scrolls them.',
    build: () => globins.scale(),
    registry: generatedRegistry,
  },
  {
    dataset: 'globins',
    id: 'scale-navigators',
    title: '10,000 rows with an overview',
    summary:
      'The 10,000 rows with a vertical overview of every row at the right (summarized by the shared reduction) and a miniature of the columns and rows above: drag either window to move the rows.',
    build: navigators.globinsNavigators,
    registry: generatedRegistry,
  },
  {
    dataset: 'globins',
    id: 'lazy',
    title: '10,000,000 rows, lazy',
    summary:
      'A generated variant of ten million rows in data order: only the rows in view are ever made, a row not yet read labeled by its number until its tile arrives.',
    build: () => globins.lazy(),
    registry: generatedRegistry,
  },
];

/** The Complex (a stand-in until its dataset is recorded): its chains and ligands as blocks, three ways. */
export const ComplexViews: readonly ExampleView[] = [
  {
    dataset: 'complex',
    id: 'blocks',
    title: 'Chains as blocks',
    summary:
      'CDK2 and cyclin A in a row of blocks, polymers first, then ATP as a one-position block; the two copies of the complex merged, each block labeled ×2.',
    build: () => complexBlocks.build(),
  },
  {
    dataset: 'complex',
    id: 'blocks-split',
    title: 'Copies split',
    summary: 'The same blocks with the copies split (`copies: each`): a block a chain and a ligand.',
    build: complexBlocks.split,
  },
  {
    dataset: 'complex',
    id: 'blocks-stacked',
    title: 'Chains stacked',
    summary: 'A block a row of the tracks, aligned at their first residue, one zoom (`arrange: stack`).',
    build: complexBlocks.stacked,
  },
  {
    dataset: 'complex',
    id: 'blocks-ring',
    title: 'Chains on a ring',
    summary: 'The blocks on a ring, a block an arc with a gap between, their labels outside with the features’.',
    build: complexBlocks.ring,
  },
  {
    dataset: 'complex',
    id: 'pae',
    title: 'PAE across the chains (synthetic)',
    summary:
      'A section a chain and a synthetic predicted aligned error between them as a grid (low across the interface, generated from a seed): its x follows CDK2’s section, its y cyclin A’s; synthetic contacts at the interface as dots on it (random positions there, not derived from coordinates). Hover a cell: both residues light in their sections; zoom or scroll it: both sections move.',
    build: () => complexPae.build(),
  },
  {
    dataset: 'complex',
    id: 'contacts',
    title: 'Contacts between the chains (synthetic)',
    summary:
      'Synthetic contacts at the interface (random positions there, not derived from coordinates) as arcs from CDK2’s block to cyclin A’s: a pair, both ends lit on hover.',
    build: complexPae.contacts,
  },
];

export const TrnaViews: readonly ExampleView[] = [
  {
    dataset: 'trna',
    id: 'arcs',
    title: 'Base pairs as arcs',
    summary:
      'tRNA-Phe (PDB 1EHZ, chain A): its 27 base pairs as arcs over the letters, the four stems nested, the six tertiary pairs (four of them crossing the stems’ pairs) over them; colored by kind (canonical, wobble, non-canonical); its 14 modified residues drawn as their parent letters, their codes in a track above. Hover an arc: both ends lit; click: the pair selected.',
    build: () => trna.build(),
  },
  {
    dataset: 'trna',
    id: 'split',
    title: 'Crossing pairs below',
    summary:
      'The largest nested set of pairs above a middle line, the pairs that cross it below (`crossing: split`), over the letters.',
    build: trna.split,
  },
  {
    dataset: 'trna',
    id: 'ring',
    title: 'Pairs as chords',
    summary: 'The tRNA on a ring, its pairs as chords through the hole (`side: below`).',
    build: trna.ring,
  },
  {
    dataset: 'trna',
    id: 'dot-plot',
    title: 'Base-pair dot plot',
    summary:
      'The arcs over the letters and the dot plot under them, its axes following their window: above the diagonal the family’s covariation on 1EHZ’s positions (mutual information over the Rfam seed’s columns, computed through 1EHZ’s place in the seed), below it the structure’s pairs mirrored. Hover a cell: a crosshair, both bases and the arcs ending there lit; drag a rectangle to select two ranges, Alt+drag to take a rectangle out.',
    build: trna.dotPlot,
  },
  {
    dataset: 'trna',
    id: 'seed-row',
    title: 'A seed row with its SS_cons pairs',
    summary:
      'The seed row 1EHZ is placed on (the row of highest identity, found by aligning 1EHZ’s parent letters to every row): its letters and the pairs of the seed’s `SS_cons` projected onto it, nested only.',
    build: trna.seedRow,
  },
  {
    dataset: 'trna',
    id: 'seed',
    title: 'Rfam RF00005 seed',
    summary:
      'The seed (954 rows of 118 columns, Stockholm as published) under its conservation logo. The `SS_cons` arc and the link from 1EHZ’s residues wait for a pairs node and a map on an alignment’s columns.',
    build: trna.seed,
    registry: () => [{ kind: 'format-runtime', runtime: AlignmentFileRuntime }],
  },
];

/**
 * The dot plot at 4,000 × 4,000: a generated grid of 16M values, so a view of its own, kept out of the lists the quick
 * tests walk (its bound is the heavy lane's).
 */
export const TrnaScaleView: ExampleView = {
  dataset: 'trna',
  id: 'dot-plot-scale',
  title: 'Dot plot at 4,000 × 4,000',
  summary:
    'A generated 4,000 × 4,000 grid on a generated RNA, the plot alone under its ruler: the frame target’s scale.',
  build: () => trna.scale(),
};

/** What the Files views need besides the default registry: the records runtime (the default viewer spec has it). */
const records = (): SeqRegistryEntry[] => [{ kind: 'format-runtime', runtime: RecordsRuntime }];

/** The table files, each drawing the same six domains: the file, its title and what is special about its form. */
export const TableFiles: readonly (readonly [file: DomainFile, id: string, title: string, form: string])[] = [
  ['domains.tsv', 'table-tsv', 'TSV', 'tab-separated text, a header row'],
  ['domains.csv', 'table-csv', 'CSV', 'comma-separated text, a header row (RFC 4180)'],
  ['domains-aos.json', 'table-json-aos', 'JSON, records', 'an array of objects, a row each'],
  ['domains-soa.json', 'table-json-soa', 'JSON, columns', 'an object of columns, one array a field'],
  ['domains-aos.msgpack', 'table-msgpack-aos', 'MessagePack, records', 'an array of maps, a row each'],
  ['domains-soa.msgpack', 'table-msgpack-soa', 'MessagePack, columns', 'a map of columns, one array a field'],
];

/** Views made for the Files stories (`Features/Files`): sources read from files. */
export const FilesViews: readonly ExampleView[] = [
  {
    dataset: 'files',
    id: 'uniprot',
    title: 'UniProt’s entry',
    summary:
      'p53 from UniProt’s own JSON (the REST entry of P04637, trimmed): the letters and the features read by the records runtime through the source’s transform, a track a feature type; R175, R248 and R273 lie inside the DNA-binding domain.',
    build: () => p53Entry.build(),
    registry: records,
  },
  ...TableFiles.map(
    ([file, id, title, form]): ExampleView => ({
      dataset: 'files',
      id,
      title: `Domains from ${title}`,
      summary: `p53’s six domains from a table file: ${form}. Every form of the table draws the same.`,
      build: () => tables.build(file),
      registry: records,
    }),
  ),
  {
    dataset: 'files',
    id: 'plasmid-gff3',
    title: 'Plasmid features from GFF3',
    summary:
      'The plasmid’s features as GFF3 (nine columns, a region row saying it is circular), as a ring and its detail: bla and a primer site end past the length in the file and are drawn as one item across the origin.',
    build: () => plasmidGff.build(PlasmidGffFile),
    registry: records,
  },
];

export const ExampleViews: readonly ExampleView[] = [
  ...P53Views,
  ...FilesViews,
  ...FeatureViews,
  ...LayerViews,
  ...PlasmidViews,
  ...ComplexViews,
  ...GlobinsViews,
  ...TrnaViews,
];
