/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqViewSpec, described as data on `@molstar/core/repr` (_spec/README.md, Language): the first slice's kinds. What
// the description cannot say (a sequence has letters or a source, a codon map links a protein to a gene) is the checks
// of `check.ts`, in the `seq.lang.*` codes.

import { patchLanguage } from '@molstar/core/repr/patch';
import { defineLanguage, type Language, type ParamType, type ValidateOptions } from '@molstar/core/repr/schema';
import { ScaleParamsDescription, ThemeLanguage, ThemeLanguageId } from './theme-language.js';

/**
 * The version of the language this description is (the package's `LanguageVersion` is the same). Pre-release (`0.x`):
 * a change moves it with no migration (_spec/README.md, Versions).
 */
export const SeqLanguageVersion = '0.4.0';

/**
 * A channel slot (_spec/README.md "Themes", Forms): a string (the `uniform` form: a color, a token, a word), a number
 * (a uniform size or opacity), or a node of `seq.theme`, read and printed in place.
 */
const channel = (doc: string) =>
  ({
    type: 'union',
    of: [{ type: 'string' }, { type: 'number' }, { type: 'node', language: ThemeLanguageId }],
    doc,
  }) as const satisfies ParamType & { doc: string };

export const SeqLanguage = defineLanguage({
  id: 'seq',
  version: SeqLanguageVersion,
  root: 'spec',
  doc: 'SeqViewSpec: what data there is (sequences, annotations, the maps between them) and how to show it (sections of tracks of layers).',
  kinds: {
    spec: {
      doc: 'The root of a spec.',
      parents: [],
      params: {
        version: {
          type: 'string',
          pattern: '^\\d+\\.\\d+\\.\\d+$',
          required: true,
          doc: 'The language version the spec is written in, `major.minor.patch`.',
        },
        title: { type: 'string', doc: 'A title for people.' },
      },
    },
    data: {
      doc: 'The container of the data of a spec; its id is `data`.',
      parents: ['spec'],
      ids: 'required',
    },
    views: {
      doc: 'The container of the views of a spec; its id is `views`.',
      parents: ['spec'],
      ids: 'required',
    },
    refs: {
      doc: 'The container of the user-defined refs of a spec (bookmarks, regions, sets); its id is `refs`.',
      parents: ['spec'],
      ids: 'required',
    },
    ref: {
      doc: 'A user-defined name for what an address names: `/refs/<id>` stands wherever an address does (a bookmark, a region, a set).',
      parents: ['refs'],
      ids: 'required',
      params: {
        to: {
          type: 'union',
          of: [
            { type: 'ref', doc: 'The address it names: a node, a node data makes, positions, items or a rect.' },
            {
              type: 'list',
              of: { type: 'ref', doc: 'An address of the set.' },
              'min-length': 1,
              doc: 'A set: what each address names, together.',
            },
          ],
          required: true,
          doc: 'What the ref names: an address, or a set of addresses (one type of result).',
        },
        label: { type: 'string', doc: 'A label for people (its id when absent).' },
        color: { type: 'string', doc: 'A color it is shown with: a hex literal or `token:<name>`.' },
      },
    },
    snapshots: {
      doc: 'The container of the story steps of a spec; its id is `snapshots`.',
      parents: ['spec'],
      ids: 'required',
    },
    snapshot: {
      doc: 'A story step (`/snapshots/<id>`): a full state, or a stored patch with no base shown over the spec (or over the step it starts `from`); entering it is navigation, never an edit.',
      parents: ['snapshots'],
      ids: 'required',
      params: {
        label: { type: 'string', doc: 'A label for people (its id when absent).' },
        patch: {
          type: 'node',
          language: 'seq.patch',
          doc: 'What the step changes: a patch of this language with no `base`, its addresses refs of the spec (one of `patch` and `state`).',
        },
        state: {
          type: 'node',
          language: 'seq',
          kinds: ['spec'],
          doc: 'The step as a full state: a spec of its own, its refs resolved inside it; it does not follow edits of the spec (one of `patch` and `state`).',
        },
        from: {
          type: 'ref',
          kinds: ['snapshot'],
          doc: 'The step a patch applies over (a delta), instead of the spec.',
        },
        keep: {
          type: 'string',
          enum: ['full', 'spec'],
          doc: 'How the story’s encoder stores the step whatever is smaller: `full` a frozen full state, `spec` a patch that follows the spec.',
        },
        view: {
          type: 'record',
          fields: {
            view: { type: 'ref', kinds: ['view'], doc: 'The view shown (the first when absent).' },
            sections: {
              type: 'map',
              of: {
                type: 'record',
                fields: {
                  start: { type: 'number', required: true, doc: 'The first position shown (a text’s anchor).' },
                  end: { type: 'number', required: true, doc: 'The end of the window, a position.' },
                  top: { type: 'number', min: 0, doc: 'The vertical scroll, px.' },
                  scale: { type: 'number', min: 0, doc: 'px a position (a text’s letter size).' },
                  bands: {
                    type: 'map',
                    of: {
                      type: 'record',
                      fields: {
                        bandTop: { type: 'number', min: 0, doc: 'The band’s top row, in display rows.' },
                        rowScale: { type: 'number', min: 0, doc: 'px a row.' },
                      },
                    },
                    doc: 'A band that scrolls its own rows, by track address.',
                  },
                  y: {
                    type: 'record',
                    fields: {
                      start: { type: 'number', required: true, doc: 'Its first position.' },
                      end: { type: 'number', required: true, doc: 'Its end, a position.' },
                    },
                    doc: 'A grid’s second axis, as positions.',
                  },
                },
              },
              doc: 'Each section’s viewport by section address, as positions (never offsets).',
            },
          },
          doc: 'What the step shows: the view and each section’s viewport; a story names it on purpose, as data of the step.',
        },
        duration: { type: 'number', min: 0, doc: 'The transition into the step, ms (`story.transition.duration`).' },
        easing: {
          type: 'string',
          enum: ['linear', 'ease-in', 'ease-out', 'ease-in-out'],
          doc: 'How the transition’s time runs (`story.transition.easing`).',
        },
        linger: { type: 'number', min: 0, doc: 'How long the step is shown, ms (`story.linger`).' },
      },
    },
    source: {
      doc: 'Where data is read from: a file or service by uri, or a value inline.',
      parents: ['data'],
      ids: 'required',
      params: {
        format: {
          type: 'string',
          pattern: '^[a-z][a-z0-9_-]*$',
          required: true,
          doc: 'The format the source is read as (`fasta`, `bed`, `bigwig`).',
        },
        uri: { type: 'string', doc: 'Where it is: a url or a path.' },
        index: { type: 'string', doc: 'Where its index is, for formats read by region.' },
        value: { type: 'any', doc: 'The data inline, for small values (instead of `uri`).' },
        fields: {
          type: 'map',
          of: { type: 'string', doc: 'The name the source gives the field.' },
          doc: 'Which of the columns of the source is which field.',
        },
        assembly: { type: 'string', doc: 'The assembly a genomic source is on (`GRCh38`).' },
      },
    },
    sequence: {
      doc: 'A sequence: letters of an alphabet, inline or from a source; a coordinate system.',
      parents: ['data'],
      ids: 'required',
      params: {
        alphabet: {
          type: 'string',
          enum: ['protein', 'dna', 'rna'],
          required: true,
          doc: 'The alphabet of its letters.',
        },
        letters: { type: 'string', doc: 'The letters inline (instead of `source`).' },
        source: { type: 'ref', kinds: ['source'], doc: 'The source that holds the letters (instead of `letters`).' },
        cyclic: { type: 'boolean', default: false, doc: 'Whether the sequence wraps around.' },
        name: {
          type: 'string',
          doc: 'The name its source knows it by (a FASTA record, a chromosome); the node’s id when absent.',
        },
      },
    },
    ranges: {
      doc: 'Ranges on a sequence, 0-based and half-open: an annotation. On a cyclic sequence a range may start after it ends: it wraps the origin.',
      parents: ['data'],
      ids: 'required',
      params: {
        target: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence the ranges are on.' },
        rows: {
          type: 'table',
          columns: {
            start: { type: 'column', of: 'int32', required: true, doc: 'The first position of each range.' },
            end: { type: 'column', of: 'int32', required: true, doc: 'The position after the last, of each range.' },
            id: { type: 'column', of: 'category', doc: 'An identifier of each range.' },
            label: { type: 'list', of: { type: 'string', doc: 'A label.' }, doc: 'A label of each range.' },
            category: { type: 'column', of: 'category', doc: 'A category of each range (a theme colors by it).' },
            strand: { type: 'column', of: 'category', doc: 'The strand of each range: `+` or `-`.' },
          },
          doc: 'The ranges inline, a column for each field (instead of `source`).',
        },
        source: { type: 'ref', kinds: ['source'], doc: 'The source that holds the ranges (instead of `rows`).' },
        map: {
          type: 'string',
          enum: ['codon'],
          doc: '`codon`: each item (its rows, grouped by `id`) is a coding sequence on the target, in genomic order with its `strand`: the item gets a coordinate system of its residues and a codon map onto the target, made when something names it.',
        },
      },
    },
    values: {
      doc: 'A number per position of a sequence: a signal.',
      parents: ['data'],
      ids: 'required',
      params: {
        target: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence the values are on.' },
        values: {
          type: 'column',
          of: 'float32',
          doc: 'One value per position, inline (instead of `source`); NaN is missing.',
        },
        source: {
          type: 'ref',
          kinds: ['source'],
          doc: 'The source that holds the values, read by region (instead of `values`).',
        },
      },
    },
    alignment: {
      doc: 'An alignment: rows (a coordinate system of its own, in data order) over columns (the section’s system), inline or from a source.',
      parents: ['data'],
      ids: 'required',
      params: {
        alphabet: {
          type: 'string',
          enum: ['protein', 'dna', 'rna'],
          required: true,
          doc: 'The alphabet of its letters.',
        },
        rows: {
          type: 'table',
          columns: {
            label: {
              type: 'list',
              of: { type: 'string', doc: 'A label.' },
              required: true,
              doc: 'Each row’s label (its key): `NAME/27-137` numbers its first residue 27 (EC-39).',
            },
            letters: {
              type: 'list',
              of: { type: 'string', doc: 'A row’s aligned letters.' },
              required: true,
              doc: 'Each row’s aligned letters, every row as wide; gaps `-` and `.`.',
            },
          },
          doc: 'The rows inline (instead of `source`).',
        },
        source: {
          type: 'ref',
          kinds: ['source'],
          doc: 'The source that holds the alignment (an aligned FASTA or Stockholm file, or one serving rows by region), instead of `rows`.',
        },
        reference: {
          type: 'string',
          doc: 'The reference row’s label: what `identity` and `matches-reference` compare with (none when absent).',
        },
      },
    },
    tree: {
      doc: 'A tree over an alignment’s rows (Newick): its leaves match rows by label, else by the label before its last `/` (EC-43).',
      parents: ['data'],
      ids: 'required',
      params: {
        of: { type: 'ref', kinds: ['alignment'], required: true, doc: 'The alignment whose rows it orders.' },
        newick: { type: 'string', required: true, doc: 'The tree as Newick text.' },
      },
    },
    map: {
      doc: 'A map from the positions of one sequence to those of another, so data on one shows on the other.',
      parents: ['data'],
      ids: 'required',
      params: {
        from: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence the map starts on.' },
        to: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence the map ends on.' },
        method: {
          type: 'string',
          enum: ['offset', 'codon'],
          required: true,
          doc: '`offset`: a position plus `by`; `codon`: a protein on its gene through `exons`.',
        },
        by: { type: 'integer', doc: 'For `offset`: what is added to a position of `from`.' },
        exons: {
          type: 'table',
          columns: {
            start: { type: 'column', of: 'int32', required: true, doc: 'The first position of each exon on `to`.' },
            end: { type: 'column', of: 'int32', required: true, doc: 'The position after the last, of each exon.' },
          },
          doc: 'For `codon`: the exons on `to`, in genomic order and not overlapping; `strand` says which way the protein runs. On a cyclic `to`, the last exon may end before it starts: it crosses the origin.',
        },
        strand: { type: 'string', enum: ['+', '-'], doc: 'For `codon`: the strand of the gene on `to`.' },
        table: {
          type: 'integer',
          min: 1,
          doc: 'For `codon`: the NCBI genetic code the codons are read with; 1 (the standard code) when absent.',
        },
      },
    },
    view: {
      doc: 'A view: sections of tracks of layers, laid out as `layout` says.',
      parents: ['views'],
      ids: 'required',
      params: {
        layout: {
          type: 'string',
          enum: ['tracks', 'text', 'rings'],
          default: 'tracks',
          doc: 'How the sections are laid out (a section’s own `layout` before it): `tracks` (rows along one axis), `text` (lines of letters) or `rings` (concentric rings, the first track outermost).',
        },
        export: {
          type: 'record',
          fields: {
            legend: {
              type: 'string',
              enum: ['none', 'below', 'right'],
              doc: 'Legends in an export: `none` (when absent), in a band `below` the data or to its `right`.',
            },
            window: {
              type: 'boolean',
              doc: 'Whether an export draws a ring’s window (false when absent: an exported map is a figure).',
            },
          },
          doc: 'How the view exports (a figure says how it exports); an export call’s options override it.',
        },
      },
    },
    scale: {
      doc: 'A shared scale: layers whose themes name it share its domain and palette (an open end: the union over their data in view).',
      parents: ['views'],
      ids: 'required',
      params: ScaleParamsDescription,
    },
    section: {
      doc: 'A section of a view: a region of one coordinate system (or its `block` children), with its tracks.',
      parents: ['view'],
      ids: 'required',
      params: {
        of: {
          type: 'ref',
          kinds: ['sequence', 'alignment'],
          doc: 'The coordinate system the section shows (an alignment: its columns); required unless the section has `block` children.',
        },
        range: {
          type: 'list',
          of: { type: 'integer', min: 0, doc: 'A position.' },
          'min-length': 2,
          'max-length': 2,
          doc: 'The region shown, `[start, end)`; all of the system when absent.',
        },
        follow: {
          type: 'ref',
          kinds: ['section'],
          doc: 'A section whose region this one shows, mapped onto its own system.',
        },
        direction: {
          type: 'union',
          of: [
            {
              type: 'string',
              enum: ['genomic', 'reversed', 'follow'],
            },
            {
              type: 'record',
              fields: {
                along: {
                  type: 'ref',
                  kinds: ['sequence', 'map'],
                  required: true,
                  doc: 'The system (or a map from it) the section reads along.',
                },
              },
            },
          ],
          doc: 'Which way positions run left to right: `genomic` (up), `reversed` (down), `follow` (the way the leader reads, through the follow path), or `{ along: <system or map> }` (the way that system reads here). `follow` for a section that follows another, else `genomic`, when absent; never inferred from what the section draws.',
        },
        pieces: {
          type: 'union',
          of: [
            { type: 'string', enum: ['whole', 'coding'] },
            { type: 'ref', kinds: ['ranges', 'map'] },
            {
              type: 'list',
              of: { type: 'ref', kinds: ['ranges', 'map'], doc: 'A ranges node or a codon map.' },
              'min-length': 1,
            },
          ],
          doc: 'The positions shown, the rest hidden as breaks: `whole` (the system, or its `range`), `coding` (the coding runs of the codon map the section follows or reads along), a ranges node or a codon map, or a list of them (their union); intersected with `range`. `whole` when absent.',
        },
        collapsed: {
          type: 'union',
          of: [
            { type: 'string', enum: ['coding'] },
            { type: 'ref', kinds: ['ranges', 'map'] },
            {
              type: 'list',
              of: { type: 'ref', kinds: ['ranges', 'map'], doc: 'A ranges node or a codon map.' },
              'min-length': 1,
            },
          ],
          doc: 'The `pieces` the `section.collapse` command switches to from `whole` (`coding` when absent).',
        },
        pad: {
          type: 'integer',
          min: 0,
          doc: 'Positions kept each side of every piece (splice sites); overlapping pads merge. 0 when absent.',
        },
        layout: {
          type: 'string',
          enum: ['tracks', 'text', 'rings'],
          doc: 'How the section is laid out: `tracks`, `text` (lines of letters, wrapped to the width) or `rings` (its tracks as concentric rings round its whole axis); the view’s `layout` when absent.',
        },
        group: {
          type: 'integer',
          min: 0,
          doc: 'Text: positions a group (a gap between groups; a line holds whole groups); 10 when absent, 0 for none.',
        },
        numbers: {
          type: 'string',
          enum: ['start', 'end', 'both', 'none'],
          doc: 'Text: each line’s first (`start`, left) or last (`end`, right) position’s number, or both, or none; `start` when absent.',
        },
        group_labels: {
          type: 'string',
          enum: ['above', 'inline', 'none'],
          doc: 'Text: a number a group, above its last position (`above`, when absent), before the group on its line (`inline`), or none.',
        },
        align: {
          type: 'string',
          enum: ['start', 'center'],
          doc: 'Text: lines at the start of the width (when absent) or centered in it.',
        },
        max_width: { type: 'number', min: 0, doc: 'Text: the most px a line takes; the whole width when absent.' },
        height: {
          type: 'number',
          min: 0,
          doc: 'Text: the section’s height in px (it scrolls inside); its content’s, up to the render config’s `text.maxHeight`, when absent. Rings: the rings’ diameter in px; the width’s, up to `rings.maxDiameter`, when absent.',
        },
        rotation: {
          type: 'number',
          doc: 'Rings: where the origin is, degrees clockwise from the top; 0 when absent.',
        },
        labels: {
          type: 'string',
          enum: ['outside', 'none'],
          doc: 'Rings: the items’ labels placed outside the rings with leaders (`outside`, when absent), or none.',
        },
        arrange: {
          type: 'string',
          enum: ['row', 'stack'],
          doc: 'Blocks in tracks: side by side along one axis (`row`, when absent) or a row of tracks each, aligned at their first position (`stack`).',
        },
        copies: {
          type: 'string',
          enum: ['merge', 'each'],
          doc: 'Blocks that are copies (the same system and range): one block labeled `×N` (`merge`, when absent) or a block each.',
        },
      },
    },
    block: {
      doc: 'A block of a section: a piece of it on a system of its own (a chain of a complex, a ligand); blocks follow one another in order.',
      parents: ['section'],
      ids: 'required',
      params: {
        of: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The coordinate system the block shows.' },
        range: {
          type: 'list',
          of: { type: 'integer', min: 0, doc: 'A position.' },
          'min-length': 2,
          'max-length': 2,
          doc: 'The region shown, `[start, end)` (over a cyclic origin when start is past end); all of the system when absent.',
        },
        label: { type: 'string', doc: 'The block’s label (its id when absent).' },
      },
    },
    track: {
      doc: 'A track of a section: a row of layers.',
      parents: ['section'],
      ids: 'required',
      params: {
        label: { type: 'string', doc: 'The label of the track (its id when absent).' },
        height: { type: 'number', min: 0, doc: 'The height of the track in px.' },
        hidden: {
          type: 'boolean',
          default: false,
          doc: 'Whether the track is hidden: kept in the spec, its row collapsed; its layers keep their own `hidden`.',
        },
        for_each: {
          type: 'record',
          fields: {
            source: {
              type: 'ref',
              kinds: ['source'],
              required: true,
              doc: 'The source whose dynamic nodes the template shows.',
            },
            node_kind: {
              type: 'string',
              enum: ['ranges', 'values'],
              required: true,
              doc: 'The kind of node a track is made for: an annotation, which its layers draw.',
            },
            depth: {
              type: 'integer',
              min: 0,
              doc: 'The nodes of the sources nested this many levels below `source` (0, its own, when absent).',
            },
            order: {
              type: 'string',
              pattern: '^-?[a-z_$][a-z0-9_]*$',
              doc: 'The order of the tracks: `key` (when absent) or a param of the nodes, descending with a leading `-`.',
            },
            filter: {
              type: 'map',
              of: { type: 'string', doc: 'The value the param must have, as written.' },
              doc: 'The params a node must have to get a track, each the value given.',
            },
          },
          doc: 'A template: a track for each dynamic node of a kind a source holds, in order, its layers drawing that node (a layer without `data` draws the node).',
        },
      },
    },
    layer: {
      doc: 'A layer of a track: data drawn as one type of mark.',
      parents: ['track'],
      ids: 'required',
      params: {
        type: {
          type: 'string',
          enum: ['letters', 'ranges', 'values', 'ruler', 'secondary', 'variants', 'frame', 'rows', 'logo'],
          required: true,
          doc: 'What it draws: the letters of the system, ranges, values, the ruler of the section’s axis, secondary structure (ranges whose `kind` field says helix, strand or turn), variants (ranges of one position each), frames (ranges outlined across the tracks `across` names, titled by their labels), an alignment’s rows (`rows`) or its conservation logo (`logo`); `secondary` and `variants` draw `ranges` where a reader lacks them.',
        },
        data: {
          type: 'ref',
          kinds: ['sequence', 'ranges', 'values', 'alignment'],
          doc: 'The data it draws (the letters of the section for `letters`, its axis for `ruler`).',
        },
        color: channel('The marks’ color: a color, a token, `light-dark(a, b)`, or a theme form.'),
        ink: channel('`letters`: the letters’ color (`auto`: follows its fill, 4.5:1), or a theme form.'),
        size: channel('The marks’ size in the layer’s unit (a share of the slot, or px), or a theme form.'),
        opacity: channel('The marks’ opacity, 0 to 1 (never a way to hide: a filter hides), or a theme form.'),
        line: channel(
          '`ranges`: how an item is drawn as a line: `solid`, `dashed`, `dotted`, `wavy`, `double`, `box` or `highlight`; blocks when absent. Any other value is drawn solid and said once (`seq.render.line-style`).',
        ),
        arrow: channel(
          '`ranges`: an arrowhead on the item’s end, start or both (`end`, `start`, `both`), on the piece holding that end, never on a cut one; `none` when absent.',
        ),
        shape: channel('`variants`: the head’s shape (`circle`, `square`, `triangle`, `diamond`).'),
        scale: {
          type: 'ref',
          kinds: ['scale'],
          doc: '`values`: the scale its y reads (shared with the layers naming it).',
        },
        domain: {
          type: 'list',
          of: { type: 'number', doc: 'An end.' },
          'min-length': 2,
          'max-length': 2,
          doc: '`values`: its y domain inline, `[a, b]` (open over the values in view when absent).',
        },
        hidden: { type: 'boolean', default: false, doc: 'Whether the layer is hidden: kept in the spec, not drawn.' },
        complement: {
          type: 'string',
          enum: ['auto', 'always', 'never'],
          doc: 'For `letters` of nucleotides: `auto` (when absent) shows a letter as its complement where it is read against its strand on screen (the axis reversed against its run); `always`, `never`. The stored letter is what search and status read.',
        },
        across: {
          type: 'list',
          of: { type: 'string', doc: 'A track’s id in the layer’s section.' },
          'min-length': 1,
          'max-length': 2,
          doc: 'For `frame`: the tracks its outlines cover, the first and the last (the bands between them too); its own track’s when absent.',
        },
        numbering: {
          type: 'string',
          doc: 'For `ruler`: the numbering it labels by (`ordinal`, `index`, or one the system declares); the section’s default numbering (`ordinal`) when absent.',
        },
        encoding: {
          type: 'string',
          enum: ['bars', 'line', 'area'],
          doc: 'For `values`: `bars` (when absent), `line` (through the cells’ centers, held over half cells at its ends and gaps, broken by missing values) or `area` (the line and the area under it).',
        },
        style: {
          type: 'string',
          enum: ['blocks'],
          doc: 'For `secondary`: the figure’s form; `blocks` (when absent): helices and strands as blocks, a strand ending in an arrowhead.',
        },
        tree: {
          type: 'ref',
          kinds: ['tree'],
          doc: 'For `rows`: the tree drawn beside the rows and read by `order: tree`.',
        },
        order: {
          type: 'union',
          of: [
            { type: 'string', enum: ['data', 'tree'] },
            {
              type: 'record',
              fields: {
                field: { type: 'string', required: true, doc: 'A field of the rows table (`label`, `first`).' },
                direction: { type: 'string', enum: ['ascending', 'descending'], doc: '`ascending` when absent.' },
              },
            },
          ],
          doc: 'For `rows`: the rows’ order: `data` (when absent), `tree` (the tree’s leaf order; clades fold only there), or a field of the rows table (missing last, ties in data order).',
        },
        collapsed: {
          type: 'list',
          of: { type: 'string', doc: 'A clade: its node’s label, else `first..last` (its first and last leaves).' },
          doc: 'For `rows`: the clades folded when the view opens (session state after: `clades.keep` writes the folds shown here).',
        },
        pinned: {
          type: 'list',
          of: { type: 'string', doc: 'A row’s label.' },
          doc: 'For `rows`: rows drawn in a band at the top (and kept in their place in the order, marked).',
        },
        hidden_rows: {
          type: 'list',
          of: { type: 'string', doc: 'A row’s label.' },
          doc: 'For `rows`: rows not shown (a count chip at the header’s foot shows them again).',
        },
        height: {
          type: 'union',
          of: [
            { type: 'string', enum: ['fill'] },
            { type: 'integer', min: 1 },
          ],
          doc: 'For `rows`: the band’s height, `fill` (when absent: the section’s room left) or a number of rows.',
        },
        summary: {
          type: 'string',
          enum: ['consensus', 'bar'],
          doc: 'For `rows`: a folded clade’s cells: `consensus` (when absent: the majority letter, opaque by its share) or `bar` (the non-gap share).',
        },
        numbers: {
          type: 'string',
          enum: ['gutters', 'none'],
          doc: 'For `rows`: each row’s first and last residue numbers in gutters (`gutters`, when absent) or none.',
        },
        tree_width: {
          type: 'number',
          min: 0,
          doc: 'For `rows`: the tree’s width in px (120 when absent: the view’s label column widens by it; 0 hides it).',
        },
        logo_scale: {
          type: 'string',
          enum: ['bits', 'probability'],
          doc: 'For `logo`: letter heights by information (`bits`, when absent: log2(k) − H) or by `probability`.',
        },
        correction: {
          type: 'string',
          enum: ['none', 'small-sample'],
          doc: 'For `logo` in `bits`: the small-sample correction (`none` when absent: exact for given frequencies).',
        },
        gaps: {
          type: 'string',
          enum: ['ignore', 'count'],
          doc: 'For `logo`: gaps left out of the frequencies (`ignore`, when absent) or counted.',
        },
        code: {
          type: 'string',
          enum: ['auto', 'one', 'three'],
          doc: 'For `letters` of residues a codon wide: `auto` (when absent) three-letter codes where the narrowest part in view has room; `one`, `three`.',
        },
      },
    },
  },
});

export type SeqLanguageType = typeof SeqLanguage;

/** The language of the spec's stored patches (a story step's `patch`): `seq.patch`, `repr/patch`'s over this one. */
export const SeqPatchLanguage: Language = patchLanguage(SeqLanguage);

/**
 * The languages the spec's slots name, by id: its stored patches' (_spec/README.md, Language, Story steps) and the
 * theme sublanguage (_spec/README.md "Themes").
 */
export function seqLanguages(id: string): Language | undefined {
  // a full-state story step is a spec of its own (`snapshot.state`)
  if (id === SeqLanguage.id) return SeqLanguage;
  return id === SeqPatchLanguage.id ? SeqPatchLanguage : id === ThemeLanguageId ? ThemeLanguage : undefined;
}

/** What validation, normalization and `applyPatch` take to read a spec's slots (`validate` in `ApplyOptions`). */
export const SeqValidateOptions: ValidateOptions = { languages: seqLanguages };
export type SeqKinds = SeqLanguageType['kinds'];
