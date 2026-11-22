/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { describe, expect, it } from 'vitest';
import { DefaultOpenLimits } from '@molstar/seqx-core/data';
import {
  droppedUri,
  DropRuntimes,
  DropSamples,
  formatOfName,
  MaxLength,
  refuseSize,
  sniff,
  specOfDropped,
} from '../dropped-file.js';

// A file dropped on a view, for each format built: named, opened by its runtime, its spec shown by a headless viewer whose
// assets port answers the dropped bytes only (nothing is read from anywhere else), laid out with no issue and rendered.

const bytes = (file: string) => new Uint8Array(readFileSync(`${PackageFolder}${file}`));

async function show(name: string, data: Uint8Array) {
  const dropped = await specOfDropped(name, data);
  if ('refused' in dropped)
    return { dropped, reads: [] as string[], issues: [] as string[], sections: [] as string[], svg: '' };
  const base = DefaultSeqViewerSpec();
  const reads: string[] = [];
  const viewer = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [
      ...(base.registry ?? []),
      ...DropRuntimes.map((runtime) => ({ kind: 'format-runtime' as const, runtime })),
    ],
    ports: {
      assets: () => ({
        reader: {
          readText: async (uri) => {
            reads.push(uri);
            return new TextDecoder().decode(data);
          },
          readBytes: async (uri) => {
            reads.push(uri);
            return data;
          },
        },
      }),
    },
  });
  viewer.resize(1000);
  await viewer.load(dropped.spec);
  await viewer.settle();
  const issues = [...viewer.modelIssues.value, ...viewer.layoutIssues.value].map((i) => `${i.code} ${i.address ?? ''}`);
  const sections = viewer.view.resolved.value?.sections.map((s) => s.address) ?? [];
  const svg = await viewer.renderSvg({});
  viewer.retire();
  return { dropped, reads, issues, sections, svg: svg.svg };
}

describe('a dropped file', () => {
  it.each([
    ['2OCJ.seq', 'fixtures/p53/2OCJ.seq', 'seq'],
    ['fab.fasta', 'fixtures/antibody/1N8Z.fasta', 'fasta'],
    ['PF00042.sto', 'fixtures/globins/PF00042-seed.sto', 'stockholm'],
    ['domains.tsv', 'fixtures/files/domains.tsv', 'tsv'],
    ['domains.csv', 'fixtures/files/domains.csv', 'csv'],
    ['domains.json', 'fixtures/files/domains-aos.json', 'json'],
    ['2ocj.cif.gz', 'fixtures/structures/2OCJ.cif.gz', 'mmcif'],
    ['pUC19.gb', 'fixtures/plasmid/pUC19c-rotated.gb', 'genbank'],
  ])('%s is read as %s: shown with no issue, read from the dropped bytes only', async (name, file, format) => {
    const r = await show(name, bytes(file));
    expect('refused' in r.dropped, JSON.stringify(r.dropped)).toBe(false);
    if ('refused' in r.dropped) return;
    expect(r.dropped.format).toBe(format);
    expect(r.issues).toEqual([]);
    expect(r.sections.length).toBeGreaterThan(0);
    expect(r.svg.startsWith('<svg')).toBe(true);
    // every location the viewer read is the dropped file's
    expect([...new Set(r.reads)]).toEqual([droppedUri(name)]);
  });

  it('an unknown extension names what is read; a file that does not parse is refused with the reason, never thrown (EC-77)', async () => {
    expect(
      'refused' in formatOfName('notes.txt') && (formatOfName('notes.txt') as { refused: string }).refused,
    ).toContain('.seq');
    const bad = await specOfDropped('broken.seq', new TextEncoder().encode('this is not a seq file\n\u0000'));
    expect('refused' in bad && bad.refused).toContain('broken.seq');
    const empty = await specOfDropped('empty.tsv', new TextEncoder().encode('a\tb\n'));
    expect('refused' in empty).toBe(true);
  });

  it.each(DropSamples.map(([name, file]) => [name, file]))(
    'the "Try %s" button’s sample loads a view: shown with no issue',
    async (name, file) => {
      const r = await show(name, bytes(file));
      expect('refused' in r.dropped, JSON.stringify(r.dropped)).toBe(false);
      expect(r.issues).toEqual([]);
      expect(r.sections.length).toBeGreaterThan(0);
    },
  );

  it('pUC19.gb: one cyclic sequence, its letters and its features (bla one item across the origin)', async () => {
    const dropped = await specOfDropped('pUC19.gb', bytes('fixtures/plasmid/pUC19c-rotated.gb'));
    if ('refused' in dropped) throw new Error(dropped.refused);
    expect(dropped.summary).toBe('genbank: 1 sequences');
    const printed = JSON.stringify(dropped.spec);
    expect(printed).toContain('"cyclic":true');
    const r = await show('pUC19.gb', bytes('fixtures/plasmid/pUC19c-rotated.gb'));
    expect(r.svg).toContain('bla (AmpR)');
  });

  it('a name with characters of a URL is one location (EC-200: names cover real names)', () => {
    expect(droppedUri('my file #1.seq')).toBe('dropped/my%20file%20%231.seq');
  });

  it('a table whose row ends at a million rows reads without a stack overflow; one ending past the cap is refused, said', async () => {
    const n = 300_000;
    const rows = ['start\tend'];
    for (let i = 0; i < n; i++) rows.push(`${i + 1}\t${i + 1}`);
    const big = await specOfDropped('big.tsv', new TextEncoder().encode(rows.join('\n')));
    expect('refused' in big, JSON.stringify('refused' in big && big.refused)).toBe(false);
    const far = await specOfDropped('far.tsv', new TextEncoder().encode(`start\tend\n1\t${MaxLength + 1}\n`));
    expect('refused' in far && far.refused).toContain(String(MaxLength));
  });

  it('a file past the store’s bound is refused before a byte is read', () => {
    expect(refuseSize('big.cif', DefaultOpenLimits.maxBytes + 1)).toContain('big.cif');
    expect(refuseSize('ok.cif', 1000)).toBeUndefined();
  });

  it.each([
    ['a.sto', 'not a stockholm file\n', '# STOCKHOLM'],
    ['a.fasta', 'ACGT\n', '>'],
    ['a.json', 'name,x\n', '{'],
    ['a.cif', 'loop_\n', 'data_'],
    ['a.gb', '>not genbank\n', 'LOCUS'],
  ])(
    '%s with the wrong content is refused by its first bytes, saying what it should start with',
    async (name, text, want) => {
      const r = await specOfDropped(name, new TextEncoder().encode(text));
      expect('refused' in r && r.refused).toContain(want);
    },
  );

  it('gzip is read only as mmCIF and GenBank', async () => {
    const gz = new Uint8Array([0x1f, 0x8b, 8, 0]);
    expect(sniff('a.tsv', 'tsv', gz)).toContain('gzipped');
    expect(sniff('a.cif.gz', 'mmcif', gz)).toBeUndefined();
    expect(sniff('a.gb.gz', 'genbank', gz)).toBeUndefined();
  });
});
