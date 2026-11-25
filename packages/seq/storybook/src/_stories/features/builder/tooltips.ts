/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// Words on rendered marks: `views.tooltip` for every view, `view.tooltip` for one. By data, a node's items speak their
// own fields; by range, a table's rows speak on every mark they cover. Markdown is the template's; values are text.

export function build(): Node {
  const b = createSeqSpec({ title: 'Tooltips' });
  const p53 = b.data.sequence('p53', {
    alphabet: 'protein',
    letters:
      'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD',
  });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: {
      start: [0, 93, 324],
      end: [61, 292, 355],
      label: ['TAD', 'DBD', 'TET'],
      uniprot: ['PF08563', 'PF00870', 'PF07710'],
    },
  });
  b.refs.ref('hotspot', { to: `${p53}/(range :from 175 :to 175)` });
  // by data: each domain's own fields, on its own marks
  b.views.tooltip('domains', { on: domains, label: 'Pfam', text: '`{uniprot}`' });
  // by range: a table of rows on the protein, a row each on every mark it covers (the letters, the domains)
  b.views.tooltip('notes', {
    target: p53,
    rows: { start: [174, 247], end: [175, 248], tooltip: ['R175H', 'R248Q'] },
    label: 'Hotspot',
    text: '**{}**, [show all](!seq-focus=/refs/hotspot)',
  });
  b.views.view('main', (v) => {
    const s = v.section('p53', { of: p53 });
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('domains', { type: 'ranges', data: domains });
    // this view's alone: an address is a one-row table
    v.tooltip('dbd', { on: `${p53}/(range :from 94 :to 292)`, text: 'binds *DNA*' });
  });
  return b.build();
}
