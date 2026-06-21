/**
 * Build script for bond-order perception validation data.
 *
 * For each entry in data/manifest.json:
 *   1. Downloads the PDB file from RCSB and extracts the first copy of the ligand.
 *   2. Downloads the CCD CIF for the comp_id and extracts non-single/aromatic bonds.
 *
 * Output files (data/*.pdb, data/*_expected.json) are gitignored — run this script
 * locally before running the 'bond-order perception (data)' Jest suite.
 *
 * Usage: node build-data.mjs
 * Requires: Node >= 18 (native fetch), and a built lib/ directory (npm run build)
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

// eslint-disable-next-line
const __dirname = path.dirname(fileURLToPath(import.meta.url));
const dataDir = path.join(__dirname, 'data');

// Import the Mol* CIF parser and CCD schema from the compiled lib/ directory.
const libRoot = path.resolve(__dirname, '../../../../../../../lib');
const { CIF } = await import(path.join(libRoot, 'mol-io/reader/cif.js'));

const manifest = JSON.parse(fs.readFileSync(path.join(dataDir, 'manifest.json'), 'utf8'));

for (const { pdbId, compId } of manifest) {
    const stem = `${pdbId}_${compId}`;
    const pdbOut = path.join(dataDir, `${stem}.pdb`);
    const expectedOut = path.join(dataDir, `${stem}_expected.json`);
    if (fs.existsSync(pdbOut) && fs.existsSync(expectedOut)) {
        console.log(`Skipping ${pdbId} / ${compId} (already built)`);
        continue;
    }
    console.log(`Processing ${pdbId} / ${compId}...`);
    await buildPdbExcerpt(pdbId, compId);
    await buildExpectedBonds(pdbId, compId);
    console.log(`  done.`);
}

// --- PDB excerpt --------------------------------------------------------------

async function buildPdbExcerpt(pdbId, compId) {
    const url = `https://files.rcsb.org/download/${pdbId}.pdb`;
    console.log(`  Fetching ${url}`);
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.status}`);
    const text = await res.text();

    const lines = text.split('\n');
    let firstSeqNum = null;
    const kept = [];

    for (const line of lines) {
        if (!line.startsWith('HETATM') && !line.startsWith('ATOM  ')) continue;
        // PDB columns (0-indexed): residue name 17-19, chain 21, seq number 22-25
        const resName = line.slice(17, 20).trim();
        if (resName !== compId) continue;
        const seqNum = line.slice(22, 26).trim();
        if (firstSeqNum === null) firstSeqNum = seqNum;
        if (seqNum !== firstSeqNum) continue;
        kept.push(line);
    }

    if (kept.length === 0) throw new Error(`No HETATM records found for ${compId} in ${pdbId}`);
    kept.push('END');

    const outPath = path.join(dataDir, `${pdbId}_${compId}.pdb`);
    fs.writeFileSync(outPath, kept.join('\n') + '\n');
    console.log(`  Wrote ${outPath} (${kept.length - 1} atoms)`);
}

// --- CCD expected bonds -------------------------------------------------------

async function buildExpectedBonds(pdbId, compId) {
    const url = `https://files.rcsb.org/ligands/download/${compId}.cif`;
    console.log(`  Fetching ${url}`);
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.status}`);
    const text = await res.text();

    const bonds = await parseCcdBonds(text, compId);
    const outPath = path.join(dataDir, `${pdbId}_${compId}_expected.json`);
    fs.writeFileSync(outPath, JSON.stringify(bonds, null, 2) + '\n');
    console.log(`  Wrote ${outPath} (${bonds.length} non-single bonds)`);
}

/**
 * Parse the _chem_comp_bond loop from a CCD CIF file using the Mol* CIF parser.
 * Returns only bonds that are non-single or aromatic (the ones perception must assign).
 * When two terminal heavy atoms of the same element are interchangeable partners for a
 * DOUB bond (e.g. nitro O, carboxylate O), the `b` field is an array of both names.
 */
async function parseCcdBonds(cifText, compId) {
    const parsed = await CIF.parseText(cifText).run();
    if (parsed.isError) throw new Error(`CIF parse error for ${compId}: ${parsed.message}`);
    const block = parsed.result.blocks[0];
    if (!block) throw new Error(`No CIF block found for ${compId}`);
    const db = CIF.schema.CCD(block);

    // --- atom elements ---
    const atomElement = new Map();
    const atomCount = db.chem_comp_atom._rowCount;
    for (let i = 0; i < atomCount; i++) {
        atomElement.set(
            db.chem_comp_atom.atom_id.value(i),
            db.chem_comp_atom.type_symbol.value(i).toUpperCase()
        );
    }
    const isH = (name) => { const el = atomElement.get(name); return el === 'H' || el === 'D'; };

    // --- bonds ---
    const bondCount = db.chem_comp_bond._rowCount;
    const { atom_id_1, atom_id_2, value_order, pdbx_aromatic_flag } = db.chem_comp_bond;

    // Heavy-atom degree (ignoring bonds to H/D) — used to detect terminal atoms.
    const heavyDegree = new Map();
    for (let i = 0; i < bondCount; i++) {
        const a = atom_id_1.value(i), b = atom_id_2.value(i);
        if (!isH(a) && !isH(b)) {
            heavyDegree.set(a, (heavyDegree.get(a) ?? 0) + 1);
            heavyDegree.set(b, (heavyDegree.get(b) ?? 0) + 1);
        }
    }
    const isTerminal = (name) => (heavyDegree.get(name) ?? 0) === 1;

    // For each non-H center, collect its terminal heavy neighbours grouped by element.
    // Used to detect interchangeable DOUB/SING pairs (e.g. nitro O, carboxylate O).
    const centerTerminals = new Map();
    for (let i = 0; i < bondCount; i++) {
        const a = atom_id_1.value(i), b = atom_id_2.value(i);
        const order = value_order.value(i), aro = pdbx_aromatic_flag.value(i);
        if (aro === 'y') continue;
        for (const [center, terminal] of [[a, b], [b, a]]) {
            if (isH(center) || isH(terminal) || !isTerminal(terminal)) continue;
            const el = atomElement.get(terminal);
            if (!centerTerminals.has(center)) centerTerminals.set(center, new Map());
            const byEl = centerTerminals.get(center);
            if (!byEl.has(el)) byEl.set(el, { doub: [], sing: [] });
            if (order === 'doub') byEl.get(el).doub.push(terminal);
            else if (order === 'sing') byEl.get(el).sing.push(terminal);
        }
    }
    // A center has interchangeable terminals when it has ≥1 DOUB and ≥1 SING terminal
    // of the same element — the double bond could have been placed on any of them.
    const interchangeable = new Map();
    for (const [, byEl] of centerTerminals) {
        for (const { doub, sing } of byEl.values()) {
            if (doub.length > 0 && sing.length > 0) {
                const all = [...doub, ...sing];
                for (const name of all) interchangeable.set(name, all);
            }
        }
    }

    const result = [];
    for (let i = 0; i < bondCount; i++) {
        const a = atom_id_1.value(i), b = atom_id_2.value(i);
        const order = value_order.value(i), aro = pdbx_aromatic_flag.value(i);
        if (aro === 'y') {
            result.push({ a, b, aromatic: true });
        } else if (order === 'doub') {
            // The interchangeable terminal can be on either side of the bond in the CCD
            // (e.g. phosphate lists the =O as atom_id_1, nitro lists it as atom_id_2).
            // Emit the array on the varying (terminal) side, keeping the center as `a`.
            const equivB = interchangeable.get(b);
            const equivA = interchangeable.get(a);
            if (equivB) result.push({ a, b: equivB, order: 2 });
            else if (equivA) result.push({ a: b, b: equivA, order: 2 });
            else result.push({ a, b, order: 2 });
        } else if (order === 'trip') {
            result.push({ a, b, order: 3 });
        }
    }
    return result;
}
