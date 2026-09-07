#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

/**
 * Build script for bond-order perception validation data (consumed by index.ts).
 *
 * For each entry in <data-dir>/manifest.json:
 *   1. Downloads the PDB file from RCSB and extracts the first copy of the ligand.
 *   2. Downloads the CCD CIF for the comp_id and extracts non-single/aromatic bonds.
 *
 * Output files (<data-dir>/*.pdb, *_expected.json) are gitignored — run this script before the
 * bond-order validation CLI.
 *
 * Usage: node lib/commonjs/cli/bond-order-validation/build-data.js [--data-dir PATH] [--force]
 * Requires: Node >= 18 (native fetch).
 */

import * as argparse from 'argparse';
import * as fs from 'fs';
import * as path from 'path';
import { CIF } from '../../mol-io/reader/cif';

interface ExpectedBond {
    a: string;
    b: string | string[];
    order?: 2 | 3;
    aromatic?: true;
}

async function run(manifestPath: string, force: boolean) {
    const dataDir = path.dirname(manifestPath);
    const manifest: { pdbId: string, compId: string }[] =
        JSON.parse(fs.readFileSync(manifestPath, 'utf8'));

    for (const { pdbId, compId } of manifest) {
        const stem = `${pdbId}_${compId}`;
        const pdbOut = path.join(dataDir, `${stem}.pdb`);
        const expectedOut = path.join(dataDir, `${stem}_expected.json`);
        if (!force && fs.existsSync(pdbOut) && fs.existsSync(expectedOut)) {
            console.log(`Skipping ${pdbId} / ${compId} (already built)`);
            continue;
        }
        console.log(`Processing ${pdbId} / ${compId}...`);
        await buildPdbExcerpt(dataDir, pdbId, compId);
        await buildExpectedBonds(dataDir, pdbId, compId);
        console.log(`  done.`);
    }
}

// --- PDB excerpt --------------------------------------------------------------

async function buildPdbExcerpt(dataDir: string, pdbId: string, compId: string) {
    const url = `https://files.rcsb.org/download/${pdbId}.pdb`;
    console.log(`  Fetching ${url}`);
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.status}`);
    const text = await res.text();

    const lines = text.split('\n');
    let firstSeqNum: string | null = null;
    let firstAltLoc: string | null = null;
    const kept: string[] = [];
    const keptSerials = new Set<string>();

    for (const line of lines) {
        if (!line.startsWith('HETATM') && !line.startsWith('ATOM  ')) continue;
        // PDB columns (0-indexed): residue name 17-19, altLoc 16, chain 21, seq number 22-25
        const resName = line.slice(17, 20).trim();
        if (resName !== compId) continue;
        const seqNum = line.slice(22, 26).trim();
        if (firstSeqNum === null) firstSeqNum = seqNum;
        if (seqNum !== firstSeqNum) continue;
        // Keep a single conformer: blank altLoc plus the first non-blank altLoc seen (PDB lists the
        // highest-occupancy conformer first). Dropping the rest avoids duplicate atoms, which would
        // otherwise be perceived as a second copy and flagged as spurious bonds (e.g. 2RNT / GPG).
        const altLoc = line[16];
        if (altLoc !== ' ') {
            if (firstAltLoc === null) firstAltLoc = altLoc;
            if (altLoc !== firstAltLoc) continue;
        }
        kept.push(line);
        // Atom serial number 6-10; kept lines retain their original serials so CONECT records match.
        keptSerials.add(line.slice(6, 11).trim());
    }

    if (kept.length === 0) throw new Error(`No HETATM records found for ${compId} in ${pdbId}`);

    // Keep CONECT records whose primary atom (serial 6-10) belongs to the excerpt. Bonds to atoms in
    // other residues reference serials absent from the excerpt and are silently ignored by Mol* when
    // the PDB is parsed, so there is no need to filter the bonded-atom columns.
    const conect: string[] = [];
    for (const line of lines) {
        if (!line.startsWith('CONECT')) continue;
        if (keptSerials.has(line.slice(6, 11).trim())) conect.push(line);
    }

    const outLines = [...kept, ...conect, 'END'];
    const outPath = path.join(dataDir, `${pdbId}_${compId}.pdb`);
    fs.writeFileSync(outPath, outLines.join('\n') + '\n');
    console.log(`  Wrote ${outPath} (${kept.length} atoms, ${conect.length} CONECT)`);
}

// --- CCD expected bonds -------------------------------------------------------

async function buildExpectedBonds(dataDir: string, pdbId: string, compId: string) {
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
async function parseCcdBonds(cifText: string, compId: string): Promise<ExpectedBond[]> {
    const parsed = await CIF.parseText(cifText).run();
    if (parsed.isError) throw new Error(`CIF parse error for ${compId}: ${parsed.message}`);
    const block = parsed.result.blocks[0];
    if (!block) throw new Error(`No CIF block found for ${compId}`);
    const db = CIF.schema.CCD(block);

    // --- atom elements ---
    const atomElement = new Map<string, string>();
    const atomCount = db.chem_comp_atom._rowCount;
    for (let i = 0; i < atomCount; i++) {
        atomElement.set(
            db.chem_comp_atom.atom_id.value(i),
            db.chem_comp_atom.type_symbol.value(i).toUpperCase()
        );
    }
    const isH = (name: string) => { const el = atomElement.get(name); return el === 'H' || el === 'D'; };

    // --- bonds ---
    const bondCount = db.chem_comp_bond._rowCount;
    const { atom_id_1, atom_id_2, value_order, pdbx_aromatic_flag } = db.chem_comp_bond;

    // Heavy-atom degree (ignoring bonds to H/D) — used to detect terminal atoms.
    const heavyDegree = new Map<string, number>();
    for (let i = 0; i < bondCount; i++) {
        const a = atom_id_1.value(i), b = atom_id_2.value(i);
        if (!isH(a) && !isH(b)) {
            heavyDegree.set(a, (heavyDegree.get(a) ?? 0) + 1);
            heavyDegree.set(b, (heavyDegree.get(b) ?? 0) + 1);
        }
    }
    const isTerminal = (name: string) => (heavyDegree.get(name) ?? 0) === 1;

    // For each non-H center, collect its terminal heavy neighbours grouped by element.
    // Used to detect interchangeable DOUB/SING pairs (e.g. nitro O, carboxylate O).
    const centerTerminals = new Map<string, Map<string, { doub: string[], sing: string[] }>>();
    for (let i = 0; i < bondCount; i++) {
        const a = atom_id_1.value(i), b = atom_id_2.value(i);
        const order = value_order.value(i), aro = pdbx_aromatic_flag.value(i);
        if (aro === 'y') continue;
        for (const [center, terminal] of [[a, b], [b, a]]) {
            if (isH(center) || isH(terminal) || !isTerminal(terminal)) continue;
            const el = atomElement.get(terminal)!;
            if (!centerTerminals.has(center)) centerTerminals.set(center, new Map());
            const byEl = centerTerminals.get(center)!;
            if (!byEl.has(el)) byEl.set(el, { doub: [], sing: [] });
            if (order === 'doub') byEl.get(el)!.doub.push(terminal);
            else if (order === 'sing') byEl.get(el)!.sing.push(terminal);
        }
    }
    // A center has interchangeable terminals when it has ≥1 DOUB and ≥1 SING terminal
    // of the same element — the double bond could have been placed on any of them.
    const interchangeable = new Map<string, string[]>();
    for (const [, byEl] of centerTerminals) {
        for (const { doub, sing } of byEl.values()) {
            if (doub.length > 0 && sing.length > 0) {
                const all = [...doub, ...sing];
                for (const name of all) interchangeable.set(name, all);
            }
        }
    }

    const result: ExpectedBond[] = [];
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

// --- driver ------------------------------------------------------------------

const parser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Download RCSB/CCD data for bond-order perception validation'
});
parser.add_argument('--manifest-path', '-m', {
    help: 'Path to the manifest JSON file; downloaded *.pdb / *_expected.json are written to the same directory',
    default: path.resolve(process.cwd(), 'src/cli/bond-order-validation/data/manifest.json')
});
parser.add_argument('--force', '-f', {
    action: 'store_true',
    help: 'Rebuild entries even if their data files already exist'
});
interface Args {
    manifest_path: string;
    force: boolean;
}
const args: Args = parser.parse_args();

run(args.manifest_path, args.force).catch(e => {
    console.error(e);
    process.exit(1);
});
