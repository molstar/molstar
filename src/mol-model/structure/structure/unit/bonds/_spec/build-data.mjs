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
 * Requires: Node >= 18 (native fetch)
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

// eslint-disable-next-line
const __dirname = path.dirname(fileURLToPath(import.meta.url));
const dataDir = path.join(__dirname, 'data');

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

    const bonds = parseCcdBonds(text, compId);
    const outPath = path.join(dataDir, `${pdbId}_${compId}_expected.json`);
    fs.writeFileSync(outPath, JSON.stringify(bonds, null, 2) + '\n');
    console.log(`  Wrote ${outPath} (${bonds.length} non-single bonds)`);
}

/**
 * Parse the _chem_comp_bond loop from a CCD CIF file.
 * Returns only bonds that are non-single or aromatic (the ones perception must assign).
 */
function parseCcdBonds(cifText, compId) {
    const lines = cifText.split('\n');

    // Find the _chem_comp_bond loop and collect column header indices
    let inBondLoop = false;
    const colIndex = {};
    let colCount = 0;
    const dataLines = [];

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();

        if (line === 'loop_') {
            // Peek ahead: is this the _chem_comp_bond loop?
            let j = i + 1;
            while (j < lines.length && lines[j].trim().startsWith('_chem_comp_bond.')) j++;
            if (j > i + 1) {
                inBondLoop = true;
                colCount = 0;
                for (let k = i + 1; k < j; k++) {
                    const col = lines[k].trim().replace('_chem_comp_bond.', '');
                    colIndex[col] = colCount++;
                }
                i = j - 1; // skip headers, outer loop will increment
            } else {
                inBondLoop = false;
            }
            continue;
        }

        if (inBondLoop) {
            if (line.startsWith('_') || line === '#' || line === '') {
                if (line.startsWith('_') || line === '#') inBondLoop = false;
                continue;
            }
            dataLines.push(line);
        }
    }

    const idxA = colIndex['atom_id_1'];
    const idxB = colIndex['atom_id_2'];
    const idxOrder = colIndex['value_order'];
    const idxAro = colIndex['pdbx_aromatic_flag'];

    if (idxA === undefined || idxB === undefined || idxOrder === undefined || idxAro === undefined) {
        throw new Error(`Missing expected _chem_comp_bond columns in CIF for ${compId}`);
    }

    const result = [];
    for (const line of dataLines) {
        const fields = line.split(/\s+/);
        const a = fields[idxA];
        const b = fields[idxB];
        const order = fields[idxOrder];
        const aromatic = fields[idxAro];

        if (aromatic === 'Y') {
            result.push({ a, b, aromatic: true });
        } else if (order === 'DOUB') {
            result.push({ a, b, order: 2 });
        } else if (order === 'TRIP') {
            result.push({ a, b, order: 3 });
        }
        // SING non-aromatic: skip
    }
    return result;
}
