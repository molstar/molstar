/**
 * Copyright (c) 2025-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */
import { ChainIdColorThemeParams } from '../../mol-theme/color/chain-id';

export type Data = Record<string, string>;

/**
 * Read a JSON file and return a promise that resolves to a Map object.
 * @param file A JSON file with the following format:
 * {pdb_chain: color} where:
 * pdb_chain is the name of the PDB chain (e.g. 4ug0_LY)
 * color is a hex color code (e.g. #FF0000)
 */
export async function readJSONFile(file: File): Promise<Data[]> {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            try {
                const text = event.target?.result as string;
                const data = JSON.parse(text);
                //let transformedData: SimpleData[];
                let transformedData: Data[];
                if (Array.isArray(data)) {
                    transformedData = data.map((item: any) => ({
                        pdb_chain: item.pdb_chain,
                        color: item.color
                    }));
                } else {
                    // Handle object format: { "4ug0_LY": "#FF0000", ... }
                    transformedData = Object.entries(data).map(
                        ([pdb_chain, color]) => ({
                        pdb_chain,
                        color: color as string
                    }));
                }
                resolve(transformedData);
            } catch (error) {
                reject(error);
            }
        };
        reader.onerror = (error) => reject(error);
        reader.readAsText(file);
    });
}

/**
 * Reads a text file and return a promise that resolves to an array of Data objects.
 * @param file A space or tab separated text file with a header line and lines of 
 * data in the following format:
 * pdb_name RP_name class color
 * pdb_name contains the label of a molecule and identifier for a PDB chain separated
 * by "_" (e.g. 4ug0_LY)
 * RP_name is the name of the ribosomal protein (e.g. RPL26)
 * class is an integer representing the color class (e.g. 1) - this is ignored for now
 * color is a hex color code (e.g. #FF0000)
 */
export async function readFile(file: File): Promise<Data[]> {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            console.log('Started parsing file:', file.name);
            const text = event.target?.result as string;
            const lines = text.split('\n');
            const data: Data[] = [];
            if (lines.length > 0) {
                const header = lines[0];
                console.log('Filename:', file.name);
                console.log('Header:', header);
            }
            const numLines = lines.length;
            console.log('Total lines in file:', numLines);
            const nprint = Math.max(1, Math.floor(numLines / 4));
            for (let i = 1; i < numLines; i++) {
                const line = lines[i].trim();
                // if (i % nprint === 0) {
                //     console.log(`Parsing line ${i}: ${line}`);
                // }
                if (!line) continue;
                const parts = line.split(/\s+/);
                if (parts.length < 4) continue; // skip malformed lines
                const [pdb_name, RP_name, , colorStr] = parts;
                const [pdb_id, pdb_chain] = pdb_name.split('_');
                if (!pdb_id || !pdb_chain) continue; // skip malformed pdb_name
                data.push({
                    pdb_name,
                    RP_name,
                    color: colorStr,
                    pdb_id,
                    pdb_chain
                });
                if (i % nprint === 0) {
                    console.log(`Parsed line ${i}: pdb_name=${pdb_name}, RP_name=${RP_name}, color=${colorStr}, pdb_id=${pdb_id}, pdb_chain=${pdb_chain}`);
                }
            }
            console.log('Finished parsing file. Total entries:', data.length);
            resolve(data);
        };
        reader.onerror = (error) => reject(error);
        reader.readAsText(file);
    });
}

/**
 * Save the given ChainIdColorThemeParams as a JSON file.
 * @param themeParams The ChainIdColorThemeParams to save.
 */
export function saveColorTheme(themeParams: ChainIdColorThemeParams) {
    const output = {
        params: themeParams
    };
    const jsonString = JSON.stringify(output, null, 2);
    const blob = new Blob([jsonString], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'color_theme.json';
    a.click();
    URL.revokeObjectURL(url);
}

/**
 * Load ChainIdColorThemeParams from a JSON file.
 * @param file The JSON file to load.
 * @return A promise that resolves to the ChainIdColorThemeParams.
 */
export function loadColorTheme(file: File): Promise<ChainIdColorThemeParams> {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            const text = event.target?.result as string;
            const json = JSON.parse(text);
            resolve(json.params);
        };
        reader.onerror = (error) => reject(error);
        reader.readAsText(file);
    });
} 