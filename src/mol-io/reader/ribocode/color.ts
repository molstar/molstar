/**
 * Copyright (c) 2025-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */

export interface SimpleData {
    pdb_chain: string;
    color: string;
}

// This function reads a JSON file and returns a promise that resolves to a Map object
// The file is expected to be a JSON file with the following format:
// {pdb_chain: color}
// Where:
// pdb_chain is the name of the PDB chain (e.g. 4ug0_LY)
// color is a hex color code (e.g. #FF0000)
export async function readJSONFile(file: File): Promise<SimpleData[]> {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            try {
                const text = event.target?.result as string;
                const data = JSON.parse(text);
                const transformedData = data.map((item: { pdb_chain: string; color: string }) => ({
                    pdb_chain: item.pdb_chain,
                    color: item.color
                }));
                resolve(transformedData);
            } catch (error) {
                reject(error);
            }
        };
        reader.onerror = (error) => reject(error);
        reader.readAsText(file);
    });
}

export interface Data {
    pdb_name: string;
    RP_name: string;
    color: string;
    pdb_id: string;
    pdb_chain: string;
}

// This function reads a file and returns a promise that resolves to an array of Data objects
// The file is expected to be a text file with the following format:
// pdb_name RP_name class color
// Where:
// pdb_name is the name of the PDB chain (e.g. 4ug0_LY)
// RP_name is the name of the ribosomal protein (e.g. RPL26)
// class is an integer representing the color class (e.g. 1) - this is ignored for now
// color is a hex color code (e.g. #FF0000)
export async function readFile(file: File): Promise<Data[]> {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            const text = event.target?.result as string;
            const lines = text.split('\n');
            const data: Data[] = [];
            if (lines.length > 0) {
                const header = lines[0];
                console.log('Filename:', file.name);
                console.log('Header:', header);
            }
            for (let i = 1; i < lines.length; i++) {
                const line = lines[i];
                if (line.trim() === '') continue;
                const [pdb_name, RP_name, , colorStr] = line.split(' ');
                const [pdb_id, pdb_chain] = pdb_name.split('_');
                data.push({
                    pdb_name,
                    RP_name,
                    //class: parseInt(classStr, 10),
                    color: colorStr, // Can later be parsed using mol-util.color.color.Color.fromHexStyle(colorStr),
                    pdb_id,
                    pdb_chain
                });
            }
            resolve(data);
        };
        reader.onerror = (error) => reject(error);
        reader.readAsText(file);
    });
}