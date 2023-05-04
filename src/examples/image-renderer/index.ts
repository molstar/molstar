/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Example command-line application generating images of PDB structures
 * Build: npm install --no-save gl jpeg-js pngjs  // these packages are not listed in dependencies for performance reasons
 *        npm run build
 * Run:   node lib/commonjs/examples/image-renderer 1cbs ../outputs_1cbs/
 */

import { ArgumentParser } from 'argparse';
import fs from 'fs';
import path from 'path';
import gl from 'gl';
import pngjs from 'pngjs';
import jpegjs from 'jpeg-js';

import { Download, ParseCif } from '../../mol-plugin-state/transforms/data';
import { ModelFromTrajectory, StructureComponent, StructureFromModel, TrajectoryFromMmCif } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { HeadlessPluginContext } from '../../mol-plugin/headless-plugin-context';
import { DefaultPluginSpec } from '../../mol-plugin/spec';
import { ExternalModules, STYLIZED_POSTPROCESSING } from '../../mol-plugin/util/headless-screenshot';
import { setFSModule } from '../../mol-util/data-source';


setFSModule(fs);

interface Args {
    pdbId: string,
    outDirectory: string
}

function parseArguments(): Args {
    const parser = new ArgumentParser({ description: 'Example command-line application generating images of PDB structures' });
    parser.add_argument('pdbId', { help: 'PDB identifier' });
    parser.add_argument('outDirectory', { help: 'Directory for outputs' });
    const args = parser.parse_args();
    return { ...args };
}

async function main() {
    const args = parseArguments();
    const url = `https://www.ebi.ac.uk/pdbe/entry-files/download/${args.pdbId}.bcif`;
    console.log('PDB ID:', args.pdbId);
    console.log('Source URL:', url);
    console.log('Outputs:', args.outDirectory);

    // Create a headless plugin
    const externalModules: ExternalModules = { gl, pngjs, 'jpeg-js': jpegjs };
    const plugin = new HeadlessPluginContext(externalModules, DefaultPluginSpec(), { width: 800, height: 800 });
    await plugin.init();

    // Download and visualize data in the plugin
    const update = plugin.build();
    const structure = update.toRoot()
        .apply(Download, { url, isBinary: true })
        .apply(ParseCif)
        .apply(TrajectoryFromMmCif)
        .apply(ModelFromTrajectory)
        .apply(StructureFromModel);
    const polymer = structure.apply(StructureComponent, { type: { name: 'static', params: 'polymer' } });
    const ligand = structure.apply(StructureComponent, { type: { name: 'static', params: 'ligand' } });
    polymer.apply(StructureRepresentation3D, {
        type: { name: 'cartoon', params: { alpha: 1 } },
        colorTheme: { name: 'sequence-id', params: {} },
    });
    ligand.apply(StructureRepresentation3D, {
        type: { name: 'ball-and-stick', params: { sizeFactor: 1 } },
        colorTheme: { name: 'element-symbol', params: { carbonColor: { name: 'element-symbol', params: {} } } },
        sizeTheme: { name: 'physical', params: {} },
    });
    await update.commit();

    // Export images
    fs.mkdirSync(args.outDirectory, { recursive: true });
    await plugin.saveImage(path.join(args.outDirectory, 'basic.png'));
    await plugin.saveImage(path.join(args.outDirectory, 'basic.jpg'));
    await plugin.saveImage(path.join(args.outDirectory, 'large.png'), { width: 1600, height: 1200 });
    await plugin.saveImage(path.join(args.outDirectory, 'large.jpg'), { width: 1600, height: 1200 });
    await plugin.saveImage(path.join(args.outDirectory, 'stylized.png'), undefined, STYLIZED_POSTPROCESSING);
    await plugin.saveImage(path.join(args.outDirectory, 'stylized.jpg'), undefined, STYLIZED_POSTPROCESSING);
    await plugin.saveImage(path.join(args.outDirectory, 'stylized-compressed-jpg.jpg'), undefined, STYLIZED_POSTPROCESSING, undefined, 10);

    // Export state loadable in Mol* Viewer
    await plugin.saveStateSnapshot(path.join(args.outDirectory, 'molstar-state.molj'));

    // Cleanup
    await plugin.clear();
    plugin.dispose();
}

main();
