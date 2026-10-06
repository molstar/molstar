/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alex Chan <smalldirkalex@gmail.com>
 *
 * Thanks to @author Adam Midlik <midlik@gmail.com> for the example code ../image-renderer and https://github.com/midlik/surface-calculator i can make reference to,
 *
 * Example command-line application generating and exporting PubChem SDF structures
 * Build: bun add -d webgpu @napi-rs/canvas jpeg-js pngjs  // these packages are not listed in dependencies for performance reasons
 *        bun run build
 * Run:   node lib/commonjs/examples/glb-export 2519 ../outputs_2519/
 */

import { ArgumentParser } from 'argparse';
import fs from 'fs';
import path from 'path';

import { Task } from '../../mol-task';
import { Download } from '../../mol-plugin-state/transforms/data';
import { GlbExporter } from '../../extensions/geo-export/glb-exporter';
import { Box3D } from '../../mol-math/geometry';
import { ModelFromTrajectory, StructureFromModel, TrajectoryFromSDF } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { HeadlessPluginContext } from '../../mol-plugin/headless-plugin-context';
import { DefaultPluginSpec } from '../../mol-plugin/spec';
import { loadHeadlessModules } from '../../mol-plugin/util/headless-modules';
import { setFSModule } from '../../mol-util/data-source';

setFSModule(fs);

// cid `2519` for Caffeine
interface Args {
    cid: string,
    outDirectory: string
}

function parseArguments(): Args {
    const parser = new ArgumentParser({ description: 'Example command-line application exporting .glb file of SDF structures from PubChem' });
    parser.add_argument('cid', { help: 'PubChem identifier' });
    parser.add_argument('outDirectory', { help: 'Directory for outputs' });
    const args = parser.parse_args();
    return { ...args };
}

async function main() {
    const args = parseArguments();
    const root = 'https://pubchem.ncbi.nlm.nih.gov/rest';
    const url = `${root}/pug/compound/cid/${args.cid}/sdf?record_type=3d`;

    console.log('PubChem CID:', args.cid);
    console.log('Source URL:', url);
    console.log('Outputs:', args.outDirectory);

    // Create a headless plugin
    const externalModules = await loadHeadlessModules();
    const plugin = await HeadlessPluginContext.create(externalModules, DefaultPluginSpec());
    await plugin.init();

    // Download and visualize data in the plugin
    const update = plugin.build();
    const structure = await update.toRoot()
        .apply(Download, { url, isBinary: false })
        .apply(TrajectoryFromSDF)
        .apply(ModelFromTrajectory)
        .apply(StructureFromModel)
        .apply(StructureRepresentation3D, {
            type: { name: 'ball-and-stick', params: { size: 'physical' } },
            colorTheme: { name: 'element-symbol', params: { carbonColor: { name: 'element-symbol', params: {} } } },
            sizeTheme: { name: 'physical', params: {} },
        })
        .commit();

    plugin.canvas3d!.commit(true);
    const renderObjects = structure.data!.repr.renderObjects;

    const boundingSphere = plugin.canvas3d?.boundingSphereVisible!;
    const boundingBox = Box3D.fromSphere3D(Box3D(), boundingSphere);

    const renderObjectExporter = new GlbExporter(boundingBox);

    // Match the image-renderer example: callers commonly pass a new output
    // directory, so create it before the asynchronous export task writes the
    // GLB file.
    await fs.promises.mkdir(args.outDirectory, { recursive: true });
    await plugin.runTask(Task.create('Export Geometry', async ctx => {
        for (let i = 0, il = renderObjects.length; i < il; ++i) {
            await renderObjectExporter.add(renderObjects[i], plugin.canvas3d?.webgl, ctx);
        }

        const blob = await renderObjectExporter.getBlob(ctx);
        const buffer = await blob.arrayBuffer();
        await fs.promises.writeFile(path.join(args.outDirectory, `${args.cid}.glb`), Buffer.from(buffer));
    }));

    // Cleanup
    await plugin.clear();
    plugin.dispose();
}

main();
