/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sukolsak Sakshuwong <sukolsak@stanford.edu>
 */

import { PluginComponent } from '../../mol-plugin-state/component';
import { PluginContext } from '../../mol-plugin/context';
import { Task } from '../../mol-task';
import { ObjExporter } from './export';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateSelection } from '../../mol-state';
import { SetUtils } from '../../mol-util/set';
import { zip } from '../../mol-util/zip/zip';

export class GeometryControls extends PluginComponent {
    getFilename() {
        const models = this.plugin.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model)).map(s => s.obj!.data);
        const uniqueIds = new Set<string>();
        models.forEach(m => uniqueIds.add(m.entryId.toUpperCase()));
        const idString = SetUtils.toArray(uniqueIds).join('-');
        return `${idString || 'molstar-model'}`;
    }

    exportObj() {
        const task = Task.create('Export OBJ', async ctx => {
            try {
                const renderObjects = this.plugin.canvas3d?.getRenderObjects()!;

                await ctx.update({ message: 'Rendering...', isIndeterminate: false, current: 0, max: renderObjects.length });

                const filename = this.getFilename();
                const objExporter = new ObjExporter(filename);
                for (let i = 0, il = renderObjects.length; i < il; ++i) {
                    objExporter.add(renderObjects[i]);
                    await ctx.update({ current: i });
                }
                const { obj, mtl } = objExporter.getData();

                const asciiWrite = (data: Uint8Array,  str: string) => {
                    for (let i = 0, il = str.length; i < il; ++i) {
                        data[i] = str.charCodeAt(i);
                    }
                };
                const objData = new Uint8Array(obj.length);
                asciiWrite(objData, obj);
                const mtlData = new Uint8Array(mtl.length);
                asciiWrite(mtlData, mtl);

                const zipDataObj = {
                    [filename + '.obj']: objData,
                    [filename + '.mtl']: mtlData
                };
                const zipData = await zip(ctx, zipDataObj);
                return {
                    zipData,
                    filename: filename + '.zip'
                };
            } catch (e) {
                this.plugin.log.error('' + e);
                throw e;
            }
        });

        return this.plugin.runTask(task, { useOverlay: true });
    }

    constructor(private plugin: PluginContext) {
        super();
    }
}