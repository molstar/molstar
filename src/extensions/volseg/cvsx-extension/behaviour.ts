import { DataFormatProvider } from "../../../mol-plugin-state/formats/provider";
import { PluginBehavior } from "../../../mol-plugin/behavior";
import { ParamDefinition } from "../../../mol-util/param-definition";
import { CVSXFormatProvider } from "./formats";

/** Collection of things that can be register/unregistered in a plugin */
interface Registrables {
    // customModelProperties?: CustomModelProperty.Provider<any, any>[],
    // customStructureProperties?: CustomStructureProperty.Provider<any, any>[],
    // representations?: StructureRepresentationProvider<any>[],
    // colorThemes?: ColorTheme.Provider[],
    // lociLabels?: LociLabelProvider[],
    // dragAndDropHandlers?: DragAndDropHandler[],
    dataFormats?: { name: string, provider: DataFormatProvider }[],
    // actions?: StateAction[],
}

export const CVSXSpec = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'CVSXSpec',
    category: 'misc',
    display: {
        name: 'CVSXSpec',
        description: 'CVSXSpec extension',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private readonly registrables: Registrables = {
            dataFormats: [
                // { name: 'MVSJ', provider: MVSJFormatProvider },
                // { name: 'MVSX', provider: MVSXFormatProvider },
                { name: 'CVSX', provider: CVSXFormatProvider}
            ],
            // actions: [
            //     LoadMvsData,
            // ]
        };

        register(): void {
            for (const format of this.registrables.dataFormats ?? []) {
                this.ctx.dataFormats.add(format.name, format.provider);
            }
            // for (const action of this.registrables.actions ?? []) {
            //     this.ctx.state.data.actions.add(action);
            // }
        }
        update(p: { autoAttach: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            // for (const prop of this.registrables.customModelProperties ?? []) {
            //     this.ctx.customModelProperties.setDefaultAutoAttach(prop.descriptor.name, this.params.autoAttach);
            // }
            // for (const prop of this.registrables.customStructureProperties ?? []) {
            //     this.ctx.customStructureProperties.setDefaultAutoAttach(prop.descriptor.name, this.params.autoAttach);
            // }
            return updated;
        }
        unregister() {
            for (const format of this.registrables.dataFormats ?? []) {
                this.ctx.dataFormats.remove(format.name);
            }
            // for (const action of this.registrables.actions ?? []) {
            //     this.ctx.state.data.actions.remove(action);
            // }
        }
    },
    params: () => ({
        autoAttach: ParamDefinition.Boolean(false),
    })
});