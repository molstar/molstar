/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginBehavior } from '../behavior';
import { Binding } from '../../../mol-util/binding';
import { ModifiersKeys } from '../../../mol-util/input/input-observer';

const M = ModifiersKeys;
const Key = Binding.TriggerKey;

const DefaultSnapshotControlsBindings = {
    next: Binding([
        Key('ArrowRight', M.create({ control: true })),
        Key('GamepadY'),
    ]),
    previous: Binding([
        Key('ArrowLeft', M.create({ control: true })),
        Key('GamepadX'),
    ]),
    first: Binding([
        Key('ArrowUp', M.create({ control: true })),
    ]),
    last: Binding([
        Key('ArrowDown', M.create({ control: true })),
    ]),
};
const SnapshotControlsParams = {
    bindings: PD.Value(DefaultSnapshotControlsBindings, { isHidden: true }),
};
type SnapshotControlsProps = PD.Values<typeof SnapshotControlsParams>

export const SnapshotControls = PluginBehavior.create<SnapshotControlsProps>({
    name: 'snapshot-controls',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<SnapshotControlsProps> {
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.interaction.keyReleased, ({ code, modifiers, key }) => {
                if (!this.ctx.canvas3d || this.ctx.isBusy) return;

                const b = this.params.bindings;
                const { snapshot } = this.ctx.managers;

                if (Binding.matchKey(b.next, code, modifiers, key)) {
                    snapshot.applyNext(1);
                }

                if (Binding.matchKey(b.previous, code, modifiers, key)) {
                    snapshot.applyNext(-1);
                }

                if (Binding.matchKey(b.first, code, modifiers, key)) {
                    const e = snapshot.state.entries.get(0)!;
                    const s = snapshot.setCurrent(e.snapshot.id);
                    if (s) return this.ctx.state.setSnapshot(s);
                }

                if (Binding.matchKey(b.last, code, modifiers, key)) {
                    const e = snapshot.state.entries.get(snapshot.state.entries.size - 1)!;
                    const s = snapshot.setCurrent(e.snapshot.id);
                    if (s) return this.ctx.state.setSnapshot(s);
                }
            });
        }
    },
    params: () => SnapshotControlsParams,
    display: { name: 'Snapshot Controls' }
});
