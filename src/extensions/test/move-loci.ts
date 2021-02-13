/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { Loci, EmptyLoci } from '../../mol-model/loci';
import { ButtonsType, ModifiersKeys } from '../../mol-util/input/input-observer';
import { Binding } from '../../mol-util/binding';
import { PluginStateTransform, PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { Vec3, Mat3 } from '../../mol-math/linear-algebra';
import { deepEqual, UUID } from '../../mol-util';
import { StateTransformer, StateSelection } from '../../mol-state';
import { StructureElement, Model, ElementIndex, Unit } from '../../mol-model/structure';
import { OrderedSet, SortedArray } from '../../mol-data/int';
import { AssignableArrayLike, TypedArray } from '../../mol-util/type-helpers';
import { ObjectControls } from '../../mol-canvas3d/controls/object';

const B = ButtonsType;
const M = ModifiersKeys;
const Trigger = Binding.Trigger;

const DefaultMoveLociBindings = {
    dragMove: Binding([
        Trigger(B.Flag.Primary, M.create({ alt: true })),
    ], 'Move Loci', 'Drag using ${triggers}'),
};
const MoveLociParams = {
    bindings: PD.Value(DefaultMoveLociBindings, { isHidden: true }),
};
type MoveLociProps = PD.Values<typeof MoveLociParams>

export const MoveLoci = PluginBehavior.create<MoveLociProps>({
    name: 'move-loci',
    category: 'interaction',
    display: { name: 'Move Loci on Canvas' },
    ctor: class extends PluginBehavior.Handler<MoveLociProps> {
        register(): void {
            this.ctx.state.data.actions.add(PatchedModel);

            let busy = false;
            let pendingDir: Vec3 | undefined;
            let currentLoci: Loci = EmptyLoci;
            const dir = Vec3();
            const m3 = Mat3();
            const ref = Vec3();

            this.subscribeObservable(this.ctx.behaviors.interaction.drag, ({ current, buttons, modifiers, pageStart, pageEnd }) => {
                if (!this.ctx.canvas3d) return;

                const binding = this.params.bindings.dragMove;
                if (!Binding.match(binding, buttons, modifiers)) return;

                if (!Loci.areEqual(currentLoci, current.loci)) {
                    currentLoci = current.loci;
                    pendingDir = undefined;
                }

                const loci = Loci.normalize(current.loci, this.ctx.managers.interactivity.props.granularity);
                if (Loci.isEmpty(loci) || !StructureElement.Loci.is(loci)) return;

                Loci.getCenter(loci, ref);
                const direction = ObjectControls.panDirection(Vec3(), pageStart, pageEnd, this.ctx.canvas3d.camera, ref);

                if (!busy) {
                    const { unit, indices } = loci.elements[0];

                    const patchedModel = this.ctx.state.data.select(StateSelection.Generators.ofTransformer(PatchedModel)).filter(r => r.obj && Model.getRoot(r.obj.data) === Model.getRoot(unit.model))[0];
                    if (!patchedModel?.obj) return;

                    const u = unit.remapModel(patchedModel.obj.data);

                    if (pendingDir) {
                        Vec3.add(dir, pendingDir, direction);
                    } else {
                        Vec3.copy(dir, direction);
                    }
                    // console.log('move', Vec3.magnitude(dir));
                    Vec3.transformMat3(dir, dir, Mat3.directionTransform(m3, u.conformation.operator.inverse));

                    const b = this.ctx.state.data.build();
                    b.to(patchedModel.transform.ref).update(old => {
                        old.coords.length = 0;
                        OrderedSet.forEach(indices, v => {
                            const index = u.elements[v];
                            const c = u.conformation.invariantPosition(index, Vec3());
                            old.coords.push({ index, position: Vec3.add(c, c, dir) });
                        });
                    });

                    busy = true;
                    pendingDir = undefined;
                    b.commit({ canUndo: true, doNotLogTiming: true }).then(() => {
                        busy = false;
                    });
                } else if (pendingDir === undefined) {
                    // console.log('pendingPos');
                    pendingDir = direction;
                } else {
                    // console.log('wait');
                    Vec3.add(pendingDir, pendingDir, direction);
                }
            });
        }

        unregister(): void {
            this.ctx.state.data.actions.remove(PatchedModel);
        }
    },
    params: () => MoveLociParams
});

type PatchedModel = typeof PatchedModel
const PatchedModel = PluginStateTransform.BuiltIn({
    name: 'patched-model',
    display: { name: 'Patched Model', description: 'Create a patched model.' },
    isDecorator: true,
    from: SO.Molecule.Model,
    to: SO.Molecule.Model,
    params(a) {
        return {
            coords: PD.Value<CoordsPatch>([])
        };
    }
})({
    apply({ a, params, cache }) {
        if (a.data.parent) throw new Error('Only the root model can be patched.');

        Model.CoordinatesHistory.set(a.data, new CoordinatesHistory());
        const { x, y, z } = a.data.atomicConformation;
        const coords = new BufferedCoordinates(x, y, z);
        (cache as any).bufferedCoords = coords;
        if (params.coords.length === 0) return new SO.Molecule.Model(a.data);

        return new SO.Molecule.Model(getPatchedModel(a.data, coords.swap(), params.coords));
    },
    update: ({ a, b, oldParams, newParams, cache }) => {
        if (a.data !== Model.getRoot(b.data)) return StateTransformer.UpdateResult.Recreate;
        if (deepEqual(oldParams, newParams)) return StateTransformer.UpdateResult.Unchanged;

        const coords = (cache as { bufferedCoords: BufferedCoordinates }).bufferedCoords;
        b.data = getPatchedModel(b.data, coords.swap(), newParams.coords);
        return StateTransformer.UpdateResult.Updated;
    }
});

type CoordsPatch = { index: ElementIndex, position: Vec3 }[]

type Coords = {
    x: AssignableArrayLike<number>
    y: AssignableArrayLike<number>
    z: AssignableArrayLike<number>
}

const CoordsChangeProp = '__CoordsChange__';
export { CoordsChange };
type CoordsChange = SortedArray<ElementIndex>
const CoordsChange = {
    set: (model: Model, coordsChange: CoordsChange) => {
        model._staticPropertyData[CoordsChangeProp] = coordsChange;
    },
    get: (model: Model): CoordsChange | undefined => {
        return model._staticPropertyData[CoordsChangeProp];
    }
};

class CoordinatesHistory implements Model.CoordinatesHistory {
    areEqual(elements: SortedArray<ElementIndex>, kind: Unit.Kind, model: Model) {
        const coordsChange = CoordsChange.get(Model.getRoot(model));
        if (coordsChange) {
            return !SortedArray.areIntersecting(elements, coordsChange);
        }
        return false;
    }
}

function getPatchedModel(model: Model, coords: Coords, patch: CoordsPatch): Model {
    const { x, y, z } = coords;
    const root = Model.getRoot(model);
    const coordsChange: ElementIndex[] = [];

    for (const { index, position } of patch) {
        x[index] = position[0];
        y[index] = position[1];
        z[index] = position[2];
        coordsChange.push(index);
    }

    CoordsChange.set(root, SortedArray.ofUnsortedArray(coordsChange));

    const atomicConformation = {
        ...root.atomicConformation,
        id: UUID.create22(),
        x, y, z
    };

    return {
        ...root,
        id: UUID.create22(),
        parent: root,
        atomicConformation
    };
}

class BufferedCoordinates {
    private x: BufferedNumberArray
    private y: BufferedNumberArray
    private z: BufferedNumberArray

    swap() {
        return {
            x: this.x.swap(),
            y: this.y.swap(),
            z: this.z.swap()
        };
    }

    constructor(x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number>) {
        this.x = new BufferedNumberArray(x, Float32Array);
        this.y = new BufferedNumberArray(y, Float32Array);
        this.z = new BufferedNumberArray(z, Float32Array);
    }
}

class BufferedNumberArray {
    private prev: TypedArray | undefined
    private next: TypedArray | undefined

    swap() {
        if (!this.prev) this.prev = new this.ctor(this.orig);
        if (!this.next) this.next = new this.ctor(this.prev);
        const ret = this.next;
        ret.set(this.prev);
        this.next = this.prev;
        this.prev = ret;
        return ret;
    }

    constructor(private orig: ArrayLike<number>, private ctor: { new(a: ArrayLike<number>): TypedArray }) { }
}