import { Transformer } from 'mol-state/transformer';
import { StateObject } from 'mol-state/object';
import { Task } from 'mol-task';
import { TransformTree } from 'mol-state/tree/tree';
import { StateTreeBuilder } from 'mol-state/tree/builder';
import { State } from 'mol-state/state';
import * as util from 'util'

export type TypeClass = 'root' | 'shape' | 'prop'
export interface ObjProps { label: string }
export interface TypeInfo { name: string, class: TypeClass }

const _obj = StateObject.factory<TypeInfo, ObjProps>()
const _transform = Transformer.factory('test');

export class Root extends _obj('root', { name: 'Root', class: 'root' }) { }
export class Square extends _obj<{ a: number }>('square', { name: 'Square', class: 'shape' }) { }
export class Circle extends _obj<{ r: number }>('circle', { name: 'Circle', class: 'shape' }) { }
export class Area extends _obj<{ volume: number }>('volume', { name: 'Volume', class: 'prop' }) { }

export const CreateSquare = _transform<Root, Square, { a: number }>({
    name: 'create-square',
    from: [Root],
    to: [Square],
    apply({ params: p }) {
        return new Square({ label: `Square a=${p.a}` }, p);
    },
    update({ b, newParams: p }) {
        b.props.label = `Square a=${p.a}`
        b.data.a = p.a;
        return Transformer.UpdateResult.Updated;
    }
});

export const CreateCircle = _transform<Root, Circle, { r: number }>({
    name: 'create-circle',
    from: [Root],
    to: [Square],
    apply({ params: p }) {
        return new Circle({ label: `Circle r=${p.r}` }, p);
    },
    update({ b, newParams: p }) {
        b.props.label = `Circle r=${p.r}`
        b.data.r = p.r;
        return Transformer.UpdateResult.Updated;
    }
});

export const CaclArea = _transform<Square | Circle, Area, {}>({
    name: 'calc-area',
    from: [Square, Circle],
    to: [Area],
    apply({ a }) {
        if (a instanceof Square) return new Area({ label: 'Area' }, { volume: a.data.a * a.data.a });
        else if (a instanceof Circle) return new Area({ label: 'Area' }, { volume: a.data.r * a.data.r * Math.PI });
        throw new Error('Unknown object type.');
    },
    update({ a, b }) {
        if (a instanceof Square) b.data.volume = a.data.a * a.data.a;
        else if (a instanceof Circle) b.data.volume = a.data.r * a.data.r * Math.PI;
        else throw new Error('Unknown object type.');
        return Transformer.UpdateResult.Updated;
    }
});

export async function runTask<A>(t: A | Task<A>): Promise<A> {
    if ((t as any).run) return await (t as Task<A>).run();
    return t as A;
}

export async function testState() {
    const state = State.create();

    const tree = state.definition.tree;
    const builder = StateTreeBuilder.create(tree)
    builder.toRoot<Root>()
        .apply(CreateSquare, { a: 10 }, { ref: 'square' })
        .apply(CaclArea);
    const tree1 = builder.getTree();

    printTTree(tree1);

    const tree2 = TransformTree.updateParams<typeof CreateSquare>(tree1, 'square', { a: 15 });
    printTTree(tree1);
    printTTree(tree2);

    const state1 = await State.update(state, tree1);
    console.log('----------------');
    console.log(util.inspect(state1.objects, true, 3, true));

    console.log('----------------');
    const jsonString = JSON.stringify(TransformTree.toJSON(tree2), null, 2);
    const jsonData = JSON.parse(jsonString);
    printTTree(tree2);
    console.log(jsonString);
    const treeFromJson = TransformTree.fromJSON(jsonData);
    printTTree(treeFromJson);

    console.log('----------------');
    const state2 = await State.update(state1, treeFromJson);
    console.log(util.inspect(state2.objects, true, 3, true));
}

testState();


//test();

export function printTTree(tree: TransformTree) {
    let lines: string[] = [];
    function print(offset: string, ref: any) {
        const t = tree.nodes.get(ref)!;
        const tr = t.value;

        const name = tr.transformer.id;
        lines.push(`${offset}|_ (${ref}) ${name} ${tr.params ? JSON.stringify(tr.params) : ''}, v${t.value.version}`);
        offset += '   ';

        t.children.forEach(c => print(offset, c!));
    }
    print('', tree.rootRef);
    console.log(lines.join('\n'));
}