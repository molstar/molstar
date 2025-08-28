import { OrderedMap } from 'immutable';
import { Observable, Subject, Unsubscribable } from 'rxjs';
import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginContext } from '../context';
import { Color } from '../../mol-util/color';


export namespace SequenceColoring {
    export interface Provider {
        readonly name: string,
        color: (loci: StructureElement.Loci) => Color | undefined, // TODO change to Location, to align with color themes?, or have Theme directly here?
        subcribeForUpdates?: (plugin: PluginContext, requestUpdate: (structure: Structure) => void) => Unsubscribable,
    }

    export class Registry {
        private readonly _providers = OrderedMap<string, Provider>().asMutable();
        private readonly _updates = new Subject<Structure>();
        private readonly _subscriptions: { [providerName: string]: Unsubscribable | undefined } = {};

        constructor(public readonly plugin: PluginContext) { }

        register(provider: Provider) {
            if (this._providers.has(provider.name)) {
                this.unregister(provider.name);
            }
            this._providers.set(provider.name, provider);
            this._subscriptions[provider.name] = provider.subcribeForUpdates?.(this.plugin, structure => this._updates.next(structure));
        }
        unregister(providerName: string) {
            if (!this._providers.has(providerName)) return;
            this._subscriptions[providerName]?.unsubscribe();
            delete this._subscriptions[providerName];
            this._providers.delete(providerName);
        }

        providers(): IterableIterator<Provider> {
            return this._providers.values();
        }

        get updates(): Observable<Structure> { return this._updates; }
    }
}
