How CIF schemas work
========

CIF representation (simplified):

```ts
type Block = (name: string) => Category | undefined
type Category = (name: string) => CIFField | undefined
type CIFField = { getNumber: (row) => number, getString: (row) => string }
```

This is obviously not strongly typed + the "fields" dont know what type they are. To solve this, we create a type to describe what a field contains and how to map it to a column (which is "typed"):

```ts
type FieldSchema<T> = { T: T /* remember the type */, createColumn: CIFField => Column<T> }
```

Category schema is just an object whose properties are all instances of "field schemas", its "shape" has the type:

```ts
type CategorySchema = { [fieldName: string]: FieldSchema<any> }
```

We can declare our first category "schema":

```ts
const my_category = {
  num_field: { T: 0 as number, createColumn: f => Column(f => f.getNumber) }
  str_field: { T: '' as string, createColumn: f => Column(f => f.getString) }
}
```

Notice that the type of ``my_category`` is not specified. Assigning it explictly would hide the actual property names which we do not want. Moreover, the names of the properties must match the names of the fields in the actual category (optionally, a field ``alias`` can be added to the field schema).

Given a category schema, we need to construct a type that defines the typed category itself:

```ts
type TypedCategory<Schema extends CategorySchema> = { [F in keyof Schema]: Column<Schema[F]['T']> }
```

In other words, the type ``TypedCategory`` has a property of type ``Column<_>`` for each property of the schema. ``Schema[F]['T']`` just says: extract type of property called ``T`` from property ``F`` in ``Schema``. ``Schema extends CategorySchema`` says that all properties of ``Schema`` must be of type ``FieldSchema<any>``.

Finally, we just define a mapping, ``toTypedCategory``:

```ts
function toTypedCategory<Schema extends CategorySchema>(schema: Schema, category: Category): TypedCategory<Shape> {
    const typedCategory: any = {};
    for (const key in Object.keys(schema)) {
        // remember a category is just a function that assigns a Field to a name
        const field = category(key);
        typedCategory[key] = field 
            ? schema[key].createFolumn(field)
            : UndefinedColumn; // a column that always returns 0 or empty string depending on type
    }
    return typedCategory;
}
```

And that's all there is to it. Extending the types to the "block" level is left as an exercise to the reader.

----------------


**Note:** To create a type alias for a category defined this way we can do:

```ts
type MyCategory = TypedCategory<typeof my_category>
function makeMyTypedCategory(c: Category): MyCategory { return toTypedCategory(my_category, c); }
```
