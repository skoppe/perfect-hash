# Perfect-hash

Generates perfect hash function for a set of keys.

Can be used to quickly determine if a key belongs in a set.

# Examples

Either supply a list of strings:

```Dlang
import std.random;
import perfecthash.generator;
auto keys = ["asdf", "fdsa", "car", "garage", "elephant", "kangeroo"];
auto rnd = Random(9123);
auto ph = createPerfectHashFunction(keys, rnd);
foreach(idx, key; keys)
  assert(ph(key) == idx);
```

Or a list of KeyValue where you yourself determine to which value each string gets hashed too. Usefull if a set of strings consists of subsets, and you quickly want to determine which subset it belongs to.

Here is a set of strings which contain either animals '0' or things '1'
```Dlang
import std.random;
import perfecthash.generator;
enum Class {
  Animal = 0,
  Thing = 1
}
auto keys = [KeyValue("gorilla",Class.Animal), KeyValue("car",Class.Thing), KeyValue("house",Class.Thing), KeyValue("kangeroo",Class.Animal), KeyValue("elephant",Class.Animal), KeyValue("door",Class.Thing)];
auto rnd = Random(9812);
auto ph = createPerfectHashFunction(keys, rnd);
assert(ph("car") == Class.Thing);
assert(ph("elephant") == Class.Animal);
```

To quickly determine whether a string belongs in a set:

```Dlang
import std.random;
import perfecthash.generator;
auto keys = ["asdf", "fdsa", "car", "garage", "elephant", "kangeroo"];
auto rnd = Random(9123);
auto ph = createPerfectHashFunction(keys, rnd);
assert(keys[ph("garage")] == "garage");
assert(keys[ph("notfound")] != "notfound");
```
