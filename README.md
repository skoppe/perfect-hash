# Perfect-hash

<img src="https://github.com/skoppe/perfect-hash/workflows/build/badge.svg"/>

Generates perfect hash function for a set of keys.

Can be used to quickly determine if a key belongs in a set.

# Examples

Either supply a list of strings:

```d
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
```d
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

```d
import std.random;
import perfecthash.generator;
auto keys = ["asdf", "fdsa", "car", "garage", "elephant", "kangeroo"];
auto rnd = Random(9123);
auto ph = createPerfectHashFunction(keys, rnd);
assert(keys[ph("garage")] == "garage");
assert(keys[ph("notfound")] != "notfound");
```

Call `toDModule("functionName", "my.module")` on a perfect hash function to get a string which contains a static version of the hash function. You can save this to a D file and include in your project.

# Our sponsors

[<img src="https://raw.githubusercontent.com/libmir/mir-algorithm/master/images/symmetry.png" height="80" />](http://symmetryinvestments.com/)
