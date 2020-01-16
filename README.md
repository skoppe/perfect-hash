# Perfect-hash

<img src="https://github.com/skoppe/perfect-hash/workflows/build/badge.svg"/>

Generates a perfect hash function for a set of strings.

Can be used to quickly determine if a string belongs in a set or a subset. Instead of comparing strings (slow) you compare its hash (fast). The set of string must be known upfront.

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

Or a list of KeyValue where you yourself determine to which value each string gets hashed to. Useful if a set of strings consists of subsets and you quickly want to determine which subset it belongs to.

Here is a set of strings which contains either animals or things.
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

To quickly determine whether a string is in a set:

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
