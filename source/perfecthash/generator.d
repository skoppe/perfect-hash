module perfecthash.generator;

// Based on: http://ilan.schnell-web.net/prog/perfect-hash/algo.html
// which is an illustration of an algorithm by Z. J. Czech, G. Havas and B.S. Majewski which is described in their 1992 paper "An optimal algorithm for generating minimal perfect hash functions" which appeared in Information Processing Letters, 43(5):257-264, 1992.

import std.range : ElementType;
import std.array : Appender;
import std.typecons : Tuple;

alias KeyValue = Tuple!(string, "key", uint, "value");

struct Vertex {
  size_t value;
  size_t edgeIndex;
}

struct Edge {
  size_t[2] vertices;
  size_t hash;
}

struct PerfectHashFunction {
  const(size_t)[] G;
  HashFunction[2] hashFunction;
  this(const(size_t)[] G, HashFunction[2] hashFunction) @safe nothrow {
    this.G = G;
    this.hashFunction = [hashFunction[0].clone, hashFunction[1].clone];
  }
  size_t opCall(string key) @safe nothrow pure {
    return (G[hashFunction[0](key)] + G[hashFunction[1](key)]) % G.length;
  }
}

auto createPerfectHashFunction(Random)(string[] keys, Random rnd, ulong maxHnsecs = 1_000_000) @safe nothrow {
  import std.range : enumerate;
  import std.algorithm : map;
  import std.array : array;
  return createPerfectHashFunction(keys.enumerate.map!(k => KeyValue(k.value, cast(uint)k.index)).array(), rnd, maxHnsecs);
}

auto createPerfectHashFunction(Range, Random)(Range keys, Random rnd, ulong maxHnsecs = 1_000_000, size_t startN = size_t.max) @safe nothrow if (is(ElementType!Range == KeyValue)) {
  import std.algorithm : map, maxElement, max;
  import std.array : array;
  import std.datetime : Clock;
  auto start = Clock.currStdTime();
  auto maxKeyLength = keys.map!(k => k.key.length).maxElement();
  auto k = keys.length;
  auto n = startN == size_t.max ? k * 2 : startN;
  assert(n > k);
  Vertex[] vertices = new Vertex[n];
  Edge[] edges = new Edge[k * 2];
  HashFunction[2] hashFunctions = [HashFunction(maxKeyLength), HashFunction(maxKeyLength)];
  PerfectHashFunction bestEffort;
  size_t[] rawBits = new size_t[vertices.length / size_t.sizeof + 1];
  Appender!(ToCheck[]) toCheck;
  for(;;) {
    if (Clock.currStdTime() - start > maxHnsecs && bestEffort.G.length != 0)
      break;
    hashFunctions[0].randomize(n, 0, rnd);
    hashFunctions[1].randomize(n, 1, rnd);
    if (!calcEdges(keys, hashFunctions, edges, vertices))
      continue;
    if (assignVertexValuesOnlyWhenACyclic(edges, vertices, n, rawBits, toCheck))
      continue;
    bestEffort = PerfectHashFunction(extractG(vertices).array(), hashFunctions);
    if (n == k)
      break;
    if (Clock.currStdTime() - start > maxHnsecs)
      break;
    n = n - (max(1, (n - k) / 20));
    vertices = vertices[0 .. n];
    vertices[] = Vertex(0,0);
  }
  return bestEffort;
}

auto recreateEdgesAndVertices(Range)(PerfectHashFunction hf, Range keys) {
  import std.exception : enforce;
  import std.typecons : tuple;
  auto n = hf.G.length;
  auto k = keys.length;
  Vertex[] vertices = new Vertex[n];
  Edge[] edges = new Edge[k * 2];
  size_t[] rawBits = new size_t[vertices.length / size_t.sizeof + 1];
  Appender!(ToCheck[]) toCheck;
  enforce(calcEdges(keys, hf.hashFunction, edges, vertices), "invalid hash function");
  enforce(!assignVertexValuesOnlyWhenACyclic(edges, vertices, n, rawBits, toCheck), "invalid hash function");
  return tuple!("edges", "vertices")(edges, vertices);
}

unittest {
  import std.random;
  import std.range : enumerate;
  import std.algorithm : map, maxElement;
  import std.array : array;
  import unit_threaded;
  import std.stdio;

  auto rnd = Random(9812);
  auto keys = [
    "aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "jj", "kk", "ll",
    "mm", "nn", "oo", "pp", "qq", "rr", "ss", "tt", "uu", "vv", "ww", "xx", "yy",
    "zz"
  ].enumerate.map!(e => KeyValue(e[1], cast(uint) e[0])).array();
  foreach(_; 0 .. 100) {
    auto ph = createPerfectHashFunction(keys, rnd);
    foreach (kv; keys)
    {
      ph(kv.key).shouldEqual(kv.value);
    }
  }
}

version (unittest) {
  string randomIdentifier(Random)(ref Random rnd, size_t len = 6) {
    import std.random : uniform;
    import std.conv : text;
    import std.algorithm : map;
    import std.range : iota;

    enum chars = "abcdefghijklmnopqrstuvwxyz";
    return iota(0, len).map!(_ => chars[uniform(0, chars.length, rnd)]).text();
  }
}

@("keys.random.5")
unittest {
  import std.random;
  import std.range : enumerate, iota;
  import std.algorithm : map, maxElement;
  import std.array : array;
  import unit_threaded;
  import std.stdio;

  auto rnd = Random(9812);
  auto len = 3365;
  auto kvs = iota(0, len).map!(_ => rnd.randomIdentifier()).enumerate.map!(e => KeyValue(e[1], cast(uint)e[0])).array();

  foreach(_; 0 .. 100) {
    auto ph = kvs.createPerfectHashFunction(rnd);
    foreach (kv; kvs) {
      ph(kv.key).shouldEqual(kv.value);
    }
  }
}

string toDModule(ref PerfectHashFunction fun, string symbolName, string moduleName = "") @safe pure {
  import std.conv : to;
  import std.string : replace;
  size_t N = fun.G.length;
  size_t coeffLen = fun.hashFunction[0].coeff.length;
  string NType = N > ushort.max ? "uint" : N > ubyte.max ? "ushort" : "ubyte";
  auto func = q{auto $symbolName(string key) @safe nothrow pure {
    static $NType hash(alias coeff)(string key) {
        size_t t;
        foreach(idx, c; key)
            t += c * coeff[idx % $coeffLen];
        return t % $N;
    }
    static immutable $NType[$N] G = $G;
    static immutable $NType[$coeffLen] coeffA = $coeffA;
    static immutable $NType[$coeffLen] coeffB = $coeffB;
    return (G[hash!(coeffA)(key)] + G[hash!(coeffB)(key)]) % $N;
}}.replace("$NType", NType).replace("$moduleName", moduleName).replace("$symbolName", symbolName).replace("$N", N.to!string).replace("$coeffLen", coeffLen.to!string).replace("$G", fun.G.to!string).replace("$coeffA", fun.hashFunction[0].coeff.to!string).replace("$coeffB", fun.hashFunction[1].coeff.to!string);
  if (moduleName.length > 0)
    return "module $moduleName;\n".replace("$moduleName", moduleName) ~ func;
  return func;
}

@safe nothrow unittest {
  import std.random;
  auto keys = ["asdf", "fdsa", "car", "garage", "elephant", "kangeroo"];
  auto rnd = Random(9123);
  foreach(i; 0..10000) {
    auto ph = createPerfectHashFunction(keys, rnd, 10);
    foreach(idx, key; keys)
      assert(ph(key) == idx);
  }
}

@safe nothrow unittest {
  import std.random;
  auto keys = [KeyValue("asdf",0), KeyValue("fdsa",1), KeyValue("car",1), KeyValue("garage",2), KeyValue("elephant",2), KeyValue("kangeroo",0)];
  auto rnd = Random(9812);
  foreach(i; 0..10000) {
    auto ph = createPerfectHashFunction(keys, rnd, 10);
    foreach(key; keys)
      assert(ph(key.key) == key.value);
  }
}

auto extractG(const ref Vertex[] vertices) {
  import std.algorithm : map;
  return vertices.map!(v => v.value);
}

bool calcEdges(Range)(Range keys, const HashFunction[2] hashFunction, ref Edge[] edges, ref Vertex[] vertices) if (is(ElementType!Range == KeyValue)) {
  import std.algorithm : multiSort;
  assert(edges.length == 2*keys.length);
  foreach(idx, key; keys) {
    auto a = hashFunction[0](key.key);
    auto b = hashFunction[1](key.key);
    if (a == b) // each edge must point to 2 distinct vertices
      return false;
    edges[idx * 2].vertices[0] = a;
    edges[idx * 2].vertices[1] = b;
    edges[idx * 2].hash = key.value;
    edges[idx * 2 + 1].vertices[0] = b;
    edges[idx * 2 + 1].vertices[1] = a;
    edges[idx * 2 + 1].hash = key.value;
  }
  edges.multiSort!("a.vertices[0] < b.vertices[0]", "a.vertices[1] < b.vertices[1]");
  size_t pos = 0;
  foreach(idx, ref v; vertices) {
    if (idx < edges[pos].vertices[0]) {
      v.edgeIndex = size_t.max;
      v.value = 0;
    } else if (idx == edges[pos].vertices[0]) {
      v.value = 0;
      v.edgeIndex = pos;
      pos++;
      for(;;) {
        if (pos < edges.length) {
          if (edges[pos].vertices[0] != idx)
            break;
          else if (edges[pos].vertices[1] == edges[pos - 1].vertices[1]) // cannot have duplicated edges
            return false;
        } else {
          if (vertices.length > idx + 1)
            vertices[idx + 1 .. $] = Vertex(0, size_t.max);
          return true;
        }
        pos++;
      }
    }
  }
  return true;
}

@("calcEdges")
unittest {
  import std.random;
  import std.range : enumerate;
  import std.algorithm : map, maxElement;
  import std.array : array;
  import unit_threaded;
  import std.stdio;

  auto rnd = Random(9812);
  auto keys = [
    "aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "jj", "kk", "ll",
    "mm", "nn", "oo", "pp", "qq", "rr", "ss", "tt", "uu", "vv", "ww", "xx", "yy",
    "zz"
               ].enumerate.map!(e => KeyValue(e[1], cast(uint) e[0])).array();
  auto k = keys.length;
  auto n = k*3;
  auto maxKeyLength = keys.maxElement!(k => k.key.length).key.length;
  HashFunction[2] hashFunctions = [
    HashFunction(maxKeyLength), HashFunction(maxKeyLength)
  ];
  Vertex[] vertices = new Vertex[n];
  Edge[] edges = new Edge[k * 2];

  for(int x = 0;x < 100000;) {
    hashFunctions[0].randomize(n, 0, rnd);
    hashFunctions[1].randomize(n, 1, rnd);
  if (calcEdges(keys, hashFunctions, edges, vertices)) {
    foreach(idx, v; vertices) {
      if (v.edgeIndex == size_t.max)
        continue;
      edges[v.edgeIndex].vertices[0].shouldEqual(idx);
    }
    x++;
  }
  }
}

// since the edges are sorted by vertex indices, we can get all outgoing edges of a vertex
// by starting at the index in the edges list for the particular vertex and extract all edges
// as long as it starts from that vertex
auto outgoing(const ref Edge[] edges, const ref Vertex vertex) @safe nothrow {
  struct Impl {
    size_t idx;
    size_t vertexId;
    const Edge[] edge;
    bool empty() @safe nothrow {
      return idx >= edge.length || edge[idx].vertices[0] != vertexId;
    }
    auto front() @safe nothrow {
      return edge[idx];
    }
    void popFront() @safe nothrow {
      idx++;
    }
  }
  assert(vertex.edgeIndex != size_t.max);
  return Impl(vertex.edgeIndex, edges[vertex.edgeIndex].vertices[0], edges);
}

struct ToCheck {
  size_t vertexId;
  size_t parentId;
}

auto assumeNoThrow(T)(lazy T t) @trusted {
  try {
    return t();
  } catch (Exception e) {
    assert(0, e.message);
  }
}

@safe nothrow unittest {
  auto edges = [Edge([0,1],0), Edge([0,2],2), Edge([1,0],0), Edge([1,2],1), Edge([2,0],2), Edge([2,1],1)];
  auto vertices = [Vertex(0, 0), Vertex(0, 2), Vertex(0, 4)];
  size_t[] rawBits = new size_t[vertices.length / size_t.sizeof + 1];
  Appender!(ToCheck[]) toCheck;
  assert(assignVertexValuesOnlyWhenACyclic(edges, vertices, 3, rawBits, toCheck));
}

// returns true when cyclic
bool assignVertexValuesOnlyWhenACyclic(const ref Edge[] edges, ref Vertex[] vertices, size_t n, size_t[] rawBits, ref Appender!(ToCheck[]) toCheck) @trusted nothrow {
  import std.bitmanip : BitArray;
  import std.array : Appender;
  rawBits[] = 0;
  auto bitArray = BitArray(rawBits, rawBits.length * size_t.sizeof * 8);
  assert(bitArray.length >= vertices.length, "need as many bits as vertices");
  toCheck.shrinkTo(0).assumeNoThrow();
  foreach(idx, startVertex; vertices) {
    if (bitArray[idx])
      continue;
    startVertex.value = 0;
    if (startVertex.edgeIndex == size_t.max)
      continue;
    bitArray[idx] = true;
    foreach(edge; edges.outgoing(startVertex)) {
      auto neighbour = edge.vertices[1];
      toCheck.put(ToCheck(neighbour, idx));
      bitArray[neighbour] = true;
      auto diff = (n - startVertex.value + edge.hash) % n;
      vertices[neighbour].value = diff;
    }
    while (toCheck.data.length > 0) {
      auto check = toCheck.data[$-1];
      auto vertex = vertices[check.vertexId];
      toCheck.shrinkTo(toCheck.data.length - 1).assumeNoThrow();
      if (vertex.edgeIndex == size_t.max)
        continue;
      foreach(edge; edges.outgoing(vertex)) {
        assert(edge.vertices[0] == check.vertexId);
        auto neighbour = edge.vertices[1];
        if (neighbour == check.parentId)
          continue;
        if (bitArray[neighbour])
          return true;
        auto diff = (n - vertex.value + edge.hash) % n;
        vertices[neighbour].value = diff;
        bitArray[neighbour] = true;
        toCheck.put(ToCheck(neighbour, check.vertexId));
      }
    }
  }
  return false;
}

struct HashFunction {
  size_t[] coeff;
  size_t n;
  size_t o;
  size_t opCall(string key) const @safe nothrow pure {
    size_t t;
    foreach(idx, c; key)
      t += (c+o) * coeff[idx % coeff.length];
    return t % n;
  }
  this(size_t k) @safe nothrow {
    coeff.length = k;
  }
  void randomize(Random)(size_t n, size_t o, auto ref Random rnd) @safe nothrow {
    import std.random : uniform;
    this.n = n;
    this.o = o;
    foreach(ref v; coeff[0..coeff.length])
      v = uniform(0, n, rnd).assumeNoThrow();
  }
  auto clone() {
    auto clone = this;
    clone.coeff = coeff.dup();
    return clone;
  }
}
