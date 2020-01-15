module perfecthash.generator;
import std.random;
import std.range : ElementType;
import std.array : Appender;

struct Input {
  string[] keys;
  size_t maxKeyLength;
  size_t n;
}

struct State {
  string[] keys;
  size_t maxKeyLength;
  size_t n;
  Vertex[] vertices;
  Edge[] edges;
  HashFunction[2] hashFunction;
}

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
  this(const(size_t)[] G, HashFunction[2] hashFunction) {
    this.G = G;
    this.hashFunction = hashFunction;
  }
  size_t opCall(string key) {
    return (G[hashFunction[0](key)] + G[hashFunction[1](key)]) % G.length;
  }
}

string toDModule(ref PerfectHashFunction fun, string symbolName, string moduleName = "") {
  import std.conv : to;
  import std.string : replace;
  size_t N = fun.G.length;
  size_t coeffLen = fun.hashFunction[0].coeff.length;
  string NType = N > ushort.max ? "uint" : N > ubyte.max ? "ushort" : "ubyte";
  auto func = q{auto $symbolName(string key) {
      static $NType hash(alias coeff)(string key) {
        size_t t;
        foreach(idx, c; key)
          t += c * coeff[idx % $coeffLen];
        return t % $N;
      }
      static $NType[$N] G = $G;
      static $NType[$coeffLen] coeffA = $coeffA;
      static $NType[$coeffLen] coeffB = $coeffB;
      return (G[hash!(coeffA)(key)] + G[hash!(coeffB)(key)]) % $N;
    }
  }.replace("$NType", NType).replace("$moduleName", moduleName).replace("$symbolName", symbolName).replace("$N", N.to!string).replace("$coeffLen", coeffLen.to!string).replace("$G", fun.G.to!string).replace("$coeffA", fun.hashFunction[0].coeff.to!string).replace("$coeffB", fun.hashFunction[1].coeff.to!string);
  if (moduleName.length > 0)
    return "module $moduleName;\n".replace("$moduleName", moduleName) ~ func;
  return func;
}

unittest {
  import std.random;
  auto keys = ["asdf", "fdsa", "car", "garage", "elephant", "kangeroo"];
  auto rnd = Random(9123);
  foreach(i; 0..10000) {
    auto ph = createPerfectHashFunction(keys, rnd, 10);
    foreach(idx, key; keys)
      assert(ph(key) == idx);
  }
}

unittest {
  import std.random;
  auto keys = [KeyValue("asdf",0), KeyValue("fdsa",1), KeyValue("car",1), KeyValue("garage",2), KeyValue("elephant",2), KeyValue("kangeroo",0)];
  auto rnd = Random(9812);
  foreach(i; 0..10000) {
    auto ph = createPerfectHashFunction(keys, rnd, 10);
    foreach(key; keys)
      assert(ph(key.key) == key.value);
  }
}

import std.typecons : Tuple;
alias KeyValue = Tuple!(string, "key", uint, "value");

auto createPerfectHashFunction(Random)(string[] keys, Random rnd, ulong maxHnsecs = 1_000_000) {
  import std.range : enumerate;
  import std.algorithm : map;
  import std.array : array;
  return createPerfectHashFunction(keys.enumerate.map!(k => KeyValue(k.value, cast(uint)k.index)).array(), rnd, maxHnsecs);
}

auto createPerfectHashFunction(Range, Random)(Range keys, Random rnd, ulong maxHnsecs = 1_000_000, size_t startN = size_t.max) if (is(ElementType!Range == KeyValue)) {
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
  HashFunction[2] hashFunctions = [HashFunction(maxKeyLength, n), HashFunction(maxKeyLength, n)];
  PerfectHashFunction bestEffort;
  size_t[] rawBits = new size_t[vertices.length / size_t.sizeof + 1];
  Appender!(ToCheck[]) toCheck;
  for(;;) {
    if (Clock.currStdTime() - start > maxHnsecs && bestEffort.G.length != 0)
      break;
    hashFunctions[0].randomize(0, n, rnd);
    hashFunctions[1].randomize(0, n, rnd);
    if (!calcEdges(keys, hashFunctions, edges, vertices))
      continue;
    if (isCyclic(edges, vertices, rawBits, toCheck))
      continue;
    assignVertexValues(edges, vertices, n, rawBits, toCheck);
    bestEffort = PerfectHashFunction(extractG(vertices).array(), hashFunctions);
    if (n == k)
      break;
    if (Clock.currStdTime() - start > maxHnsecs)
      break;
    n = n - (max(1, (n - k) / 20));
    vertices = vertices[0 .. n];
    vertices[] = Vertex(0,0);
    hashFunctions = [HashFunction(maxKeyLength, n), HashFunction(maxKeyLength, n)];
  }
  return bestEffort;
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
    if (a == b)
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
    if (idx < edges[pos].vertices[0])
      v.edgeIndex = size_t.max;
    else if (idx == edges[pos].vertices[0]) {
      v.edgeIndex = pos;
      auto neighbour = edges[pos].vertices[1];
      pos++;
      for(;;) {
        if (pos < edges.length) {
          if (edges[pos].vertices[0] != idx)
            break;
          else if (edges[pos].vertices[1] == neighbour)
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

auto outgoing(const ref Edge[] edges, const ref Vertex vertex) {
  struct Impl {
    size_t idx;
    size_t vertexId;
    const Edge[] edge;
    bool empty() {
      return idx >= edge.length || edge[idx].vertices[0] != vertexId;
    }
    auto front() {
      return edge[idx];
    }
    void popFront() {
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

void assignVertexValues(const ref Edge[] edges, ref Vertex[] vertices, size_t n, size_t[] rawBits, ref Appender!(ToCheck[]) toCheck) {
  import std.bitmanip : BitArray;
  import std.array : Appender;
  rawBits[] = 0;
  auto bitArray = BitArray(rawBits, rawBits.length * size_t.sizeof * 8);
  assert(bitArray.length >= vertices.length);
  toCheck.shrinkTo(0);
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
      toCheck.shrinkTo(toCheck.data.length - 1);
      if (vertex.edgeIndex == size_t.max)
        continue;
      foreach(edge; edges.outgoing(vertex)) {
        auto neighbour = edge.vertices[1];
        if (neighbour == check.parentId)
          continue;
        auto diff = (n - vertex.value + edge.hash) % n;
        vertices[neighbour].value = diff;
        bitArray[neighbour] = true;
        toCheck.put(ToCheck(neighbour, check.vertexId));
      }
    }
  }
}

unittest {
  auto edges = [Edge([0,1],0), Edge([0,2],2), Edge([1,0],0), Edge([1,2],1), Edge([2,1],1), Edge([2,0],2)];
  auto vertices = [Vertex(0, 0), Vertex(0, 2), Vertex(0, 4)];
  size_t[] rawBits = new size_t[vertices.length / size_t.sizeof + 1];
  Appender!(ToCheck[]) toCheck;
  assert(isCyclic(edges, vertices, rawBits, toCheck));
}

bool isCyclic(const ref Edge[] edges, const ref Vertex[] vertices, size_t[] rawBits, ref Appender!(ToCheck[]) toCheck) {
  import std.bitmanip : BitArray;
  import std.array : Appender;
  rawBits[] = 0;
  auto bitArray = BitArray(rawBits, rawBits.length * size_t.sizeof * 8);
  assert(bitArray.length >= vertices.length);
  toCheck.shrinkTo(0);
  foreach(idx, startVertex; vertices) {
    if (bitArray[idx])
      continue;
    if (startVertex.edgeIndex == size_t.max)
      continue;
    bitArray[idx] = true;
    foreach(edge; edges.outgoing(startVertex)) {
      auto neighbour = edge.vertices[1];
      toCheck.put(ToCheck(neighbour, idx));
      bitArray[neighbour] = true;
    }
    while (toCheck.data.length > 0) {
      auto check = toCheck.data[$-1];
      auto vertex = vertices[check.vertexId];
      toCheck.shrinkTo(toCheck.data.length - 1);
      if (vertex.edgeIndex == size_t.max)
        continue;
      foreach(edge; edges.outgoing(vertex)) {
        auto neighbour = edge.vertices[1];
        if (neighbour == check.parentId)
          continue;
        if (bitArray[neighbour])
          return true;
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
  size_t opCall(string key) const @safe nothrow {
    size_t t;
    foreach(idx, c; key)
      t += c * coeff[idx % coeff.length];
    return t % n;
  }
  this(size_t k, size_t n) {
    coeff.length = k;
    this.n = n;
  }
  void randomize(Random)(size_t lower, size_t upper, auto ref Random rnd) @safe {
    import std.random : uniform;
    foreach(ref v; coeff[0..coeff.length])
      v = uniform(lower, upper, rnd);
  }
}
