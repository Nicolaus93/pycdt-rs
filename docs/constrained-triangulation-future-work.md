# Future constrained-triangulation work

The optimized inserter intentionally supports non-intersecting, non-overlapping
2D constraints over an existing point set. The following features require
changes to the data model and public API; they should not be approximated with
epsilon offsets or untracked topology edits.

## Intersecting constraints

A proper crossing must be constructed as a new stable point at the exact
intersection of the two supporting lines. Both physical constrained edges must
be removed and replaced by four edges incident to that point. The affected
triangles then need a local cavity retriangulation, reciprocal-neighbor repair,
incidence-map updates, and constrained Delaunay restoration outside the new
constraint edges.

One logical constraint can consequently become a chain of physical edges. The
triangulation must retain logical constraint IDs and an ordered chain per ID so
splitting an old constraint updates every owner instead of merely changing a
set of vertex pairs. Tests need multiple crossings, endpoint-on-edge
T-junctions, crossing order independence, exact rationally representable and
ill-conditioned intersections, neighbor reciprocity, chain ownership, and
Delaunay legality of every unconstrained edge.

## Overlapping constraints

Collinear overlap cannot be represented faithfully by the current
`HashSet<(usize, usize)>`: one physical edge may belong to several logical
constraints, and a partial overlap introduces split points at both overlap
boundaries. Supporting it requires edge-to-constraint ownership sets (or
reference counts), ordered chains for each logical constraint, exact projection
ordering, and atomic splitting of all owners. Removing or querying one logical
constraint must not discard an edge still owned by another.

Tests should cover equal reversed segments, containment, partial overlap on
either end, several coincident constraints, overlap combined with an existing
vertex, ownership removal, and insertion-order independence.

## Stable point identities and API changes

Intersection construction adds points after triangulation, so callers need a
stable point/vertex identity that survives internal storage changes. A future
API should return a constraint ID and expose its physical vertex chain, define
whether generated vertices are visible in `points`, and report unsupported or
invalid input with a structured error rather than a boolean. Serialization and
Python bindings must preserve generated-point and logical-constraint IDs.

An atomic insertion transaction is also needed: validate and construct all
split points and chains first, then commit topology and ownership together, or
restore the prior triangulation on any failure.

## Symbolic perturbation

Exact zero predicates currently select documented degenerate fallbacks.
Full symbolic perturbation needs a consistent simulation-of-simplicity ordering
based on stable point identities, applied uniformly to orientation, incircle,
fan selection, intersection ordering, and cavity retriangulation. Applying a
local tie-break in only one predicate can produce contradictory topology.

Tests must permute input and insertion order over collinear/cocircular fixtures,
verify deterministic physical topology and chains, and exercise generated
intersection points whose algebraic coordinates compare equal.

## Predicate caching

General orientation/incircle caching is deferred pending profiling. The current
one-orientation flip classifier carries only the combinatorial side labels
needed by a single insertion; it is not a global predicate cache. A future cache
must define keys for stable point identities (including generated points),
invalidation rules, memory bounds, and measured hit rates before adding its
complexity.
