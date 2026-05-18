"""
Formula:
    int_{M_{0,K}} prod_i psi_i^{k_i}
        = sum_{K-admissible partitions P = {P_1,...,P_r} of {1,...,n}}
            (-1)^{n + l(P)} * int_{M_{0,l(P)}} prod_j psi_j^{k_{P_j} - |P_j| + 1}
Notes:
    A K-admissible partition is a set-partition of {1,...,n} where every block is a simplex of K
    k_{P_j} = sum_{i in P_j} k_i
    l(P) = r
"""

from fractions import Fraction # using this so that g > 0 could be brought in later
from math import factorial


# Build the simplicial complex
def _downward_closure(facets):
    """ Return the set of all simplices (as frozensets) given a list of facets """
    simplices = set()
    for facet in facets:
        facet = list(facet)
        # Enumerate all nonempty subsets of `facet` via bitmask
        k = len(facet)
        for mask in range(1, 1 << k):
            subset = frozenset(facet[i] for i in range(k) if mask & (1 << i))
            simplices.add(subset)
    return simplices


# Standard M_{0,n} intersection number
def _standard_intersection(exponents):
    """ int_{M_{0,n}} prod psi_j^{a_j} where exponents = (a_1, ..., a_n) """
    n = len(exponents)
    if n < 3:
        return Fraction(0)
    if any(a < 0 for a in exponents):
        return Fraction(0)
    if sum(exponents) != n - 3:
        return Fraction(0)

    # Multinomial coefficient (n - 3)! / prod a_j!
    num = factorial(n - 3)
    den = 1
    for a in exponents:
        den *= factorial(a)
    return Fraction(num, den)


# Enumerate K-admissible partitions of {1,...,n}
def _admissible_partitions(n, simplices):
    """
    Generate each K-admissible partition of {1,...,n} as a tuple of frozensets
    """
    universe = frozenset(range(1, n + 1))

    simplices_by_min = {}
    for s in simplices:
        m = min(s)
        simplices_by_min.setdefault(m, []).append(s)

    def recurse(remaining, blocks_so_far):
        if not remaining:
            yield tuple(blocks_so_far)
            return
        smallest = min(remaining)
        # Block must contain `smallest` and be a subset of `remaining`
        for block in simplices_by_min.get(smallest, ()):
            if block <= remaining:
                yield from recurse(remaining - block, blocks_so_far + [block])

    yield from recurse(universe, [])



# Primary function
def calc_psi_num(simplicial_complex, powers_vector):
    """
    simplicial_complex : list of lists
        The top-dimensional simplices (facets) of a simplicial complex K.
        Vertices should be labeled 1, ..., n. Order doesn't matter.
    powers_vector : sequence of int
        (k_1, ..., k_n)
    Returns a Fraction
    """
    powers_vector = list(powers_vector)
    n = len(powers_vector)

    # Validation checks:
    # Vertex set of K equals {1,...,n}.
    vertices_in_K = set()
    for facet in simplicial_complex:
        vertices_in_K.update(facet)
    if vertices_in_K != set(range(1, n + 1)):
        raise ValueError(
            f"Vertex set of K is {sorted(vertices_in_K)}, "
            f"expected {{1,...,{n}}} to match len(powers_vector)={n}."
        )

    # Dimension check: sum k_i must equal n - 3
    if sum(powers_vector) != n - 3:
        raise ValueError(
            f"sum(powers_vector) = {sum(powers_vector)} but n - 3 = {n - 3}; "
            f"the integrand has the wrong degree."
        )

    # Non-negative powers
    if any(k < 0 for k in powers_vector):
        raise ValueError("powers_vector must have non-negative integer entries.")

    # Check that there are no facets F1, F2 in K such that
    # F1 union F2 = {1,...,n}
    universe = set(range(1, n + 1))
    facets_as_sets = [set(f) for f in simplicial_complex]
    bad_pair = None
    for i, F1 in enumerate(facets_as_sets):
        # universe \ F1 is what F2 must (super-)cover.
        need = universe - F1
        for F2 in facets_as_sets[i:]:
            if need <= F2:
                bad_pair = (F1, F2)
                break
        if bad_pair is not None:
            break
    if bad_pair is not None:
        F1, F2 = bad_pair
        raise ValueError(
            f"For g = 0, K must not contain two facets whose union is "
            f"{{1,...,{n}}}; found facets {sorted(F1)} and {sorted(F2)}. "
        )

    # Build the complex
    simplices = _downward_closure(simplicial_complex)

    # Sum over K-admissible partitions
    total = Fraction(0)
    for partition in _admissible_partitions(n, simplices):
        ell = len(partition)
        sign = 1 if (n + ell) % 2 == 0 else -1
        exps = []
        for block in partition:
            k_block = sum(powers_vector[i - 1] for i in block)
            exps.append(k_block - len(block) + 1)
        total += sign * _standard_intersection(exps)

    return total


# Example calculation that gives a negative number
K = [[1],[2,3],[2,5],[4]]
powers_vector = [0,2,0,0,0]
print(calc_psi_num(K, powers_vector))

