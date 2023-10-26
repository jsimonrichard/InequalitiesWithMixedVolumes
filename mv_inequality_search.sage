from math import prod
from collections import defaultdict

def mixed_volume(A: Polyhedron, B: Polyhedron):
    """
    Compute the mixed volume (times 2) of two polyhedra A and B in R^2.
    """
    return ((A+B).volume() - A.volume() - B.volume())

def generate_inequality_terms(n: int, dimension=2, term_factors=3):
    """
    Generate all possible terms for {n} bodies in R^{dimension} with {term_factors} factors
    """
    if n != dimension * term_factors:
        raise ValueError("n must be equal to dimension*term_action so that the polynomial remains multi-linear")
    
    bodies = Set(range(1,n+1))

    term_stack = [[Set(), bodies]]
    final_terms = []

    while len(term_stack) > 0:
        term, remaining = term_stack.pop()
        if len(term) == term_factors:
            final_terms.append(term)
        else:
            for mv_group in Subsets(remaining, dimension):
                term_stack.append([term + Set([mv_group]), remaining - mv_group])
    
    return Set(final_terms)

def compute_inequality_terms(terms: Set, bodies: Set, mixed_volume=mixed_volume):
    """
    Compute the values of polynomial inequality terms
    """
    return tuple(
        prod(mixed_volume(*map(lambda i: bodies[i-1], mv_group)) for mv_group in term)
        for term in terms
    )

def generate_vector_combinations(n: int, term_factors=3):
    """
    Generate all required vector multiset combinations with appropriate multiplicity
    """
    max_multiplicity = term_factors # Otherwise all terms collapse to 0 (pigeonhole principle)
    min_multiplicity = 2 # Otherwise the single line segment can be transformed with linearity
    combinations = []
    for k in range(2, floor(n/2)+1):
        # k: number of ordered vectors
        # All vectors will have at least 2, so we'll add 2 after generating the combinations
        m = n - min_multiplicity*k
        multiset = list(range(k)) * (max_multiplicity-min_multiplicity)
        base = list(range(k)) * min_multiplicity
        combinations.extend(map(
            lambda c: base + c,
            Combinations(multiset, m).list()
        ))
    return combinations

def compute_rays(n: int, dimension=2, term_factors=3, mixed_volume=mixed_volume):
    """
    Find the rays in term value space produced by n bodies when the terms follow
    the form V(a_1,b_1)*...*V(c_{term_factors},d_{term_factors})
    """
    terms = generate_inequality_terms(n, dimension=dimension, term_factors=term_factors)
    vector_combinations = generate_vector_combinations(n, term_factors=term_factors)

    rays = []
    configurations = []
    for combination in vector_combinations:
        for arrangement in Arrangements(combination, n):
            rays.append(compute_inequality_terms(
                terms,
                list(map(lambda a: Polyhedron([(0,0), (1,a)]), arrangement)),
                mixed_volume=mixed_volume
            ))
            configurations.append(arrangement)
    return rays, configurations

def compute_convex_hull(n: int, mixed_volume=mixed_volume):
    """
    Find convex hull of rays in term value space produced by n bodies when the terms follow
    the form V(a,b)V(c,d) where a neq b neq c neq d.
    """
    rays, _ = compute_rays(n, mixed_volume=mixed_volume)
    pmv = Polyhedron(rays=rays)
    return pmv

# def generate_plucker_vectors(n: int, terms: Set):
#     """
#     Generate vectors of the coefficients for the corresponding 2D Plucker relations
#     over a set of n bodies. Coefficients are ordered according to the terms argument.
#     """
#     term_list = list(terms)
#     bodies = Set(range(1,n+1))
#     plucker_vectors = []
#     for i,j,k,l in Subsets(bodies, 4):
#         v = vector([0]*len(terms))
#         # Set the coefficients corresponding to the plucker relation
#         # with the indices i,j,k,l
#         v[term_list.index(Set([Set([i,j]), Set([k,l])]))] = -1
#         v[term_list.index(Set([Set([i,k]), Set([j,l])]))] = 1
#         v[term_list.index(Set([Set([i,l]), Set([j,k])]))] = -1
#         plucker_vectors.append(v)
#     return plucker_vectors

def permute_bodies(p: Permutation, term: Set):
    """
    Permute the bodies in a term according to the permutation p
    """
    return Set([Set(map(p, mv_group)) for mv_group in term])

def term_permutation_from_body_permutation(p: Permutation, term_list: list[Set]):
    """
    Compute the permutation of terms induced by a permutation of bodies
    """
    return Permutation([term_list.index(permute_bodies(p, term))+1 for term in term_list])

def get_orbit_of_relation(rel: vector, term_group: PermutationGroup):
    orbit = []
    for p in term_group:
        result = p.matrix() * rel
        if result not in orbit:
            orbit.append(result)
    return orbit
    
def get_all_orbits(rels: list[vector], term_group: PermutationGroup):
    orbits = []
    for rel in rels:
        not_in_orbit = True
        for orbit in orbits:
            if rel in orbit:
                not_in_orbit = False
                break
        if not_in_orbit:
            orbits.append(get_orbit_of_relation(rel, term_group))
    return orbits