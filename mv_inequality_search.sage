from math import prod
from collections import defaultdict

def mixed_volume(A: Polyhedron, B: Polyhedron):
    """
    Compute the mixed volume of two polyhedra A and B in R^2.
    """
    return 1/2 * ((A+B).volume() - A.volume() - B.volume())

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

body_types = Set([
    Polyhedron([(0,0), (1,0)]),
    Polyhedron([(0,0), (0,1)]) 
])

def compute_rays(n: int, mixed_volume=mixed_volume, body_types=body_types):
    """
    Find the rays in term value space produced by n bodies when the terms follow
    the form V(a,b)V(c,d) where a neq b neq c neq d.
    """
    terms = generate_inequality_terms(n)
    body_space = cartesian_product([body_types] * n)
    return list(filter(
        lambda x: sum(x) != 0, 
        map(
            lambda s: compute_inequality_terms(terms, s, mixed_volume=mixed_volume),
            body_space
        )
    ))

def compute_convex_hull(n: int, mixed_volume=mixed_volume, body_types=body_types):
    """
    Find convex hull of rays in term value space produced by n bodies when the terms follow
    the form V(a,b)V(c,d) where a neq b neq c neq d.
    """
    rays = compute_rays(n, mixed_volume=mixed_volume, body_types=body_types)
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

def get_orbit_of_relation(rel: Expression, term_group: PermutationGroup, x_vec: vector):
    """
    Compute the orbit of a relation under the action of the term group
    """

    # Encode the relation as a tuple of sets (functionally tuples of tuples)
    v_rel = list(rel.vector()[1:])
    d_rel = defaultdict(list)
    for i, coefficient in enumerate(v_rel):
        d_rel[coefficient].append(i+1)
    coefficient_values = list(d_rel.keys())
    encoded_rel = tuple(tuple(d_rel[key]) for key in d_rel)

    # Use GAP to calculate orbit
    o = term_group.orbit(encoded_rel, action="OnTuplesSets")

    # Decode the orbit elements back into relations
    return Set([
        sum(
            sum(coefficient_values[i] * x_vec[term_index-1] for term_index in term_group)
            for i, term_group in enumerate(orbit_element)
        )
        for orbit_element in o
    ])