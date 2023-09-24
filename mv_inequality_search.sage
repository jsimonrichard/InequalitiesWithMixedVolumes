def mixed_volume(A: Polyhedron, B: Polyhedron):
    """
    Compute the mixed volume of two polyhedra A and B in R^2.
    """
    return 1/2 * ((A+B).volume() - A.volume() - B.volume())

# Generate polynomial inequality terms
def generate_inequality_terms(n: int):
    """
    Generate all possible terms of polynomial inequality terms
    of the form V(a,b)V(c,d) where a neq b neq c \neq d are elements
    of an n-body set.
    """
    bodies = Set(range(1,n+1))
    term_list = []
    for mv_1 in Subsets(bodies, 2):
        if n - mv_1[0] < 2:
            # Not enough space for another mixed volume factor
            continue
        
        # Get subset of remaining bodies so that terms aren't repeated
        remaining = Set(range(mv_1[0]+1, n+1)) - [mv_1[1]]

        for mv_2 in Subsets(remaining, 2):
            term_list.append(Set([mv_1, mv_2]))
    
    assert len(term_list) == binomial(n, 2) * binomial(n-2, 2) / 2

    return Set(term_list)

def compute_inequality_terms(terms: Set, bodies: Set, mixed_volume=mixed_volume):
    """
    Compute the values of polynomial inequality terms of the form V(a,b)V(c,d)
    where a neq b neq c neq d provided a set of bodies.
    """
    return tuple([
        mixed_volume(bodies[term[0][0]-1], bodies[term[0][1]-1]) * mixed_volume(bodies[term[1][0]-1], bodies[term[1][1]-1])
        for term in terms
    ])

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

def generate_plucker_vectors(n: int, terms: Set):
    """
    Generate vectors of the coefficients for the corresponding 2D Plucker relations
    over a set of n bodies. Coefficients are ordered according to the terms argument.
    """
    term_list = list(terms)
    bodies = Set(range(1,n+1))
    plucker_vectors = []
    for i,j,k,l in Subsets(bodies, 4):
        v = vector([0]*len(terms))
        # Set the coefficients corresponding to the plucker relation
        # with the indices i,j,k,l
        v[term_list.index(Set([Set([i,j]), Set([k,l])]))] = -1
        v[term_list.index(Set([Set([i,k]), Set([j,l])]))] = 1
        v[term_list.index(Set([Set([i,l]), Set([j,k])]))] = -1
        plucker_vectors.append(v)
    return plucker_vectors

def term_action(p: Permutation, i: int, term_list):
    term = term_list[i-1]
    i,j,k,l = term[0][0], term[0][1], term[1][0], term[1][1]
    new_term = Set([Set([p(i), p(j)]), Set([p(k), p(l)])])
    return term_list.index(new_term)+1

def term_permutation(term_list, p: Permutation):
    return Permutation([term_action(p, i, term_list) for i in range(1,len(terms)+1)])